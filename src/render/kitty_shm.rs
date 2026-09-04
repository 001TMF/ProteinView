//! Shared-memory transport for the Kitty graphics protocol.
//!
//! The escape-code transport has to base64 every pixel of every frame and push
//! it down the PTY, and because the protocol has no delta mechanism that is the
//! *whole* frame, every frame.  At a Retina-resolution viewport that is tens of
//! megabytes of pixel data per frame; compressing it first trades a large chunk
//! of the frame budget for a smaller — but still substantial — write, and the
//! terminal then has to decompress it again.
//!
//! Kitty's `t=s` transmission medium removes the transfer entirely.  The client
//! puts the pixels in a POSIX shared memory object and sends only its *name*;
//! the terminal maps the same pages and reads them directly.  The escape
//! sequence shrinks from hundreds of kilobytes to about sixty bytes, no
//! compression happens on either side, and the pixels are written exactly once
//! — straight into the shared mapping.
//!
//! This obviously only works when the terminal is on the same machine, so the
//! escape-code path in [`super::kitty_png`] remains for SSH sessions and for
//! terminals that speak the Kitty protocol without supporting `t=s`.  Support is
//! established once at startup by [`probe`] rather than assumed.
//!
//! # Object lifetime
//!
//! Per the protocol, *the terminal* unlinks and closes the object once it has
//! read it, so a fresh object is needed for every frame.  Names are therefore
//! cycled through a small ring of slots, and each slot is unlinked before being
//! recreated: a frame the terminal never consumed leaks one object until its
//! slot comes round again, and [`unlink_all`] clears whatever is left at exit.
//! Our own mapping can be dropped as soon as the pixels are written — the
//! object outlives it, and that is what the terminal opens.

use std::ffi::CString;

/// Number of shared-memory names cycled through, so the terminal is never
/// still reading the object we are about to replace.  At 30 fps this leaves
/// over a tenth of a second of slack, orders of magnitude more than a terminal
/// needs to drain an escape sequence.
pub const SLOTS: u32 = 4;

/// Slot reserved for the startup capability probe, outside the display ring.
const PROBE_SLOT: u32 = SLOTS;

/// Image id used for the probe.  Distinct from the display image id so a
/// query can never disturb the image currently on screen.
const PROBE_IMAGE_ID: u32 = 7317;

/// POSIX shared memory object names are limited to 31 characters on macOS, so
/// this stays deliberately terse: `/pv<pid>-<slot>`.
fn slot_name(slot: u32) -> String {
    format!("/pv{}-{}", std::process::id(), slot)
}

/// A freshly created shared memory object, mapped writable.
///
/// Dropping this unmaps our view but deliberately leaves the object in place:
/// the terminal opens it by name after it reads the escape sequence.
pub struct ShmFrame {
    name: String,
    ptr: *mut u8,
    len: usize,
}

impl ShmFrame {
    /// Create slot `slot`'s object at `len` bytes and map it writable.
    ///
    /// Any object left in the slot by an unconsumed earlier frame is unlinked
    /// first. Returns `None` if the platform refuses at any step, which the
    /// caller should treat as "fall back to the escape-code transport".
    pub fn create(slot: u32, len: usize) -> Option<Self> {
        if len == 0 {
            return None;
        }
        let name = slot_name(slot);
        let cname = CString::new(name.as_str()).ok()?;

        // SAFETY: `cname` is a valid NUL-terminated C string that outlives all
        // of these calls, `len` is non-zero, and every failure path releases
        // whatever the previous step acquired.
        unsafe {
            // Drop a stale object from a frame the terminal never read.
            libc::shm_unlink(cname.as_ptr());

            let fd = libc::shm_open(
                cname.as_ptr(),
                libc::O_CREAT | libc::O_EXCL | libc::O_RDWR,
                0o600 as libc::c_uint,
            );
            if fd < 0 {
                return None;
            }

            if libc::ftruncate(fd, len as libc::off_t) != 0 {
                libc::close(fd);
                libc::shm_unlink(cname.as_ptr());
                return None;
            }

            let ptr = libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ | libc::PROT_WRITE,
                libc::MAP_SHARED,
                fd,
                0,
            );
            // The mapping keeps the object alive; the descriptor is not needed.
            libc::close(fd);

            if ptr == libc::MAP_FAILED {
                libc::shm_unlink(cname.as_ptr());
                return None;
            }

            Some(Self {
                name,
                ptr: ptr.cast::<u8>(),
                len,
            })
        }
    }

    /// The object's name, to be sent as the escape sequence's payload.
    pub fn name(&self) -> &str {
        &self.name
    }

    /// The mapped bytes, for the caller to write pixels into.
    pub fn as_mut_slice(&mut self) -> &mut [u8] {
        // SAFETY: `ptr` came from a successful `mmap` of exactly `len` writable
        // bytes and is unmapped only in `Drop`, so the slice cannot outlive it.
        unsafe { std::slice::from_raw_parts_mut(self.ptr, self.len) }
    }
}

impl Drop for ShmFrame {
    fn drop(&mut self) {
        // SAFETY: `ptr`/`len` are exactly what `mmap` returned, and this runs
        // once because `Drop` runs once.  The object itself is intentionally
        // left linked for the terminal to open.
        unsafe {
            libc::munmap(self.ptr.cast::<libc::c_void>(), self.len);
        }
    }
}

/// Unlink every name this process may have created.
///
/// The terminal normally does this itself, so in the common case every call
/// here is a no-op; it exists so a frame that was dropped mid-flight — on
/// quit, or on a switch out of FullHD — cannot outlive the process.
pub fn unlink_all() {
    for slot in 0..=PROBE_SLOT {
        if let Ok(cname) = CString::new(slot_name(slot)) {
            // SAFETY: valid NUL-terminated string; unlinking an absent name
            // just returns an error we do not care about.
            unsafe {
                libc::shm_unlink(cname.as_ptr());
            }
        }
    }
}

/// Ask the terminal whether it can actually read pixels out of shared memory.
///
/// Uses the handshake the protocol documents for exactly this: a *query* action
/// — which validates without storing or displaying anything — followed by a
/// primary device attributes request.  A terminal that supports the graphics
/// protocol must answer the query before it answers the DA request, so a DA
/// response arriving alone is a definitive "no".  The DA response also bounds
/// the wait, with `timeout` as a backstop for terminals that answer neither.
///
/// Both responses are consumed here so they cannot surface later as spurious
/// key events.
pub fn probe(timeout: std::time::Duration) -> bool {
    let Some(mut frame) = ShmFrame::create(PROBE_SLOT, 4) else {
        return false;
    };
    frame.as_mut_slice().copy_from_slice(&[0, 0, 0, 255]);
    let name = frame.name().to_string();
    drop(frame);

    let supported = query_terminal(&name, timeout);

    // The terminal unlinks the object only if it read it; on the "no" path it
    // is still there.
    if let Ok(cname) = CString::new(name) {
        // SAFETY: valid NUL-terminated string.
        unsafe {
            libc::shm_unlink(cname.as_ptr());
        }
    }

    supported
}

/// Send the query + DA pair and classify the replies.
fn query_terminal(shm_name: &str, timeout: std::time::Duration) -> bool {
    use std::io::Write;

    let payload = base64_encode(shm_name.as_bytes());
    let query = format!("\x1b_Gi={PROBE_IMAGE_ID},a=q,f=32,t=s,s=1,v=1,S=4;{payload}\x1b\\\x1b[c");

    let mut stdout = std::io::stdout();
    if stdout.write_all(query.as_bytes()).is_err() || stdout.flush().is_err() {
        return false;
    }

    let deadline = std::time::Instant::now() + timeout;
    let mut buf = Vec::with_capacity(256);
    let mut chunk = [0u8; 128];

    loop {
        let remaining = deadline.saturating_duration_since(std::time::Instant::now());
        if remaining.is_zero() {
            return false;
        }
        match read_stdin_timeout(&mut chunk, remaining) {
            Some(0) | None => return false,
            Some(n) => buf.extend_from_slice(&chunk[..n]),
        }

        if let Some(ok) = classify(&buf) {
            return ok;
        }
        // Guard against a terminal that streams unrelated input at us.
        if buf.len() > 4096 {
            return false;
        }
    }
}

/// `Some(true)` once a graphics `OK` has arrived, `Some(false)` once the device
/// attributes reply has arrived without one, `None` while still undecided.
fn classify(buf: &[u8]) -> Option<bool> {
    if let Some(start) = find(buf, b"\x1b_G") {
        if let Some(end) = find(&buf[start..], b"\x1b\\") {
            let body = &buf[start..start + end];
            return Some(find(body, b";OK").is_some());
        }
        // Graphics reply started but has not finished; keep reading.
        return None;
    }
    // A DA reply is `ESC [ ? ... c`.
    let da = find(buf, b"\x1b[?")?;
    buf[da..].iter().position(|b| *b == b'c').map(|_| false)
}

fn find(haystack: &[u8], needle: &[u8]) -> Option<usize> {
    haystack
        .windows(needle.len())
        .position(|window| window == needle)
}

/// Read from stdin, waiting at most `timeout`.  `None` on error or timeout.
fn read_stdin_timeout(buf: &mut [u8], timeout: std::time::Duration) -> Option<usize> {
    let mut fds = libc::pollfd {
        fd: libc::STDIN_FILENO,
        events: libc::POLLIN,
        revents: 0,
    };
    let millis = timeout.as_millis().min(i32::MAX as u128) as libc::c_int;

    // SAFETY: a single well-formed `pollfd` describing stdin.
    let ready = unsafe { libc::poll(&mut fds, 1, millis) };
    if ready <= 0 {
        return None;
    }

    // SAFETY: `buf` is a valid writable slice of the length passed in.
    let n = unsafe {
        libc::read(
            libc::STDIN_FILENO,
            buf.as_mut_ptr().cast::<libc::c_void>(),
            buf.len(),
        )
    };
    if n < 0 { None } else { Some(n as usize) }
}

/// Standard base64, for the handful of bytes in a shared memory object name.
fn base64_encode(data: &[u8]) -> String {
    use base64::Engine;
    base64::engine::general_purpose::STANDARD.encode(data)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn slot_names_fit_the_posix_limit() {
        // macOS caps shared memory object names at 31 characters.
        for slot in 0..=PROBE_SLOT {
            let name = slot_name(slot);
            assert!(name.starts_with('/'), "{name} must be an absolute name");
            assert!(name.len() <= 31, "{name} is {} chars", name.len());
        }
    }

    #[test]
    fn create_maps_writable_memory_and_reads_back() {
        let len = 64 * 1024;
        let mut frame = ShmFrame::create(0, len).expect("shm unavailable");
        assert_eq!(frame.as_mut_slice().len(), len);
        frame.as_mut_slice()[0] = 0xAB;
        frame.as_mut_slice()[len - 1] = 0xCD;
        assert_eq!(frame.as_mut_slice()[0], 0xAB);
        assert_eq!(frame.as_mut_slice()[len - 1], 0xCD);
        drop(frame);
        unlink_all();
    }

    #[test]
    fn recreating_a_slot_replaces_a_stale_object() {
        // Simulates a frame the terminal never consumed: the slot is still
        // linked when the ring comes back round to it.
        let first = ShmFrame::create(1, 4096).expect("shm unavailable");
        let name = first.name().to_string();
        drop(first);
        let second = ShmFrame::create(1, 8192).expect("stale slot blocked reuse");
        assert_eq!(second.name(), name);
        drop(second);
        unlink_all();
    }

    /// The whole transport rests on this: after we unmap, a *separate* opener
    /// -- the terminal -- can still find the object by name and read the
    /// pixels we wrote.
    #[test]
    fn another_opener_reads_the_pixels_after_we_unmap() {
        let len = 4096;
        let mut frame = ShmFrame::create(3, len).expect("shm unavailable");
        let name = frame.name().to_string();
        for (i, byte) in frame.as_mut_slice().iter_mut().enumerate() {
            *byte = (i % 256) as u8;
        }
        drop(frame); // exactly what the render path does before the terminal reads

        let cname = CString::new(name.as_str()).unwrap();
        // SAFETY: mirrors what the terminal does -- open the named object
        // read-only, map it, read it, then unlink and unmap.
        unsafe {
            let fd = libc::shm_open(cname.as_ptr(), libc::O_RDONLY, 0 as libc::c_uint);
            assert!(fd >= 0, "reopening {name} failed");
            let ptr = libc::mmap(
                std::ptr::null_mut(),
                len,
                libc::PROT_READ,
                libc::MAP_SHARED,
                fd,
                0,
            );
            libc::close(fd);
            assert_ne!(ptr, libc::MAP_FAILED, "mapping {name} failed");

            let seen = std::slice::from_raw_parts(ptr.cast::<u8>(), len);
            let expected: Vec<u8> = (0..len).map(|i| (i % 256) as u8).collect();
            assert_eq!(
                seen,
                expected.as_slice(),
                "pixels did not survive unmapping"
            );

            libc::munmap(ptr, len);
            libc::shm_unlink(cname.as_ptr());
        }
    }

    #[test]
    fn zero_length_is_rejected() {
        assert!(ShmFrame::create(2, 0).is_none());
    }

    #[test]
    fn classify_detects_ok_error_and_da() {
        assert_eq!(classify(b"\x1b_Gi=7317;OK\x1b\\"), Some(true));
        assert_eq!(classify(b"\x1b_Gi=7317;EBADF:no shm\x1b\\"), Some(false));
        // Graphics reply still arriving.
        assert_eq!(classify(b"\x1b_Gi=7317;O"), None);
        // DA alone means the graphics protocol went unanswered.
        assert_eq!(classify(b"\x1b[?62;c"), Some(false));
        // DA still arriving.
        assert_eq!(classify(b"\x1b[?62;"), None);
        assert_eq!(classify(b""), None);
    }

    #[test]
    fn classify_prefers_graphics_reply_over_da() {
        assert_eq!(classify(b"\x1b_Gi=7317;OK\x1b\\\x1b[?62;c"), Some(true));
    }
}
