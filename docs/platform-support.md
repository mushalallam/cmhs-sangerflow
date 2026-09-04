# Platform support and first-launch security

CMHS SangerFlow processes chromatograms locally. Its browser interface binds only to
`127.0.0.1` and uses a random session token. Platform security warnings concern the downloaded
executable's publisher identity; they do not mean that sequence data are being uploaded.

## Version 0.3 support matrix

| Platform | Package | Native build and test environment | Current limitation |
| --- | --- | --- | --- |
| Apple Silicon macOS | `macOS-Apple-Silicon.zip` | macOS 15 arm64 | Not Developer ID signed or notarized |
| Intel macOS | `macOS-Intel.zip` | macOS 15 x86-64 | Not Developer ID signed or notarized |
| Windows x86-64 | `Windows-x86_64.zip` | Windows Server 2025 x86-64 | Not Authenticode signed |
| Linux x86-64 | `Linux-x86_64.tar.gz` | Ubuntu 22.04 x86-64 | Requires glibc 2.35 or newer |

Windows ARM64, Linux ARM64, and musl-based distributions such as Alpine are not currently
packaged. The Python wheel remains an alternative on any supported platform with Python 3.10+
and compatible dependencies.

## Verify a release download

Every native archive has a matching `.sha256` sidecar. Both are GitHub release assets. The
digest GitHub displays beside the `.sha256` file is the digest of that sidecar itself; open or
download the sidecar to obtain the expected digest of the application archive.

```bash
# macOS
shasum -a 256 CMHS-SangerFlow-*.zip

# Linux
sha256sum CMHS-SangerFlow-*.tar.gz
```

```powershell
# Windows PowerShell
Get-FileHash .\CMHS-SangerFlow-*.zip -Algorithm SHA256
```

The value must exactly match the archive digest published on the release page or inside its
sidecar. Do not bypass an operating-system warning if it does not match.

## macOS Gatekeeper

Version 0.3 is ad-hoc signed by the packaging tool, not signed with an Apple Developer ID. A
downloaded copy may be blocked or may appear to freeze while Gatekeeper performs online checks.

First try **System Settings → Privacy & Security → Open Anyway**. If a verified download remains
stalled, stop it with **Control-C**, then remove quarantine only from the extracted SangerFlow
folder:

```bash
xattr -dr com.apple.quarantine "/path/to/extracted/CMHS-SangerFlow-folder"
```

This is an interim testing procedure, not the final distribution design. The permanent fix is a
Developer ID-signed application with hardened runtime, Apple notarization, and a stapled ticket.
Apple documents these requirements in [Signing your apps for Gatekeeper](https://developer.apple.com/developer-id/)
and [Notarizing macOS software before distribution](https://developer.apple.com/documentation/security/notarizing-macos-software-before-distribution).

## Windows SmartScreen and Smart App Control

The current executable can display **Windows protected your PC** or be blocked by managed-device
policy because it has no trusted publisher signature. After checking the SHA-256 value, a user
may select **More info → Run anyway** where organizational policy permits it.

The permanent direct-download improvement is consistent Authenticode/Artifact Signing with an
RFC 3161 timestamp. Microsoft notes that a newly signed application may still show SmartScreen
warnings while publisher or file reputation develops. Microsoft Store distribution is the most
reliable warning-free path. See Microsoft's [SmartScreen reputation guidance](https://learn.microsoft.com/windows/apps/package-and-deploy/smartscreen-reputation)
and [SignTool documentation](https://learn.microsoft.com/windows/win32/seccrypto/signtool).

## Linux

Linux does not have one universal Gatekeeper-equivalent download prompt. The tar archive keeps
the executable bit, but some extraction tools can remove it. Restore it with:

```bash
chmod +x sangerflow Start-CMHS-SangerFlow.sh
```

The current executable is built on Ubuntu 22.04 to avoid depending on a newer glibc. It is not a
universal Linux binary. Planned broader distribution includes x86-64 and ARM64 AppImages plus a
Debian/Ubuntu package, with the portable archive retained as a fallback.

## Release acceptance criteria

A future package may be described as warning-free only after the release workflow verifies:

- macOS Developer ID signature, hardened runtime, accepted notarization, and stapled ticket;
- Windows trusted Authenticode signature and RFC 3161 timestamp;
- native GUI startup checks on every advertised architecture;
- archive SHA-256 publication and installation testing from a fresh downloaded copy.
