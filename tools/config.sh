#!/bin/sh
# Rust build configuration (plan/06 §6.9), run by ./configure ("unix") and
# ./configure.win ("windows") from the package root:
#   - finds cargo on PATH, then in ~/.cargo/bin (Rust is never installed here);
#   - reports the cargo and rustc versions in the install log (CRAN policy);
#   - stops if rustc is older than the minimum version (R-D7);
#   - writes src/Makevars (src/Makevars.win) from the .in template.
#
# PSP_CARGO_TARGET_DIR, if set, is cargo's target directory (default
# src/rust/target). Set it during development to keep cargo's output off
# synced folders and to build incrementally; cargo then also uses the user's
# own CARGO_HOME instead of one inside the build tree.

PLATFORM="$1"
RUST_MIN_MAJOR=1
RUST_MIN_MINOR=84

fail() {
  echo "-------------------------- ERROR: CONFIGURATION FAILED --------------------------"
  echo ""
  echo "$1"
  echo ""
  echo "panelspatialprobit needs Rust (cargo and rustc >= ${RUST_MIN_MAJOR}.${RUST_MIN_MINOR}) to build"
  echo "from source. Install it with rustup (https://rustup.rs) or the system package"
  echo "manager, for example 'apt install cargo' or 'brew install rust'."
  echo "----------------------------------------------------------------------------------"
  exit 1
}

CARGO="$(command -v cargo 2>/dev/null)"
if [ -z "${CARGO}" ] && [ -x "${HOME}/.cargo/bin/cargo" ]; then
  CARGO="${HOME}/.cargo/bin/cargo"
fi
if [ -z "${CARGO}" ]; then
  CARGO_EXE="$(command -v cargo.exe 2>/dev/null)"
  [ -n "${CARGO_EXE}" ] && CARGO="${CARGO_EXE}"
fi
[ -n "${CARGO}" ] || fail "cargo was not found on PATH or in ~/.cargo/bin."

# rustc is looked up next to cargo first, as the build does.
CARGO_DIR="$(dirname "${CARGO}")"
PATH="${CARGO_DIR}:${PATH}"
export PATH

CARGO_VERSION="$(cargo --version 2>/dev/null)"
RUSTC_VERSION="$(rustc --version 2>/dev/null)"
[ -n "${RUSTC_VERSION}" ] || fail "rustc was not found (cargo is at ${CARGO})."
echo "using Rust package manager: '${CARGO_VERSION}'"
echo "using Rust compiler: '${RUSTC_VERSION}'"

VERSION="$(echo "${RUSTC_VERSION}" | sed -n 's/^rustc \([0-9][0-9]*\)\.\([0-9][0-9]*\).*/\1 \2/p')"
[ -n "${VERSION}" ] || fail "Could not read the rustc version from '${RUSTC_VERSION}'."
set -- ${VERSION}
if [ "$1" -lt "${RUST_MIN_MAJOR}" ] ||
   { [ "$1" -eq "${RUST_MIN_MAJOR}" ] && [ "$2" -lt "${RUST_MIN_MINOR}" ]; }; then
  fail "rustc $1.$2 is too old."
fi

if [ -n "${PSP_CARGO_TARGET_DIR}" ]; then
  TARGET_DIR="${PSP_CARGO_TARGET_DIR}"
  CARGO_HOME_SETUP=""
  echo "using cargo target directory: '${TARGET_DIR}'"
else
  # Paths are relative to src/, where make runs. Packages must not write to
  # the home directory while installing, so cargo's home is in the build tree.
  TARGET_DIR="rust/target"
  CARGO_HOME_SETUP='CARGO_HOME=.cargo; export CARGO_HOME;'
fi

# A tarball made after tools/vendor.sh carries the crates: build offline from
# them, with cargo's home in the build tree, and remove them after linking.
if [ -f src/rust/vendor.tar.xz ]; then
  echo "using vendored Rust crates (offline build)"
  VENDOR_STEP='"$(TAR)" -xf rust/vendor.tar.xz -C rust; mkdir -p .cargo; cp rust/vendor-config.toml .cargo/config.toml;'
  CARGO_HOME_SETUP='CARGO_HOME=.cargo; export CARGO_HOME;'
  OFFLINE="--offline"
  CLEAN_VENDOR="rm -rf rust/vendor .cargo"
else
  VENDOR_STEP=""
  OFFLINE=""
  CLEAN_VENDOR=""
fi

BEFORE_CARGO_BUILD=""
if [ "$(uname)" = "Darwin" ]; then
  # The cc crate (savvy's C helper) otherwise targets the SDK version, which
  # can be newer than the deployment target of R's compilers.
  MACOSX_DEPLOYMENT_TARGET="${MACOSX_DEPLOYMENT_TARGET:-11.0}"
  echo "using macOS deployment target: ${MACOSX_DEPLOYMENT_TARGET}"
  BEFORE_CARGO_BUILD="MACOSX_DEPLOYMENT_TARGET=${MACOSX_DEPLOYMENT_TARGET}; export MACOSX_DEPLOYMENT_TARGET;"
fi

if [ "${PLATFORM}" = "windows" ]; then
  CC="$("${R_HOME}/bin/R" CMD config CC)"
  case "$(${CC} -dumpmachine)" in
    aarch64*|arm64*) TARGET="aarch64-pc-windows-gnullvm" ;;
    *) TARGET="x86_64-pc-windows-gnu" ;;
  esac
  echo "using Rust target: '${TARGET}'"
  IN=src/Makevars.win.in
  OUT=src/Makevars.win
else
  TARGET=""
  IN=src/Makevars.in
  OUT=src/Makevars
fi

sed \
  -e "s|@CARGO_DIR@|${CARGO_DIR}|" \
  -e "s|@TARGET@|${TARGET}|" \
  -e "s|@TARGET_DIR@|${TARGET_DIR}|" \
  -e "s|@CARGO_HOME_SETUP@|${CARGO_HOME_SETUP}|" \
  -e "s|@BEFORE_CARGO_BUILD@|${BEFORE_CARGO_BUILD}|" \
  -e "s|@VENDOR_STEP@|${VENDOR_STEP}|" \
  -e "s|@OFFLINE@|${OFFLINE}|" \
  -e "s|@CLEAN_VENDOR@|${CLEAN_VENDOR}|" \
  "${IN}" > "${OUT}" || fail "Could not write ${OUT}."
