#!/bin/sh
# Vendors the Rust dependencies into src/rust/vendor.tar.xz (plan/06 §6.9),
# so the package builds offline. R CMD check warns when an install downloads
# crates, so run this before `R CMD build` (CI does). Run from the package
# root; needs cargo and python3. The tarball is not committed.
#
# Crates that no target compiles are replaced by empty stubs to keep the
# tarball small. windows-sys and windows-link are reachable only through
# sysctl, a dependency of faer's gemm-common on Apple targets, and then only
# on Windows (winapi-util), so they are never built.
set -e

STUBS="windows-sys windows-link"

if [ -d "${HOME}/.cargo/bin" ]; then
  PATH="${PATH}:${HOME}/.cargo/bin"
  export PATH
fi

cd src/rust
rm -rf vendor vendor.tar.xz
cargo vendor --quiet --locked --versioned-dirs vendor > /dev/null

for crate in ${STUBS}; do
  for dir in vendor/${crate}-[0-9]*; do
    [ -d "${dir}" ] || continue
    # Keep Cargo.toml (name, version, features) and an empty library; the
    # checksum file keeps the package hash that Cargo.lock records but lists
    # no files, so cargo does not check the removed ones.
    python3 - "${dir}" <<'EOF'
import json, os, shutil, sys
d = sys.argv[1]
ck = json.load(open(os.path.join(d, ".cargo-checksum.json")))
for name in os.listdir(d):
    if name not in ("Cargo.toml", ".cargo-checksum.json"):
        p = os.path.join(d, name)
        shutil.rmtree(p) if os.path.isdir(p) else os.remove(p)
toml = open(os.path.join(d, "Cargo.toml")).read()
# Drop sections that point at files that no longer exist.
out, skip = [], False
for line in toml.splitlines():
    if line.startswith("["):
        # [lib] stays: published manifests set autolib = false, and its path
        # is src/lib.rs, where the stub goes.
        skip = line.strip() == "[build-dependencies]" or line.startswith(("[[bench]]", "[[example]]", "[[test]]", "[[bin]]"))
    if line.startswith(("build =", "readme =", "include =", "exclude =")):
        continue
    if not skip:
        out.append(line)
open(os.path.join(d, "Cargo.toml"), "w").write("\n".join(out) + "\n")
os.makedirs(os.path.join(d, "src"))
open(os.path.join(d, "src", "lib.rs"), "w").write("#![no_std]\n// Stub: never compiled on any target the package is built for.\n")
json.dump({"files": {}, "package": ck["package"]}, open(os.path.join(d, ".cargo-checksum.json"), "w"))
EOF
    echo "stubbed ${dir}"
  done
done

tar -cJf vendor.tar.xz vendor
rm -rf vendor
ls -l vendor.tar.xz
