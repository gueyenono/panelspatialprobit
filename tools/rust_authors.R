# Writes inst/AUTHORS: the Rust crates compiled into the package, with their
# versions, authors and licences, from `cargo metadata` (plan/06 §6.9). CRAN
# asks that the authorship and copyright of bundled code be stated; DESCRIPTION
# refers to this file. Rerun after any change to src/rust/Cargo.lock.
# Needs cargo and jsonlite. Run from the package root:
#   Rscript tools/rust_authors.R
cargo <- Sys.which("cargo")
if (!nzchar(cargo)) cargo <- file.path(Sys.getenv("HOME"), ".cargo", "bin", "cargo")
# Normal dependencies of the pspr library on the platforms the package is built
# for; the stubbed windows-* crates (tools/vendor.sh) are never compiled.
meta <- jsonlite::fromJSON(system2(cargo, c("metadata", "--format-version", "1", "--locked",
                                            "--manifest-path", "src/rust/Cargo.toml"),
                                   stdout = TRUE), simplifyVector = FALSE)
own <- c("pspcore", "pspr")
stub <- c("windows-sys", "windows-link")
pkgs <- Filter(function(p) !(p$name %in% c(own, stub)), meta$packages)
pkgs <- pkgs[order(vapply(pkgs, function(p) tolower(p$name), ""))]
authors <- function(p) {
  a <- unlist(p$authors)
  if (length(a)) paste(sub(" <[^>]*>$", "", a), collapse = ", ") else "the crate's authors"
}
lines <- c(
  "Rust crates compiled into panelspatialprobit",
  "",
  "The package's Rust code (src/rust/pspcore, src/rust/pspr) links the crates",
  "below, which are distributed with the source package in",
  "src/rust/vendor.tar.xz. Each is copyright its authors, as listed, and used",
  "under the licence shown (an SPDX expression; where it offers a choice, the",
  "package uses the crate under the MIT licence). The windows-sys and",
  "windows-link crates are included only as empty stubs: no supported",
  "platform compiles them.",
  "",
  unlist(lapply(pkgs, function(p) {
    c(sprintf("%s %s", p$name, p$version),
      sprintf("  Authors: %s", authors(p)),
      sprintf("  Licence: %s", if (is.null(p$license)) "see the crate" else p$license),
      if (!is.null(p$repository)) sprintf("  Source: %s", p$repository),
      "")
  })))
writeLines(lines, "inst/AUTHORS")
message("Wrote inst/AUTHORS (", length(pkgs), " crates)")
