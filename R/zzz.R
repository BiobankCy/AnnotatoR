.onAttach <- function(libname, pkgname) {
  version <- utils::packageVersion(pkgname)
  packageStartupMessage(
    sprintf("Welcome to %s, an R package for the creation of general-purpose,\nIonReporter-compatible annotation files starting from a list of gene names.", pkgname),
    sprintf("\n%s v.%s | Copyright (C) 2026 Dionysios Fanidis | GPL-3", pkgname, version),
    sprintf("\nType citation(\"%s\") for how to cite.", pkgname)
  )
}
