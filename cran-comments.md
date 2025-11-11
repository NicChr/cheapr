Updated to version 1.4.0

Regarding the CRAN check note:

This package makes use of the non-API call to DATAPTR in limited areas
where performance is crucial and its usage is deemed safe. The macros 
'UNSAFE_STRING_PTR' and 'UNSAFE_VECTOR_PTR' are aliases for DATAPTR and
are named as such to signal that these pointers are to be used with care.
I hope that this is not an issue for the purposes of publishing to CRAN.

We checked 3 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

* Checked and passed using rhub v2.0.1 in the following environments:

## Test environments
- R-hubv2 windows (R-devel)
- R-hubv2 linux (R-devel)
- R-hubv2 macos-arm64 (R-devel)
- R-hubv2 clang-asan
- R-hubv2 ubuntu-gcc12
- R-hubv2 ubuntu-clang

Additionally checked on win-builder.r-project.org:

- windows (R-devel)
- windows (R-release)
- windows (R-old release)

and R-CMD-check github actions:

- macos-latest (release)
- ubuntu-latest (devel)
- ubuntu-latest (oldrel-1)
- ubuntu-latest (release)
- windows-latest (release)

