# BSA 1.1.0

BSA3 runs from data on my scratch -- this is intended to show how to modify 
BSA6 to run on a different experiment design. However, BSA3 does not run from 
package data (and isn't intended to).

BSA6 is confimed to work from package data

## significant changes

the analyzer function and all that it wraps has been updated to `na.rm` -- the 
ramifications of this haven't been well thought through, but they should be.

## minor changes

`picker2` was updated so that input is an already sample-filtered dataframe

## new deps

futile.logger was added

# BSA 1.0.0

The documentation at this point should be minimally helpful, and much of the 
template has been removed from the README. This serves as the first 'official' 
release, though there will be updates to at least the documentation shortly.

# BSA 0.0.0

NEW FEATURES

* Added a `NEWS.md` file to track changes to the package.

SIGNIFICANT USER-VISIBLE CHANGES

* Your main changes to a function `foo()` or parameter `param`.

BUG FIXES

* Your bug fixes. See more details at
<http://bioconductor.org/developers/package-guidelines/#news>.
