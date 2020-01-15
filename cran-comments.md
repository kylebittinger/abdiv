## Resubmission notes

Thank you to Martina Schmirl for the helpful comments on our initial
submission. We have improved the documentation and are resubmitting.

"Please always write package names, software names and API names in
single quotes in title and description. e.g: --> 'vegan'"

Done.

"If there are references describing the methods in your package, please
add these in the description field of your DESCRIPTION file in the form..."

Thank you, we have no such references.

"Please add \value to .Rd files regarding exported methods and explain
the functions results in the documentation.
e.g. kullback_leibler_divergence.Rd"

We have added \value{} to all exported functions.

"In some functions the return value is NaN, however, it would be better
to use stop() as it "stops execution of the current expression and
executes an error action.". e.g. kullback_leibler_divergence.Rd
Please change the code (or mention it in \value{} in the documentation)."

We have carefully considered whether to use stop() in various edge cases.
We call stop() when the type of input is not correct or when the input size is
wrong, for example in match_to_tree() and minkowski(). However, for many of
the functions here, the result is undefined even though the input is
considered valid. In this case, we return NaN and have documented it in
\value{}.

## Test environments
* local OS X install, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.1
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

## Downstream dependencies

There are currently no downstream dependencies for this package.
