library(devtools)
library(tidyverse)
library(earth)
`%nin%` <- purrr::negate(`%in%`)

check() # then...
# commit and push to github



# write a function, use roxygen skeleton to document it
# COMMON ROXYGEN TAGS
# @description @examples @examplesIf
# @export @family @inheritParams @param
# @returns @seealso @rdname
# To wrap blocks to be less then 80 characters wide: 
#  use Ctrl/Cmd + Shift + /
#  Or menu: code -> re-flow comment

# To other documentation:
#   \code{\link{function}}: function in this package.
# \code{\link[MASS]{abbey}}: function in another package.
# \link[=dest]{name}: link to dest, but show name.
# \code{\link[MASS:abbey]{name}}: link to function in another package, but show name.
# \linkS4class{abc}: link to an S4 class.
# To the web:
#   \url{http://rstudio.com}: a url.
# \href{http://rstudio.com}{Rstudio}:, a url with custom link text.
# \email{hadley@@rstudio.com} (note the doubled @): an email address.

# formatting with LATEX-style...
# \emph{will give italics}
# \strong{will give bold}
# \code{inline code format}
# \preformatted{for multiline code}
# \itemize{
    # \item First thing
    # \item Second thing
# }
# \tabular, etc.

# Special characters:
#   Use @@ to get @
#   Use \% to get %
#   Use \\ to get \

# use @inheritParams AnotherFunctionWithASharedArgument
# to copy parameters from that other function (not recursively tho)

# ... then...
devtools::document()
check()

################################################################## commit and push to GH again


##### TRY INSTALLING FROM GITHUB
library(devtools)
devtools::install_github(repo = "https://github.com/joshua-nugent/tmle4rcts")
library(tmle4rcts)
#library(tidyverse)
?test_me_out
?get.inference
?get.IC.variance
?Stage2
?do.TMLE
?do.data.adapt



################ OR...
# skip it and just install freshest version
install()
library(tmle4rcts)
dat <- simulate_clustered_data()
Stage2(data.input = dat, goal = "RD",
       psi = dat$unit_mean_diff[1], remove.pscore = F,
       do.data.adapt = T, verbose = F, sample.effect = F,
       cand.Qform = c("glm"),
       cand.gform = c("glm"),
       cand.QAdj = c("covar0", "covar1", "covar2"),
       cand.gAdj = c("covar0", "covar1", "covar2"))
Stage2(data.input = dat, goal = "aRR",
       psi = dat$unit_mean_diff[1], remove.pscore = T,
       do.data.adapt = T, verbose = F, sample.effect = F,
       cand.Qform = c("glm", "lasso"),
       cand.gform = c("glm","step","mars"),
       cand.QAdj = c("covar0", "covar1", "covar2"),
       cand.gAdj = c("covar0", "covar1", "covar2"))
head(dat)




#combs <- do.call(c, lapply(seq_along(ccs), combn, x = ccs, simplify = FALSE))
#combs
(combs <- do.call(c, lapply(seq_along(1:length(ccs)), combn, x = 1:length(ccs), simplify = FALSE)))
length(combs)
dim(expand.grid(forms = forms, combs = combs))
expand.grid(forms = forms, combs = combs)
forms[expand.grid(forms = forms, combs = combs)[[14,"forms"]]]
ccs[expand.grid(forms = forms, combs = combs)[[14,"combs"]]]


devtools::install_github(repo = "https://github.com/joshua-nugent/tmle4rcts")
library(tmle4rcts)
forms <- c("glm","step", "lasso")
ccs <- c("a", "b", "c")
get.cand.adj(cand.vars = ccs, cand.algos = forms, data.adapt.complexity = "high")
get.cand.adj(cand.vars = ccs, cand.algos = forms, data.adapt.complexity = "med")
get.cand.adj(cand.vars = ccs, cand.algos = forms, data.adapt.complexity = "low")




install()
dat <- simulate_clustered_data(treatment_arm_clusters = 10,
                               control_arm_clusters = 10)
library(tmle4rcts)
Stage2(data.input = dat, goal = "RD",
       psi = dat$unit_mean_diff[1], remove.pscore = F,
       do.data.adapt = T, verbose = F, sample.effect = F, do.cv.variance = T,
       cand.Qform = c("glm"), cand.gform = c("glm"),
       cand.QAdj = c("covar0", "covar2"),
       cand.gAdj = c("covar0", "covar2"))
Stage2(data.input = dat, goal = "RD",
       psi = dat$unit_mean_diff[1], remove.pscore = F,
       do.data.adapt = T, verbose = F, sample.effect = F, do.cv.variance = F,
       cand.Qform = c("glm"), cand.gform = c("glm"),
       cand.QAdj = c("covar0", "covar2"),
       cand.gAdj = c("covar0", "covar2"))

dat <- simulate_clustered_data(treatment_arm_clusters = 30,
                               control_arm_clusters = 30)
library(tmle4rcts)
Stage2(data.input = dat, goal = "RD",
       psi = dat$unit_mean_diff[1], remove.pscore = F,
       do.data.adapt = T, verbose = F, sample.effect = F, do.cv.variance = T,
       cand.Qform = c("glm"), cand.gform = c("glm"),
       cand.QAdj = c("covar0", "covar2"),
       cand.gAdj = c("covar0", "covar2"))

















