##  ====================================================================================================================================================
##  |  LinkageAnalyzer Copyright                                                                                                                       |
##  |  ----------------------------------------------------------------------------------------------------------------------------------------------  |
##  |  a.   Copyright ©2014, The University of Texas Southwestern Medical Center.  All rights reserved; and                                            |
##  |  b.   This software and any related documentation constitutes published and/or unpublished works and may contain valuable trade secrets and      |
##  |       proprietary information belonging to The University of Texas Southwestern Medical Center (UT SOUTHWESTERN).  None of the foregoing         |
##  |       material may be copied, duplicated or disclosed without the express written permission of UT SOUTHWESTERN.  IN NO EVENT SHALL UT           |
##  |       SOUTHWESTERN BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING   |
##  |       OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF UT SOUTHWESTERN HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.         |
##  |       UT SOUTHWESTERN SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND        |
##  |       FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". UT          |
##  |       SOUTHWESTERN HAS NO OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.                                   |
##  |  c.   This software contains copyrighted materials from R-package, ggplot2, gplots, gridExtra, lme4, logging, mclust, plyr and stringr.          |
##  |       Corresponding terms and conditions apply.                                                                                                  |
##  ====================================================================================================================================================    

##  ====================================================================================================================================================
##  |  This file is part of LinkageAnalyzer.													       |
##  |																		       |
##  |  LinkageAnalyzer is free software: you can redistribute it and/or modify									       |
##  |  it under the terms of the GNU General Public License as published by									       |
##  |  the Free Software Foundation, either version 3 of the License, or									       |
##  |  (at your option) any later version.													       |
##  |																		       |
##  |  LinkageAnalyzer is distributed in the hope that it will be useful,									       |
##  |  but WITHOUT ANY WARRANTY; without even the implied warranty of										       |
##  |  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the										       |
##  |  GNU General Public License for more details.												       |
##  |																		       |
##  |  You should have received a copy of the GNU General Public License									       |
##  |  along with LinkageAnalyzer.  If not, see <http://www.gnu.org/licenses/>.									       |
##  ====================================================================================================================================================


#' this function calculates the p val of the genotype variable according to anova
#' of two glmer fits or two glm fits
#' @keywords internal
anova_test <- function(data, bin, test, silent = T, log_file, tail,
                       null.model = NULL, fit.null = FALSE) {
  ## require(lme4)
  # check options
  if (test != "wG2" && test != "woG2") {
    report("e", "Unrecognized option for test parameter!", log_file)
  }
  if (FALSE) {
    null.model = NULL;
    fit.null = TRUE
  }
  # the link family depends on whether the response is binary or continuous
  if (bin == T) {
    fam <- "binomial"
  } else {
    fam <- "gaussian"
  }

  ## check if we fit null or fit alt
  if (!fit.null) {
    bFitAlt <- TRUE
  } else {
    bFitAlt <- FALSE
  }
  if (fit.null || is.null(null.model)) {
    bFitNull <- TRUE
  } else {
    bFitNull <- FALSE
  }

  # wrap testing commands into a dummy function glmer/glm test, capture warnings
  # issued when glmer does not converge or when a continous response variable is
  # being used (glmer only)
  test_commands <- function() {
    p.val <- NA
    direction <- NA
    fit0 <- null.model

    if (test == "wG2") {
      if (length(unique(data$sex)) == 1) {
        if (bFitAlt ) {
          # full model
          if (bin == TRUE) {
            fit <- glmer(pt ~ gt + (1 | mother), family = fam, nAGQ = 10, data = data)
          } else {
            fit <- lmer(pt ~ gt + (1 | mother), data = data, REML = FALSE)
          }
        }
        if (bFitNull) {
          # reduced model
          if (bin == TRUE) {
            fit0 <- glmer(pt ~ (1 | mother), family = fam, nAGQ = 10, data = data)
          } else {
            fit0 <- lmer(pt ~ (1 | mother), data = data, REML = FALSE)
          }
        }
      } else {
        if (bFitAlt ) {
          # full model
          if (bin == TRUE) {
            fit <- glmer(pt ~ gt + sex + (1 | mother), family = fam, nAGQ = 10, data = data)
          } else {
            fit <- lmer(pt ~ gt + sex + (1 | mother), data = data, REML = FALSE)
          }
        }
        if (bFitNull) {
          # reduced model
          if (bin == TRUE) {
            fit0 <- glmer(pt ~ sex + (1 | mother), family = fam, nAGQ = 10, data = data)
          } else {
            fit0 <- lmer(pt ~ sex + (1 | mother), data = data, REML = FALSE)
          }
        }
      }

      if (bFitAlt) {
        p.val <- anova(fit, fit0)$"Pr(>Chisq)"[2]
        if (is.na(p.val)) {
          p.val <- 1
        }
        direction <- fixef(fit)["gt"] > 0  # TRUE means protective effect, FALSE means harmful effect (desired)
      }
    } else { ## without G2
      if (length(unique(data$sex)) == 1) {
        ## not adjust for sex
        if (bFitAlt ) {
          # full model
          if (bin == TRUE) {
            fit <- glm(pt ~ gt, family = fam, data = data)
          } else {
            fit <- lm(pt ~ gt, data = data)
          }
        }
        if (bFitNull) {
          # reduced model
          if (bin == TRUE) {
            fit0 <- glm(pt ~ 1, family = fam, data = data)
          } else {
            fit0 <- lm(pt ~ 1, data = data)
          }
        }
      } else {
        ## adjust for sex
        if (bFitAlt ) {
          # full model
          if (bin == TRUE) {
            fit <- glm(pt ~ gt + sex, family = fam, data = data)
          } else {
            fit <- lm(pt ~ gt + sex, data = data)
          }
        }
        if (bFitNull) {
          # reduced model
          if (bin == TRUE) {
            fit0 <- glm(pt ~ sex, family = fam, data = data)
          } else {
            fit0 <- lm(pt ~ sex, data = data)
          }
        }
      }
      if (bFitAlt) {
        p.val <- anova(fit, fit0, test = "Chisq")[["Pr(>Chi)"]][2]
        if (is.na(p.val)) {
          p.val <- 1
        }
        direction <- coef(fit)["gt"] > 0
      }
    }

    return(list(pvalue = p.val, direction = direction, null.model = fit0))
  } ## end woG2

  # put the testing commands in one block
  assign("last.warning", NULL, envir = baseenv())
  pval_1tail <- tryCatch({
    tmp <- test_commands()
    if (fit.null) {
      ## fit null model, so pval is just empty
      pval <- NA
    } else {
      pval <- convert_tail(tmp$direction, tmp$pvalue, tail)  # if everything is ok
    }
    return(list(pvalue = pval, direction = tmp$direction, null.model = tmp$null.model))
  }, error = function(err) {
    if (silent == F) {
      if ("message" %in% err$message ) {
        m <- err$message
      } else {
        m <- "unknown"
      }
      msg <- sprintf("Error occured in ANOVA: %s", m)
      report("m", msg, log_file)
      reportError(err)
      ## cat("err=\n")
      ## print(err)
    }
    pval <- 1
    return(list(pvalue = pval, direction = NA, null.model = NA, error.occured = TRUE))
  })

  ## warning has occured during fitting, but error not occured
  if ( exists("last.warning", envir = baseenv()) &&
      !is.null(get("last.warning", baseenv())) &&
      ! "error.occured" %in% names(pval_1tail) ) {
    ## suppressWarnings((tmp <- test_commands()))  # run again, suppress warnings

    if (bin == F) {
      # fisher exact test

      fisher <- list(p.val = 0)
      # this will effectively avoid fisher's test when the response is continuous the
      # result will always be the one from glmer
    } else {
      fisher <- fisher.test(factor(data$gt), factor(data$pt))
    }

    # trust fisher or glmer?
    if (convert_tail(tmp$direction, fisher$p.val, tail)/
        convert_tail(tmp$direction, tmp$pvalue, tail) > 100) {
      report("m", "Fisher's exact test is used as a fallback", log_file)
      tmp$pvalue <- fisher$p.val
    }

    pval <- convert_tail(tmp$direction, tmp$pvalue, tail)
    return(list(pvalue = pval, direction = tmp$direction, null.model = tmp$fit0))
  }

  ## return(list(pvalue = pval_1tail))
  return(pval_1tail)
}
