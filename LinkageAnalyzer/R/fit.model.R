run.mixed.effect.alt.model <- function(isBinary, alt.model, pheno) {
  assign("last.warning", NULL, envir = baseenv())
  alt <- tryCatch(
      {
        if (isBinary) {
          alt <- glmer(as.formula(alt.model), data = pheno, family = "binomial")
        } else {
          alt <- lmer(as.formula(alt.model), data = pheno, REML = FALSE)
        }
      },
      warning = function(warn) {
        logwarn(paste0("Fit [ ", alt.model, " binary =", isBinary, " ] has warnings ", str(warn)))
        msg <- as.character(get("last.warning", baseenv()))
        return(list(returncode = 1, message = msg, warning = warn, isBinary = isBinary, alt.model = alt.model))
      },
      error = function(err) {
        logwarn(paste0("Fit [ ", alt.model, " binary =", isBinary, " ] failed ", str(err)))
        reportError(err)
        msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
        return(list(returncode = 1, message = msg, error = err, isBinary = isBinary, alt.model = alt.model))
      })
  alt
}

run.fixed.effect.alt.model <- function(isBinary, alt.model, pheno, pheno.name, tail) {
  fid <- numGT <- NULL ## bypass CRAN check
  tmp <- ddply(pheno, .(fid), function(x) {
    c(numGT = length(unique(x$gt)))})
  tmp <- subset(tmp, numGT == 1)$fid
  reduced.pheno <- subset (pheno, !fid %in% tmp )
  reduced.model <- str_replace(alt.model, "\\(1\\|fid\\)", "fid")
  if (length(unique(reduced.pheno$fid)) == 1) {
    reduced.model <- str_replace(reduced.model, "\\+ fid", "")
  }
  reduced.model <- str_replace(reduced.model, "\\(1\\|mother\\)", "mother")
  if (length(unique(reduced.pheno$mother)) == 1) {
    reduced.model <- str_replace(reduced.model, "\\+ mother", "")
  }
  reduced.pheno$mother <- factor(reduced.pheno$mother)

  ## check if the model can be fit
  model.fittable <- TRUE
  pval <- NA
  while (model.fittable) {
    ## check sample size
    if (nrow(reduced.pheno) == 0) {
      loginfo("Insufficient sample size to fit models, set pvalue as one")
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## check factor variables
    model.var <- str_split(reduced.model, " ")[[1]]
    model.var <- model.var [!model.var %in% c("", "~", "1", "+") ]
    model.var <- reduced.pheno[, model.var]
    model.var <- model.var[complete.cases(model.var), ]
    useless.var <- apply(model.var, 2, function(x) {length(unique(x)) == 1})

    if (any(useless.var)) {
      loginfo(paste0("Reduced model has useless variable:", names(model.var)[useless.var], ", set pvalue to one."))
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## check rank
    reduced.model.matrix <- model.matrix(as.formula(reduced.model), reduced.pheno)
    if ( qr(reduced.model.matrix)$rank < ncol(reduced.model.matrix)) {
      loginfo(paste0("Insufficient sample size to fit models, set pvalue to one."))
      pval <- 1
      model.fittable <- FALSE
      break
    }

    ## try alt model
    if (isBinary) {
      if (is.factor(reduced.pheno[, pheno.name])) {
        ## convert AFFECTED/UNAFFECTED to 0/1
        reduced.pheno[, pheno.name] <- as.numeric(reduced.pheno[, pheno.name]) - 1
      } else {
        reduced.pheno[, pheno.name] <- as.numeric(reduced.pheno[, pheno.name])
      }
      alt <- tryCatch({
        alt <- glm(as.formula(reduced.model), data = reduced.pheno, family = "binomial")
      }, warning = function(x) { x }, error = function(x) { x })

      ## when there are separation problems, use firth regression
      if ((inherits(alt, "warning") &&
             !is.null(alt$message) &&
               alt$message == "glm.fit: fitted probabilities numerically 0 or 1 occurred") ||
          (inherits(alt, "glm") && vcov(alt)["gt", "gt"] > 1000)) {
        logwarn("Refit using Firth regression")
        alt <- logistf(as.formula(reduced.model), data = reduced.pheno, family = "binomial")
        alt <- with(alt, cbind(coefficients, prob))
        stopifnot(colnames(alt)[2] == "prob")
        pval <- alt["gt", "prob"]
        direction <- alt["gt", "coefficients"] > 0
        if (!is.na(direction)) {
          pval <- convert_tail(direction, pval, tail)
          return(pval)
        } else {
          logwarn("Refit using logistf() failed, set pvalue to one.\n")
          pval <- 1
          return(pval)
        }
      }

      ## check error
      if (inherits(alt, "error") && !is.null(alt$message)) {
        logerror(alt$message)
      }
      if (! inherits(alt, "glm")) {
        pval <- 1
        return(pval)
      }
    } else {
      alt <- lm(as.formula(reduced.model), data = reduced.pheno)
      ## check error
      if (inherits(alt, "error") && !is.null(alt$message)) {
        logerror(alt$message)
      }
      if (! inherits(alt, "lm")) {
        pval <- 1
        return(pval)
      }
    }

    idx <- which(colnames(summary(alt)$coefficients) == "Pr(>|t|)" |
                 colnames(summary(alt)$coefficients) == "Pr(>|z|)")
    coef <- summary(alt)$coefficients
    if ("gt" %in% rownames(coef) ) {
      pval <- coef["gt", idx]
      direction <- coef(alt)["gt"] > 0
      if (!is.na(direction)) {
        pval <- convert_tail(direction, pval, tail)
      }
    } else{
      ## fitting failed, so set pval to one
      logwarn("Refit using lm()/glm()/logsitf() failed, set pvalue to one.\n")
      pval <- 1
    }
    break
  }
  pval
}

create.null.model <- function(pheno, pheno.name, test) {
  null.model <- sprintf("%s ~ 1 ", pheno.name)
  has.random.effect = FALSE
  isBinary = is.factor(pheno[,pheno.name])

  nFam <- length(unique(pheno$fid))
  if (nFam > 1) {
    null.model <- paste0(null.model, " + (1|fid) ")
    has.random.effect = TRUE
  }

  nSex <- length(unique(pheno$sex)) 
  if ( nSex > 1) {
    null.model <- paste0(null.model, " + sex")
  }

  if (test == "wG2") {
    if (all(is.na(pheno$mother))) {
      msg <- "Does not have mother info at all (wG2 model cannot work)!!"
      logerror(msg)
      return(list(returncode = 1, message = msg))
    }
    if (length(unique(pheno$mother))==1) {
      msg <- "Cannot model G2 effect as there is only one G2 mouse."
      logerror(msg)
      return(list(returncode = 1, message = msg))
    }
    null.model <- paste0(null.model, " + (1|mother)")
    has.random.effect = TRUE
  }

  msg <- sprintf("%d family, %d sex detected, %s model used",
                 nFam, nSex, ifelse(has.random.effect, "random-effect", "fixed-effect"))
  loginfo(msg)
    
  ret <- list(formula = null.model,
              has.random.effect = has.random.effect,
              isBinary = isBinary,
              pheno.name = pheno.name)
  class(ret) <- c(class(ret), "null.model")
  return(ret)
}

isResponseConstant <- function(model, pheno) {
  if (inherits(model, "null.model") && !is.null(model$pheno.name)) {
    pheno.name <- model$pheno.name
  } else {
    pheno.name <- str_extract(model, "^[^~ ]+")
  }
  
  y <- pheno[, pheno.name]
  if (length(y) == 0 || length(unique(y)) == 1) {
    return(TRUE)
  }
  return(FALSE)
}

countModelParam <- function(model, pheno) {
  formula <- as.formula(model$formula)
  has.random.effect <- model$has.random.effect
  stopifnot(!is.null(formula), !is.null(has.random.effect))
  
  if (has.random.effect) {
    num.fixed.eff <- ncol(model.matrix(nobars(formula), data = pheno))
    num.random.eff <- do.call(sum, lapply(findbars(formula),
                             function(x) {
                               eff <- (sub("1 \\| ", "", deparse(x)))
                               length(unique(pheno[, eff]))
                             }))
    return(num.fixed.eff + num.random.eff)
  }

  num.fixed.eff <- ncol(model.matrix(formula, data = pheno))
  return(num.fixed.eff)
}

#' if success, will return append "result" field to @param null.model
#' otherwise, will return a list including error informatino
fit.null.model <- function(null.model, pheno) {
  stopifnot(inherits(null.model, "null.model"))

  isBinary              <- null.model$isBinary
  has.random.effect     <- null.model$has.random.effect
  formula               <- null.model$formula
  
  logdebug("Fit null model: %s", formula)
  ## snapshot("fit.null.model", "fit.null.model.Rdata")
  if (isResponseConstant(null.model, pheno)) {
    return(list(returncode = 1, message = "Response is constant - cannot fit the model"))
  }
  if (grepl("mother", formula) && length(unique(pheno$mother)) == 1)  {
    return(list(returncode = 1, message = "Cannot access mother effect when input only contains one mother mouse"))
  }
  nSample <- nrow(pheno)
  nParam <- countModelParam(null.model, pheno)
  if (nSample <= nParam) {
    return(list(returncode = 1, message = "Sample size is smaller than free model parameters"))
  }
  if (isBinary) {
    loginfo("Phenotypes are treated as binary\n")
    null <- tryCatch(
        {
          if (has.random.effect) {
            glmer(as.formula(formula), data = pheno, family = "binomial")
          } else {
            glm(as.formula(formula), data = pheno, family = "binomial")
          }
        },
        error = function(err) {
          reportError(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          return(list(returncode = 1, message = msg, error = err))
        })
    if (!isSuccess(null)) {
      null.model$formula <- formula <- str_replace(formula, "\\+ sex", "")
      loginfo("Refit null model using less covariates: %s\n", formula)
      null <- tryCatch(
          {
            if (has.random.effect) {
              null <- glmer(as.formula(formula), data = pheno, family = "binomial")
            } else {
              null <- glm(as.formula(formula), data = pheno, family = "binomial")
            }
            null
          },
          error = function(err) {
            reportError(err)
            msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
            return(list(returncode = 1, message = msg, error = err))
          })
    }
  } else { ## qtl
    loginfo("Phenotypes are treated as continuous\n")
    null <- tryCatch(
        {
          if (has.random.effect) {
            null <- lmer(formula, data = pheno, REML = FALSE)
          } else {
            null <- lm(formula, data = pheno)
          }
          null
        },
        error = function(err) {
          reportError(err)
          msg <- ifelse(is.null(err[["message"]]), "UnknownError", err$message)
          return(list(returncode = 1, message = msg, error = err))
        })
  }
  
  if (!isSuccess(null)) {
    return(null)
  } else {
    null.model$result <- null
    return(null.model)
  }
}

check.mixed.effect.alt.model.converge <- function(alt) {
  if (.hasSlot(alt, "optinfo")) {
    converge =  is.null(alt@optinfo$conv$lme4$messages[[1]])
    if (!converge) {
      msg <- alt@optinfo$conv$lme4$messages[[1]]
      err <- list(message = msg)
      class(err) <- "error"
      logwarn(paste0("Model may suffer from convergence: ", msg))
      alt <- list(returncode = 1, message = msg, error = err)
      return(alt)
    }
  }
  return(alt)
}
