# this function calculates the p val of the genotype variable according to anova
# of two glmer fits or two glm fits
anova_test <- function(data, bin, test, silent = T, log_file, tail) {
  # check options
  if (test != "wG2" && test != "woG2") {
    report("e", "Unrecognized option for test parameter!", log_file)
  }
  
  # the link family depends on whether the response is binary or continuous
  if (bin == T) {
    fam <- "binomial"
  } else {
    fam <- "gaussian"
  }
  
  # wrap testing commands into a dummy function glmer/glm test, capture warnings
  # issued when glmer does not converge or when a continous response variable is
  # being used (glmer only)
  test_commands <- function() {
    if (test == "wG2") {
      if (length(unique(data$sex)) == 1) {
        fit <- glmer(pt ~ gt + (1 | mother), family = fam, nAGQ = 10, data = data)  # full model
        fit0 <- glmer(pt ~ (1 | mother), family = fam, nAGQ = 10, data = data)  # reduced model
      } else {
        fit <- glmer(pt ~ gt + sex + (1 | mother), family = fam, nAGQ = 10, 
          data = data)  # full model
        fit0 <- glmer(pt ~ sex + (1 | mother), family = fam, nAGQ = 10, data = data)  # reduced model
      }
      
      p.val <- anova(fit, fit0)$"Pr(>Chisq)"[2]
      if (is.na(p.val)) {
        p.val <- 1
      }
      direction <- fixef(fit)["gt"] > 0  # TRUE means protective effect, FALSE means harmful effect (desired)
    } else {
      if (length(unique(data$sex)) == 1) {
        fit <- glm(pt ~ gt, family = fam, data = data)  # full model
        fit0 <- glm(pt ~ 1, family = fam, data = data)  # reduced model
      } else {
        fit <- glm(pt ~ gt + sex, family = fam, data = data)  # full model
        fit0 <- glm(pt ~ sex, family = fam, data = data)  # reduced model
      }
      
      p.val <- anova(fit, fit0, test = "Chisq")[["Pr(>Chi)"]][2]
      if (is.na(p.val)) {
        p.val <- 1
      }
      direction <- coef(fit)["gt"] > 0
    }
    return(list(p.val = p.val, direction = direction))
  }
  
  # put the testing commands in one block
  pval_1tail <- tryCatch({
    tmp <- test_commands()
    convert_tail(tmp$direction, tmp$p.val, tail)  # if everything is ok
  }, warning = function(warn) {
    suppressWarnings((tmp <- test_commands()))  # run again, suppress warnings
    
    if (bin == F) {
      # fisher exact test
      
      fisher <- list(p.val = 0)
      # this will effectively avoid fisher's test when the response is continuous the
      # result will always be the one from glmer
    } else {
      fisher <- fisher.test(factor(data$gt), factor(data$pt))
    }
    
    # trust fisher or glmer?
    if (convert_tail(tmp$direction, fisher$p.val, tail)/convert_tail(tmp$direction, 
      tmp$p.val, tail) > 100) {
      report("m", "Fisher's exact test is used as a fallback", log_file)
      tmp$p.val <- fisher$p.val
    }
    
    return(convert_tail(tmp$direction, tmp$p.val, tail))
  }, error = function(err) {
    if (silent == F) {
      report("m", "Error occured!", log_file)
    }
    return(1)
  })
  
  return(pval_1tail)
} 
