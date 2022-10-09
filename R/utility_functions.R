`%nin%` <- purrr::negate(`%in%`)

#' Harmonic mean
#'
#' @description Calculates the harmonic mean of a vector of numbers.
#'
#' @param a Vector of values
#'
#' @return The harmonic mean of the values: 1 / mean(1/a).
#' 
#' @export
mean.harmonic <- function(a){
  1/mean(1/a) #compute the harmonic mean
}



make_adj_var_df <- function(vars, Q_or_g = "Q"){
  if(is.null(vars)){
    stop("No adjustment variables selected. This should be 'U' in this case... something wrong? See utility_functions.R")
  }
  adj_row <- data.frame(matrix(NA, nrow = 1, ncol = length(vars)))
  adj_row[1,] <- vars
  colnames(adj_row) <- paste0(Q_or_g,"Adj",1:length(vars))
  return(adj_row)
}


#' Generate set of working models for adaptive pre-specification
#'
#' @description Function to get candidate adjustment strategies (variables +
#'   algorithms) for estimating the outcome regression and the propensity score.
#' 
#' @param cand.vars Character vector of variable names that may be used in the
#'   adjustment set. Default is NULL, in which case it returns a single
#'   unadjusted working model.
#' @param cand.algos Character vector of algorithm possibilities. See 'Qform' in
#'   Stage2 documentation for supported options. Default is NULL, in which case
#'   only GLMs are returned.
#' @param data.adapt.complexity Either 'low' (default), 'med', or 'high'. 'low'
#'   returns GLMs with a single adjustment variable, while 'med' adds all
#'   cand.algos with the full cand.vars specification. 'high' includes all
#'   non-redundant subsets of algorithms and variables, and may lead to
#'   overfitting; use carefully.
#'
#' @return A data frame where each row has a candidate model: The first column
#'   is the algorithm and the second column is a list of variables to use with
#'   it.
#' @export
#'
#' @examples
get.cand.adj <- function(cand.vars, cand.algos = NULL, data.adapt.complexity){
  unadjusted <- cbind.data.frame(cand.algos = 'glm', cand.varset = 'U')
  nvars <- length(cand.vars)
  if(nvars == 1){
    if(cand.vars == "U"){
      return(unadjusted)
    } else {
      data.adapt.complexity <- "low"
    }
  }
  
  cand.vars <- cand.vars[cand.vars %nin% 'U'] # We manually add U, so ignore if they add it as well.
  
  if(is.null(cand.algos)) {
    ###### simple Adaptive Prespec - GLMs with one adjustment variable each.
    sets <- expand.grid(cand.algos = 'glm', cand.varset = union('U', cand.vars))
  } else {                  
    ###### fancy adaptive prespec with expanded algorithms
    simple_sets <- cbind.data.frame(cand.algos = "glm", cand.varset = cand.vars)
    if(nvars > 1){
      full_set <- cbind.data.frame(cand.algos = cand.algos, cand.varset = NA) %>%
        mutate(cand.varset = list(cand.vars))
      
      # Find all of the proper subsets of length > 1
      combos <- do.call(c, lapply(seq_along(cand.vars), combn, x = cand.vars, simplify = FALSE))
      not_listed_yet <- intersect(which(lapply(combos, length) != 1), which(lapply(combos, length) != nvars))
      combos <- combos[not_listed_yet]
      middle_sets <- expand.grid(cand.algos = cand.algos, cand.varset = combos)
    }
    
    # Merge relevant combos (removing some unneeded duplicates) depending on desired level of complexity    
    if(data.adapt.complexity == "high"){
      message("High complexity; APS considering all unique subsets of cand.adj with each cand.form")
      sets <- rbind(simple_sets, full_set, middle_sets,
                    unadjusted)
    } else if (data.adapt.complexity == "med") {
      message("Medium complexity; APS considering GLMs with one adjustment variable, plus the full set of cand.adj with each cand.form")
      sets <- rbind(simple_sets, full_set, unadjusted)
    } else {
      if(length(cand.algos > 1)){
        message("Low complexity being used to prevent overfitting; APS considering GLMs with one adjustment variable only")
      }
      sets <- rbind(simple_sets, unadjusted)
    }
  }
  return(sets)
}




#' Check data inputs for Stage2()
#'
#' @description Performs several checks and throws warnings before running any
#'   statistical analysis.
#'
#' @param data.input Data frame input.
#' @param verbose Logical indicating if you want a lot of details to be printed
#'   to the console.
#'
#' @return A data freame with the necessary structure for Stage2().
#' @export
#'
#' @examples
preprocess_and_check_data <- function(data.input, verbose){
  if("U" %nin% colnames(data.input)){ # Will now add dummy column "U" if not in dataset already
    data.input <- data.input %>% dplyr::mutate(U = 1)
  }
  if("Y" %nin% colnames(data.input)){
    stop("No Y column in data.input")
  }
  if("id" %nin% colnames(data.input)){
    warning("data.input does not have a cluster 'id' column, assuming each row is an independent observation",
            call. = F)
    data.input <- cbind.data.frame(data.input, id = 1:nrow(data.input))
  }
  if("alpha" %nin% colnames(data.input)){
    if(verbose) print("data.input does not have an 'alpha' column, assuming each row is weighted equally.")
    data.input <- data.input %>% dplyr::mutate(alpha = 1)
  }
  if("U" %in% colnames(data.input))
    if(sum(data.input$U == 1) != nrow(data.input)){
      warning("You have a variable 'U' in your dataset that is not a column of 1s. Replacing with a column of 1s. If you have a variable named U, please rename.")
      data.input <- data.input %>% dplyr::mutate(U = 1)
    }
  return(data.input)
}






#' Aggregate influence curves to the cluster level
#'
#' Aggregates individual-level influence curve (IC) estimates to the cluster
#' level. This allows for appropriate statistical inference that respects the
#' cluster as the independent unit, whether the target parameter is an
#' individual-level or cluster-level effect.
#'
#' @param recordIC The individual-level IC estimates
#' @param id Vector of cluster membership IDs
#'
#' @return A vector of aggregated IC values with length equal to the number of
#'   unique cluster IDs.
#' @export
#'
#' @examples
aggregate_IC <- function(recordIC, id) {
  if (is.null(id)) return(recordIC)
  aggregatedIC <- as.matrix(stats::aggregate(recordIC, list(id=id), sum)[, -1, drop = FALSE])
  num.records <- nrow(recordIC)
  num.clusters <- nrow(aggregatedIC)
  aggregatedIC <- aggregatedIC * num.clusters / num.records
  return(aggregatedIC)
}



#' Calculate confidence intervals and p-values on relative or absolute scale
#'
#' @param goal String specifying the scale of the target parameter. Default is
#'   \code{RD}, risk/rate difference. Any other values assume that input values
#'   are given on the log scale, and the function will exponentiate the
#'   estimated target parameter and confidence interval bounds to output an
#'   artihmetic risk/rate ratio.
#' @param psi True value (if known) of the target parameter, for example, in a
#'   simulation study.
#' @param psi.hat Estimated value of the target parameter.
#' @param se Standard error of the estimated target parameter.
#' @param df Degrees of freedom for the Student's \emph{t} distribution as an
#'   approximation of the asymptotic normal distribution. If \code{df > 40}, the
#'   value is ignored and a normal distribution is used for inference.
#' @param sig.level Desired significance (alpha) level. Defaults to 0.05.
#' @param one.sided Logical indicating that a one-sided test is desired.
#'   Defaults to \code{FALSE}.
#' @param alt.smaller If one-sided test is desired, is the alternative
#'   hypothesis that the intervention arm will have a smaller value that the
#'   control arm? For example, if you expect a public health intervention to
#'   reduce the mean of a disease outcome, use alt.smaller = TRUE (this
#'   corresponds to a null hypothesis that the intervention did not reduce the
#'   mean disease outcome).
#'   
#' @return A one-row data frome with the estimated target parameter value
#'   (\code{est}), the (two-sided) confidence interval \code{CI.lo},
#'   \code{CI/hi}, the standard error of the estimate, the (possibly one-sided)
#'   p-value, and bias/coverage/rejection indicators (if true value or target
#'   parameter is supplied). NOTE: If \code{goal != "RD"}, the output standard
#'   error will be on the log scale.
#' @export
#'
#' @examples
get.inference <- function(goal = 'RD', psi = NA, psi.hat, se, df = 99, sig.level = 0.05, 
                          one.sided = F, alt.smaller = NULL){
  
  # test statistic
  # (on the log-transformed scale if goal is arithmetic RR or odds ratio)
  tstat <- psi.hat / se
  
  if(df > 40){
    # assume normal distribution
    cutoff <- stats::qnorm(sig.level/2, lower.tail=F)
    # one.sided hypothesis test 
    if(one.sided){
      pval <- stats::pnorm(tstat, lower.tail=alt.smaller) 
    } else{
      pval<- 2*stats::pnorm(abs(tstat), lower.tail=F) 
    }
  } else {
    # use Student's t-distribution
    # print('Using t-distribution')
    cutoff <- stats::qt(sig.level/2, df=df, lower.tail=F)
    # one.sided hypothesis test if specified
    if(one.sided){
      pval <- stats::pt(tstat, df=df, lower.tail= alt.smaller ) 
    } else{
      pval <- 2*stats::pt(abs(tstat), df=df, lower.tail=F)
    }
  }
  # 95% confidence interval 
  CI.lo <- (psi.hat - cutoff*se)
  CI.hi <- (psi.hat + cutoff*se)
  
  # If on relative (log) scale, transform back 
  if(goal != 'RD'){
    psi.hat <- exp(psi.hat)
    CI.lo <- exp(CI.lo)
    CI.hi <- exp(CI.hi)
  }  
  
  # bias
  bias <- (psi.hat - psi)
  # confidence interval coverage?
  cover <-  CI.lo <= psi & psi <= CI.hi
  # reject the null?
  reject <- as.numeric(pval < sig.level)
  return(data.frame(est = psi.hat, CI.lo, CI.hi, se = se, pval, bias, cover, reject))
}








#' Calculate the variance of the influence curve
#'
#' @description NOTE: UPDATE OF UNSCALING HAPPENS HERE
#'
#' @param goal One of 'aRR'  (arithmetic risk ratio), 'RD' (risk difference), or
#'   'OR' (odds ratio)
#' @param target Target of inference, either cluster-level ('clust') or pooled
#'   individivual-level effect ('indv')
#' @param Vdata data set
#' @param R1
#' @param R0
#' @param sample.effect Logical indicating whether to perform inference for the
#'   sample average treatment effect (if TRUE) or the population average
#'   treatment effect (if FALSE).
#' @param scale_value maximum value for outcome scaling
#' @param scale_value_min minimum value for outcome scaling
#' @param doing.CV
#'
#' @return on log scale for if goal='aRR' or 'OR' ... estimated IC & variance -
#'   preserving/breaking the match
#' @export
#'
#' @examples
get.IC.variance <- function(goal, target, Vdata, R1 = NA, R0 = NA,
                            scale_value = 1, scale_value_min = 0,
                            doing.CV = F,
                            verbose = F, sample.effect){
  
  # number of randomized units
  J <- length(unique(Vdata$id))
  
  # calculate the relevant components of the IC 
  if(sample.effect){
    # default - assume interest is in the sample effect
    DY1 <- Vdata$alpha*Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star)
    DY0 <- Vdata$alpha*Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star)
  } else {
    # calculate the IC for population effect (extra term for DW)
    DY1 <- Vdata$alpha*( Vdata$H.1W*(Vdata$Y - Vdata$Qbar1W.star) + Vdata$Qbar1W.star - R1 )
    DY0 <- Vdata$alpha*( Vdata$H.0W*(Vdata$Y - Vdata$Qbar0W.star) + Vdata$Qbar0W.star - R0 )	
  }
  
  # unscale 
  DY1 <- DY1*(scale_value - scale_value_min) + scale_value_min
  DY0 <- DY0*(scale_value - scale_value_min) + scale_value_min
  
  # if individual-level data, then need to aggregate the IC to the cluster-level 
  # approach for aggregation depends on the target effect
  if(length(DY1) > J) {
    if(target=='clust'){
      # Data are indv-level; target is cluster-level 
      if(!doing.CV & verbose) print('data = indv; target = clust')
      DY1 <- aggregate(DY1, by = list(Vdata$id), sum)[,-1]
      DY0 <- aggregate(DY0, by = list(Vdata$id), sum)[,-1]
    } else {
      # Data are indv-level; target is indv-level 
      if(!doing.CV & verbose) print('data = indv; target = indv')
      DY1 <- c(aggregate_IC(as.matrix(DY1), id = Vdata$id))
      DY0 <- c(aggregate_IC(as.matrix(DY0), id = Vdata$id))
    }
    # for the pair-matched IC also need to aggregate to the cluster-level
    # Vdata <- aggregate(Vdata, by=list(Vdata$id), mean)[,-1]
  } 
  
  # INFLUCENCE CURVES ARE NOW AT THE LEVEL OF THE RANDOMIZED UNIT
  if(goal=='RD'){
    # going after RD, easy IC
    DY <-  DY1 - DY0
  } else if (goal=='aRR'){ 
    # going after aRR, then get IC estimate on log scale
    #	i.e. Delta method for log(aRR) = log(R1) - log(R0)
    DY <- 1/R1*DY1 - 1/R0*DY0
  } else if(goal=='OR'){
    # Delta method for log(OR)
    DY <- 1/R1*DY1 + 1/(1-R1)*DY1 - 1/(1-R0)*DY0 - 1/R0*DY0
  }
  
  if(!doing.CV & verbose) print(paste0('Mean DY (Solve EIF): ', mean(DY) ))
  
  
  # estimated variance for txt specific means or if break the match	
  var.R1 <- stats::var(DY1) / J ## Doesn't work for LOOCV
  var.R0 <- stats::var(DY0) / J ## Doesn't work for LOOCV
  var.break <- stats::var(DY) / J
  #print(var.R1)
  
  if( 'pair' %in% colnames(Vdata) ){
    # estimated variance if preserve the match
    pairC <- aggregate(Vdata, by=list(Vdata$id), mean)[,'pair']
    pairs <- unique(pairC)
    n.pairs <- length(pairs)
    DY.paired <-  rep(NA, n.pairs)
    for(i in 1:n.pairs){		
      these<- pairC %in% pairs[i] 
      DY.paired[i]<- 0.5*sum(DY[ these] )			
    }
    
    var.pair <- stats::var(DY.paired) / n.pairs
  } else {
    DY.paired <- var.pair <- NA
  }
  
  return(list(R1=R1, R0=R0, DY1=DY1, DY0=DY0,var.R1 = var.R1, var.R0 = var.R0,
              DY=DY, var.break=var.break, 
              DY.paired=DY.paired, var.pair=var.pair))
}












