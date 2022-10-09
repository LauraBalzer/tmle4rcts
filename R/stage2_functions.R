#' Two-stage or single-stage RCT estimation
#'
#' Takes data frame input.data and calculates arm-specific and treatment effect estimates.
#' 
#' Outcomes can be bounded continuous (in which case they are scaled to [0,1]
#' before estimation and inference, then unscaled - See Ch7 TLB) or binary. Does
#' not work for multicategory outcomes.
#' 
#' The input.data object must include treatment indicator column 'A' and outcome
#' 'Y'. If the trial pair-matched and want to keep pairs (break.match = F), the
#' column with pair IDs must be labeled as 'pair'.
#' 
#' For cluster randomized trials, should also include weights 'alpha', and cluster id 'id'. 
#' 
#' column for weights (alpha)
#'   Set = 1 for individually randomized trials.
#'   For cluster randomized trials: 
#'     value of the weights depends on the target of inference and data level 
#'     Details in Benitez et al. https://arxiv.org/abs/2110.09633v2
#'     let J=number of clusters, N_j = cluster-specific sample size, N_tot = total # participants= sum_j N_j
#'       if target='clust' with cluster-level data, alpha=1 
#'       if target='clust' with indv-level data, alpha = 1/N_j
#'       if target='indv' with cluster-level data, alpha = J/N_tot*N_j
#'       if target='indv' with indv-level data, then alpha = 1
#'     for demonstration, see sim2.R in https://github.com/LauraBalzer/Comparing_CRT_Methods
#'   weights must sum to the total # of randomized units
#' 
#' For testing, or for a prespecified (not-adaptive) estimation approach,
#' specify conditional mean outcome adjustment variables (QAdj) and method
#' (Qform), and propensity score adjustment variables (gAdj) and method (gform).
#' Note that this is NOT recommended in general--- instead use Adaptive
#' Prespecification (do.data.adapt = T).
#'
#'
#' @param goal aRR = arithmetic risk ratio;  RD=risk difference; OR= odds ratio
#'   (not recommended)
#' @param target target of inference: cluster-level ("clust") or pooled-indv
#'   effect ("indv"). With appropriate weights, this can recover cluster-level
#'   effects from individual data, or individual-level effects from clustered
#'   data (see note above re: alpha weights).
#' @param data.input Observed data with treatment column A, outcome Y. See
#'   further specifications above.
#' @param QAdj Character string of variable names to force in to the adjustment
#'   set for the outcome regression. Default is NULL. If non-null, will override
#'   'cand.QAdj'. For an unadjusted outcome regression, specify QAdj = 'U'.
#' @param Qform String pre-specifying which algorithm to use to estimate the
#'   outcome regression. Currently implemented options: 'glm', 'lasso', 'mars'
#'   (multivariate adaptive regression splines), 'mars.corP' (MARS that
#'   pre-screens in variables correlated with the outcome), 'step' (stepwise
#'   selection with main terms), and 'step.interaction' (stepwise selection with
#'   interaction terms).
#' @param gAdj Character string of variable names to force in to the adjustment
#'   set for the propensity score. Default is NULL. If non-null, will override
#'   'cand.gAdj'. For an unadjusted outcome regression, specify gAdj = 'U'.
#' @param gform String pre-specifying which algorithm to use to estimate the
#'   propensity score. Same options as 'Qform'.
#' @param do.data.adapt Should adaptive pre-specification (Balzer 2016) be used?
#'   Defaults to FALSE.
#' @param data.adapt.complexity
#' @param cand.QAdj Character vector of candidate adjustment variable(s) for
#'   conditional mean outcome
#' @param cand.Qform Character vector with candidate adjustment approache(s) for
#'   conditional mean outcome. See 'Qform' for options.
#' @param cand.gAdj Character vector of candidate adjustment variable(s) for
#'   propensity score
#' @param cand.gform Character vector with candidate adjustment approache(s) for
#'   propensity score. See 'Qform' for options.
#' @param V Number of folds for cross-validation steps. Default is 5, once the
#'   number of independent units exceeds 40. For fewer than 40 independent
#'   units, defaults to LOOCV.
#' @param remove.pscore Relevant when do.data.adapt = T. Should adjustment
#'   variables that are selected for the outcome regression be removed from the
#'   candidate set for the propensity score? See 'Collaborative TMLE'. Recommend
#'   remove.pscore = TRUE when there are few independent units; defaults to TRUE
#'   if the number of independent units is under 40.
#' @param do.cv.variance Should cross-validated variance estimate be used for
#'   inference? Default FALSE.
#' @param sample.effect Is the target the sample effect among the observed
#'   units? If FALSE, target is population effect. Default is TRUE.
#' @param break.match If pair-matched trial, should pairs be ignored during
#'   inference?
#' @param one.sided Indicator specifying whether a one-sided p-value should be
#'   returned. Defaults to FALSE.
#' @param alt.smaller Only needed if one.sided = T. Specifies whether the
#'   alternative hypothesis is that the treatment arm has a lower level of the
#'   outcome. If the outcome is 'bad' (say, incidence of an infectious disease),
#'   and the treatment is hypothesized to be 'good' and reduce incidence,
#'   alt.smaller = T.
#' @param verbose Prints more information about what's happening under the
#'   hood.
#' @param psi In (say) a simulation study, if the true value of the treatment
#'   effect is known, it can be input here and the outcome will include
#'   an indicator for rejection of the null hypothesis and coverage. Default is
#'   NA.
#' @param return.IC indicator of whether to return the influence curve
#'
#' @return Point estimate and inference
#' @export
#'
#' @examples
Stage2 <- function(goal = 'aRR', target = 'indv', data.input, 
                   QAdj = NULL, Qform = 'glm',
                   gAdj = NULL, gform = 'glm',
                   do.data.adapt = F, data.adapt.complexity = NULL,
                   cand.QAdj = NULL, cand.Qform = 'glm',
                   cand.gAdj = NULL, cand.gform = 'glm',
                   V = 5, remove.pscore = NULL, do.cv.variance = F,
                   sample.effect = T,
                   break.match = T, one.sided = F, alt.smaller = NULL,
                   verbose = F, psi = NA,
                   return.IC = F){
  
  #=====================================================
  # INPUT CHECKS
  
  # if doing a one-sided test, need to specify the alternative
  # alt.smaller = T if intervention reduces mean outcome
  # alt.smaller = F if intervention increases mean outcome
  if(one.sided & is.null(alt.smaller)){
    stop('For one-sided test, need to specify the direction of the alternative hypothesis')
  }
  if(goal %nin% c("aRR", "RD", "OR")){
    stop("Please specify goal = 'aRR' for risk/rate ratio, 'RD' for risk difference, or 'OR' for odds ratio")
  }
  if(!is.null(QAdj) & !is.null(cand.QAdj) & do.data.adapt){
    stop("You have specified both forced-in adjustment variables (QAdj) for the outcome regression, as well as a candidate set to choose from (cand.QAdj). Please specify one only and leave the other as NULL.")
  }
  if(!is.null(gAdj) & !is.null(cand.gAdj) & do.data.adapt){
    stop("You have specified both forced-in adjustment variables (gAdj) for the propensity score, as well as a candidate set to choose from (cand.gAdj). Please specify one only and leave the other as NULL.")
  }
  # Make sure data has correct structure  
  data.input <- preprocess_and_check_data(data.input, verbose = verbose)

  
  
  
  #=====================================================
  # TRANSFORM the outcome as in Chpt7 of TLB if not bounded in [0,1]
  if(max(data.input[,'Y']) > 1){
    scale_value <- max(data.input[,'Y'])
  } else {
    scale_value <- 1
  }
  if(min(data.input[,'Y']) < 0){
    scale_value_min <- min(data.input[,'Y'])
  } else {
    scale_value_min <- 0
  }
  data.input[,'Y'] <- (data.input[,'Y'] - scale_value_min) / (scale_value - scale_value_min)

  
  

  #=====================================================
  # ADAPTIVE PRESPECIFICATION
  # update: flexibility in CV-scheme and candidate prediction algorithms
  if(do.data.adapt){
    select <- do.adaptive.prespec(goal = goal, target = target, break.match = break.match, 
                                  Ldata = data.input, V = V,
                                  cand.QAdj = cand.QAdj, cand.Qform = cand.Qform,
                                  cand.gAdj = cand.gAdj, cand.gform = cand.gform,
                                  remove.pscore = remove.pscore,
                                  QAdj = QAdj, gAdj = gAdj,
                                  scale_value = scale_value, scale_value_min = scale_value_min,
                                  verbose = verbose, sample.effect = sample.effect, data.adapt.complexity)
    
    Q.index <- select$Q.index
    QAdj <- select$QAdj
    Qform <- select$Qform
    g.index <- select$g.index
    gAdj <- select$gAdj	
    gform <- select$gform
    #if(verbose){print(paste(QAdj, gAdj))}
  } else {
    # QAdj <- gAdj <- 'U' # THIS WAS DUMB - OVERRODE ANY USER-SPECIFIED NON-NULL QAdj or gAdj
    Q.index <- g.index <- 1
  }
  
  # RUN FULL TMLE WITH ADJUSTMENT SET 
  # runs all code for point estimation on scaled outcome
  est <- do.TMLE(goal = goal, target = target, train =  data.input,
                 QAdj = QAdj, Qform = Qform, 
                 gAdj = gAdj, gform = gform,
                 scale_value = scale_value, scale_value_min = scale_value_min,
                 doing.CV = F, verbose = verbose, sample.effect = sample.effect)  

  n.clust <- length(unique(data.input$id)) 
  # Get point estimates of the treatment-specific mean
  R1 <- est$R1
  R0 <- est$R0

  
  # Now: for the intervention effect 
  #  the point estimate on the relevant scale for getting inference
  if(goal == 'aRR') {
    psi.hat <- log(R1/R0)
  } else if (goal == 'RD'){
    psi.hat <- R1 - R0
  } else if (goal == 'OR'){
    psi.hat <- log( R1/(1-R1) * (1-R0)/R0 )
  }
  
  if(break.match){
    # if breaking the match, set df to (#clusters - 2)
    df <- n.clust - 2
    var.hat <- est$var.break
  } else{
    # if preserving the match, set df to (#pairs-1)
    df <- length(unique(data.input$pair)) -1 
    var.hat <- est$var.pair
  }
  
  ############################## GET INFERENCE for arms and   
  if(do.cv.variance){
    #print(select)
    Txt.CV <- get.inference(psi.hat = R1, se = sqrt(select$var.CV.1), # need to change SE based on 'select' object
                            df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
    Con.CV <- get.inference(psi.hat = R0, se = sqrt(select$var.CV.0),
                            df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  }
  
  Txt <- get.inference(psi.hat = R1, se = sqrt(est$var.R1), df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  Con <- get.inference(psi.hat = R0, se = sqrt(est$var.R0), df = (n.clust-2))[, c('est','CI.lo','CI.hi','se')]
  
  inference <- get.inference(goal = goal, psi=psi, psi.hat = psi.hat,
                             se = sqrt(var.hat), df = df,
                             one.sided = one.sided, alt.smaller = alt.smaller)

  Q_adj_and_form_final <- make_adj_var_df(vars = QAdj, Q_or_g = "Q")
  g_adj_and_form_final <- make_adj_var_df(vars = gAdj, Q_or_g = "g")
  
  if(do.cv.variance){
    # if getting cross-validated inference
    inference.CV <- get.inference(goal=goal, psi=psi, psi.hat=psi.hat,
                                  se=sqrt(select$var.CV), df=df,
                                  one.sided=one.sided, alt.smaller = alt.smaller)
    
    est.df <- data.frame(Txt=Txt, TxtCV = Txt.CV, Con=Con, ConCV=Con.CV,  psi=psi, inference, CV=inference.CV,
                         Qform = est$Qform, #Q_adj_and_form_final,
                         gform = est$gform)#, #g_adj_and_form_final)
  } else {
    est.df <- data.frame(Txt = Txt, Con = Con, psi = psi, inference, 
                         Qform = est$Qform, #Q_adj_and_form_final,
                         gform = est$gform)#, g_adj_and_form_final)
  }
  
  if(return.IC){
    Stage2_output <- list(IC=est, est.df=est.df)
  } else{
    Stage2_output <- est.df
  }
  # Making the return a little weirder for the purposes of easy sim studies?
  return(list(output = Stage2_output,
              QAdj = QAdj, gAdj = gAdj,
              Qvars = paste0("Adjustment variable(s) for outcome regression: ",
                             paste(QAdj, sep = "", collapse = ", ")),
              Qform = paste0("Adjustment algorithm for outcome regression: ",
                             paste(est$Qform, sep = "", collapse = ", ")),
              gvars = paste0("Adjustment variable(s) for propensity score: ",
                             paste(gAdj, sep = "", collapse = ", ")),
              gform = paste0("Adjustment algorithm for propensity score: ",
                             paste(est$gform, sep = "", collapse = ", "))))
}










