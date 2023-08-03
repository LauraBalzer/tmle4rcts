



#' Implement adaptive pre-specification
#'
#' @description function to implement adaptive prespecification as described in
#'   Balzer et al. "Adaptive pre-specification in randomized trials with and
#'   without pair-matching"
#'
#' @param goal goal of analysis ('aRR', 'OR', 'RD')
#' @param target Target of estimation and inference: cluster-level ("clust") or
#'   pooled-indv effect ("indv")
#' @param break.match Logical indicating whether to break matched pairs (if they
#'   exist)
#' @param Ldata Input data set
#' @param V
#' @param cand.QAdj Character vector of names of candidate adjustment variables
#'   for the outcome regression
#' @param cand.Qform Character vector of names of candidate adjustment
#'   algorithms for the outcome regression
#' @param cand.gAdj Character vector of names of candidate adjustment variables
#'   for the propensity score
#' @param cand.gform Character vector of names of candidate adjustment
#'   algorithms for the propensity score
#' @param remove.pscore if T, remove the variable(s) selected for adjustment in
#'   the outcome regression from candidates for the pscore... should only be
#'   used if doing adaptive prespec in RCT with few indpt units
#' @param QAdj selected adjustment variable for the outcome regression, if
#'   specifying in advannce
#' @param gAdj selected adjustment variable for the pscore regression, if
#'   specifying in advance
#' @param scale_value maximum value for outcome scaling
#' @param scale_value_min minimum value for outcome scaling
#' @param verbose
#' @param sample.effect If the target of inference the sample effect or the
#'   population effect?
#' @param data.adapt.complexity Set to 'low', 'med', or 'high'. See
#'   ?get.cand.adj for details.
#'
#' @return Working model selection for candidate TMLE
#' @export
#' 
#' @examples
do.adaptive.prespec <- function(goal, target='indv', break.match=T, Ldata, V=5,
                                cand.QAdj, cand.Qform,
                                cand.gAdj, cand.gform,
                                remove.pscore, QAdj = NULL, gAdj = NULL,
                                scale_value, scale_value_min, verbose = F,
                                sample.effect, data.adapt.complexity){
  
  # get the indpt units (will be each observation in and individual-level RCT)
  if( !break.match ){
    Ldata$indpt.unit <- Ldata$pair
  } else {
    Ldata$indpt.unit <- Ldata$id
  }
  unique.ids <- unique(Ldata$indpt.unit)
  
  # get folds
  if(length(unique.ids) > 40){
    # V-fold CV
    folds <- get.folds(V=V, Y=Ldata$Y, ids=unique.ids)
  } else {
    # leave-one-out CV
    folds <- vector("list", length(unique.ids))
    for(v in seq(length(unique.ids))){
      folds[[v]] <- unique.ids[v]
    }
  }
  
  if(is.null(remove.pscore)){
    if(length(unique.ids) > 40){
      remove.pscore <- F
    } else {
      remove.pscore <- T
    }
    message(paste0("remove.pscore not speficied; setting to ", remove.pscore, " based on # of independent units"))
  }

  
  #=====================================================
  # GENERATE COVARIATE SET / WORKING MODEL COMBINATIONS FOR OUTCOME REGRESSION
  # Figure out complexity level first
  if(is.null(data.adapt.complexity)){
    nvars <- length(union(cand.gAdj, cand.QAdj))
    ratio <- length(unique.ids) / nvars
    low_cutoff <- 10
    med_cutoff <- 50
    if(ratio < low_cutoff){
      message(paste0("Fewer than ", low_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'low'"))
      data.adapt.complexity <- "low"
    } else if (ratio < med_cutoff){
      message(paste0("Between ", low_cutoff, " and ", med_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'med'"))
      data.adapt.complexity <- "med"
    } else {
      message(paste0("More than ", med_cutoff, " independent units per covariate. Setting data.adapt.complexity = 'high'"))
      data.adapt.complexity <- "high"
    }
  }
  cand.Q <- get.cand.adj(cand.vars = cand.QAdj, cand.algos = cand.Qform,
                         data.adapt.complexity = data.adapt.complexity)
  cand.QAdj <- cand.Q$cand.varset
  cand.Qform <- cand.Q$cand.algos


  
  
  #=====================================================
  # FIGURE OUT WHICH COMBINATION FOR OUTCOME REGRESSION IS BEST
  if( is.null(QAdj) ){
    
    if(verbose) print("Examining sets for outcome regression")
    
    # do adaptive pre-specification to select from candidate approaches for Qbar
    select.Q <- suppressWarnings(CV.selector(goal=goal, target=target, break.match=break.match, Ldata=Ldata,
                                              CAND.ADJ = cand.QAdj, CAND.FORM=cand.Qform, forQ=T, verbose = verbose,
                                              scale_value=scale_value, scale_value_min=scale_value_min,
                                              folds=folds, sample.effect = sample.effect
                                              ))
    #if(verbose) print(select.Q)
    Q.index <- select.Q$adj.index
    QAdj <- select.Q$Adj
    Qform <- select.Q$form
    #if(verbose){print("Selected for Q:") print(QAdj); print(Qform)}
    
    # if select unadjusted estimator for QbarAW=E(Y|A,W), then stop
    if(sum('U' %in% QAdj) == 1){ 
      g.index <- -99; gAdj <- 'U'; gform <- 'glm'
      var.CV <- select.Q$var.CV
    }
    
    if((sum(QAdj %in% 'U')) == 0 & remove.pscore){
        if(verbose) print('remove.pscore = T; removing selected QAdj variable(s) from candidates for gAdj')
        cand.gAdj <- cand.gAdj[cand.gAdj %nin% QAdj]
        if(length(cand.gAdj) == 0){
          if(verbose) print("all covars used for outcome regression. setting cand.gAdj = U")
          cand.gAdj <- "U"
        }
    }
  }

  
  #=====================================================
  # GENERATE COVARIATE SET / WORKING MODEL COMBINATIONS FOR PROPENSITY SCORE
  cand.g <- get.cand.adj(cand.vars = cand.gAdj, cand.algos = cand.gform,
                         data.adapt.complexity = data.adapt.complexity)
  cand.gAdj <- cand.g$cand.varset
  cand.gform <- cand.g$cand.algos


  
  #=====================================================
  # FIGURE OUT WHICH COMBINATION FOR PROPENSITY SCORE IS BEST
  if( is.null(gAdj) ){ 		
    if(verbose) print("Examining sets for propensity score")
    
    select.G <- suppressWarnings(CV.selector(goal = goal, target = target,
                                              break.match = break.match, Ldata = Ldata,
                                              CAND.ADJ = cand.gAdj, CAND.FORM = cand.gform,
                                              forQ = F,  verbose = verbose,
                                              # input selected variables/form of the outcome regression
                                              QAdj = QAdj, Qform = Qform,
                                              scale_value = scale_value, scale_value_min = scale_value_min,
                                              folds = folds,
                                              sample.effect = sample.effect))
    
    g.index <- select.G$adj.index
    gAdj <- select.G$Adj
    gform <- select.G$form
    var.CV <- select.G$var.CV		
    var.CV.1 <- select.G$var.CV.1		
    var.CV.0 <- select.G$var.CV.0		
    #if(verbose){print(select.G);    print("Selected for g:");  print(gAdj);  print(gform)}
  }		
  output_aps <- list(Q.index = Q.index, QAdj = QAdj, Qform = Qform, 
                     g.index = g.index, gAdj = gAdj, gform = gform,
                     var.CV = var.CV, var.CV.1 = var.CV.1, var.CV.0 = var.CV.0)
  #if(verbose) print(output_aps)
  return(output_aps)
}



#-----------------------------------------------------#-----------------------------------------------------
# CV.selector: function to estimate the cross-validated risk
#		Loss function is the squared-IC; Risk is then the variance of the TMLE
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD), 
#   target of inference 
#		indicator to break the match (break.match)
#		dataset (Ldata)
#		candidate adjustment variables; they do not have to be at the cluster-level
#		indicator if for the conditional mean outcome (forQ)
#		selected adjustment variable for the outcome regression (QAdj)
#	output: selection for adjustment variable (corresponding to a TMLE)
#-----------------------------------------------------#-----------------------------------------------------


CV.selector <- function(goal, target, break.match, Ldata, CAND.ADJ, CAND.FORM, 
                        forQ, QAdj=NULL, Qform=NULL, verbose,
                        scale_value, scale_value_min, folds, sample.effect){
  # if(length(CAND.FORM) == 1){
  #   # if exploring only one estimation algorithm (usually GLM) then need to replicate the number forms
  #   CAND.FORM <- rep(CAND.FORM, length(CAND.ADJ))
  # }
  
  if(length(CAND.FORM) != length(CAND.ADJ)){ # After changes to get.cand.adj, this should not happen.
    stop('PROBLEM: MISMATCH LENGTHS OF ADJ VAR SETS AND MODELS')
  }

  # Number of candidate estimators is given by length Qform//gform
  num.tmles <- length(CAND.FORM)
  CV.risk <- var.CV <- var.CV.1 <- var.CV.0 <- rep(NA, num.tmles)
  
  for(k in 1: num.tmles){	
    
    if(forQ){
      # print("evaluating a Q option:")
      # print(CAND.ADJ[[k]])
      # print(CAND.FORM[k])
      # if selecting the adjustment approach for the outcome regression
      IC.temp <- get.IC.CV(goal = goal, target = target, break.match = break.match, Ldata = Ldata,
                           QAdj = CAND.ADJ[[k]], Qform = CAND.FORM[k],
                           gAdj = NULL, gform = 'glm', verbose = verbose,
                           scale_value = scale_value, scale_value_min = scale_value_min, 
                           folds = folds, sample.effect = sample.effect)
    } else {
      # print("evaluating g option:")
      # print(CAND.ADJ[[k]])
      # print(CAND.FORM[k])
      # if collaboratively selecting the adjustment approach for the pscore
      IC.temp <- get.IC.CV(goal = goal, target = target, break.match = break.match, Ldata = Ldata, 
                           QAdj = QAdj, Qform = Qform, verbose = verbose,
                           gAdj = CAND.ADJ[[k]], gform = CAND.FORM[k],
                           scale_value = scale_value, scale_value_min = scale_value_min, 
                           folds = folds, sample.effect = sample.effect)
    }
    # if(verbose) print(IC.temp)
    # estimating the CV risk for each candidate
    CV.risk[k]<- IC.temp$CV.risk
    # estimating the CV variance for that TMLE
    var.CV[k] <- IC.temp$var.CV
    var.CV.1[k] <- IC.temp$var.CV.1
    var.CV.0[k] <- IC.temp$var.CV.0
  }
  #print(CV.risk)

  # select the candidate estimator resulting in the smallest CV-risk
  adj.index <- which.min(CV.risk)
  output_best <- list(CV.risk = CV.risk[adj.index],
                      adj.index = adj.index, Adj = CAND.ADJ[[adj.index]], form = CAND.FORM[adj.index],
                      var.CV = var.CV[adj.index],
                      var.CV.1 = var.CV.1[adj.index],
                      var.CV.0 = var.CV.0[adj.index])
  return(output_best)
}




#' Calculate a cross-validated estimate of the influence curve 
#'
#' UPDATES
#' previous version only did leave-one-out (unit or pair)
#' this version generalizes to V-fold CV if V>=40
#' - can input the number of folds V (default=10)
#' - folds created stratified on binary outcomes (by default)
#' - if stratify=T and # observations in a given class is <V, 
#' then sets V=min observations in that fold
#'
#' See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching."
#'
#' @param goal goal of analysis ('aRR' or 'RD)
#' @param target cluster/individual effect (JN: Not implemented yet?)
#' @param break.match indicator to break the match
#' @param Ldata input dataset
#' @param QAdj adjustment variable for the outcome regression
#' @param Qform adjustment approach for outcome regression
#' @param gAdj adjustment variable for the pscore
#' @param gform adjustment approach for pscore regression
#' @param scale_value For outcomes not bounded within [0,1], the maximum value of the outcome.
#' @param scale_value_min For outcomes not bounded within [0,1], the minimum value of the outcome.
#' @param folds Number of folds in cross-validation.
#'
#' @return cross-validated estimate of the IC for pair
#' @export
#'
#' @examples
get.IC.CV <- function(goal, target, break.match, verbose,
                      Ldata, QAdj, Qform, gAdj=NULL, gform='glm', 
                      scale_value, scale_value_min, folds, sample.effect){
  
  nFolds <- length(folds)
  DY.CV <- CV.risk <- DY1.CV <- DY0.CV <- NULL
  
  # doing a cross-validated estimate
  for(i in 1:nFolds) {
    
    these <- Ldata$indpt.unit %in% folds[[i]]  ########  IMPORTANT!!!!
    valid <- Ldata[these, ]
    train <- Ldata[!these,]
    
    # run full TMLE algorithm on the training set
    train.out <- do.TMLE(goal=goal, target=target, train=train, QAdj=QAdj, Qform=Qform, 
                         gAdj=gAdj, gform=gform,
                         scale_value=scale_value, scale_value_min=scale_value_min,
                         doing.CV=T, verbose=F, sample.effect = sample.effect)	
    
    # get the relevant components of the IC for the validation set, 
    # using fits based on the training set
    valid.out <- do.TMLE.validset(goal=goal, target=target, valid=valid, train.out=train.out,
                                  scale_value=scale_value, scale_value_min=scale_value_min,
                                  sample.effect = sample.effect)	

    # estimating the CV risk for each candidate
    # risk = Expectation of loss with loss as IC-sq
    # risk = variance of TMLE
    if(break.match){
      DY.CV <- c(DY.CV, valid.out$DY)
      CV.risk <- c(CV.risk, mean(valid.out$DY^2))
    } else {
      DY.CV <- c(DY.CV, valid.out$DY.paired)
      CV.risk <- c(CV.risk, mean(valid.out$DY.paired^2))
    }
    DY1.CV <- c(DY1.CV, valid.out$DY1)
    DY0.CV <- c(DY0.CV, valid.out$DY0)
  }
  #if(verbose){print(valid.out)}
  
  # average across folds
  CV.risk <- mean(CV.risk)
  
  # estimating the CV variance for that TMLE
  var.CV <- stats::var(DY.CV) / length(DY.CV)
  var.CV.1 <- stats::var(DY1.CV) / length(DY1.CV)
  var.CV.0 <- stats::var(DY0.CV) / length(DY0.CV)
  
  return(list(CV.risk = CV.risk, var.CV = var.CV,
              var.CV.0 = var.CV.0, var.CV.1 = var.CV.1))
}



#-----------------------------------------------------#-----------------------------------------------------
# do.TMLE.for.valid: function to obtain a cross-validated estimate of the influence curve
#	for observations in the validation set
#		See Balzer et al. "Adaptive pre-specification in randomized trials with and without pair-matching"
#
#	input: goal of analysis ('aRR' or 'RD'),
#		validation dataset ('valid') 
#		TMLE-fits from training set (train.out)
#	output: cross-validated estimate of the IC,
#		cross-validated risk estimate (loss=IC^2)
#-----------------------------------------------------#-----------------------------------------------------

do.TMLE.validset <- function(goal, target, valid,
                             train.out, scale_value, scale_value_min,
                             sample.effect){
	
	# J <- length(unique(valid$id) )

	#=============================================
	# Step1 - initial estimation of E(Y|A,W)= Qbar(A,W)
	#=============================================
	valid <- do.Init.Qbar(train=valid, QAdj=train.out$QAdj, Qform=train.out$Qform, 
	                     glm.out=train.out$Q.out)$train
	
	#=============================================
	# Step2: Calculate the clever covariate
	#=============================================
	valid <- get.clever.cov(train=valid, gAdj=train.out$gAdj, gform=train.out$gform, 
	                        p.out=train.out$p.out)$train

	#=============================================
	# Step3: Targeting - 			
	#=============================================
	valid <- do.targeting(train=valid, eps=train.out$eps, goal=goal)
	
	#=============================================
	# Step5: Variance estimation using treatment-specific means from training set
	#=============================================
	 get.IC.variance(goal=goal, target=target, Vdata=valid, R1=train.out$R1, R0=train.out$R0, 
	                 scale_value=scale_value, scale_value_min=scale_value_min,
	                 doing.CV=T, sample.effect = sample.effect)
}








#-----------------------------------------------------#-----------------------------------------------------

# adapted from .cvFolds from cvAUC package: https://CRAN.R-project.org/package=cvAUC
# by Erin LeDell 
# **** WARNING - stratify=T option is currently broken for cluster randomizzed trials! *****
get.folds <- function(V, Y, ids, stratify=F){
  
  if(stratify & length(unique(Y))==2){
    # stratify on the outcome
    classes <- tapply(1:length(Y), INDEX=Y, FUN=split, 1)
    ids.Y1 <- ids[classes$`1`]
    ids.noY1 <- ids[classes$`0`]
    if(length(ids.Y1) < V | length(ids.noY1) < V) {
      V <- min( length(ids.Y1), length(ids.noY1))
    }
    ids.Y1.split <- split(sample(length(ids.Y1)), rep(1:V, length=length(ids.Y1)))
    ids.noY1.split <- split(sample(length(ids.noY1)), rep(1:V, length=length(ids.noY1)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- c(ids.Y1[ids.Y1.split[[v]]], ids.noY1[ids.noY1.split[[v]]])
    }
    
  } else {
    # dont stratify on the outcome
    ids.split <- split(sample(length(ids)), rep(1:V, length=length(ids)))
    folds <- vector("list", V)
    for (v in seq(V)){
      folds[[v]] <- ids[ids.split[[v]]]
    }
    
  }
  folds
}
