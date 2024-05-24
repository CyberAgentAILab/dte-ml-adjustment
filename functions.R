##------------------------------------------------------------------------------
## Data generating process
##------------------------------------------------------------------------------
dgp = function(n, rho){
  # Generate binary treatment
  d  = rbinom(n, 1, rho)
  # Number of covariates (fixed)
  num_covariates = 100
  # Generate covariate matrix
  covariates = matrix(runif(n*num_covariates), nrow=n, ncol=num_covariates)
  # Define coefficients
  beta_main = c(rep(1, num_covariates/2), rep(0, num_covariates/2))
  # Error term
  u  = rnorm(n, mean=0, sd=1) 
  # Generate outcome
  y  = covariates %*% beta_main + (covariates^2) %*% beta_main + d + u 
  # Create a data frame
  df.sim = data.frame(cbind(y,d,covariates))
  # Rename variables
  colnames(df.sim) = c("y", "d", paste0("x", seq(1,num_covariates,1)))
  # Return generated dataframe
  return(df.sim)
}

##------------------------------------------------------------------------------
## Data generating process (with decaying sequence of coefficients)
##------------------------------------------------------------------------------
dgp_sequence = function(n, rho, dgp_number){
  # Generate binary treatment
  d  = rbinom(n, 1, rho)
  # Number of covariates
  num_covariates = 100
  # Generate covariate matrix
  covariates = matrix(runif(n*num_covariates), nrow=n, ncol=num_covariates)
  # Define coefficients
  beta_main = c(rep(2/dgp_number, num_covariates/2), rep(0, num_covariates/2))
  u  = rnorm(n, mean=0, sd=1) 
  y  = covariates %*% beta_main + (covariates^2) %*% beta_main + d + u 
  ## data 
  df.sim = data.frame(cbind(y,d,covariates))
  colnames(df.sim) = c("y", "d", paste0("x", seq(1,num_covariates,1)))
  return(df.sim)
}

##------------------------------------------------------------------------------
## Empirical CDF by treatment group
##------------------------------------------------------------------------------
estimate_cdf = function(d.treat, vec.y, vec.d, mat.x, vec.loc){
  ##--------------
  ## setup   
  ##--------------
  n.loc   = length(vec.loc)    ## # of locations 
  col.sub = which(vec.d == d.treat)
  sub.y   = vec.y[col.sub]    ## y - subsample
  sub.x   = mat.x[col.sub,]   ## x - subsample
  sub.folds = folds[col.sub] ## indicator for folds - subsample
  
  ##--------------------------------
  ## empirical distribution function
  ##--------------------------------
  mat.d.y = 1 * outer(sub.y, vec.loc, "<=")  ## n x n.loc
  vec.cdf = as.matrix( colMeans(mat.d.y) )   ## n.loc x 1  
  
  return(vec.cdf)
}

##------------------------------------------------------------------------------
## Regression-adjusted distribution functions
##------------------------------------------------------------------------------
cdf.ML.adj = function(d.treat, vec.y, vec.d, mat.x, vec.loc, model, cross_fitting){
  ##--------------
  ## setup   
  ##--------------
  n.loc   = length(vec.loc)    ## # of locations 
  ind.pte = seq(1,length(vec.loc), by = h.pte)
  n.pte   = length(ind.pte) 
  col.sub = which(vec.d == d.treat)
  sub.y   = vec.y[col.sub]    ## y - subsample
  sub.x   = mat.x[col.sub,]   ## x - subsample
  sub.folds = folds[col.sub] ## indicator for folds - subsample
  
  ##--------------------------------
  ## empirical distribution function
  ##--------------------------------
  mat.d.y = 1 * outer(sub.y, vec.loc, "<=")  ## n x n.loc
  vec.cdf = as.matrix( colMeans(mat.d.y) )   ## n.loc x 1  
  vec.sum = as.matrix( colSums(mat.d.y) ) 
  
  ##--------------
  ## pdf estimation     
  ##--------------
  temp.cdf = vec.cdf[ind.pte]
  vec.pdf  = temp.cdf - c(0, temp.cdf[1:(n.pte-1)])    ## initial prob + difference 
  
  ##--------------
  ## variance      
  ##--------------
  vec.cdf.var = vec.cdf * (1 - vec.cdf)
  vec.pdf.var = vec.pdf * (1 - vec.pdf)
  
  ##--------------------------------
  ## Nuisance estimation
  ##--------------------------------
  vec.cdf.ra = matrix(0, n.loc, 1)
  mu.sub = matrix(NA, length(sub.y), n.loc) # Predicted means in sub sample (treatment or control)
  mu.all = matrix(NA, length(vec.y), n.loc) # Predicted means in full sample
  
  if (cross_fitting==1){
    for(j in 1:n.loc){
      y = mat.d.y[,j] ## jth column - outcome
      # Conduct cross-fitting
      for(f in 1:F){
        if (model=="gradient_boosting"){
          # Fit ML model using data in treatment group d, but not in fold f
          bst = xgboost(data=sub.x[sub.folds!=f,], label=y[sub.folds!=f],
                        max_depth =3, eta=0.1, nrounds = 300, verbose = 0,
                        objective="binary:logistic")
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(bst, sub.x[sub.folds==f,])
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(bst, mat.x[folds==f,])
        }
        
        else if (model=="logistic_regression"){
          # Fit ML model using data in treatment group d, but not in fold f
          logit = fastglm(x=cbind(1, sub.x[sub.folds!=f,]),
                          y= y[sub.folds!=f],
                          family = binomial(link = "logit"))
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(logit, cbind(1,sub.x[sub.folds==f,]), type="response")
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(logit, cbind(1, mat.x[folds==f,]), type="response")
        }
        
        else if (model == "lasso"){
          cvlasso = cv.glmnet(sub.x[sub.folds!=f,], y[sub.folds!=f], alpha=1, family = "binomial", parallel = TRUE)
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(cvlasso, sub.x[sub.folds==f,], s=cvlasso$lambda.min, type="response")
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(cvlasso, mat.x[folds==f,], s=cvlasso$lambda.min, type="response")
          
        }
        else if (model == "random_forest"){
          rf = randomForest(sub.x[sub.folds!=f,], as.factor(y[sub.folds!=f]), ntree=500)
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(rf, sub.x[sub.folds==f,])
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(rf, mat.x[folds==f,])
          
        }
        
        else if (model == "regression_forest"){
          regf = regression_forest(sub.x[sub.folds!=f,], y[sub.folds!=f], tune.parameters = "all")
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(regf, sub.x[sub.folds==f,])
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(regf, mat.x[folds==f,])
          
        }
        
        else if (model=="ols"){
          data_train = as.data.frame(cbind(y[sub.folds!=f], sub.x[sub.folds!=f,]))
          colnames(data_train) = c("y", colnames(sub.x))
          new_x_sub = as.data.frame(sub.x[sub.folds==f,])
          new_x_all = as.data.frame(mat.x[folds==f,])
          ols = lm(y~., data=data_train)
          # Predict on fold f for sub sample
          mu.sub[sub.folds==f,j] =  predict(ols, new_x_sub)
          # Predict on fold f for full sample
          mu.all[folds==f,j] =  predict(ols, new_x_all)
        }
      }
    }
  } 
  
  else if (cross_fitting==0){
    for(j in 1:n.loc){
      y = mat.d.y[,j] ## jth column - outcome
      
      if (model=="gradient_boosting"){
        # Fit ML model using data in treatment group d
        bst = xgboost(data=sub.x, label=y,
                      max_depth =3, eta=0.1, nrounds = 300, verbose = 0,
                      objective="binary:logistic")
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(bst, sub.x)
        # Predict on fold f for full sample
        mu.all[,j] =  predict(bst, mat.x)
      }
      
      else if (model=="logistic_regression"){
        # Fit ML model using data in treatment group d, but not in fold f
        logit = fastglm(x=cbind(1, sub.x),
                        y= y,
                        family = binomial(link = "logit"))
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(logit, cbind(1,sub.x), type="response")
        # Predict on fold f for full sample
        mu.all[,j] =  predict(logit, cbind(1, mat.x), type="response")
      }
      
      else if (model == "lasso"){
        cvlasso = cv.glmnet(sub.x, y, alpha=1, family = "binomial", parallel = TRUE)
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(cvlasso, sub.x, s=cvlasso$lambda.min, type="response")
        # Predict on fold f for full sample
        mu.all[,j] =  predict(cvlasso, mat.x, s=cvlasso$lambda.min, type="response")
        
      }
      else if (model == "random_forest"){
        rf = randomForest(sub.x, as.factor(y), ntree=100)
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(rf, sub.x)
        # Predict on fold f for full sample
        mu.all[,j] =  predict(rf, mat.x)
        
      }
      
      else if (model == "regression_forest"){
        regf = regression_forest(sub.x, y, tune.parameters = "all")
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(regf, sub.x)
        # Predict on fold f for full sample
        mu.all[,j] =  predict(regf, mat.x)
        
      }
      else if (model=="ols"){
        data_train = as.data.frame(cbind(y, sub.x))
        colnames(data_train) = c("y", colnames(sub.x))
        new_x_sub = as.data.frame(sub.x)
        new_x_all = as.data.frame(mat.x)
        ols = lm(y~., data=data_train)
        # Predict on fold f for sub sample
        mu.sub[,j] =  predict(ols, new_x_sub)
        # Predict on fold f for full sample
        mu.all[,j] =  predict(ols, new_x_all)
      }
    }
  }
  
  ##---------------------------------
  ## ML regression adjustment for CDF
  ##---------------------------------
  vec.cdf.ra = vec.cdf + colMeans(mu.all) - colMeans(mu.sub) ## n.loc x 1
  
  ##---------------
  ## pdf estimation     
  ##---------------
  temp.cdf.ra = vec.cdf.ra[ind.pte]
  vec.pdf.ra  = temp.cdf.ra - c(0, temp.cdf.ra[1:(n.pte-1)])    ## initial prob + difference 
  
  ##--------------
  ## variance      
  ##--------------
  vec.cdf.var.ra = vec.cdf.ra * (1 - vec.cdf.ra)
  vec.pdf.var.ra = vec.pdf.ra * (1 - vec.pdf.ra)
  
  ## return 
  return(list(mu.all = mu.all,
              cdf = vec.cdf, 
              cdf.ra = vec.cdf.ra,
              pdf = vec.pdf,
              pdf.ra = vec.pdf.ra,
              cdf.var = vec.cdf.var,
              pdf.var = vec.pdf.var,
              cdf.var.ra = vec.cdf.var.ra,
              pdf.var.ra = vec.pdf.var.ra,
              sum = vec.sum,
              num_obs = length(sub.y)))
}

##------------------------------------------------------------------------------
## Regression-adjusted DTE and PTE estimation
##------------------------------------------------------------------------------
DTE.ML.estimation = function(vec.y, vec.d, mat.x, vec.loc, h.pte, model, cross_fitting, B.size=500){
  
  num_obs = length(vec.y)
  
  ## locations   
  n.loc   = length(vec.loc)
  ind.pte = seq(1,length(vec.loc), by = h.pte)
  n.pte   = length(ind.pte) 
  
  ## containers 
  mat.sum    = matrix(NA, n.loc, 2)
  mat.cdf    = matrix(NA, n.loc, 2)
  mat.cdf.var = matrix(NA, n.loc, 2)
  mat.cdf.ra = matrix(NA, n.loc, 2)
  mat.cdf.var.ra = matrix(NA, n.loc, 2)
  mat.pdf    = matrix(NA, n.pte, 2)
  mat.pdf.var    = matrix(NA, n.pte, 2)
  mat.pdf.ra = matrix(NA, n.pte, 2)
  mat.pdf.var.ra = matrix(NA, n.pte, 2)
  mat.mu.all    = array(NA, c(2, num_obs, n.loc)) # Predictions for all observables
  treat_prob = matrix(NA, 1, 2 ) # Number of observations in each group
  
  ##------------------------------------------
  ## estimation for each group     
  ##------------------------------------------
  for(j in 1:2){  ## treatment/control group 
    
    ## data: subset: first, estimate "treatment" group and then "control" group     
    d.treat = 2 - j ## j=1 <-> 1:treatment & j=2 <-> 0:control  
    
    ## CDF estimation   
    res.cdf = cdf.ML.adj(d.treat, vec.y, vec.d, mat.x, vec.loc, model, cross_fitting)
    
    ## save 
    mat.sum[,j] = res.cdf$sum
    
    mat.cdf[,j]     = res.cdf$cdf
    mat.cdf.var[,j] = res.cdf$cdf.var
    
    mat.cdf.ra[,j]  = res.cdf$cdf.ra
    mat.cdf.var.ra[,j] = res.cdf$cdf.var.ra
    
    mat.pdf[,j] = res.cdf$pdf
    mat.pdf.var[,j] = res.cdf$pdf.var
    
    mat.pdf.ra[,j] = res.cdf$pdf.ra
    mat.pdf.var.ra[,j] = res.cdf$pdf.var.ra
    
    mat.mu.all[j,,] = res.cdf$mu.all
    treat_prob[,j] = res.cdf$num_obs
    
  }
  
  ##------------------------------------------
  ## treatment effect estimation (DTE and PTE)      
  ##------------------------------------------
  vec.dte     = mat.cdf[,1]    - mat.cdf[,2]   
  vec.pte     = mat.pdf[,1]    - mat.pdf[,2]   
  vec.dte.ra  = mat.cdf.ra[,1] - mat.cdf.ra[,2]   
  vec.pte.ra  = mat.pdf.ra[,1] - mat.pdf.ra[,2] 
  
  ##--------------------------------------------------------------- 
  ## variance for adjusted-DTE: uniform and pointwise
  ##--------------------------------------------------------------
  # Predictions in treatment and control group
  mat.mu.1 = as.matrix(mat.mu.all[1,,])
  mat.mu.0 = as.matrix(mat.mu.all[2,,])
  
  # Number of observations in each group
  num_1 = treat_prob[,1]
  num_0 = treat_prob[,2]
  
  # Indicator that y is below threshold u
  mat.y.u = 1 * outer(vec.y, vec.loc, "<=")       ## n x n.loc
  mat.d = replicate(n.loc, vec.d)                 ## n x n.loc 
  mat.dte = t(matrix((rep(vec.dte, num_obs)), nrow=n.loc, ncol=num_obs))
  mat.dte.ra = t(matrix((rep(vec.dte.ra, num_obs)), nrow=n.loc, ncol=num_obs))
  
  # Sample analog of asymptotic variance
  omega_moment = (num_obs/num_1*(mat.d*(mat.y.u-mat.mu.1)) + mat.mu.1 - num_obs/num_0*((1-mat.d)*(mat.y.u-mat.mu.0)) - mat.mu.0 - mat.dte.ra)^2  ## n x n.loc
  # Take average over all observations to obtain an estimate of variance
  omega = colMeans(omega_moment) ## 1 x n.loc
  
  ## Pointwise confidence intervals (moment conditions)
  vec.dte.ra.cv025.moment = vec.dte.ra - 1.96 * sqrt(omega/num_obs)
  vec.dte.ra.cv975.moment = vec.dte.ra + 1.96 * sqrt(omega/num_obs)
  
  # Influence function
  influence_function = num_obs/num_1*(mat.d*(mat.y.u-mat.mu.1)) + mat.mu.1 - num_obs/num_0*((1-mat.d)*(mat.y.u-mat.mu.0)) - mat.mu.0 - mat.dte.ra  ## n x n.loc
  
  # An empty matrix to store t-stats
  tstats = matrix(0, nrow = B.size, ncol = n.loc)
  boot.draw = matrix(0, nrow = B.size, ncol = n.loc)
  boot.dte.ra = t(matrix((rep(vec.dte.ra, B.size)), nrow=n.loc, ncol=B.size))
  
  for (b in 1:B.size){
    # Mammen multiplier 
    eta1 = rnorm(num_obs,0,1)
    eta2 = rnorm(num_obs,0,1)
    xi = eta1/sqrt(2) + (eta2^2-1)/2
    
    # Get bootstrap draws
    boot.draw[b,] = vec.dte.ra + 1/num_obs*(colSums(xi*influence_function))
  }
  
  # Bootstrap standard errors
  q75 = apply(boot.draw, 2, quantile, probs = 0.75)
  q25 = apply(boot.draw, 2, quantile, probs = 0.25)
  boot.se = (q75-q25)/(qnorm(0.75)-qnorm(0.25))
  
  # Compute t-stats
  tstats = abs(boot.draw-vec.dte.ra)/boot.se
  
  # Take max over t-stats
  max_tstats = apply(tstats, 1, max)
  # Take p quantile over bootstrap draws
  quantile_max_tstats = quantile(max_tstats, 0.05)
  
  ## Uniform confidence band 
  vec.dte.ra.cv025.boot = vec.dte.ra - quantile_max_tstats*boot.se
  vec.dte.ra.cv975.boot = vec.dte.ra + quantile_max_tstats*boot.se
  
  ##----------------------------------------------------------------- 
  ## variance for simple DTE: uniform and pointwise
  ##-----------------------------------------------------------------
  # Influence function
  influence_function.simple = num_obs/num_1*mat.d*(mat.y.u)-num_obs/num_0*(1-mat.d)*(mat.y.u)- mat.dte
  
  # An empty matrix to store t-stats
  tstats.simple = matrix(0, nrow = B.size, ncol = n.loc)
  boot.draw.simple = matrix(0, nrow = B.size, ncol = n.loc)
  boot.dte = t(matrix((rep(vec.dte, B.size)), nrow=n.loc, ncol=B.size))
  
  for (b in 1:B.size){
    # Mammen multiplier 
    eta1 = rnorm(num_obs,0,1)
    eta2 = rnorm(num_obs,0,1)
    xi = eta1/sqrt(2) + (eta2^2-1)/2
    
    # Get bootstrap draw
    boot.draw.simple[b,] = 1/num_obs*colSums(xi*(influence_function.simple))
  }
  
  # Sample analog of asymptotic variance
  omega_moment.simple = (num_obs/num_1*mat.d*mat.y.u-num_obs/num_0*(1-mat.d)*mat.y.u- mat.dte)^2  ## n x n.loc
  # Take average over all observations to obtain an estimate of variance
  omega.simple = colMeans(omega_moment.simple) ## 1 x n.loc
  
  # Compute t-stats
  tstats.simple = abs(boot.draw.simple)/sqrt(omega.simple/num_obs)
  
  # Take max over t-stats
  max_tstats.simple = apply(tstats.simple, 1, max)
  # Take p quantile over bootstrap draws
  quantile_max_tstats.simple = quantile(max_tstats.simple, 0.05)
  
  ## Uniform confidence band 
  vec.dte.cv025.boot = vec.dte - quantile_max_tstats.simple* sqrt(omega.simple/num_obs)
  vec.dte.cv975.boot = vec.dte + quantile_max_tstats.simple* sqrt(omega.simple/num_obs)
  
  ## Point wise confidence interval (moment conditions)
  vec.dte.cv025.moment = vec.dte - 1.96* sqrt(omega.simple/num_obs)
  vec.dte.cv975.moment = vec.dte + 1.96* sqrt(omega.simple/num_obs)
  
  ##------------------------------------  
  ## other analytic standard errors
  ##------------------------------------
  w1        = 1 / mean(vec.d)      ## inverse 
  w0        = 1/ (1 - mean(vec.d)) ## inverse
  vec.dte.var = (w1 * mat.cdf.var[,1]) + (w0 * mat.cdf.var[,2])  ## "+"
  vec.dte.var.ra = (w1 * mat.cdf.var.ra[,1]) + (w0 * mat.cdf.var.ra[,2])  ## "+"
  vec.pte.var = (w1 * mat.pdf.var[,1]) + (w0 * mat.pdf.var[,2])  ## "+"
  vec.pte.var.ra = (w1 * mat.pdf.var.ra[,1]) + (w0 * mat.pdf.var.ra[,2])  ## "+"
  
  ## confidence interval: two-side test with equal cv
  vec.cv.dte = (1.96 / sqrt(num_obs))  * sqrt(vec.dte.var) 
  mat.ci.dte = cbind(vec.dte, vec.dte) +  cbind(-vec.cv.dte, vec.cv.dte)
  
  vec.cv.dte.ra = (1.96 / sqrt(num_obs))  * sqrt(vec.dte.var.ra)
  mat.ci.dte.ra = cbind(vec.dte.ra, vec.dte.ra) +  cbind(-vec.cv.dte.ra, vec.cv.dte.ra)
  
  vec.cv.pte = (1.96 / sqrt(num_obs))  * sqrt(vec.pte.var)
  mat.ci.pte = cbind(vec.pte, vec.pte) +  cbind(-vec.cv.pte, vec.cv.pte)
  
  vec.cv.pte.ra = (1.96 / sqrt(num_obs))  * sqrt(vec.pte.var.ra)
  mat.ci.pte.ra = cbind(vec.pte.ra, vec.pte.ra) +  cbind(-vec.cv.pte.ra, vec.cv.pte.ra)
  
  ##---------------------  
  ## stack up results 
  ##---------------------
  est.dte    = data.frame( cbind(vec.loc, mat.cdf,    vec.dte, mat.ci.dte,
                                 vec.dte.cv025.moment, vec.dte.cv975.moment, 
                                 vec.dte.cv025.boot, vec.dte.cv975.boot) )
  est.dte.ra = data.frame( cbind(vec.loc, mat.cdf.ra, vec.dte.ra, mat.ci.dte.ra,
                                 vec.dte.ra.cv025.moment, vec.dte.ra.cv975.moment,
                                 vec.dte.ra.cv025.boot, vec.dte.ra.cv975.boot) )
  est.pte    = data.frame( cbind(vec.loc[ind.pte[1:n.pte]], mat.pdf,    vec.pte, mat.ci.pte) )
  est.pte.ra = data.frame( cbind(vec.loc[ind.pte[1:n.pte]], mat.pdf.ra, vec.pte.ra, mat.ci.pte.ra) )
  
  colnames(est.dte)    = c("vec.loc", "cdf1", "cdf0", "dte","cv025", "cv975", "cv025_moment", "cv975_moment", "cv025_boot", "cv975_boot")
  colnames(est.dte.ra) = c("vec.loc", "cdf1", "cdf0", "dte", "cv025", "cv975", "cv025_moment", "cv975_moment", "cv025_boot", "cv975_boot")
  colnames(est.pte)    = c("vec.loc", "pdf1", "pdf0", "pte", "cv025", "cv975")
  colnames(est.pte.ra) = c("vec.loc", "pdf1", "pdf0", "pte", "cv025", "cv975")
  
  return( list(dte     = est.dte, 
               pte     = est.pte,
               dte.ra  = est.dte.ra, 
               pte.ra  = est.pte.ra))
  
}  

##------------------------------------------------------------------------------
## Regression-adjusted QTE estimation
##------------------------------------------------------------------------------
# Simple quantile estimation
estimate_quantile = function(treatment_group, vec.y, vec.d, mat.x, probability_vector) {
  lower_bound = min(vec.y)  # Choose a suitable lower bound for searching
  upper_bound = max(vec.y) # Choose a suitable upper bound for searching
  tolerance = 1e-2   # Choose a suitable tolerance
  
  n.prob = length(probability_vector)
  vec.quantile = matrix(NA, n.prob, 1)
  
  for(j in 1:n.prob){
    probability = probability_vector[j]
    
    lower_bound_temp = lower_bound
    upper_bound_temp = upper_bound
    
    while (upper_bound_temp - lower_bound_temp > tolerance) {
      midpoint = (lower_bound_temp + upper_bound_temp) / 2
      if (estimate_cdf(treatment_group, vec.y, vec.d, mat.x, midpoint) < probability) {
        lower_bound_temp = midpoint
      } else {
        upper_bound_temp = midpoint
      }
    }
    vec.quantile[j] = midpoint
  }
  return(vec.quantile)
}

# Regression adjusted quantile estimation
estimate_reg_adj_quantile = function(treatment_group, vec.y, vec.d, mat.x, model, cross_fitting, probability_vector) {
  lower_bound = min(vec.y)  # Choose a suitable lower bound for searching
  upper_bound = max(vec.y)  # Choose a suitable upper bound for searching
  tolerance = 1e-2   # Choose a suitable tolerance
  
  n.prob = length(probability_vector)
  vec.quantile = matrix(NA, n.prob, 1)
  
  for(j in 1:n.prob){
    probability = probability_vector[j]
    
    lower_bound_temp = lower_bound
    upper_bound_temp = upper_bound
    
    while (upper_bound_temp - lower_bound_temp > tolerance) {
      midpoint = (lower_bound_temp + upper_bound_temp) / 2
      if (cdf.ML.adj(treatment_group, vec.y, vec.d, mat.x, midpoint, model, cross_fitting)$cdf.ra < probability) {
        lower_bound_temp = midpoint
      } else {
        upper_bound_temp = midpoint
      }
    }
    vec.quantile[j] = midpoint
  }
  return(vec.quantile)
}

# QTE estimation
QTE.ML.estimation = function(vec.y, vec.d, mat.x, model, cross_fitting, vec.loc.qte){
  ## locations   
  n.loc   = length(vec.loc.qte)
  
  ## container 
  mat.quantile    = matrix(NA, n.loc, 2)
  mat.quantile.ra    = matrix(NA, n.loc, 2)
  
  ## estimation for each group 
  for(j in 1:2){  ## treatment/control group 
    
    ## data: subset: first, estimate "treatment" group and then "control" group     
    d.treat = 2 - j ## j=1 <-> 1:treatment & j=2 <-> 0:control  
    
    ## Quantile estimation 
    res.quantile = estimate_quantile(d.treat, vec.y, vec.d, mat.x, vec.loc.qte)
    
    ## Regression-adjusted quantile estimation   
    res.quantile.ra = estimate_reg_adj_quantile(d.treat, vec.y, vec.d, mat.x, model, cross_fitting, vec.loc.qte)
    
    ## save 
    mat.quantile[,j] = res.quantile
    mat.quantile.ra[,j] = res.quantile.ra
  }
  
  ## quantile treatment effect 
  vec.qte     = mat.quantile[,1] - mat.quantile[,2] 
  vec.qte.ra     = mat.quantile.ra[,1] - mat.quantile.ra[,2]   
  
  est.qte = data.frame( cbind(vec.loc.qte, mat.quantile, vec.qte) )
  est.qte.ra = data.frame( cbind(vec.loc.qte, mat.quantile.ra, vec.qte.ra) )
  
  colnames(est.qte)    = c("vec.loc", "quantile1", "quantile0", "qte")
  colnames(est.qte.ra) =  c("vec.loc", "quantile1", "quantile0", "qte")
  
  return( list(qte     = est.qte, 
               qte.ra     = est.qte.ra))
  
}

##------------------------------------------------------------------------------
## Run simulations for DTE
##------------------------------------------------------------------------------
run_simulation_dte = function(sim){
  ## DGP 
  eval(parse(text=paste0("df=dgp(n, rho)")))
  
  ## data 
  vec.y = df$y
  vec.d = df$d
  mat.x = as.matrix(df[, 3:ncol(df)])
  num_obs = n ## sample size (for whole dataset)
  
  ## model
  res.ols = DTE.ML.estimation(vec.y, vec.d, mat.x, vec.loc, h.pte, "ols", 1)
  res.lasso = DTE.ML.estimation(vec.y, vec.d, mat.x, vec.loc, h.pte, "lasso", 1)
  
  results = list(res.ols = res.ols, res.lasso=res.lasso)
  v.info = paste0("dgp", "_n", n, "_round", sim)
  saveRDS(results, file=paste0("./result/dte_result_n", n, "/", v.info, ".rds"))
  
  ## return
  return (list(res.ols = res.ols, res.lasso=res.lasso))
}

##------------------------------------------------------------------------------
## Run simulations for QTE
##------------------------------------------------------------------------------
run_simulation_qte = function(sim){
  ## DGP 
  eval(parse(text=paste0("df=dgp(n, rho)")))
  
  ## data 
  vec.y = df$y
  vec.d = df$d
  mat.x = as.matrix(df[, 3:ncol(df)])
  num_obs = n ## sample size (for whole dataset)
  
  ## model
  res.ols = QTE.ML.estimation(vec.y, vec.d, mat.x, "ols", 1, quantiles)
  res.lasso = QTE.ML.estimation(vec.y, vec.d, mat.x, "lasso", 1, quantiles)
  
  results = list(res.ols = res.ols, res.lasso=res.lasso)
  v.info = paste0("dgp", "_n", n, "_round", sim)
  saveRDS(results, file=paste0("./result/qte_result_n", n, "/", v.info, ".rds"))
  
  return (list(res.ols = res.ols, res.lasso=res.lasso))
}

##------------------------------------------------------------------------------
## Run simulations for sequence of DGPs
##------------------------------------------------------------------------------

run_simulation_sequence = function(sim){
  ## DGP 
  eval(parse(text=paste0("df=dgp(n, rho, dgp_number)")))
  
  ## data 
  vec.y = df$y
  vec.d = df$d
  mat.x = as.matrix(df[, 3:ncol(df)])
  num_obs = n ## sample size (for whole dataset)
  
  ## model
  res.ols = DTE.ML.estimation(vec.y, vec.d, mat.x, vec.loc, h.pte, "ols", 1)
  res.lasso = DTE.ML.estimation(vec.y, vec.d, mat.x, vec.loc, h.pte, "lasso", 1)
  
  results = list(res.ols = res.ols, res.lasso=res.lasso)
  v.info = paste0("dgp_seq", dgp_number, "_round", sim)
  saveRDS(results, file=paste0("./result/dgp_seq", dgp_number, "/", v.info, ".rds"))
  
  ## return
  return (list(res.ols = res.ols, res.lasso=res.lasso))
}



