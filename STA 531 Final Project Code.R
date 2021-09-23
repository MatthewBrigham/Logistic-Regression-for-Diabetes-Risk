#Matthew Brigham
#STA 531
#Final Project
#
# Diabetes With Generalized Linear Mixed Models and Multiple Logistic Regression.
# MICE missing value imputation was used and one data set was used to beuild the models.
# At the end of the code, the final regression model was used on the remaining imputed data sets
# and it is shown that all of the coefficients are similar.
#
# Note: glyhb > 7 is diagnosed diabetes
#
# fit_08 is reduced GLM model


#######################################
############# Import Data #############
#######################################

orig_dat = read.csv("diabetes.csv")
head(orig_dat)

dat0 = orig_dat

#######################################
########## Data Manipulations #########
#######################################


# Identify and Impute Missing Values
    
    library(mice)
    library(VIM)
    
    # Look for Pattern of Missingness
    md.pattern(dat0)  #tells us number of observations with that missing data structure
    agg_plot = aggr(dat0, combined = F) #histogram showing proportion of missing data by variable
    
    
    # Impute NA values
    dat1=dat0
    dat1 = dat0[,-c(15, 16)] #remove categorical variables and 2nd measurements of BP
    imp = mice(dat1, maxit = 0) # 0 iterations
    
    pred_mat = imp$predictorMatrix  #get predictor matrix
    method = imp$method  #get methods, mostly pmm: predictive mean matching
    
    imp2 = mice(dat1, m=5, maxit = 50, predictorMatrix = pred_mat, method = method, seed = 500, print = F) #creates 5 datasets of imputed values
    
    # Analyze the distribution with and without imputed values for each of 5 imputed datasets
    newdat_imp_1 = complete(imp2,1)
    newdat_imp_2 = complete(imp2,2)
    newdat_imp_3 = complete(imp2,3)
    newdat_imp_4 = complete(imp2,4)
    newdat_imp_5 = complete(imp2,5)
    
    #Density Plot of Imputed Data
    par(mfrow=c(2,7))
    densityplot(imp2)

#Add new variables
    
    dat = newdat_imp_1
    
    # Rescale and square nonlinear terms (determined later)
    dat$stab.glu.2 = ((dat$stab.glu - median(dat$stab.glu))/sd(dat$stab.glu))^2  #center before squaring
    dat$age.2 =  ((dat$age - median(dat$age))/sd(dat$age))^2
    dat$bp.1s.2 = ((dat$bp.1s - median(dat$bp.1s))/sd(dat$bp.1s))^2
    
    # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
    
    for (i in 1:length(dat$glyhb)){
      if (dat$glyhb[i] >= 7){
        dat$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat$glyhb.cat[i] = 0
      }
    }
    
    
    #calculate waist to hip ratio
    dat$w.h.ratio = (dat$waist)/dat$hip
    
    #Add new variable BMI
    dat$bmi = dat$weight/((dat$height)^2)*703
    
    # Change "frame" to a factor small, medium, large
    dat$frame = factor(dat$frame, levels = c("small", "medium", "large"))

###############################  
######## Data Exploration #####
###############################  

    dat = newdat_imp_1 #first imputed data set

#summary statistics of original and imputed data sets
    
    summary(dat) #contains imputed values and new values
    summary(dat0) #contains NA values

#Linearity - restricted cubic splines for continuous variables.
    
    fit_chol = lrm(glyhb.cat ~ rcs(chol,3), x=T, y=T, data =dat)
    anova(fit_chol)  #linear p-value = 0.9129
    plot(dat$chol,dat$glyhb, main = "Chol")
    
    fit_stab.glu = lrm(glyhb.cat ~ rcs(stab.glu,3), x=T, y=T, data =dat)
    anova(fit_stab.glu)  #NONLINEAR p-value = 0.0014
    plot(dat$stab.glu,dat$glyhb, main = "Stab. Glu")
    
    fit_hdl = lrm(glyhb.cat ~ rcs(hdl,3), x=T, y=T, data =dat)
    anova(fit_hdl)  #LINEAR or NONLINEAR p-value = 0.0502
    plot(dat$hdl,dat$glyhb, main = "HDL")
    
    fit_ratio = lrm(glyhb.cat ~ rcs(ratio,3), x=T, y=T, data =dat)
    anova(fit_ratio)  #LINEAR p-value = 0.5695
    plot(dat$ratio,dat$glyhb, main = "ratio")
    
    fit_age = lrm(glyhb.cat ~ rcs(age,3), x=T, y=T, data =dat)
    anova(fit_age)  #NONLINEAR p-value = 0.0166
    plot(dat$age,dat$glyhb, main = "Age")
    
    fit_height = lrm(glyhb.cat ~ rcs(height,3), x=T, y=T, data =dat)
    anova(fit_height)  #LINEAR p-value = 0.9470
    plot(dat$height,dat$glyhb, main = "Height")
    
    fit_weight = lrm(glyhb.cat ~ rcs(weight,3), x=T, y=T, data =dat)
    anova(fit_weight)  #LINEAR or NONLINEAR p-value = 0.1086
    plot(dat$weight,dat$glyhb, main = "Weight")
    
    fit_bp.1s = lrm(glyhb.cat ~ rcs(bp.1s,3), x=T, y=T, data =dat)
    anova(fit_bp.1s)  #NONLINEAR p-value = 0.0182
    plot(dat$bp.1s,dat$glyhb, main = "Systolic BP")
    
    fit_bp.1d = lrm(glyhb.cat ~ rcs(bp.1d,3), x=T, y=T, data =dat)
    anova(fit_bp.1d)  #LINEAR p-value = 0.5562
    plot(dat$bp.1d,dat$glyhb, main = "Diastolic BP")
    
    fit_waist = lrm(glyhb.cat ~ rcs(waist,3), x=T, y=T, data =dat)
    anova(fit_waist)  #LINEAR p-value = 0.1611
    plot(dat$waist,dat$glyhb, mian = "Waist")
    
    fit_hip = lrm(glyhb.cat ~ rcs(hip,3), x=T, y=T, data =dat)
    anova(fit_hip)  #LINEAR or NONLINEAR p-value = 0.0846
    plot(dat$hip,dat$glyhb, main = "Hip")
    
    fit_time.ppn = lrm(glyhb ~ rcs(time.ppn,3), x=T, y=T, data =dat)
    anova(fit_time.ppn)  #LINEAR p-value = 0.9678
    plot(dat$time.ppn,dat$glyhb, main = "Time.ppn")
    
    fit_bmi = lrm(glyhb ~ rcs(bmi,3), x=T, y=T, data =dat)
    anova(fit_bmi)  #LINEAR or NONLINEAR p-value = 0.0706
    plot(dat$bmi,dat$glyhb, main = "BMI")
    
#Scatterplots
    
    x1 = dat[ ,c("glyhb", "chol", "stab.glu", "stab.glu.2", "hdl", "ratio", "age",       
              "bp.1s", "bmi", "age.2", "w.h.ratio" )]
    pairs(x1)
    
    par(mfrow= c(2,4))
    plot(dat$chol,dat$glyhb, main = "Chol")
    plot(dat$stab.glu,dat$glyhb, main = "Stab. Glu")
    plot(dat$hdl,dat$glyhb, main = "HDL")
    plot(dat$ratio,dat$glyhb, main = "ratio")
    plot(dat$age,dat$glyhb, main = "Age")
    plot(dat$bp.1s,dat$glyhb, main = "Systolic BP")
    plot(dat$bmi,dat$glyhb, main = "BMI")

######################################################
######## GLMM W/ and W/O Random Location and Restr. Cubic Splines
######################################################  

#Fit glmer with random effect location    
        
    #Full Model
    fit_1 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + weight + frame + bp.1s + rcs(bp.1s, 3) + bp.1d + waist*hip + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_1 = summary(fit_1); 
    sumfit_1
    
    #remove bp.1d
    fit_2 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + weight + frame + bp.1s + rcs(bp.1s, 3) + waist*hip + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_2 = summary(fit_2); 
    sumfit_2
    
    #remove frame
    fit_3 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + weight + bp.1s + rcs(bp.1s, 3) + waist*hip + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_3 = summary(fit_3); 
    sumfit_3
    
    #remove waist*hip
    fit_4 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + weight + bp.1s + rcs(bp.1s, 3) + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_4 = summary(fit_4); 
    sumfit_4
    
    #remove waist*hip
    fit_4 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + weight + bp.1s + rcs(bp.1s, 3) + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_4 = summary(fit_4); 
    sumfit_4
    
    #remove weight
    fit_5 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + bp.1s + rcs(bp.1s, 3) + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_5 = summary(fit_5); 
    sumfit_5
    
    #remove ratio
    fit_6 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age + rcs(age, 3) + height 
                  + bp.1s + rcs(bp.1s, 3) + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_6 = summary(fit_6); 
    sumfit_6
    
    #remove rcs(bp.1s, 3)
    fit_7 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age + rcs(age, 3) + height 
                  + bp.1s + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_7 = summary(fit_7); 
    sumfit_7
    
    #remove rcs(age, 3)
    fit_8 = glmer(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age + height
                  + bp.1s + gender*bmi + w.h.ratio
                  + (1 | location), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 10)
    sumfit_8 = summary(fit_8); 
    sumfit_8
    
    #likelihood ratio test for full and reduced GLMM data sets
    res_dev_glmm_loc_reduced = 155.3
    res_dev_glmm_full = 
    df_glmm_loc_reduced = 390
    df_glmm_loc_full = 
    AIC_glmm_loc_reduced = 181.3
    AIC_glmm_loc_full = 
    pchisq(res_dev_glmm_loc_reduced - res_dev_glmm_full, df_glmm_loc_reduced - df_glmm_loc_full, 
           lower.tail = F) #p-value is 0.990 - therefore can use either
    
#####################################
######## Multiple Logistic Regression Assessment
#####################################

#Compare to glmm random location with likelihood ratio test
#Restricted cubic splines were used to model terms identified as non-linear
    
    md.pattern(dat)  #tells us number of observations with that missing data structure
    agg_plot = aggr(dat, combined = F) #histogram showing proportion of missing data by variable
    
    #full model - removed frame, rcs(bp.1s, 3), rcs(stab.glu, 3), rcs(age, 3), w.h.ratio
    fit_01 = glm(glyhb.cat ~ chol + stab.glu + hdl + ratio + age +  height 
                 + weight + bp.1s + bp.1d + waist*hip + gender*bmi,
                 data = dat)
    sumfit_01 = summary(fit_01)
    sumfit_01
    
    #remove bp.1d
    fit_02 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                 + weight + frame + bp.1s + rcs(bp.1s, 3)  + waist*hip + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_02 = summary(fit_02)
    sumfit_02
    
    #remove waist*hip
    fit_02 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                 + weight + frame + bp.1s + rcs(bp.1s, 3) + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_02 = summary(fit_02)
    sumfit_02
    
    #remove rcs(bp.1s, 3)
    fit_03 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                 + weight + frame + bp.1s + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_03 = summary(fit_03)
    sumfit_03

    #remove frame
    fit_04 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                 + weight + bp.1s + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_04 = summary(fit_04)
    sumfit_04
    
    #remove weight
    fit_05 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + ratio + age + rcs(age, 3) + height 
                  + bp.1s + gender*bmi + w.h.ratio,
                  data = dat, family = binomial() )
    sumfit_05 = summary(fit_05)
    sumfit_05
    
    #remove ratio
    fit_06 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age + rcs(age, 3) + height 
                 + bp.1s + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_06 = summary(fit_06)
    sumfit_06

    #remove rcs(age, 3)    
    fit_07 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age + height 
                 + bp.1s + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_07 = summary(fit_07)
    sumfit_07

    #looks good, try removing height (other research says height is important)
    
    #remove height
    fit_08 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                 + bp.1s + gender*bmi + w.h.ratio,
                 data = dat, family = binomial() )
    sumfit_08 = summary(fit_08)
    sumfit_08
    
    confint(fit_08)
    
    #likelihood ratio test to compare GLM full and reduced Mult Log Regression
    res_dev_log_full = 147.84
    res_dev_log_red = 155.28
    df_log_full = 369
    df_log_red = 392
    AIC_log_full = 191.84
    AIC_log_red = 177.28
    pchisq(res_dev_log_red - res_dev_log_full, df_log_red - df_log_full, 
           lower.tail = F) #p-value is 0.999 - therefore can use reduced
    
    
    #likelihood ratio test for GLMM random location and Mult Log. Regression
    res_dev_glmm_loc = 155.3
    res_dev_logistic = 155.28
    df_glmm_loc = 390
    df_logistic = 392
    AIC_glmm_loc = 181.3
    AIC_logistic = 177.28
    pchisq(res_dev_glmm_loc - res_dev_logistic, df_glmm_loc - df_logistic, 
           lower.tail = F) #p-value is 0.990 - therefore can use either
    
###########################################################################
##########   Assessment of Reduced Multiple Logistic Regression Model #####    
###########################################################################

    
    #fit a logistic regression model using terms from glm (which gave me p values)
    require(rms)
    dd = datadist(dat)
    options(datadist = 'dd')
    lrm.reduced = lrm(glyhb.cat ~ chol + stab.glu + hdl + age  
                      + bp.1s + gender*bmi, data = dat)
    
    
    #Calibration using bootstraps and Calibration Plot - DOES NOT WORK (can't use glm)
    cal_log = calibrate(lrm.reduced, B=100) #does not work with glm or lrm models
    boot_strap = boot(dat,
                      predict(fit_08),
                      R = 1000,
                      cor.type = 's')
    boots = bootstrap(dat, 100)
    probs = predict(fit_08, newdata = dat[boots[1,],])
    predict_matrix = data.matrix(dat[,])
    
    #discrimination - ROC Curve
    library(pROC)
    pred = predict(fit_08, type = c("response"))
    roccurve = roc(dat$glyhb.cat ~ pred)
    plot(roccurve, main = "ROC Curve")
    auc(roccurve) #Area Under Curve = 0.9419 - good discrimination 
    
    #validation
    library(purrr) #for map()
    validate(lrm.reduced, method = "boot", B= 1000, data = dat, x = TRUE, y = TRUE)
    require(lme4)
    require(languageR)
    somers.mer(fit_08)
    
    x = bootstrap(dat, 100)
    map(x, lrm.reduced)
    somers2()
    
    
############################################
############# Assess Model GLM #############
############################################
    
    
    #Residuals
    res = residuals(fit_08, "pearson")
    hist(res)
    plot(res) #resid vs index
    min(res)
    max(res)
    which(res>25) #observation 334 has extremely large residual comparatively
    
    #try removing observation 334 - Does not improve fit
    fit_09 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                 + bp.1s + gender*bmi + w.h.ratio,
                 data = dat[-334,], family = binomial() )
    sumfit_09 = summary(fit_09)
    sumfit_09 
    plot(fit_09)
    
    #Overdispersion - no consensus on how to calc, cite source bookmarked
    r_df = df.residual(fit_08)
    r_pearson = residuals(fit_08, "pearson")
    test_stat_over = sum(r_pearson^2)
    pearson_ratio = test_stat_over/r_df #indicates possibility of underdispersion. usually
    p = pchisq(test_stat_over, df = r_df, lower.tail = F)
    p #p=0 reject null, dispersion is significantly different from zero
    
    
##########################################################
############# Fit Remaining Imputed Data Sets #############
##########################################################
#All models were ran on first imputed data set.
#Show that imputed data sets are not significantly different due to small 
#percentage of missing values.
    
    #Final Fit for remaining imputed data sets
   
    #assume same final features as calculated from 1st imputed data set
    dat2 = complete(imp2,2)
    dat3 = complete(imp2,3)
    dat4 = complete(imp2,4)
    dat5 = complete(imp2,5)
    
#Add new variables to data sets
    #dat2
    dat2$stab.glu.2 = ((dat2$stab.glu - median(dat2$stab.glu))/sd(dat2$stab.glu))^2 
    dat2$age.2 =  ((dat2$age - median(dat2$age))/sd(dat2$age))^2
    dat2$bp.1s.2 = ((dat2$bp.1s - median(dat2$bp.1s))/sd(dat2$bp.1s))^2
    for (i in 1:length(dat2$glyhb)){ # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
      if (dat2$glyhb[i] >= 7){
        dat2$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat2$glyhb.cat[i] = 0
      }
    }
    dat2$w.h.ratio = (dat2$waist)/dat2$hip #waist-hip ratio
    dat2$bmi = dat2$weight/((dat2$height)^2)*703 #Add new variable BMI
    dat2$frame = factor(dat2$frame, levels = c("small", "medium", "large"))  # Change "frame" to a factor small, medium, large
    
    #dat3
    dat3$stab.glu.2 = ((dat3$stab.glu - median(dat3$stab.glu))/sd(dat3$stab.glu))^2 
    dat3$age.2 =  ((dat3$age - median(dat3$age))/sd(dat3$age))^2
    dat3$bp.1s.2 = ((dat3$bp.1s - median(dat3$bp.1s))/sd(dat3$bp.1s))^2
    for (i in 1:length(dat3$glyhb)){ # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
      if (dat3$glyhb[i] >= 7){
        dat3$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat3$glyhb.cat[i] = 0
      }
    }
    dat3$w.h.ratio = (dat3$waist)/dat3$hip #waist-hip ratio
    dat3$bmi = dat3$weight/((dat3$height)^2)*703 #Add new variable BMI
    dat3$frame = factor(dat3$frame, levels = c("small", "medium", "large"))  # Change "frame" to a factor small, medium, large
    
    #dat4
    dat4$stab.glu.2 = ((dat4$stab.glu - median(dat4$stab.glu))/sd(dat4$stab.glu))^2 
    dat4$age.2 =  ((dat4$age - median(dat4$age))/sd(dat4$age))^2
    dat4$bp.1s.2 = ((dat4$bp.1s - median(dat4$bp.1s))/sd(dat4$bp.1s))^2
    for (i in 1:length(dat4$glyhb)){ # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
      if (dat4$glyhb[i] >= 7){
        dat4$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat4$glyhb.cat[i] = 0
      }
    }
    dat4$w.h.ratio = (dat4$waist)/dat4$hip #waist-hip ratio
    dat4$bmi = dat4$weight/((dat4$height)^2)*703 #Add new variable BMI
    dat4$frame = factor(dat4$frame, levels = c("small", "medium", "large"))  # Change "frame" to a factor small, medium, large
    
    #dat5
    dat5$stab.glu.2 = ((dat5$stab.glu - median(dat5$stab.glu))/sd(dat5$stab.glu))^2 
    dat5$age.2 =  ((dat5$age - median(dat5$age))/sd(dat5$age))^2
    dat5$bp.1s.2 = ((dat5$bp.1s - median(dat5$bp.1s))/sd(dat5$bp.1s))^2
    for (i in 1:length(dat5$glyhb)){ # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
      if (dat5$glyhb[i] >= 7){
        dat5$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat5$glyhb.cat[i] = 0
      }
    }
    dat5$w.h.ratio = (dat5$waist)/dat5$hip #waist-hip ratio
    dat5$bmi = dat5$weight/((dat5$height)^2)*703 #Add new variable BMI
    dat5$frame = factor(dat5$frame, levels = c("small", "medium", "large"))  # Change "frame" to a factor small, medium, large
    
    
    fitimp2.2 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                    + bp.1s + gender*bmi + w.h.ratio,
                    data = dat2, family = binomial() )
    sumfitimp2.2 = summary(fitimp2.2); 
    sumfitimp2.2
    plot(fitimp2.2)
    
    fitimp2.3 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                    + bp.1s + gender*bmi + w.h.ratio,
                    data = dat3, family = binomial() )
    sumfitimp2.3 = summary(fitimp2.3); 
    sumfitimp2.3
    plot(fitimp2.3)
    
    fitimp2.4 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                    + bp.1s + gender*bmi + w.h.ratio,
                    data = dat4, family = binomial() )
    sumfitimp2.4 = summary(fitimp2.4); 
    sumfitimp2.4
    plot(fitimp2.4)
    
    
    fitimp2.5 = glm(glyhb.cat ~ chol + stab.glu + rcs(stab.glu, 3) + hdl + age  
                    + bp.1s + gender*bmi + w.h.ratio,
                    data = dat5, family = binomial() )
    sumfitimp2.5 = summary(fitimp2.5); 
    sumfitimp2.5
    plot(fitimp2.5)
    
    #Pool Paramters for 5 datasets
    paramdat1 = fit_08$coefficients 
    paramdat2 = fitimp2.2$coefficients 
    paramdat3 = fitimp2.3$coefficients 
    paramdat4 = fitimp2.4$coefficients 
    paramdat5 = fitimp2.5$coefficients 
    
    paramdat1
    paramdat2
    paramdat3
    paramdat4
    paramdat5
    
    all.equal(paramdat5, paramdat4, paramdat3, paramdat2, paramdat1) #returns TRUE, all equal
    
   
    