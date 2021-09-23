# Diabetes

# Note: glyhb > 7 is diagnosed diabetes
#fit 9 is final model, not all coefficients are signif (just like in notes)
#fit 12 is final model without 2 outliers and no ratio variable using dat11

#need to check assumptions and validate and review checklist

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

    #dat1 = complete(imp2,1)  #contains 1 datasets containing imputed values, choose 1st dataset
    
    # Analyze the distribution with and without imputed values for each of 5 imputed datasets
    newdat_imp_1 = complete(imp2,1)
    newdat_imp_2 = complete(imp2,2)
    newdat_imp_3 = complete(imp2,3)
    newdat_imp_4 = complete(imp2,4)
    newdat_imp_5 = complete(imp2,5)
    
    
    par(mfrow=c(2,7))
    
    # for (i in ncol(newdat_imp_1)){
    #   if is.integer()
    #   densityplot(newdat_imp_1[, i])
    # }
    
    par(mfrow=c(2,7))
    densityplot(imp2)
    
    #xyplot(newdat_imp_1,glyhb ~ chol+stab.glu+hdl+ratio+age+height+weight+bp.1s+bp.1d+waist+hip+time.ppn,pch=18,cex=1)
    #plot(dat1[, 13], )
    
# Change categorical variables to 0 or 1
    
    dat = newdat_imp_1

    # Change glyhb to categorical: glyhb>7 = 1 , glyhb<7 = 0
    
    for (i in 1:length(dat$glyhb)){
      if (dat$glyhb[i] >= 7){
        dat$glyhb.cat[i] = 1
      } else if (dat$glyhb[i] < 7){
        dat$glyhb.cat[i] = 0
      }
    }

    #Add new variable BMI
    dat$bmi = dat$weight/((dat$height)^2)*703
    
    # Change "frame" to a factor small, medium, large
    dat$frame = factor(dat$frame, levels = c("small", "medium", "large"))
    
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
    
###################################################    
######## Fit a model with random effect for ID
###################################################    
    # nonlinear stab.glu, age, bp.1s
    # borderline non-linear: hip, weight, hdl, bmi
    library(lme4)
    
    #rescale
    dat$stab.glu.2 = ((dat$stab.glu - median(dat$stab.glu))/sd(dat$stab.glu))^2  #center before squaring
    dat$age.2 =  ((dat$age - median(dat$age))/sd(dat$age))^2
    dat$bp.1s.2 = ((dat$bp.1s - median(dat$bp.1s))/sd(dat$bp.1s))^2
    #dat$bmi = dat$bmi - median(dat$bmi)
    
    #calculate waist to hip ratio
    dat$w.h.ratio = (dat$waist)/dat$hip
    
################################################################################
    # new model with id as random intercept
    # compare results to previous model fit 9
    #fit models - backwards elimination
    fit00 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age + age.2 + height 
                  + weight + frame + bp.1s + bp.1s.2 + bp.1d + waist*hip + gender*bmi + location + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit00 = summary(fit00); 
    sumfit00
    
    #remove bp.1d
    fit01 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age + age.2 + height 
                  + weight + frame + bp.1s + bp.1s.2 + waist*hip + gender*bmi + location + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit01 = summary(fit01); 
    sumfit01
    
    #remove location
    fit02 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age + age.2 + height 
                  + weight + frame + bp.1s + bp.1s.2 + waist*hip + gender*bmi + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit02 = summary(fit02); 
    sumfit02
    
    #remove frame
    fit03 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age + age.2 + height 
                  + weight + bp.1s + bp.1s.2 + waist*hip + gender*bmi + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit03 = summary(fit03); 
    sumfit03
    
    #remove waist*hip
    fit04 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age + age.2 + height 
                  + weight + bp.1s + bp.1s.2  + gender*bmi + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit04 = summary(fit04); 
    sumfit04
    
    #remove age.2
    fit05 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age 
                  + height + weight + bp.1s + bp.1s.2  + gender*bmi + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit05 = summary(fit05); 
    sumfit05
    
    #remove weight
    fit06 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio 
                  + height + age + bp.1s + bp.1s.2  + gender*bmi + w.h.ratio
                  + (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit06 = summary(fit06); 
    sumfit06
    
    #remove height
    fit07 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age
                  + bp.1s + bp.1s.2 + gender*bmi + w.h.ratio
                  +  (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit07 = summary(fit07); 
    sumfit07
    
    
    #remove bp.1s.2
    fit08 = glmer(glyhb.cat ~ chol + stab.glu + stab.glu.2 + hdl + ratio + age
                  + bp.1s + gender*bmi + w.h.ratio
                  +  (1 | id), data = dat, family = binomial, 
                  control = glmerControl(optimizer = "bobyqa"),
                  nAGQ = 1)
    sumfit08 = summary(fit08); 
    sumfit08
    
    #Likelihood ratio test to compare models
    resid.dev_full = 136.3
    df_full = 367
    resid.dev_red = 145.9
    df_red = 390
    test_stat = resid.dev_red - resid.dev_full
    df = df_red - df_full
    pchisq(test_stat, df, lower.tail = F) #p-value=0.9934 use reduced model
    
    
################################################################################
    # #remove chol and hdl & height and weight for bmi
    # fit2 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2 + weight + gender + frame + bp.1s.2 
    #              + bp.1d + waist + hip + bmi + (1 |location), 
    #              data = dat, family = binomial, 
    #             control = glmerControl(optimizer = "bobyqa"),
    #             nAGQ = 1)
    # summary(fit2)
    # 
    # #remove gender
    # fit3 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  + height 
    #              + weight + frame + bp.1s.2 + bp.1d + waist + hip + bmi +
    #                (1 |location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit3)
    # 
    # #remove hip
    # fit4 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  + height 
    #              + weight + frame + bp.1s.2 + bp.1d + waist + bmi +
    #                (1|location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit4)
    # 
    # #remove height
    # fit5 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  + frame 
    #              + weight  + bp.1s.2 + bp.1d + waist + bmi +
    #                (1 |location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit5)
    # 
    # #remove frame
    # fit6 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  
    #              + weight + bp.1d + waist + bmi + bp.1s.2 + 
    #                (1 | location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit6)
    # 
    # #remove bp.1s.2
    # fit7 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  
    #              + weight + bp.1d + waist + bmi + 
    #                (1 | location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit7)
    # 
    # #remove bp.1d.2
    # fit8 = glmer(glyhb.cat ~ stab.glu.2 + ratio + age.2  
    #              + weight  + waist + bmi + 
    #                (1 | location), data = dat, family = binomial, 
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit8)
    # 
    # #remove age.2, and ratio, waist, weight
    # fit9 = glmer(glyhb.cat ~ stab.glu.2 + ratio + waist + weight + bmi +
    #                (1 | location), data = dat, family = binomial, #remove outliers dat[-c(134, 195),]
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit9)
    # plot(fit9)
    # 
    # inf=influence(fit9,obs=T) #decide to remove 
    # plot(inf,which="cook")
    # 
    # fit9 = glmer(glyhb.cat ~ stab.glu.2 + ratio + waist + weight + bmi +
    #                (1 | location), data = dat, family = binomial, #remove outliers dat[-c(134, 195),]
    #              control = glmerControl(optimizer = "bobyqa"),
    #              nAGQ = 1)
    # summary(fit9)
    # plot(fit9) #residual plot shows linearity assumption violated
    # 
    inf=influence(fit9,obs=T) #decide to remove observation 195
    plot(inf,which="cook")
    
    
########################################
################## Assess Model
########################################
    
      #Likelihood ratio test to compare models
      resid.dev_full = 211.2
      df_full = 385
      resid.dev_red = 214.7
      df_red = 396
      test_stat = resid.dev_red - resid.dev_full
      df = df_red - df_full
      pchisq(test_stat, df, lower.tail = F) #pvalue=0.9823 use reduced model
      
      #Residuals
      res = residuals(fit9, "pearson")
      min(res)
      
    
      #Outliers
      inf=influence(fit9,obs=T) #decide to remove observation 195
      plot(inf,which="cook")
      
      dat10 = dat[-195,]
      fit10 = glmer(glyhb.cat ~ stab.glu.2 + ratio + waist + weight + bmi +
                     (1 | location), data = dat10, family = binomial, #remove outliers dat[-c(134, 195),]
                   control = glmerControl(optimizer = "bobyqa"),
                   nAGQ = 1)
      summary(fit10)
      plot(fit10) #decide to remove observation 134 based on residuals
      
      inf=influence(fit10,obs=T) 
      plot(inf,which="cook")
    
      dat11 = dat[-c(134, 195),]
      fit11 = glmer(glyhb.cat ~ stab.glu.2 + ratio + waist + weight + bmi +
                      (1 | location), data = dat11, family = binomial, #remove outliers dat[-c(134, 195),]
                    control = glmerControl(optimizer = "bobyqa"),
                    nAGQ = 1)
      summary(fit11)
      plot(fit11)
      
      inf=influence(fit11,obs=T) 
      plot(inf,which="cook")
      
      #remove ratio
      fit12 = glmer(glyhb.cat ~ stab.glu.2 + waist + weight + bmi +
                      (1 | location), data = dat11, family = binomial, #remove outliers dat[-c(134, 195),]
                    control = glmerControl(optimizer = "bobyqa"),
                    nAGQ = 1)
      summary(fit12)
      plot(fit12)
      
      #Normality QQ Plot
      par(mfrow = c(1,1))
      qqnorm(resid(fit12))
      qqline(resid(fit12))
      
      #Overdispersion - no consensus on how to calc, cite source bookmarked
      r_df = df.residual(fit12)
      r_pearson = residuals(fit12, "pearson")
      test_stat_over = sum(r_pearson^2)
      pearson_ratio = test_stat_over/r_df #indicates possibility of underdispersion. usually
      p = pchisq(test_stat_over, df = r_df, lower.tail = F)
      p #p=0.994 reject null, dispersion is not significantly different from zero
      
      #Independence Tests
      chisq.test(newdat_imp_1[, -c(1, 7, 9, 12)]) #p<<0.001 statistically significant, dependence
      
      d_reduced = dat11[,c(3,6,7,11,15,18)]
      chisq.test(d_reduced)
      
      
