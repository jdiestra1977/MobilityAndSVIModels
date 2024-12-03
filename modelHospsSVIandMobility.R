library(tidyverse)
library(lme4)
library(plm)
library(cowplot)
library(MASS)
library(leaps)
library(caret)
library(mice)
library(glmmTMB)

allCounties <-read_csv("~/Documents/GitHub/MobilityAndSVIModels/dataModelSVI_MR_Hops.csv")
allCounties %>%
  filter(week==20) %>% arrange((SVI_State))
forMob<-allCounties %>%
  mutate(SVIGroup=cut(SVI_State,quantile(SVI_State),include.lowest=TRUE,labels=FALSE) %>% as.factor()) %>%
  dplyr::select(-mobRatioIntra,-SVI_State,-Hosp_least,-Zip,-pop2019)

forHosp<-allCounties %>%
  mutate(SVIGroup=cut(SVI_State,quantile(SVI_State),include.lowest=TRUE,labels=FALSE) %>% as.factor()) %>%
  dplyr::select(-mobRatioIntra,-SVI_State,-mobRatioOut,-Zip,-pop2019)

forSVI<-allCounties %>%
  mutate(SVIGroup=cut(SVI_State,quantile(SVI_State),include.lowest=TRUE,labels=FALSE) %>% as.factor()) %>%
  dplyr::select(-mobRatioIntra,-Hosp_least,-mobRatioOut,-Zip,-pop2019,-week)

forPopu<-allCounties %>%
  mutate(SVIGroup=cut(SVI_State,quantile(SVI_State),include.lowest=TRUE,labels=FALSE) %>% as.factor()) %>%
  dplyr::select(-mobRatioIntra,-SVI_State,-mobRatioOut,-Hosp_least,-week)

dataBySVI_groups<-forMob %>% group_by(week,SVIGroup) %>% summarise_each(mean) %>%
  left_join(forHosp %>% group_by(week,SVIGroup) %>% summarise_each(sum)) %>%
  left_join(forSVI %>% unique() %>% group_by(SVIGroup) %>% summarise_each(mean)) %>%
  left_join(forPopu %>% unique() %>% dplyr::select(-Zip) %>% group_by(SVIGroup) %>% summarise_each(sum)) %>%
  rename_at(c("mobRatioOut","Hosp_least","SVI_State","pop2019"),
            ~c("mean_MR_InSVIGroup","TotalHosp_inSVIGroup","mean_SVI_inSVIGroup","Total_popu_inSVIGroup"))

dataBySVI_groups %>%
  ggplot(aes(x=week,y=TotalHosp_inSVIGroup,color=SVIGroup,group=SVIGroup))+ theme_bw() +
  geom_line() + theme(legend.position = c(0.2,0.8))

dataBySVI_groups %>% ungroup() %>% dplyr::select(-week,-SVIGroup,-mean_SVI_inSVIGroup,-Total_popu_inSVIGroup) %>%
  pivot_longer(cols=c("mean_MR_InSVIGroup","TotalHosp_inSVIGroup")) %>%
  ggplot(aes(x=value))+ geom_histogram() + facet_wrap(~name,scales = "free") + scale_x_log10()

# Convert data to a panel data frame
data_panel_groups <- pdata.frame(dataBySVI_groups, index = c("SVIGroup", "week"))
#Since we know that SVI impacts the effect of mobility on hospitalizations and potentially
#hospitalizations, we include an interaction term between SVI and mobility ratio
# Create interaction term in the dataset
data_panel_groups <- data_panel_groups %>% mutate(week=as.numeric(week))

data_panel_groups

#I standardize mobility data, but still has some skewness
data_panel_groups_std<-data_panel_groups %>% group_by(SVIGroup) %>% 
  mutate(mean_mobRatio_std = (mean_MR_InSVIGroup - mean(mean_MR_InSVIGroup, na.rm = TRUE)) / sd(mean_MR_InSVIGroup, na.rm = TRUE))

data_panel_groups_std %>%
  ggplot(aes(x=(mean_mobRatio_std))) + geom_histogram()

## I include lags in mobility ratio, since I expect that mobility at time t influences hospitaliations on
## time t+l

data_panel_groups_std_lags<-data_panel_groups_std %>%
  arrange(SVIGroup,week) %>% group_by(SVIGroup) %>%
  mutate(mean_mobRatio_std_lag1=dplyr::lag(mean_mobRatio_std, 1),
         mean_mobRatio_std_lag2=dplyr::lag(mean_mobRatio_std, 2),
         mean_mobRatio_std_lag3=dplyr::lag(mean_mobRatio_std, 3),
         mean_mobRatio_std_lag4=dplyr::lag(mean_mobRatio_std, 4),
         mean_mobRatio_std_lag5=dplyr::lag(mean_mobRatio_std, 5),
         mean_mobRatio_std_lag6=dplyr::lag(mean_mobRatio_std, 6),
         mean_mobRatio_std_lag7=dplyr::lag(mean_mobRatio_std, 7),
         mean_mobRatio_std_lag8=dplyr::lag(mean_mobRatio_std, 8)) %>% ungroup()
names(data_panel_groups_std_lags)
#I will use stepwise regression to find the best model that could be created with
#these variables
# Drop NAs selectively for lagged terms

#I need to deal with NAs in the lags. I do not want to remove these lines. Hence,
#I will use the mice function to add extrapolated data in those NAs. Since mice
#has trouble when there is high correlation between different variables (lags), 
#I will define two data frames with different lagged terms and then I will see which
#gives me the best models.

data_panel_groups_std_lags_impar<-data_panel_groups_std_lags %>% 
  dplyr::select(-mean_MR_InSVIGroup,-mean_mobRatio_std,-mean_mobRatio_std_lag2,-mean_mobRatio_std_lag4,-mean_mobRatio_std_lag6,-mean_mobRatio_std_lag8)
data_panel_groups_std_lags_par<-data_panel_groups_std_lags %>% 
  dplyr::select(-mean_MR_InSVIGroup,-mean_mobRatio_std_lag1,-mean_mobRatio_std_lag3,-mean_mobRatio_std_lag5,-mean_mobRatio_std_lag7)

imputed_data_impar <- mice(data_panel_groups_std_lags_impar, m = 1, method = "pmm", seed = 123)
imputed_data_par <- mice(data_panel_groups_std_lags_par, m = 1, method = "pmm", seed = 123)

# cor_matrix <- cor(data_panel_groups_std_lags_par %>% dplyr::select(contains("mean_")),use = "pairwise.complete.obs")
# print(cor_matrix)

data_imputed_impar <- complete(imputed_data_impar)
data_imputed_par <- complete(imputed_data_par)

data_std_lags_complete_impar<-data_panel_groups_std_lags_impar %>% dplyr::select(week,SVIGroup,TotalHosp_inSVIGroup) %>%
  left_join(data_imputed_impar)

data_std_lags_complete_par<-data_panel_groups_std_lags_par %>% dplyr::select(week,SVIGroup,TotalHosp_inSVIGroup) %>%
  left_join(data_imputed_par)

data_std_lags_complete_both <-data_std_lags_complete_impar %>% left_join(data_std_lags_complete_par)

#These are tests with Random effects for different SVI groups
prueba_impar <- glmmTMB(TotalHosp_inSVIGroup ~ week+mean_SVI_inSVIGroup + 
                                  mean_mobRatio_std_lag1+mean_mobRatio_std_lag3+
                                  mean_mobRatio_std_lag5 + mean_mobRatio_std_lag7 +
                                  offset(log(Total_popu_inSVIGroup)) + (1 | SVIGroup), 
                                data = data_std_lags_complete_impar, 
                                family = nbinom2)

prueba_par <- glmmTMB(TotalHosp_inSVIGroup ~ week+mean_SVI_inSVIGroup + mean_mobRatio_std+
                          mean_mobRatio_std_lag2++mean_mobRatio_std_lag4+
                          mean_mobRatio_std_lag6 +  + mean_mobRatio_std_lag8 +
                          offset(log(Total_popu_inSVIGroup)) + (1 | SVIGroup), 
                        data = data_std_lags_complete_par, 
                        family = nbinom2)

prueba_both <- glmmTMB(TotalHosp_inSVIGroup ~ week+mean_SVI_inSVIGroup + mean_mobRatio_std+
                         mean_mobRatio_std_lag1+mean_mobRatio_std_lag2+mean_mobRatio_std_lag3+
                         mean_mobRatio_std_lag4+mean_mobRatio_std_lag5+
                        mean_mobRatio_std_lag6 +mean_mobRatio_std_lag7+ mean_mobRatio_std_lag8 +
                        offset(log(Total_popu_inSVIGroup)) + (1 | SVIGroup), 
                      data = data_std_lags_complete_both, 
                      family = nbinom2)

summary(prueba_impar)
summary(prueba_par)
summary(prueba_both)

step.model_impar <- stepAIC(prueba_impar, direction = "both",trace = FALSE)
step.model_par <- stepAIC(prueba_par, direction = "both",trace = FALSE)
step.model_both <- stepAIC(prueba_both, direction = "both",trace = FALSE)

AIC(step.model_impar)
AIC(step.model_par)
AIC(step.model_both)

summary(step.model_impar)
summary(step.model_par)
summary(step.model_both)

#These are tests for fixed effects on SVI groups

prueba_impar_fe <- glmmTMB(TotalHosp_inSVIGroup ~ week+ 
                          mean_mobRatio_std_lag1+mean_mobRatio_std_lag3+
                          mean_mobRatio_std_lag5 + mean_mobRatio_std_lag7 +
                          offset(log(Total_popu_inSVIGroup)) + SVIGroup, 
                        data = data_std_lags_complete_impar, 
                        family = nbinom2)

prueba_par_fe <- glmmTMB(TotalHosp_inSVIGroup ~ week+ mean_mobRatio_std+
                        mean_mobRatio_std_lag2++mean_mobRatio_std_lag4+
                        mean_mobRatio_std_lag6 +  + mean_mobRatio_std_lag8 +
                        offset(log(Total_popu_inSVIGroup)) + SVIGroup, 
                      data = data_std_lags_complete_par, 
                      family = nbinom2)

prueba_both_fe <- glmmTMB(TotalHosp_inSVIGroup ~ week+ mean_mobRatio_std+
                         mean_mobRatio_std_lag1+mean_mobRatio_std_lag2+mean_mobRatio_std_lag3+
                         mean_mobRatio_std_lag4+mean_mobRatio_std_lag5+
                         mean_mobRatio_std_lag6 +mean_mobRatio_std_lag7+ mean_mobRatio_std_lag8 +
                         offset(log(Total_popu_inSVIGroup)) + SVIGroup, 
                       data = data_std_lags_complete_both, 
                       family = nbinom2)

summary(prueba_impar_fe)
summary(prueba_par_fe)
summary(prueba_both_fe)

step.model_impar <- stepAIC(prueba_impar_fe, direction = "both",trace = FALSE)
step.model_par <- stepAIC(prueba_par_fe, direction = "both",trace = FALSE)
step.model_both <- stepAIC(prueba_both_fe, direction = "both",trace = FALSE)

AIC(step.model_impar)
AIC(step.model_par)
AIC(step.model_both)

summary(step.model_impar)
summary(step.model_par)
summary(step.model_both)

best_fixed_both<-glmmTMB(TotalHosp_inSVIGroup ~ week + mean_mobRatio_std + mean_mobRatio_std_lag2 +
                           mean_mobRatio_std_lag8 + SVIGroup + offset(log(Total_popu_inSVIGroup)), 
                         data = data_std_lags_complete_both, 
                         family = nbinom2)

summary(best_fixed_both)
#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(best_fixed_both, newdata = data_std_lags_complete_both, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- data_std_lags_complete_both %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
#    hosps_transformed1 = hosps_transformed + 1,
#    backTrans_observed = exp(hosps_transformed) - 1,
#    backTrans_fitted = exp(fitted) - 1,
#    backTrans_lower_ci = exp(lower_ci) - 1,
#    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, TotalHosp_inSVIGroup, fitted, lower_ci, upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% 
           str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = TotalHosp_inSVIGroup, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=lower_ci,ymax=upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(best_fixed_both)
residuals <- resid(best_fixed_both)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")

### Now, I add non-linear terms (exponents and interactions) and evaluate with stepwise

prueba_both_fe_big <- glmmTMB(TotalHosp_inSVIGroup ~ week + I(week^2)+ mean_mobRatio_std+I(mean_mobRatio_std^2)+
                                mean_mobRatio_std_lag1 + I(mean_mobRatio_std_lag1^2)+mean_mobRatio_std_lag2+I(mean_mobRatio_std_lag2^2)+
                              mean_mobRatio_std_lag3+I(mean_mobRatio_std_lag3^2) + mean_mobRatio_std_lag4+I(mean_mobRatio_std_lag4^2)+
                              mean_mobRatio_std_lag5+I(mean_mobRatio_std_lag5^2)+mean_mobRatio_std_lag6 +I(mean_mobRatio_std_lag6^2)+
                              mean_mobRatio_std_lag7+I(mean_mobRatio_std_lag7^2)+ mean_mobRatio_std_lag8 +I(mean_mobRatio_std_lag8^2)+
                            offset(log(Total_popu_inSVIGroup)) + SVIGroup+
                              SVIGroup*mean_mobRatio_std+SVIGroup*mean_mobRatio_std_lag1+SVIGroup*mean_mobRatio_std_lag2+
                              SVIGroup*mean_mobRatio_std_lag3+SVIGroup*mean_mobRatio_std_lag4+SVIGroup*mean_mobRatio_std_lag5+
                              SVIGroup*mean_mobRatio_std_lag6+SVIGroup*mean_mobRatio_std_lag7+SVIGroup*mean_mobRatio_std_lag8+
                              week*mean_mobRatio_std+week*mean_mobRatio_std_lag1+week*mean_mobRatio_std_lag2+
                              week*mean_mobRatio_std_lag3+week*mean_mobRatio_std_lag4+week*mean_mobRatio_std_lag5+
                              week*mean_mobRatio_std_lag6+week*mean_mobRatio_std_lag7+week*mean_mobRatio_std_lag8, 
                          data = data_std_lags_complete_both, 
                          family = nbinom2)

summary(prueba_both_fe_big)

step.model_both_big <- stepAIC(prueba_both_fe_big, direction = "both",trace = FALSE)
summary(step.model_both_big)
AIC(step.model_both_big)

prueba_both_fe_big_best <- glmmTMB(TotalHosp_inSVIGroup ~ week + I(week^2) + mean_mobRatio_std +  
                                I(mean_mobRatio_std^2) + I(mean_mobRatio_std_lag1^2) + mean_mobRatio_std_lag2 +  
                                mean_mobRatio_std_lag3 + I(mean_mobRatio_std_lag3^2) + mean_mobRatio_std_lag4 +  
                                I(mean_mobRatio_std_lag4^2) + mean_mobRatio_std_lag5 + I(mean_mobRatio_std_lag5^2) +  
                                mean_mobRatio_std_lag6 + I(mean_mobRatio_std_lag6^2) + mean_mobRatio_std_lag7 +  
                                I(mean_mobRatio_std_lag7^2) + mean_mobRatio_std_lag8 + I(mean_mobRatio_std_lag8^2) +  
                                SVIGroup + mean_mobRatio_std:SVIGroup + mean_mobRatio_std_lag2:SVIGroup +  
                                mean_mobRatio_std_lag4:SVIGroup + mean_mobRatio_std_lag5:SVIGroup +mean_mobRatio_std_lag8:SVIGroup + 
                                week:mean_mobRatio_std + week:mean_mobRatio_std_lag2 + week:mean_mobRatio_std_lag3 +
                                week:mean_mobRatio_std_lag5 + week:mean_mobRatio_std_lag6 + week:mean_mobRatio_std_lag7 + 
                                week:mean_mobRatio_std_lag8 + offset(log(Total_popu_inSVIGroup)), 
                              data = data_std_lags_complete_both, 
                              family = nbinom2)

summary(prueba_both_fe_big_best)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(prueba_both_fe_big_best, newdata = data_std_lags_complete_both, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- data_std_lags_complete_both %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
    #    hosps_transformed1 = hosps_transformed + 1,
    #    backTrans_observed = exp(hosps_transformed) - 1,
    #    backTrans_fitted = exp(fitted) - 1,
    #    backTrans_lower_ci = exp(lower_ci) - 1,
    #    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, TotalHosp_inSVIGroup, fitted, lower_ci, upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = TotalHosp_inSVIGroup, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=lower_ci,ymax=upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(prueba_both_fe_big_best)
residuals <- resid(prueba_both_fe_big_best)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")

#I will try with functions that account for the heteroscedasticity as recommended

#Now, chatGPT recommends to Avoid Transforming the Response Variable with Negative Binomial Models
#Instead of log-transforming the response variable, you can model heteroscedasticity directly.
#In glmmTMB, use the dispformula argument to allow the dispersion (variance) to depend on one or more predictors:

model_4_1<-glmmTMB(TotalHosp_inSVIGroup ~ week + I(week^2) + mean_mobRatio_std +  
                     I(mean_mobRatio_std^2) + I(mean_mobRatio_std_lag1^2) + mean_mobRatio_std_lag2 +  
                     mean_mobRatio_std_lag3 + I(mean_mobRatio_std_lag3^2) + mean_mobRatio_std_lag4 +  
                     I(mean_mobRatio_std_lag4^2) + mean_mobRatio_std_lag5 + I(mean_mobRatio_std_lag5^2) +  
                     mean_mobRatio_std_lag6 + I(mean_mobRatio_std_lag6^2) + mean_mobRatio_std_lag7 +  
                     I(mean_mobRatio_std_lag7^2) + mean_mobRatio_std_lag8 + I(mean_mobRatio_std_lag8^2) +  
                     SVIGroup + mean_mobRatio_std:SVIGroup + mean_mobRatio_std_lag2:SVIGroup +  
                     mean_mobRatio_std_lag4:SVIGroup + mean_mobRatio_std_lag5:SVIGroup +mean_mobRatio_std_lag8:SVIGroup + 
                     week:mean_mobRatio_std + week:mean_mobRatio_std_lag2 + week:mean_mobRatio_std_lag3 +
                     week:mean_mobRatio_std_lag5 + week:mean_mobRatio_std_lag6 + week:mean_mobRatio_std_lag7 + 
                     week:mean_mobRatio_std_lag8 + offset(log(Total_popu_inSVIGroup)), 
                   data = data_std_lags_complete_both, 
                   family = nbinom2,
                   dispformula = ~ mean_SVI_inSVIGroup  # Specify predictors for heteroscedasticity
)

step.model_4_1 <- stepAIC(model_4_1, direction = "both",trace = FALSE)
summary(step.model_4_1)
AIC(step.model_4_1)

model_4_1_best<-glmmTMB(TotalHosp_inSVIGroup ~ week + I(week^2) + mean_mobRatio_std +  
                          I(mean_mobRatio_std^2) + I(mean_mobRatio_std_lag1^2) + mean_mobRatio_std_lag2 +  
                          mean_mobRatio_std_lag3 + I(mean_mobRatio_std_lag3^2) + mean_mobRatio_std_lag4 +  
                          I(mean_mobRatio_std_lag4^2) + mean_mobRatio_std_lag5 + I(mean_mobRatio_std_lag5^2) +  
                          mean_mobRatio_std_lag6 + I(mean_mobRatio_std_lag6^2) + mean_mobRatio_std_lag7 +  
                          I(mean_mobRatio_std_lag7^2) + mean_mobRatio_std_lag8 + I(mean_mobRatio_std_lag8^2) +  
                          SVIGroup + mean_mobRatio_std:SVIGroup + mean_mobRatio_std_lag2:SVIGroup +  
                          mean_mobRatio_std_lag4:SVIGroup + mean_mobRatio_std_lag5:SVIGroup + 
                          mean_mobRatio_std_lag8:SVIGroup + week:mean_mobRatio_std + week:mean_mobRatio_std_lag2 + 
                          week:mean_mobRatio_std_lag3 + week:mean_mobRatio_std_lag6 + week:mean_mobRatio_std_lag7 +  
                          week:mean_mobRatio_std_lag8 + offset(log(Total_popu_inSVIGroup)), 
                   data = data_std_lags_complete_both, 
                   family = nbinom2,
                   dispformula = ~ week  # Specify predictors for heteroscedasticity
)

summary(model_4_1_best)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(model_4_1_best, newdata = data_std_lags_complete_both, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- data_std_lags_complete_both %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
    #    hosps_transformed1 = hosps_transformed + 1,
    #    backTrans_observed = exp(hosps_transformed) - 1,
    #    backTrans_fitted = exp(fitted) - 1,
    #    backTrans_lower_ci = exp(lower_ci) - 1,
    #    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, TotalHosp_inSVIGroup, fitted, lower_ci, upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = TotalHosp_inSVIGroup, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=lower_ci,ymax=upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(model_4_1_best)
residuals <- resid(model_4_1_best)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")

#Other thing we can do is, following Dongah's recommendation plus chatGPT's switch to another 
#distribution, such as gaussian or gamma

with_response_log<-data_std_lags_complete_both %>%
  mutate(hosps_transformed=log(TotalHosp_inSVIGroup+1))

log_model <- glmmTMB(
  hosps_transformed ~ week + I(week^2)+ mean_mobRatio_std+I(mean_mobRatio_std^2)+
    mean_mobRatio_std_lag1 + I(mean_mobRatio_std_lag1^2)+#mean_mobRatio_std_lag2+I(mean_mobRatio_std_lag2^2)+
    mean_mobRatio_std_lag3+I(mean_mobRatio_std_lag3^2) + mean_mobRatio_std_lag4+I(mean_mobRatio_std_lag4^2)+
    mean_mobRatio_std_lag5+I(mean_mobRatio_std_lag5^2)+mean_mobRatio_std_lag6 +I(mean_mobRatio_std_lag6^2)+
    mean_mobRatio_std_lag7+I(mean_mobRatio_std_lag7^2)+ mean_mobRatio_std_lag8 +I(mean_mobRatio_std_lag8^2)+
    offset(log(Total_popu_inSVIGroup)) + SVIGroup+
#    SVIGroup*mean_mobRatio_std+SVIGroup*mean_mobRatio_std_lag1+SVIGroup*mean_mobRatio_std_lag2+
#    SVIGroup*mean_mobRatio_std_lag3+SVIGroup*mean_mobRatio_std_lag4+SVIGroup*mean_mobRatio_std_lag5+
#    SVIGroup*mean_mobRatio_std_lag6+ #SVIGroup*mean_mobRatio_std_lag7+SVIGroup*mean_mobRatio_std_lag8+
    week*mean_mobRatio_std+week*mean_mobRatio_std_lag1+#week*mean_mobRatio_std_lag2+
    week*mean_mobRatio_std_lag3+#week*mean_mobRatio_std_lag4+week*mean_mobRatio_std_lag5+
    week*mean_mobRatio_std_lag6,#+week*mean_mobRatio_std_lag7+week*mean_mobRatio_std_lag8, 
  data = with_response_log,
  family = gaussian()  # Change family to Gaussian
)

step.model_logResponse <- stepAIC(log_model, direction = "both",trace = FALSE)
summary(step.model_logResponse)
AIC(step.model_logResponse)

log_model_best <- glmmTMB(
  hosps_transformed ~ week + I(week^2) + mean_mobRatio_std + I(mean_mobRatio_std^2) +  
    mean_mobRatio_std_lag1 + I(mean_mobRatio_std_lag1^2) + mean_mobRatio_std_lag3 +  
    I(mean_mobRatio_std_lag3^2) + I(mean_mobRatio_std_lag4^2) +  
    mean_mobRatio_std_lag5 + I(mean_mobRatio_std_lag5^2) + mean_mobRatio_std_lag6 +  
    I(mean_mobRatio_std_lag6^2) + mean_mobRatio_std_lag7 + I(mean_mobRatio_std_lag7^2) +  
    mean_mobRatio_std_lag8 + SVIGroup + mean_mobRatio_std_lag3:SVIGroup +  
    mean_mobRatio_std_lag5:SVIGroup + mean_mobRatio_std_lag6:SVIGroup +  
    week:mean_mobRatio_std + week:mean_mobRatio_std_lag1 + offset(log(Total_popu_inSVIGroup)), 
  data = with_response_log,
  family = gaussian()  # Change family to Gaussian
)

summary(log_model_best)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(log_model_best, newdata = with_response_log, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- with_response_log %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
    hosps_transformed1 = hosps_transformed + 1,
    backTrans_observed = exp(hosps_transformed) - 1,
    backTrans_fitted = exp(fitted) - 1,
    backTrans_lower_ci = exp(lower_ci) - 1,
    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, TotalHosp_inSVIGroup, fitted, lower_ci, upper_ci,backTrans_observed,
                backTrans_fitted,backTrans_lower_ci,backTrans_upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = TotalHosp_inSVIGroup, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=backTrans_lower_ci,ymax=backTrans_upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = backTrans_fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(log_model_best)
residuals <- resid(log_model_best)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")

# I want to reduce the number of variables to avoid over-fitting -------

#Then, I will statrt with varibles that I think are more relevant and their interactions,
#as well as their lagged versions

log_model_reduced <- glmmTMB(
  hosps_transformed ~ week + I(week^2)+ mean_mobRatio_std+I(mean_mobRatio_std^2)+
    mean_mobRatio_std_lag1 + I(mean_mobRatio_std_lag1^2)+mean_mobRatio_std_lag2+
    I(mean_mobRatio_std_lag2^2)+mean_mobRatio_std_lag3+I(mean_mobRatio_std_lag3^2) + 
    mean_mobRatio_std_lag4+I(mean_mobRatio_std_lag4^2)+offset(log(Total_popu_inSVIGroup)) + SVIGroup+
    SVIGroup*mean_mobRatio_std+SVIGroup*mean_mobRatio_std_lag1+SVIGroup*mean_mobRatio_std_lag2+
    SVIGroup*mean_mobRatio_std_lag3+SVIGroup*mean_mobRatio_std_lag4+
    week*mean_mobRatio_std+week*mean_mobRatio_std_lag1+week*mean_mobRatio_std_lag2+
    week*mean_mobRatio_std_lag3+week*mean_mobRatio_std_lag4, 
  data = with_response_log,
  family = gaussian()  # Change family to Gaussian
)

step.model_logResponse <- stepAIC(log_model_reduced, direction = "both",trace = FALSE)
summary(step.model_logResponse)
AIC(step.model_logResponse)

### Time restricted ------
# I am using only mobility data from the week of March 23, since we are interested in 
# what happened after restrictions, not before, where patterns are different.

dataForModel<-read_csv("~/Documents/GitHub/MobilityAndSVIModels/allCounties_for_models.csv")
labelsDataAndWeek<-dataForModel %>% dplyr::select(Date) %>% unique() %>%
  rownames_to_column(var="week")

#This is the clean data, without adding or averaging
dataForModel_clean<-dataForModel %>% left_join(labelsDataAndWeek) %>%
  dplyr::select(-Zip,-SVI,-Date)
#Mean for intra and out mobility ratio
forMobs<-dataForModel_clean %>% ungroup() %>% dplyr::select(-Hospitalized,-pop2019) %>%
  group_by(SVIGroup,week) %>% summarise_each(mean)
#Sum for population and hospitalizations
forOtherVars<-dataForModel_clean %>% ungroup() %>% dplyr::select(-mobRatioOut,-mobRatioIntra) %>%
  group_by(SVIGroup,week) %>% summarise_each(sum)

data_clean<-forMobs %>% left_join(forOtherVars) %>% 
  mutate(week=as.numeric(week),SVIGroup=as.character(SVIGroup)) %>% arrange(SVIGroup,week)

#Data clean and with standardized mobility variables
data_clean_std<-data_clean %>%
  group_by(SVIGroup) %>% 
  mutate(mean_mobRatioOut_std = (mobRatioOut - mean(mobRatioOut, na.rm = TRUE)) / sd(mobRatioOut, na.rm = TRUE)) %>%
  mutate(mean_mobRatioIntra_std = (mobRatioIntra - mean(mobRatioIntra, na.rm = TRUE)) / sd(mobRatioIntra, na.rm = TRUE)) %>%
  dplyr::select(-mobRatioOut,-mobRatioIntra) 

data_clean_std %>%
#  ggplot(aes(x=week,y=Hospitalized,group=SVIGroup,color=SVIGroup)) + 
#  ggplot(aes(x=week,y=mean_mobRatioOut_std,group=SVIGroup,color=SVIGroup)) + 
  ggplot(aes(x=week,y=mean_mobRatioIntra_std,group=SVIGroup,color=SVIGroup)) + 
  geom_line()

# Convert data to a panel data frame
data_panel_groups <- pdata.frame(data_clean_std, index = c("SVIGroup", "week"))

data_panel_groups_std_lags<-data_panel_groups %>%
  arrange(SVIGroup,week) %>% group_by(SVIGroup) %>%
  mutate(mean_mobRatioOut_std_lag1=dplyr::lag(mean_mobRatioOut_std, 1),
         mean_mobRatioOut_std_lag2=dplyr::lag(mean_mobRatioOut_std, 2),
         mean_mobRatioOut_std_lag3=dplyr::lag(mean_mobRatioOut_std, 3),
         mean_mobRatioOut_std_lag4=dplyr::lag(mean_mobRatioOut_std, 4),
         mean_mobRatioIntra_std_lag1=dplyr::lag(mean_mobRatioIntra_std, 1),
         mean_mobRatioIntra_std_lag2=dplyr::lag(mean_mobRatioIntra_std, 2),
         mean_mobRatioIntra_std_lag3=dplyr::lag(mean_mobRatioIntra_std, 3),
         mean_mobRatioIntra_std_lag4=dplyr::lag(mean_mobRatioIntra_std, 4)) %>% ungroup()

names(data_panel_groups_std_lags)

data_panel_groups_std_lags_impar<-data_panel_groups_std_lags %>% 
  dplyr::select(-mean_mobRatioOut_std,-mean_mobRatioIntra_std,-mean_mobRatioOut_std_lag2,-mean_mobRatioIntra_std_lag2,
                -mean_mobRatioOut_std_lag4,-mean_mobRatioIntra_std_lag4)
data_panel_groups_std_lags_par<-data_panel_groups_std_lags %>% 
  dplyr::select(-mean_mobRatioOut_std_lag1,-mean_mobRatioIntra_std_lag1,-mean_mobRatioOut_std_lag3,-mean_mobRatioIntra_std_lag3)

#I will use stepwise regression to find the best model that could be created with
#these variables
# Drop NAs selectively for lagged terms

#I need to deal with NAs in the lags. I do not want to remove these lines. Hence,
#I will use the mice function to add extrapolated data in those NAs. Since mice
#has trouble when there is high correlation between different variables (lags), 
#I will define two data frames with different lagged terms and then I will see which
#gives me the best models.

data_panel_groups_std_lags <- data_panel_groups_std_lags %>% mutate(week=as.factor(week))

forMobOut<-data_panel_groups_std_lags %>% 
  dplyr::select(SVIGroup,week,mean_mobRatioOut_std,mean_mobRatioOut_std_lag1,
                mean_mobRatioOut_std_lag2,mean_mobRatioOut_std_lag3,mean_mobRatioOut_std_lag4)

forMobOut_imputed <- mice(forMobOut,method ="pmm",m=1,seed=123)

forMobIntra<-data_panel_groups_std_lags %>% dplyr::select(SVIGroup,week,mean_mobRatioIntra_std,mean_mobRatioIntra_std_lag1,mean_mobRatioIntra_std_lag2,
                                                   mean_mobRatioIntra_std_lag3,mean_mobRatioIntra_std_lag4)
forMobIntra_imputed <- mice(forMobIntra,method ="pmm",m=1,seed=123)

data_imputed_forMobOut <- complete(forMobOut_imputed) %>% as_tibble()
data_imputed_forMobIntra <- complete(forMobIntra_imputed) %>% as_tibble()

#These is the data I will use for models. It has been standardized, lagged, imputed to 
#calculate values where NA.
data_std_lags_complete<-data_panel_groups_std_lags %>% dplyr::select(SVIGroup,week,Hospitalized,pop2019) %>%
  left_join(data_imputed_forMobOut) %>% left_join(data_imputed_forMobIntra) %>% mutate(week=as.numeric(week))
names(data_std_lags_complete)
#These are tests with Random effects for different SVI groups
prueba_complete <- glmmTMB(Hospitalized ~ week + mean_mobRatioOut_std +
                          mean_mobRatioOut_std_lag1+mean_mobRatioOut_std_lag2+mean_mobRatioOut_std_lag3+
                          mean_mobRatioOut_std_lag4+mean_mobRatioIntra_std+mean_mobRatioIntra_std_lag1+
                          mean_mobRatioIntra_std_lag2+mean_mobRatioIntra_std_lag3+mean_mobRatioIntra_std_lag4+
                          offset(log(pop2019)) + (1 | SVIGroup), 
                        data = data_std_lags_complete, 
                        family = nbinom2)

summary(prueba_complete)

step.model_complete <- stepAIC(prueba_complete, direction = "both",trace = FALSE)

AIC(step.model_complete)

summary(step.model_complete)

#This one is with fixed effects for SVIGroup
prueba_complete_fixed <- glmmTMB(Hospitalized ~ week + I(week^2) + 
                                   mean_mobRatioOut_std + I(mean_mobRatioOut_std^2)+mean_mobRatioIntra_std+I(mean_mobRatioIntra_std^2)+
                                   SVIGroup +
                                   mean_mobRatioOut_std_lag1+mean_mobRatioOut_std_lag2+mean_mobRatioOut_std_lag3+mean_mobRatioOut_std_lag4+
                                   SVIGroup*mean_mobRatioOut_std+SVIGroup*mean_mobRatioIntra_std+
                                   mean_mobRatioIntra_std_lag1+mean_mobRatioIntra_std_lag2+mean_mobRatioIntra_std_lag3+mean_mobRatioIntra_std_lag4+
                             offset(log(pop2019)), 
                           data = data_std_lags_complete, 
                           family = nbinom2)

summary(prueba_complete_fixed)

step.model_complete_fixed <- stepAIC(prueba_complete_fixed, direction = "both",trace = FALSE)

AIC(step.model_complete_fixed)
summary(step.model_complete_fixed)

#Now, following Dongah's recommendation I will log-transform the response variable
#and due to this, I have to switch from a negative binomial (needs discrete values)
#to a Gaussian family distribution

library(tidyverse)
library(lme4)
library(plm)
library(cowplot)
library(MASS)
library(leaps)
library(caret)
library(mice)
library(glmmTMB)

# data_std_lags_complete_logResponse<-data_std_lags_complete %>%
#   mutate(hosps_transformed=log(Hospitalized+1))
#write_csv(data_std_lags_complete_logResponse,file="~/Documents/GitHub/MobilityAndSVIModels/data_std_lags_complete_logResponse.csv")

data_std_lags_complete_logResponse<-read_csv("~/Documents/GitHub/MobilityAndSVIModels/data_std_lags_complete_logResponse.csv")
data_std_lags_complete_logResponse <- data_std_lags_complete_logResponse %>%
  mutate(SVIGroup=as.factor(SVIGroup))

log_model_complete_fixed <- glmmTMB(
  hosps_transformed ~ week + I(week^2) + 
    mean_mobRatioOut_std + #I(mean_mobRatioOut_std^2)+
    mean_mobRatioIntra_std+#I(mean_mobRatioIntra_std^2)+
    SVIGroup + #I(mean_mobRatioOut_std_lag2^2)+ #I(mean_mobRatioOut_std_lag2^3)+#I(mean_mobRatioOut_std_lag1^2)+
    #I(mean_mobRatioIntra_std_lag2^2)+I(mean_mobRatioIntra_std_lag2^3)+
    #I(mean_mobRatioOut_std_lag3^3) +
    I(mean_mobRatioOut_std_lag3^2)+I(mean_mobRatioIntra_std_lag1^2)+
    #I(mean_mobRatioOut_std_lag4^3) +
    mean_mobRatioOut_std_lag1+mean_mobRatioOut_std_lag2+mean_mobRatioOut_std_lag3+#mean_mobRatioOut_std_lag4+
#    SVIGroup*mean_mobRatioOut_std+SVIGroup*mean_mobRatioOut_std_lag1+SVIGroup*mean_mobRatioOut_std_lag2+
    #SVIGroup*mean_mobRatioIntra_std+
    SVIGroup*mean_mobRatioOut_std_lag3+week*mean_mobRatioOut_std_lag3+
#  I(mean_mobRatioOut_std_lag1^4)+I(mean_mobRatioOut_std_lag2^4)+I(mean_mobRatioOut_std_lag3^4) +I(mean_mobRatioOut_std_lag4^4)+
#  I(mean_mobRatioOut_std_lag1^5)+I(mean_mobRatioOut_std_lag2^5)+I(mean_mobRatioOut_std_lag3^5) +I(mean_mobRatioOut_std_lag4^5)+
  mean_mobRatioIntra_std_lag1+mean_mobRatioIntra_std_lag2+mean_mobRatioIntra_std_lag3+#mean_mobRatioIntra_std_lag4+
    offset(log(pop2019)), 
  data = data_std_lags_complete_logResponse,
  family = gaussian()  # Change family to Gaussian
)

step.model_logResponse <- stepAIC(log_model_complete_fixed, direction = "both",trace = FALSE)
summary(step.model_logResponse)
AIC(step.model_logResponse)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(step.model_logResponse, newdata = data_std_lags_complete_logResponse, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- data_std_lags_complete_logResponse %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
    hosps_transformed1 = hosps_transformed + 1,
    backTrans_observed = exp(hosps_transformed) - 1,
    backTrans_fitted = exp(fitted) - 1,
    backTrans_lower_ci = exp(lower_ci) - 1,
    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, Hospitalized, fitted, lower_ci, upper_ci,backTrans_observed,
                backTrans_fitted,backTrans_lower_ci,backTrans_upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = Hospitalized, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=backTrans_lower_ci,ymax=backTrans_upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = backTrans_fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(step.model_logResponse)
residuals <- resid(step.model_logResponse)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")

#### Last model I tested -------

log_model_complete_fixed <- glmmTMB(
  hosps_transformed ~ week + I(week^2) + SVIGroup + 
    mean_mobRatioOut_std + mean_mobRatioIntra_std+
    mean_mobRatioOut_std_lag1+mean_mobRatioOut_std_lag2+mean_mobRatioOut_std_lag3+mean_mobRatioOut_std_lag4+
    mean_mobRatioIntra_std_lag1+mean_mobRatioIntra_std_lag2+mean_mobRatioIntra_std_lag3+mean_mobRatioIntra_std_lag4+
    I(mean_mobRatioOut_std_lag2^2)+ I(mean_mobRatioOut_std_lag2^4)+#I(mean_mobRatioOut_std_lag1^2)+
    #I(mean_mobRatioIntra_std_lag2^2)+I(mean_mobRatioIntra_std_lag2^3)+
    I(mean_mobRatioOut_std_lag3^5) +
 #   I(mean_mobRatioOut_std_lag3^2)+I(mean_mobRatioIntra_std_lag1^2)+
    #I(mean_mobRatioOut_std_lag4^3) +
#    SVIGroup*mean_mobRatioOut_std+SVIGroup*mean_mobRatioOut_std_lag1+SVIGroup*mean_mobRatioOut_std_lag2+
    #SVIGroup*mean_mobRatioIntra_std+
    SVIGroup*mean_mobRatioOut_std_lag3+#week*mean_mobRatioOut_std_lag3+
    #  I(mean_mobRatioOut_std_lag1^4)+I(mean_mobRatioOut_std_lag2^4)+I(mean_mobRatioOut_std_lag3^4) +I(mean_mobRatioOut_std_lag4^4)+
    #  I(mean_mobRatioOut_std_lag1^5)+I(mean_mobRatioOut_std_lag2^5)+I(mean_mobRatioOut_std_lag3^5) +I(mean_mobRatioOut_std_lag4^5)+
    offset(log(pop2019)), 
  data = data_std_lags_complete_logResponse,
  family = gaussian()  # Change family to Gaussian
)

step.model_logResponse <- stepAIC(log_model_complete_fixed, direction = "both",trace = FALSE)
summary(step.model_logResponse)
AIC(step.model_logResponse)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(step.model_logResponse, newdata = data_std_lags_complete_logResponse, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- data_std_lags_complete_logResponse %>%
  mutate(
    fitted = predictions$fit,
    se.fit = predictions$se.fit,
    lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
    upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
    hosps_transformed1 = hosps_transformed + 1,
    backTrans_observed = exp(hosps_transformed) - 1,
    backTrans_fitted = exp(fitted) - 1,
    backTrans_lower_ci = exp(lower_ci) - 1,
    backTrans_upper_ci = exp(upper_ci) - 1
  ) %>%
  dplyr::select(week, SVIGroup, Hospitalized, fitted, lower_ci, upper_ci,backTrans_observed,
                backTrans_fitted,backTrans_lower_ci,backTrans_upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_line(aes(y = Hospitalized, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=backTrans_lower_ci,ymax=backTrans_upper_ci),alpha=0.2,fill="red")+
  geom_line(aes(y = backTrans_fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(step.model_logResponse)
residuals <- resid(step.model_logResponse)
#confint(simpler_model_with_SVIGroup)

#Residuals vs. Fitted Plot
fig1<-ggplot(data = data.frame(fitted = fitted_values, residuals = residuals), aes(x = fitted, y = residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  labs(title = "Residuals vs Fitted", x = "Fitted Values", y = "Residuals") +
  theme_minimal()
#Q-Q Plot
# Create a data frame for the plot
fig2<-ggplot(data = data.frame(sample = residuals), aes(sample = sample)) +
  stat_qq() +
  stat_qq_line(color = "red", linetype = "dashed") +
  labs(title = "Q-Q Plot of Residuals", x = "Theoretical Quantiles", y = "Sample Quantiles") +
  theme_minimal()
#Residual Density Plot
fig3<-ggplot(data = data.frame(residuals = residuals), aes(x = residuals)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "skyblue", color = "black") +
  geom_density(color = "red", size = 1) +
  labs(title = "Residual Density Plot", x = "Residuals", y = "Density") +
  theme_minimal()
#Scale-Location Plot
std_residuals <- residuals / sd(residuals)  # Standardized residuals
fig4<-ggplot(data = data.frame(fitted = fitted_values, std_residuals = abs(std_residuals)), aes(x = fitted, y = std_residuals)) +
  geom_point() +
  geom_smooth(method = "loess", color = "blue", se = FALSE) +
  labs(title = "Scale-Location Plot", x = "Fitted Values", y = "Standardized Residuals") +
  theme_minimal()

plot_grid(fig1,fig2,fig3,fig4,ncol=2)

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/diagnoseModelBig.png")



