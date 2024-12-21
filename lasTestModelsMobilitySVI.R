library(tidyverse)
library(lme4)
library(plm)
library(cowplot)
library(MASS)
library(leaps)
library(caret)
library(mice)
library(glmmTMB)

#allCounties <-read_csv("~/Documents/GitHub/MobilityAndSVIModels/dataModelSVI_MR_Hops.csv")

dataForModel<-read_csv("~/Documents/GitHub/MobilityAndSVIModels/allCounties_for_models.csv")
labelsDataAndWeek<-dataForModel %>% dplyr::select(Date) %>% unique() %>%
  rownames_to_column(var="week")

dataWithSVIValues<-dataForModel %>% dplyr::select(-Zip) %>%
  left_join(labelsDataAndWeek) %>% dplyr::select(-Date) %>%
  mutate(week=as.numeric(week))

# dataWithSVIValues<-dataForModel %>% dplyr::select(-Zip,-SVIGroup) %>%
#   left_join(labelsDataAndWeek) %>% dplyr::select(-Date) %>%
#   mutate(week=as.numeric(week))

dataWithSVIValues %>%
  pivot_longer(cols = c("SVI","mobRatioOut" ,"mobRatioIntra")) %>%
  mutate(hosps_log=log(Hospitalized+1)) %>%
  ggplot(aes(x=value,y=Hospitalized))+
  geom_point() + facet_wrap(~name,scales = "free_x") +
  scale_x_log10() + geom_smooth(method = "lm") 

dataWithSVIValues %>% 
  dplyr::select(contains("mobRa")) %>%
  pivot_longer(cols = contains("mobRa")) %>%
  ggplot(aes(x=value,fill=name))+
  geom_density(alpha=0.2) + scale_x_log10()

dataWithSVIValues_std<-dataWithSVIValues %>% 
  mutate(mobRatioOut_std = (mobRatioOut - mean(mobRatioOut, na.rm = TRUE)) / sd(mobRatioOut, na.rm = TRUE),
         mobRatioIntra_std = (mobRatioIntra - mean(mobRatioIntra, na.rm = TRUE)) / sd(mobRatioIntra, na.rm = TRUE),
         mobRatioOut_log=log(mobRatioOut),mobRatioIntra_log=log(mobRatioIntra),
         SVIGroup=as.factor(SVIGroup))

dataWithSVIValues_std %>% 
  dplyr::select(contains("_std")) %>%
  pivot_longer(cols = contains("_std")) %>%
  ggplot(aes(x=value,fill=name))+
  geom_density()

dataWithSVIValues_std %>% 
  dplyr::select(contains("_log")) %>%
  pivot_longer(cols = contains("_log")) %>%
  ggplot(aes(x=value,fill=name))+
  geom_density()

modelWithSVIValues <- glmmTMB(Hospitalized ~ week+SVI+SVIGroup+mobRatioOut_std+mobRatioIntra_std +
                         offset(log(pop2019)), 
                       data = dataWithSVIValues_std, 
                       family = nbinom2)

summary(modelWithSVIValues)

modelWithSVIValues <- glmmTMB(Hospitalized ~ week+SVI+#SVIGroup+
                                mobRatioOut_log+mobRatioIntra_log +
                                offset(log(pop2019)), 
                              data = dataWithSVIValues_std, 
                              family = nbinom2)

summary(modelWithSVIValues)

dataWithSVIValues_std_logs_lags<-dataWithSVIValues_std %>%
  arrange(SVI,week) %>% group_by(SVI) %>%
  mutate(mobRatioOut_std_lag1=dplyr::lag(mobRatioOut_std, 1),mobRatioOut_std_lag2=dplyr::lag(mobRatioOut_std, 2),
         mobRatioOut_std_lag3=dplyr::lag(mobRatioOut_std, 3),
         mobRatioIntra_std_lag1=dplyr::lag(mobRatioIntra_std, 1),mobRatioIntra_std_lag2=dplyr::lag(mobRatioIntra_std, 2),
         mobRatioIntra_std_lag3=dplyr::lag(mobRatioIntra_std, 3),
         mobRatioOut_log_lag1=dplyr::lag(mobRatioOut_log, 1),mobRatioOut_log_lag2=dplyr::lag(mobRatioOut_log, 2),
         mobRatioOut_log_lag3=dplyr::lag(mobRatioOut_log, 3),
         mobRatioIntra_log_lag1=dplyr::lag(mobRatioIntra_log, 1),mobRatioIntra_log_lag2=dplyr::lag(mobRatioIntra_log, 2),
         mobRatioIntra_log_lag3=dplyr::lag(mobRatioIntra_log, 3)) %>% ungroup()

dataWithSVIValues_std_logs_lags_prueba<-dataWithSVIValues_std_logs_lags %>%
  mutate(Hospitalized_log=log(Hospitalized+1))%>% drop_na()

dataWithSVIValues_std_logs_lags_prueba %>% glimpse()

# imputed_data_impar <- mice(data_panel_groups_std_lags_impar, m = 1, method = "pmm", seed = 123)
# imputed_data_par <- mice(data_panel_groups_std_lags_par, m = 1, method = "pmm", seed = 123)
# 
# # cor_matrix <- cor(data_panel_groups_std_lags_par %>% dplyr::select(contains("mean_")),use = "pairwise.complete.obs")
# # print(cor_matrix)
# 
# data_imputed_impar <- complete(imputed_data_impar)
# data_imputed_par <- complete(imputed_data_par)
# 
# data_std_lags_complete_impar<-data_panel_groups_std_lags_impar %>% dplyr::select(week,SVIGroup,TotalHosp_inSVIGroup) %>%
#   left_join(data_imputed_impar)
# 
# names(dataWithSVIValues_std_logs_lags)

modelWithSVIValues <- glmmTMB(Hospitalized ~ week+SVI+#SVIGroup+#mobRatioOut+mobRatioIntra+
                                mobRatioOut_std+mobRatioIntra_std+mobRatioOut_log+mobRatioIntra_log+
                                mobRatioOut_std_lag1+mobRatioOut_std_lag2+mobRatioOut_std_lag3+
                                mobRatioIntra_std_lag1+mobRatioIntra_std_lag2+mobRatioIntra_std_lag3+
                                mobRatioOut_log_lag1+mobRatioOut_log_lag2+mobRatioOut_log_lag3+
                                mobRatioIntra_log_lag1+mobRatioIntra_log_lag2+mobRatioIntra_log_lag3+
                                offset(log(pop2019)), 
                              data = dataWithSVIValues_std_logs_lags_prueba, 
                              family = nbinom2)

step.model_rolling <- stepAIC(modelWithSVIValues, direction = "both",trace = FALSE)
summary(step.model_rolling)
AIC(step.model_rolling)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(step.model_rolling, newdata = dataWithSVIValues_std_logs_lags_prueba, type = "response", se.fit = TRUE)

# Add predictions and confidence intervals to the data
plot_data_fixed_effects <- dataWithSVIValues_std_logs_lags_prueba %>%
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
  dplyr::select(week, SVIGroup, Hospitalized, fitted, lower_ci, upper_ci) %>%#,backTrans_observed),
#                backTrans_fitted,backTrans_lower_ci,backTrans_upper_ci) %>%
  mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))

blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")

plot_data_fixed_effects %>% dplyr::select(week,Hospitalized,fitted) %>%
  pivot_longer(cols=-week) %>%
  ggplot(aes(x=week,y=value,color=name)) +
  geom_point() + geom_smooth()

ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
  geom_point(aes(y = Hospitalized, color = "Observed"),linetype=2,linewidth=1.5) +
  geom_ribbon(aes(ymin=lower_ci,ymax=upper_ci),alpha=0.2,fill="red")+
  geom_point(aes(y = fitted, color = "Fitted"),linewidth=1.5) +
  facet_wrap(~ SVIGroup1, scales = "free_y") +
  labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
  theme(legend.position = "bottom",text=element_text(size=25))+
  scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))

ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")

# Fitted values and residuals for random effects model
fitted_values <- fitted(step.model_rolling)
residuals <- resid(step.model_rolling)
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

#Transforming y-variable to avoid heteroscedisticity 
library(zoo)

dataWithSVIValues_std_logs_lags_prueba<-dataWithSVIValues_std_logs_lags_prueba %>%
  mutate(mobRatioOut_log_rollMean=rollmean(mobRatioOut_log,k=3,fill=NA),
         mobRatioOut_log_rollMean_lag1=rollmean(mobRatioOut_log_lag1,k=3,fill=NA),
         mobRatioOut_log_rollMean_lag2=rollmean(mobRatioOut_log_lag2,k=3,fill=NA),
         mobRatioOut_log_rollMean_lag3=rollmean(mobRatioOut_log_lag3,k=3,fill=NA),
         mobRatioIntra_log_rollMean=rollmean(mobRatioIntra_log,k=3,fill=NA),
         mobRatioIntra_log_rollMean_lag1=rollmean(mobRatioIntra_log_lag1,k=3,fill=NA),
         mobRatioIntra_log_rollMean_lag2=rollmean(mobRatioIntra_log_lag2,k=3,fill=NA),
         mobRatioIntra_log_rollMean_lag3=rollmean(mobRatioIntra_log_lag3,k=3,fill=NA)) %>% drop_na()

modelWithSVIValues <- glmmTMB(Hospitalized ~ week+SVI+SVIGroup+mobRatioOut+mobRatioIntra+
                                #mobRatioOut_std+mobRatioIntra_std+
                                mobRatioOut_log+mobRatioIntra_log+
                                mobRatioOut_std_lag1+mobRatioOut_std_lag2+mobRatioOut_std_lag3+
                                mobRatioIntra_std_lag1+mobRatioIntra_std_lag2+mobRatioIntra_std_lag3+
                                mobRatioOut_log_lag1+mobRatioOut_log_lag2+mobRatioOut_log_lag3+
                                mobRatioIntra_log_lag1+mobRatioIntra_log_lag2+mobRatioIntra_log_lag3+
                                mobRatioOut_log_rollMean+mobRatioOut_log_rollMean_lag1+mobRatioOut_log_rollMean_lag2+mobRatioOut_log_rollMean_lag3+
                                mobRatioIntra_log_rollMean+mobRatioIntra_log_rollMean_lag1+mobRatioIntra_log_rollMean_lag2+mobRatioIntra_log_rollMean_lag3+
                                offset(log(pop2019)), 
                              data = dataWithSVIValues_std_logs_lags_prueba, 
                              family = nbinom2)

step.model_rolling <- stepAIC(modelWithSVIValues, direction = "both",trace = FALSE)
summary(step.model_rolling)
AIC(step.model_rolling)

log_model_rolling <- glmmTMB(
  Hospitalized_log ~ week+SVI+#SVIGroup+mobRatioOut+mobRatioIntra+
    #mobRatioOut_std+mobRatioIntra_std+
    mobRatioOut_log+mobRatioIntra_log+
    mobRatioOut_std_lag1+mobRatioOut_std_lag2+mobRatioOut_std_lag3+
    mobRatioIntra_std_lag1+mobRatioIntra_std_lag2+mobRatioIntra_std_lag3+
    mobRatioOut_log_lag1+mobRatioOut_log_lag2+mobRatioOut_log_lag3+
    mobRatioIntra_log_lag1+mobRatioIntra_log_lag2+mobRatioIntra_log_lag3+
    mobRatioOut_log_rollMean+mobRatioOut_log_rollMean_lag1+mobRatioOut_log_rollMean_lag2+mobRatioOut_log_rollMean_lag3+
    mobRatioIntra_log_rollMean+mobRatioIntra_log_rollMean_lag1+mobRatioIntra_log_rollMean_lag2+mobRatioIntra_log_rollMean_lag3+
    I(mobRatioOut_log_lag1^2)+I(mobRatioOut_log_lag1^3)+
    I(mobRatioOut_log_lag2^2)+I(mobRatioOut_log_lag2^3)+
    I(mobRatioOut_log_lag3^2)+I(mobRatioOut_log_lag3^3)+
#    I(mobRatioIntra_log_lag1^2)+I(mobRatioIntra_log_lag1^3)+
#    I(mobRatioIntra_log_lag2^2)+I(mobRatioIntra_log_lag2^3)+
#    I(mobRatioIntra_log_lag3^2)+I(mobRatioIntra_log_lag3^3)+
    offset(log(pop2019)),#family = Gamma(link = "log"), 
  family = gaussian(),  # Change family to Gaussian
  data = dataWithSVIValues_std_logs_lags_prueba
)

step.model_rolling_log <- stepAIC(log_model_rolling, direction = "both",trace = FALSE)
summary(step.model_rolling_log)

library(splines)
library(car)
log_model_rolling <- glmmTMB(
  Hospitalized_log ~ week+SVI+#SVIGroup+#mobRatioOut+mobRatioIntra+
    #mobRatioOut_std+mobRatioIntra_std+
    #mobRatioOut_log+
#    mobRatioIntra_log+
    mobRatioOut_log_rollMean+mobRatioOut_log_rollMean_lag1+mobRatioOut_log_rollMean_lag2+mobRatioOut_log_rollMean_lag3+
    mobRatioOut_log_lag1+mobRatioOut_log_lag2+mobRatioOut_log_lag3+
    I(mobRatioOut_log_lag1^2)+I(mobRatioOut_log_lag1^3)+
    I(mobRatioOut_log_lag2^2)+I(mobRatioOut_log_lag2^3)+
#    I(mobRatioOut_log_lag3^2)+I(mobRatioOut_log_lag3^3)+
    mobRatioIntra_log_lag1+mobRatioIntra_log_lag2+mobRatioIntra_log_lag3+
#    bs(mobRatioOut_log, degree = 3)+
    bs(mobRatioOut_log_lag2, degree = 2)+
#    I(mobRatioIntra_log_lag1^2)+I(mobRatioIntra_log_lag1^3)+
#    I(mobRatioIntra_log_lag2^2)+I(mobRatioIntra_log_lag2^3)+
#    I(mobRatioIntra_log_lag3^2)+I(mobRatioIntra_log_lag3^3)+
    offset(log(pop2019)),#family = Gamma(link = "log"), 
family = gaussian(),  # Change family to Gaussian
data = dataWithSVIValues_std_logs_lags_prueba
)

#model <- glmmTMB(Hospitalized_log ~ predictors, family = Gamma(link = "log"), data = ...)
#vif(log_model_rolling)

step.model_rolling <- stepAIC(log_model_rolling, direction = "both",trace = FALSE)
summary(step.model_rolling)
AIC(step.model_rolling)

#This is to see the predicted and observed values
# Generate predictions with standard errors separately
predictions <- predict(step.model_rolling, newdata = dataWithSVIValues_std_logs_lags_prueba, type = "response", se.fit = TRUE)

# # Add predictions and confidence intervals to the data
# plot_data_fixed_effects <- dataWithSVIValues_std_logs_lags_prueba %>%
#   mutate(
#     fitted = predictions$fit,
#     se.fit = predictions$se.fit,
#     lower_ci = fitted - 1.96 * se.fit,  # 95% confidence interval lower bound
#     upper_ci = fitted + 1.96 * se.fit,  # 95% confidence interval upper bound
#     hosps_transformed1 = Hospitalized_log + 1,
#     backTrans_observed = exp(Hospitalized_log) - 1,
#     backTrans_fitted = exp(fitted) - 1,
#     backTrans_lower_ci = exp(lower_ci) - 1,
#     backTrans_upper_ci = exp(upper_ci) - 1
#   ) %>%
#   dplyr::select(week, SVIGroup, Hospitalized, fitted, lower_ci, upper_ci,backTrans_observed,
#                 backTrans_fitted,backTrans_lower_ci,backTrans_upper_ci) %>%
#   mutate(SVIGroup1=SVIGroup %>% str_replace_all(c("1"="SVI group 1","2"="SVI group 2","3"="SVI group 3","4"="SVI group 4")))
# 
# blindEstos<-c("#E69F00","#56B4E9","#009E73","#0072B2","#000000","#D55E00","#CC79A7")
# 
# ggplot(plot_data_fixed_effects, aes(x = week)) + theme_bw() +
#   geom_point(aes(y = Hospitalized, color = "Observed"),linetype=2,linewidth=1.5) +
#   geom_ribbon(aes(ymin=backTrans_lower_ci,ymax=backTrans_upper_ci),alpha=0.2,fill="red")+
#   geom_point(aes(y = backTrans_fitted, color = "Fitted"),linewidth=1.5) +
#   facet_wrap(~ SVIGroup1, scales = "free_y") +
#   labs(title = "Observed vs. Fitted Values by SVI group",y = "Hospitalizations",color = "Legend") +
#   theme(legend.position = "bottom",text=element_text(size=25))+
#   scale_color_manual(values = c("Observed"="#009E73","Fitted"="#D55E00"))
# 
# ggsave(last_plot(),file="~/Documents/GitHub/MobilityAndSVIModels/predictionGaussianModelBig.png")
# 
# Fitted values and residuals for random effects model
fitted_values <- fitted(step.model_rolling_log)
residuals <- resid(step.model_rolling_log)
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

## Prueba with GAM
library(mgcv)
library(gam)
packageVersion("mgcv")

gam_model <- gam(
  Hospitalized_log ~ 
    week + 
    random(SVIGroup) +  # Alternative random effect specification
    s(mobRatioOut_log_rollMean) + 
    s(mobRatioOut_log_rollMean_lag1) + 
    s(mobRatioOut_log_rollMean_lag2) + 
    s(mobRatioOut_log_rollMean_lag3) + 
    s(mobRatioOut_log_lag1) + 
    s(mobRatioOut_log_lag2) + 
    s(mobRatioOut_log_lag3) + 
    s(mobRatioIntra_log_lag1) + 
    s(mobRatioIntra_log_lag2) + 
    s(mobRatioIntra_log_lag3) + 
    offset(log(pop2019)),
  family = gaussian(),
  data = dataWithSVIValues_std_logs_lags_prueba,
  select=TRUE
)

summary(gam_model)
# Fit an initial GAM model
initial_gam <- gam(
  Hospitalized_log ~ 
    week + 
#    s(SVIGroup, bs = "re") + 
    s(mobRatioOut_log_rollMean) + 
    s(mobRatioOut_log_rollMean_lag1) + 
    s(mobRatioOut_log_lag1, k = 4) + 
    offset(log(pop2019)),
  family = gaussian(),
  data = dataWithSVIValues_std_logs_lags_prueba
)

# Perform stepwise selection
final_gam <- step.gam(initial_gam, scope = list(
  lower = ~ week + s(SVIGroup, bs = "re") + offset(log(pop2019)),
  upper = ~ week + s(SVIGroup, bs = "re") +
    s(mobRatioOut_log_rollMean) + 
    s(mobRatioOut_log_rollMean_lag1) + 
    s(mobRatioOut_log_rollMean_lag2) + 
    s(mobRatioOut_log_rollMean_lag3) + 
    s(mobRatioOut_log_lag1, k = 4) + 
    s(mobRatioOut_log_lag2, k = 4) + 
    s(mobRatioOut_log_lag3, k = 4) + 
    s(mobRatioIntra_log_lag1, k = 4) + 
    s(mobRatioIntra_log_lag2, k = 4) + 
    s(mobRatioIntra_log_lag3, k = 4)
))


gam_model <- gam(Hospitalized_log ~ week + s(mobRatioOut_log_rollMean) + 
                   s(mobRatioOut_log_rollMean_lag1) + 
                   s(mobRatioOut_log_rollMean_lag2) + 
                   s(mobRatioOut_log_rollMean_lag3) + 
                   s(mobRatioOut_log_lag1) + 
                   s(mobRatioOut_log_lag2) + 
                   offset(log(pop2019)), 
                 family = gaussian(), data = dataWithSVIValues_std_logs_lags_prueba)

summary(gam_model)




