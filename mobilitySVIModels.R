library(tidyverse)

#Correlation and best fit for Mobility ratio
load("statsMobilityRatioAndSVI.RData")
statsMobilityRatioAndSVI %>% head()

statsMobilityRatioAndSVI %>% 
  ggplot(aes(x=as.Date(Days),y=Average))+
  geom_line() + facet_wrap(~Measure)
