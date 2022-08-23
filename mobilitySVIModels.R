library(tidyverse)

#Personal access token: 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 
#system("git remote set-url origin https://github.com/jdiestra1977/MobilityAndSVIModels.git")
#system("git branch -M main")
#system("git push -u origin main")
#Correlation and best fit for Mobility ratio
load("statsMobilityRatioAndSVI.RData")
statsMobilityRatioAndSVI %>% head()

statsMobilityRatioAndSVI %>% 
  ggplot(aes(x=as.Date(Days),y=Average))+
  geom_line() + facet_wrap(~Measure)

load("totalVisitsByZipVsTime.RData")
totalVisitsByZipVsTime %>%
  ggplot(aes(x=as.Date(Date),y=TotalVisits,group=Zip))+
  geom_line(size=1,alpha=07,color="gray")
