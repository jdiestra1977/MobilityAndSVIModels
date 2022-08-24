library(tidyverse)

#Personal access token: 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 
#system("git remote set-url origin https://github.com/jdiestra1977/MobilityAndSVIModels.git")
#system("git branch -M main")
#system("git push -u origin main")

load("totalVisitsByZipVsTime.RData")
totalVisitsByZipVsTime %>%
  ggplot(aes(x=as.Date(Date),y=InterVisitsOut+InterVisitsIn,group=Zip))+
  theme_bw()+
  geom_line(size=1,alpha=07,color="gray") +
  ylab("Total visits between zip codes") + xlab("")+
  theme(text=element_text(size=20))

#Correlation and best fit for Mobility ratio
load("statsMobilityRatioAndSVI.RData")
statsMobilityRatioAndSVI %>% 
  ggplot(aes(x=as.Date(Days),y=Average))+
  geom_line() + facet_wrap(~Measure)

load("statsDistanceRatioAndSVI.RData")
statsDistanceRatioAndSVI %>%
  ggplot(aes(x=as.Date(Days),y=Average,group=Measure))+
  geom_line()+
  facet_wrap(~Measure,scales = "free")

