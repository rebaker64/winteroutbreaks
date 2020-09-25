#code to run epiestim was provided by Alex Becker and adapted here

## Based on
## https://cmmid.github.io/topics/covid19/current-patterns-transmission/global-time-varying-transmission.html
## further adapted from
## https://github.com/aperaltasantos/covid_pt && https://cran.r-project.org/web/packages/EpiEstim/vignettes/demo.html
## Methods as described above
## Time-varying effective reproduction estimates were made with a 7-day sliding window using EpiEstim
## assuming a serial interval of 4.7 days 
## and a standard deviation of 2.9 days 


require(EpiEstim)
require(dplyr)
require(ggplot2)
require(RCurl)
require(reshape2)
require(purrr)
require(lubridate)

theme_set(theme_classic(base_size = 16))

data <- read.csv(text=getURL("https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"),skip=1) #read from github

names(data) <- c('time','county','state','fips','cases','deaths')

data %>%
  group_by(state,county) %>%
  summarize() -> test


data %>%
  #subset(state == 'Missouri') %>%
  mutate(county = as.character(county))  -> data

data$time <- as.Date(as.character(data$time),format='%Y-%m-%d')

data %>%
  mutate(Date = time) %>%
  select(time,county,state,cases,deaths,Date) %>%
  subset(state == 'New York') %>%
  subset(county == 'New York City') -> dat

dat %>%
  select(time,cases,deaths) %>%
  melt(id.vars = 'time') %>%
  ggplot(aes(time,value,color=variable))+geom_line(size=2)+
  facet_wrap(~variable,nrow=2,scales='free')


## which data source do you want to use?
## confirmed = deaths or cases
dat %>%
  mutate(Confirmed = cases) %>%
  select(Date,Confirmed) -> covid_pt


covid_pt<-covid_pt  %>%
  #subset(Date >= '2020-03-05') %>%
  mutate(epiweek = epiweek(Date))

first.date <- head(covid_pt$Date,1)

covid_pt %>%
  mutate(
    Confirmed_lag = lag(x = Confirmed,
                        n = 1,
                        order_by = Date),
    Confirmed_var=Confirmed-Confirmed_lag,
    Confirmed_sign=if_else(Confirmed_var>=0,"+","-")
  ) %>%
  subset(Date >  first.date) -> covid_pt

covid_pt  %>%
  select(
    Date,Confirmed_var
  )  %>%
  dplyr::mutate(
    t_start = dplyr::row_number() %>% as.numeric(),
    t_end = t_start + 6
  ) -> covid_r

## set negative daily counts to zero
covid_r$Confirmed_var[ covid_r$Confirmed_var < 0 ] <- 0

## parametric estimate given 'known' SI (no CIs)
res_parametric_si <-
  estimate_R(
    covid_r$Confirmed_var,
    method ="parametric_si",
    config = make_config(
      list(
        mean_si = 4.7,
        std_si = 2.9
      )
    )
  )

plot(res_parametric_si, legend = FALSE)

r_prt <- as.data.frame(res_parametric_si$R)
NYestimR <- r_prt

save(NYestimR , file="NYestimR.RData")

