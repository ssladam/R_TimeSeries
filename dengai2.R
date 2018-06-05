rm(list=ls())
pkgs <- c('tidyverse','tidyimpute', 'tibbletime', 'corrplot', 'magrittr', 'zoo', 'RColorBrewer', 'gridExtra','MASS','mice', 'forecast')
invisible(lapply(pkgs, require, character.only = T))
rm(pkgs)
select <- dplyr::select

setwd("/temp/NU/413/dengai")
#X <- read.csv("dengue_features_train.csv") #reads in as a list...
train.X <- read_csv("dengue_features_train.csv")
train.y <- read_csv("dengue_labels_train.csv")

test.X <- read_csv('dengue_features_test.csv')

#train.X %<>% left_join(train.y, by=c('city', 'year', 'weekofyear')) #brings the response into same frame

#=====================================================================================
#             Deaths on 53-week calendar, rest of data is on 52. Let's combine
#=====================================================================================
train.X %<>% filter(weekofyear!=53) #It was all blank data previously, safe to just drop

#create a holder variable to mark where we need to forecast deaths in the test set...
test.X.week53 <- test.X %>% filter(weekofyear==53)
test.X %<>% filter(weekofyear!=53)

joinyear <- train.y %>% filter(weekofyear==53) %>% select(year) %>% as.list()
joinyear <- joinyear$year
incity <- train.y %>% filter(weekofyear==53) %>% select(city) %>% as.list()
incity <- incity$city

mergecases <- train.y %>% filter(year %in% joinyear & weekofyear==1) %>% select(total_cases) +
  train.y %>% filter(weekofyear==53) %>% select(total_cases) %>% as.list()
mergecases <- mergecases$total_cases

for (i in (1:length(joinyear))) {
  #curr_deaths <- train.y[(train.y$weekofyear==1 & train.y$year==(joinyear$year[i] + 1) & train.y$city==incity[i]),'total_cases']
  train.y[(train.y$weekofyear==1 & train.y$year==joinyear[i] & train.y$city==incity[i]),'total_cases'] <-
    as.integer(mergecases[i])
}
#now get rid of the week-53 in the training set...
train.y %<>% filter(weekofyear!=53)

rm(i, incity, joinyear, mergecases)


#=====================================================================================
#             Missing data? Perform imputation
#=====================================================================================
full.X <- rbind(train.X, test.X)

is_missing <- function(input.data){
  input.data %>% group_by(city) %>%
    summarise_all(funs(round(100*mean(is.na(.)),2))) %>% t %>% as.data.frame()
}

# count missing values (as percent)
full.X %>% is_missing()
full.X.imp <- mice(data=full.X, method = 'rf', ntree=3)
full.X.comp <- as.tibble(complete(full.X.imp, 1))
full.X.comp %>% is_missing()
#reanalysis_sat_precip_amt_mm could not be fully imputed, let's not use that predictor
full.X.comp[is.na(full.X.comp$reanalysis_sat_precip_amt_mm),] %>% glimpse()
rm(full.X.imp)

# plot correlation matrix
full.X.comp %>% 
  dplyr::select(-city, -year, -weekofyear, -week_start_date) %>%
  cor(use = 'pairwise.complete.obs') -> M1
corrplot(M1, type="lower", method="color",
         col=brewer.pal(n=8, name="RdBu"),diag=FALSE)
#reanalysis_sat_precip_amt_mm is perfectly correlated with precipitation_amt_mm, just use that

#=====================================================================================
#             Create features
#=====================================================================================
#BUG! you must split by city first, or else you roll over between cities
num <- 3
full.X.comp$precip1_ms3 <- rollapplyr(full.X.comp$precipitation_amt_mm,list(-(num:1)),sum,fill=NA)
full.X.comp$precip2_ms3 <- rollapplyr(full.X.comp$reanalysis_precip_amt_kg_per_m2,list(-(num:1)),sum,fill=NA)
full.X.comp$avgTemp1_ms3 <- rollapplyr(full.X.comp$station_avg_temp_c,list(-(num:1)),sum,fill=NA)
full.X.comp$avgTemp2_ms3 <- rollapplyr(full.X.comp$reanalysis_avg_temp_k,list(-(num:1)),sum,fill=NA)
full.X.comp$maxTemp1_ms3 <- rollapplyr(full.X.comp$station_max_temp_c,list(-(num:1)),sum,fill=NA)
full.X.comp$maxTemp2_ms3 <- rollapplyr(full.X.comp$reanalysis_max_air_temp_k,list(-(num:1)),sum,fill=NA)
full.X.comp$humid_ms3 <- rollapplyr(full.X.comp$reanalysis_specific_humidity_g_per_kg,list(-(num:1)),sum,fill=NA)
full.X.comp$ndvi_all <- full.X.comp$ndvi_ne + full.X.comp$ndvi_nw + full.X.comp$ndvi_se + full.X.comp$ndvi_sw
full.X.comp$ndvi_ms3 <- rollapplyr(full.X.comp$ndvi_all,list(-(num:1)),sum,fill=NA)

mydel <- full.X.comp %>% mutate(season= case_when(
  weekofyear > (52-4) ~ as.integer(0),
  weekofyear > (52-4-12) ~ as.integer(3),
  weekofyear > (52-4-12-12) ~ as.integer(2),
  weekofyear > (52-4-12-12-12) ~ as.integer(1),
  TRUE ~ as.integer(0)
))

mydel %<>% mutate(year = as.integer(format(as.Date(week_start_date),'%Y')))

mydel %<>% group_by(city, year, season) %>% mutate(seas_avg_temp = mean(reanalysis_avg_temp_k))
mydel %<>% group_by(city, year, season) %>% mutate(seas_min_temp = mean(reanalysis_min_air_temp_k))


boxplot(full.X.comp$reanalysis_air_temp_k)

mydel %>% filter(city=='sj') %>% select(seas_avg_temp, seas_min_temp)

full.sj <- full.X.comp %>% filter(city=='sj')
full.iq <- full.X.comp %>% filter(city=='iq')

#Now let's create a new batch of features... this time focus on longer periods!

full.sj$percip.lagw1 <- rollapplyr(full.sj$precipitation_amt_mm,list(-(7:1)),mean,fill=NA)
full.sj$percip.lagw2 <- rollapplyr(full.sj$precipitation_amt_mm,list(-(14:8)),mean,fill=NA)
full.sj$percip.lagw3 <- rollapplyr(full.sj$precipitation_amt_mm,list(-(21:15)),mean,fill=NA)
full.iq$percip.lagw1 <- rollapplyr(full.iq$precipitation_amt_mm,list(-(7:1)),mean,fill=NA)
full.iq$percip.lagw2 <- rollapplyr(full.iq$precipitation_amt_mm,list(-(14:8)),mean,fill=NA)
full.iq$percip.lagw3 <- rollapplyr(full.iq$precipitation_amt_mm,list(-(21:15)),mean,fill=NA)

full.sj$temp.lagw1 <- rollapplyr(full.sj$reanalysis_avg_temp_k,list(-(7:1)),mean,fill=NA)
full.sj$temp.lagw2 <- rollapplyr(full.sj$reanalysis_avg_temp_k,list(-(14:8)),mean,fill=NA)
full.sj$temp.lagw3 <- rollapplyr(full.sj$reanalysis_avg_temp_k,list(-(21:15)),mean,fill=NA)
full.iq$temp.lagw1 <- rollapplyr(full.iq$reanalysis_avg_temp_k,list(-(7:1)),mean,fill=NA)
full.iq$temp.lagw2 <- rollapplyr(full.iq$reanalysis_avg_temp_k,list(-(14:8)),mean,fill=NA)
full.iq$temp.lagw3 <- rollapplyr(full.iq$reanalysis_avg_temp_k,list(-(21:15)),mean,fill=NA)

full.sj$humid.lagw1 <- rollapplyr(full.sj$reanalysis_specific_humidity_g_per_kg,list(-(7:1)),mean,fill=NA)
full.sj$humid.lagw2 <- rollapplyr(full.sj$reanalysis_specific_humidity_g_per_kg,list(-(14:8)),mean,fill=NA)
full.sj$humid.lagw3 <- rollapplyr(full.sj$reanalysis_specific_humidity_g_per_kg,list(-(21:15)),mean,fill=NA)
full.iq$humid.lagw1 <- rollapplyr(full.iq$reanalysis_specific_humidity_g_per_kg,list(-(7:1)),mean,fill=NA)
full.iq$humid.lagw2 <- rollapplyr(full.iq$reanalysis_specific_humidity_g_per_kg,list(-(14:8)),mean,fill=NA)
full.iq$humid.lagw3 <- rollapplyr(full.iq$reanalysis_specific_humidity_g_per_kg,list(-(21:15)),mean,fill=NA)

full.sj$ndvi.lagw1 <- rollapplyr(full.sj$ndvi_all,list(-(7:1)),mean,fill=NA)
full.sj$ndvi.lagw2 <- rollapplyr(full.sj$ndvi_all,list(-(14:8)),mean,fill=NA)
full.sj$ndvi.lagw3 <- rollapplyr(full.sj$ndvi_all,list(-(21:15)),mean,fill=NA)
full.iq$ndvi.lagw1 <- rollapplyr(full.iq$ndvi_all,list(-(7:1)),mean,fill=NA)
full.iq$ndvi.lagw2 <- rollapplyr(full.iq$ndvi_all,list(-(14:8)),mean,fill=NA)
full.iq$ndvi.lagw3 <- rollapplyr(full.iq$ndvi_all,list(-(21:15)),mean,fill=NA)


nrow.sj.train <- train.X %>% filter(city=='sj') %>% nrow()
nrow.iq.train <- train.X %>% filter(city=='iq') %>% nrow()
nrow.sj.test <- test.X %>% filter(city=='sj') %>% nrow()
nrow.iq.test <- test.X %>% filter(city=='iq') %>% nrow()
sj.train <- head(full.sj, n=nrow.sj.train)
sj.test <- tail(full.sj, n=nrow.sj.test)
iq.train <- head(full.iq, n=nrow.iq.train)
iq.test <- tail(full.iq, n=nrow.iq.test)

full.train <- rbind(sj.train, iq.train)
full.test <- rbind(sj.test, iq.test)
full.X <- rbind(full.train, full.test)

#impute missing values...
full.X %>% is_missing()
full.X.imp <- mice(data=full.X, method = 'rf', ntree=3)
full.X.comp <- as.tibble(complete(full.X.imp, 1))
full.X.comp %>% is_missing()
#reanalysis_sat_precip_amt_mm could not be fully imputed, let's not use that predictor
full.X.comp[is.na(full.X.comp$reanalysis_sat_precip_amt_mm),] %>% glimpse()
rm(full.X.imp)



#Now break them back into train / test set, now that imputation has finished...
train.X <- head(full.X.comp, n=nrow(train.X))
#train.X %>% is_missing()
test.X <- tail(full.X.comp, n=nrow(test.X))
#test.X %>% is_missing()







#=====================================================================================
#             Split data set into two cities
#=====================================================================================
train.X.sj <- train.X %>% filter(city=='sj')
train.y.sj <- train.y %>% filter(city=='sj')
train.X.iq <- train.X %>% filter(city=='iq')
train.y.iq <- train.y %>% filter(city=='iq')

test.X.sj <- test.X %>% filter(city=='sj')
test.X.iq <- test.X %>% filter(city=='iq')


train.X.sj %<>% mutate('total_cases' = train.y.sj$total_cases)
train.X.sj$case_lag1 <- lag(train.X.sj$total_cases, 1)
train.X.sj$case_lag2 <- lag(train.X.sj$total_cases, 2)
train.X.sj$case_lag3 <- lag(train.X.sj$total_cases, 3)
train.X.iq %<>% mutate('total_cases' = train.y.iq$total_cases)
train.X.iq$case_lag1 <- lag(train.X.iq$total_cases, 1)
train.X.iq$case_lag2 <- lag(train.X.iq$total_cases, 2)
train.X.iq$case_lag3 <- lag(train.X.iq$total_cases, 3)


test.X.sj$predicted <- 0
test.X.sj$case_lag1 <- 0
test.X.sj$case_lag2 <- 0
test.X.sj$case_lag3 <- 0
test.X.sj$case_lag1[1] <- train.X.sj$total_cases[nrow(train.X.sj)]
test.X.sj$case_lag2[1] <- train.X.sj$total_cases[(nrow(train.X.sj)-1)]
test.X.sj$case_lag2[2] <- train.X.sj$total_cases[(nrow(train.X.sj))]
test.X.sj$case_lag3[1] <- train.X.sj$total_cases[(nrow(train.X.sj)-2)]
test.X.sj$case_lag3[2] <- train.X.sj$total_cases[(nrow(train.X.sj)-1)]
test.X.sj$case_lag3[3] <- train.X.sj$total_cases[(nrow(train.X.sj))]

test.X.iq$predicted <- 0
test.X.iq$case_lag1 <- 0
test.X.iq$case_lag2 <- 0
test.X.iq$case_lag3 <- 0
test.X.iq$case_lag1[1] <- train.X.iq$total_cases[nrow(train.X.iq)]
test.X.iq$case_lag2[1] <- train.X.iq$total_cases[(nrow(train.X.iq)-1)]
test.X.iq$case_lag2[2] <- train.X.iq$total_cases[(nrow(train.X.iq))]
test.X.iq$case_lag3[1] <- train.X.iq$total_cases[(nrow(train.X.iq)-2)]
test.X.iq$case_lag3[2] <- train.X.iq$total_cases[(nrow(train.X.iq)-1)]
test.X.iq$case_lag3[3] <- train.X.iq$total_cases[(nrow(train.X.iq))]

#=====================================================================================
#             final imputation for missing values
#=====================================================================================
train.X.sj.imp <- mice(data=train.X.sj, method = 'rf', ntree=3)
train.X.sj <- as.tibble(complete(train.X.sj.imp, 1))
train.X.sj %>% is_missing()
rm(train.X.sj.imp)
train.X.iq.imp <- mice(data=train.X.iq, method = 'rf', ntree=3)
train.X.iq <- as.tibble(complete(train.X.iq.imp, 1))
train.X.iq %>% is_missing()
rm(train.X.iq.imp)





cust.xreg <- c('station_max_temp_c', 'station_avg_temp_c', 'precip1_ms3',
               'reanalysis_specific_humidity_g_per_kg', 'case_lag1',
               'ndvi_se', 'ndvi_sw','ndvi_ne','ndvi_nw')
cust.xreg1 <- c('precip1_ms3', 'avgTemp1_ms3', 'maxTemp1_ms3', 'humid_ms3', 'ndvi_all', 'ndvi_ms3', 'case_lag1')
cust.xreg2 <- c('precip2_ms3', 'avgTemp2_ms3', 'maxTemp2_ms3', 'humid_ms3', 'ndvi_all', 'ndvi_ms3', 'case_lag1')
cust.xreg3 <- c('precip1_ms3', 'avgTemp2_ms3', 'ndvi_se', 'ndvi_sw','ndvi_ne','ndvi_nw', 'ndvi_ms3', 'case_lag1')

cust.xreg4 <- c('case_lag1', 'percip.lagw1', 'percip.lagw2', 'percip.lagw3', 'temp.lagw1', 
                'temp.lagw2', 'temp.lagw3', 'humid.lagw1', 'humid.lagw2', 'humid.lagw3')
#,'ndvi.lagw1', 'ndvi.lagw2', 'ndvi.lagw3')

cust.xreg.lag1 <- c('case_lag1', 'percip.lagw1', 'temp.lagw1', 'humid.lagw1', 'ndvi.lagw1')
cust.xreg.lag2 <- c('case_lag1', 'percip.lagw2', 'temp.lagw2', 'humid.lagw2', 'ndvi.lagw2')
cust.xreg.lag3 <- c('case_lag1', 'percip.lagw3', 'temp.lagw3', 'humid.lagw3', 'ndvi.lagw3')

cust.xreg.caselag <- c('case_lag1', 'case_lag2','case_lag3','percip.lagw2', 'temp.lagw2', 'humid.lagw2', 'ndvi.lagw2')


#=====================================================================================
#             Create the model
#=====================================================================================

get_best_model <- function(input.data, input.xreg) {
  best_mae <- 1000
  true.mae = vector()
  split <- floor(nrow(input.data)*.2)
  X.train <- ts(input.data$total_cases[1:(nrow(input.data)-split)], frequency = 52)
  X.test <- ts(input.data$total_cases[(nrow(input.data)-split):nrow(input.data)], frequency = 52, start=(end(X.train) + c(0,1)))
  xreg.data = input.data %>% select(input.xreg)
  xreg.train <- xreg.data[1:(nrow(input.data)-split),]
  xreg.test <- xreg.data[(nrow(input.data)-split):nrow(input.data),]
  grid = seq(1:3)
  
  
  for (i in grid) {
    X.fit <- nnetar(X.train, p=1, repeats = 10, maxit=256, scale.inputs = TRUE, xreg=xreg.train)
    X.fc <- forecast(X.fit, h=split, xreg=xreg.test)
    model.mae <- accuracy(X.fc)[3]
    
    #for (j in seq((nrow(input.data)-split),(nrow(input.data)-1) )) {
    for (j in 1:nrow(xreg.test) ) {
      #print(j)
      point.fc <- forecast(X.fit, h=1, xreg=xreg.test[j,])
      point.fc <- as.numeric(point.fc$mean)
      if(point.fc < 0) {point.fc <- 0}
      point.mae <- abs(point.fc - X.test[j])
      true.mae <- append(true.mae, point.mae)
      #paste(j, ":", true.mae)
      # test.X.iq$predicted[i] <- point.fc
      # xreg.iq$case_lag1[i+1] <- point.fc
    }
    
    #print(mean(true.mae))
    
    mae <- mean(true.mae)
    if (mae < best_mae) {
      best_fit <- X.fit
      best_mae <- mae
    }
    cat('Iteration:', i, 'model.mae:', model.mae, ' true mae:', mae, '\n')
  }
  X.fit <- nnetar(X.train, xreg=xreg.train, model = best_fit)
  X.fc <- forecast(X.fit, h=split, xreg=xreg.test)
  plot(X.fc, type='l')
  lines(X.test,col='red')
  X.fit <- nnetar(ts(input.data$total_cases, frequency=52), xreg=xreg.data, model = best_fit)
  return(X.fit)
}

best_sj <- get_best_model(train.X.sj, cust.xreg) #5.8, 5.9
best_sj <- get_best_model(train.X.sj, cust.xreg1) #6.3, 6.3
best_sj <- get_best_model(train.X.sj, c(cust.xreg, cust.xreg1)) #4.8, 4.9
best_sj <- get_best_model(train.X.sj, cust.xreg2) #6.4, 6.3
best_sj <- get_best_model(train.X.sj, c(cust.xreg, cust.xreg2)) #4.4, 4.6
best_sj <- get_best_model(train.X.sj, c(cust.xreg, cust.xreg2, cust.xreg1)) #3.9,4.2 
best_sj <- get_best_model(train.X.sj, c(cust.xreg, cust.xreg3))

#remove first year and a half since apparently disease not common before then
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , cust.xreg) #2.2
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , cust.xreg1) #2.6, but big spikes!
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , c(cust.xreg, cust.xreg1)) #1.3
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , cust.xreg2) #2.6
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , c(cust.xreg, cust.xreg2)) #1.0
best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , c(cust.xreg, cust.xreg2, cust.xreg1)) #0.66





#=====================================================================================
#             Now that we have best models, let's fit it to the test data
#=====================================================================================

test.xreg <- c(cust.xreg2, cust.xreg4) #11.3
test.xreg <- c(cust.xreg.lag1) #9.2
test.xreg <- c(cust.xreg.lag2) #7.9
test.xreg <- c(cust.xreg.lag3) #10.5
test.xreg <- c(cust.xreg.caselag, cust.xreg)

xreg.sj = test.X.sj %>% select(test.xreg)
best_sj <- get_best_model(train.X.sj, test.xreg)
for (i in (1:(nrow(test.X.sj)))) {
  point.fc <- forecast(best_sj, h=1, xreg=xreg.sj[i,])
  point.fc <- as.numeric(point.fc$mean)
  if(point.fc < 0) {point.fc <- 0}
  test.X.sj$predicted[i] <- point.fc
  if (i<(nrow(test.X.sj)-1)) {
    xreg.sj$case_lag1[i+1] <- point.fc
    if (i<(nrow(test.X.sj)-2) & 'case_lag2' %in% colnames(xreg.sj)) {
      xreg.sj$case_lag2[i+2] <- point.fc
      if (i<(nrow(test.X.sj)-3) & 'case_lag3' %in% colnames(xreg.sj)) {
        xreg.sj$case_lag3[i+3] <- point.fc
      }
    }
  }
}
# we didn't forecast the last points, to avoid writing a lag1 that goes beyond the index. Let's do it now
# point.fc <- forecast(best_sj, h=1, xreg=xreg.sj[(i+1),])
# point.fc <- as.numeric(point.fc$mean)
# if(point.fc < 0) {point.fc <- 0}
# test.X.sj$predicted[(i+1)] <- point.fc

plot(test.X.sj$predicted, type='l', main='SJ predicted cases')


best_iq <- get_best_model(tail(train.X.iq,n=(520-52-26)) , test.xreg) #1.0
xreg.iq = test.X.iq %>% select(test.xreg)
for (i in (1:(nrow(test.X.iq)))) {
  point.fc <- forecast(best_sj, h=1, xreg=xreg.iq[i,])
  point.fc <- as.numeric(point.fc$mean)
  if(point.fc < 0) {point.fc <- 0}
  test.X.iq$predicted[i] <- point.fc
  if (i<(nrow(test.X.iq)-1)) {
    xreg.iq$case_lag1[i+1] <- point.fc
    if (i<(nrow(test.X.iq)-2) & 'case_lag2' %in% colnames(xreg.iq)) {
      xreg.iq$case_lag2[i+2] <- point.fc
      if (i<(nrow(test.X.iq)-3) & 'case_lag3' %in% colnames(xreg.iq)) {
        xreg.iq$case_lag3[i+3] <- point.fc
      }
    }
  }
}
# we didn't forecast the last point, to avoid writing a lag1 that goes beyond the index. Let's do it now
# point.fc <- forecast(best_iq, h=1, xreg=xreg.iq[(i+1),])
# point.fc <- as.numeric(point.fc$mean)
# if(point.fc < 0) {point.fc <- 0}
# test.X.iq$predicted[(i+1)] <- point.fc

plot(test.X.iq$predicted, type='l', main='IQ predicted cases')







#====================================================================
#  Gather all the data for submission
#====================================================================
submission = read.csv('submission_format.csv')
inner_join(submission, rbind(test.X.sj,test.X.iq), by=c('city', 'year', 'weekofyear')) %>%
  dplyr::select(city, year, weekofyear, total_cases=predicted) ->
  prediction.df

prediction.df$total_cases <- round(prediction.df$total_cases)

write.csv(prediction.df, file='AdamsStewartSubmission.csv', row.names = FALSE)

