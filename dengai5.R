rm(list=ls())
pkgs <- c('tidyverse','tidyimpute', 'tibbletime', 'data.table', 'corrplot', 'magrittr', 'zoo', 'RColorBrewer', 'gridExtra','MASS','mice', 'forecast')
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
#             Cases on 53-week calendar, rest of data is on 52. Let's combine
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
full.X.imp <- mice(data=full.X, method = 'rf', ntree=3, m=1)
full.X.comp <- as.tibble(complete(full.X.imp, 1))
full.X.comp %>% is_missing()
#reanalysis_sat_precip_amt_mm could not be fully imputed, let's not use that predictor
full.X.comp[is.na(full.X.comp$reanalysis_sat_precip_amt_mm),] %>% glimpse()
rm(full.X.imp)

full.X.comp %<>% select(-reanalysis_sat_precip_amt_mm)

# my.na.locf <- function(x) {
#   v <- !is.na(x)
#   c(NA, x[v])[cumsum(v)+1]
# }
# 
# full.X[]=lapply(full.X, my.na.locf)
# 
# full.X.comp <- full.X

#=====================================================================================
#             Create features
#=====================================================================================
add_group_summaries <- function(d, 
                                groupingVars, 
                                ...) {
  groupingSyms <- rlang::syms(groupingVars)
  d <- ungroup(d) # just in case
  dg <- group_by(d, !!!groupingSyms)
  ds <- summarize(dg, ...)
  ds <- ungroup(ds)
  left_join(d, ds, by= groupingVars)
}

#seasonal variables
full.X.comp %<>% mutate(season= case_when(
  weekofyear > (52-4) ~ as.integer(0),
  weekofyear > (52-4-12) ~ as.integer(3),
  weekofyear > (52-4-12-12) ~ as.integer(2),
  weekofyear > (52-4-12-12-12) ~ as.integer(1),
  TRUE ~ as.integer(0)
))

full.X.comp %<>% add_group_summaries(c('city', 'year', 'season'),
                                     seas_min_temp = min(reanalysis_min_air_temp_k),
                                     seas_max_temp = mean(reanalysis_max_air_temp_k),
                                     seas_avg_temp = max(reanalysis_avg_temp_k))


full.sj <- full.X.comp %>% filter(city=='sj')
full.iq <- full.X.comp %>% filter(city=='iq')

zz<-full.sj %>% select(-c(city, year,weekofyear, week_start_date, season))
zz.colnames <- zz%>%colnames()
setDT(zz)[,paste0('ndvi_ne', 1:4) := shift(ndvi_ne, 1:4)][]
setDT(zz)[,paste0('ndvi_nw', 1:4) := shift(ndvi_nw, 1:4)][]
setDT(zz)[,paste0('ndvi_se', 1:4) := shift(ndvi_se, 1:4)][]
setDT(zz)[,paste0('ndvi_sw', 1:4) := shift(ndvi_sw, 1:4)][]
setDT(zz)[,paste0('precipitation_amt_mm', 1:4) := shift(precipitation_amt_mm, 1:4)][]
setDT(zz)[,paste0('reanalysis_air_temp_k', 1:4) := shift(reanalysis_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_avg_temp_k', 1:4) := shift(reanalysis_avg_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_dew_point_temp_k', 1:4) := shift(reanalysis_dew_point_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_max_air_temp_k', 1:4) := shift(reanalysis_max_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_min_air_temp_k', 1:4) := shift(reanalysis_min_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_precip_amt_kg_per_m2', 1:4) := shift(reanalysis_precip_amt_kg_per_m2, 1:4)][]
setDT(zz)[,paste0('reanalysis_relative_humidity_percent', 1:4) := shift(reanalysis_relative_humidity_percent, 1:4)][]
setDT(zz)[,paste0('reanalysis_tdtr_k', 1:4) := shift(reanalysis_tdtr_k, 1:4)][]
setDT(zz)[,paste0('station_avg_temp_c', 1:4) := shift(station_avg_temp_c, 1:4)][]
setDT(zz)[,paste0('station_max_temp_c', 1:4) := shift(station_max_temp_c, 1:4)][]
setDT(zz)[,paste0('station_min_temp_c', 1:4) := shift(station_min_temp_c, 1:4)][]
setDT(zz)[,paste0('station_precip_mm', 1:4) := shift(station_precip_mm, 1:4)][]
setDT(zz)[,paste0('seas_min_temp', 1:4) := shift(seas_min_temp, 1:4)][]
setDT(zz)[,paste0('seas_max_temp', 1:4) := shift(seas_max_temp, 1:4)][]
setDT(zz)[,paste0('seas_avg_temp', 1:4) := shift(seas_avg_temp, 1:4)][]
full.sj <- inner_join(full.sj, zz)

zz<-full.iq %>% select(-c(city, year,weekofyear, week_start_date, season))
zz.colnames <- zz%>%colnames()
setDT(zz)[,paste0('ndvi_ne', 1:4) := shift(ndvi_ne, 1:4)][]
setDT(zz)[,paste0('ndvi_nw', 1:4) := shift(ndvi_nw, 1:4)][]
setDT(zz)[,paste0('ndvi_se', 1:4) := shift(ndvi_se, 1:4)][]
setDT(zz)[,paste0('ndvi_sw', 1:4) := shift(ndvi_sw, 1:4)][]
setDT(zz)[,paste0('precipitation_amt_mm', 1:4) := shift(precipitation_amt_mm, 1:4)][]
setDT(zz)[,paste0('reanalysis_air_temp_k', 1:4) := shift(reanalysis_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_avg_temp_k', 1:4) := shift(reanalysis_avg_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_dew_point_temp_k', 1:4) := shift(reanalysis_dew_point_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_max_air_temp_k', 1:4) := shift(reanalysis_max_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_min_air_temp_k', 1:4) := shift(reanalysis_min_air_temp_k, 1:4)][]
setDT(zz)[,paste0('reanalysis_precip_amt_kg_per_m2', 1:4) := shift(reanalysis_precip_amt_kg_per_m2, 1:4)][]
setDT(zz)[,paste0('reanalysis_relative_humidity_percent', 1:4) := shift(reanalysis_relative_humidity_percent, 1:4)][]
setDT(zz)[,paste0('reanalysis_tdtr_k', 1:4) := shift(reanalysis_tdtr_k, 1:4)][]
setDT(zz)[,paste0('station_avg_temp_c', 1:4) := shift(station_avg_temp_c, 1:4)][]
setDT(zz)[,paste0('station_max_temp_c', 1:4) := shift(station_max_temp_c, 1:4)][]
setDT(zz)[,paste0('station_min_temp_c', 1:4) := shift(station_min_temp_c, 1:4)][]
setDT(zz)[,paste0('station_precip_mm', 1:4) := shift(station_precip_mm, 1:4)][]
setDT(zz)[,paste0('seas_min_temp', 1:4) := shift(seas_min_temp, 1:4)][]
setDT(zz)[,paste0('seas_max_temp', 1:4) := shift(seas_max_temp, 1:4)][]
setDT(zz)[,paste0('seas_avg_temp', 1:4) := shift(seas_avg_temp, 1:4)][]
full.iq <- inner_join(full.iq, zz)




nrow.sj.train <- train.X %>% filter(city=='sj') %>% nrow()
nrow.sj.test <- test.X %>% filter(city=='sj') %>% nrow()
nrow.iq.train <- train.X %>% filter(city=='iq') %>% nrow()
nrow.iq.test <- test.X %>% filter(city=='iq') %>% nrow()
sj.train <- head(full.sj, n=nrow.sj.train)
sj.test <- tail(full.sj, n=nrow.sj.test)
sj.cases <- train.y %>% filter(city=='sj') %>% select(total_cases)
iq.train <- head(full.iq, n=nrow.iq.train)
iq.test <- tail(full.iq, n=nrow.iq.test)
iq.cases <- train.y %>% filter(city=='iq') %>% select(total_cases)


xreg.full <- full.sj %>% select(-c(city,year,week_start_date)) %>% colnames()


sj.train <- head(full.sj, n=nrow.sj.train)
sj.test <- tail(full.sj, n=nrow.sj.test)
iq.train <- head(full.iq, n=nrow.iq.train)
iq.test <- tail(full.iq, n=nrow.iq.test)

full.train <- rbind(sj.train, iq.train)
full.test <- rbind(sj.test, iq.test)
full.X <- rbind(full.train, full.test)


#Now break them back into original structure train / test set, now that feature creation has finished
train.X <- head(full.X, n=nrow(train.X))
test.X <- tail(full.X, n=nrow(test.X))

#=====================================================================================
#             Split data set into two cities, setup case lags
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
#             Prepare for modelling
#=====================================================================================


get_best_model <- function(input.data, input.xreg, grid.len=3, weights=1000) {
  nMaxNWts <- weights
  best_mae <- 1000
  true.mae = vector()
  split <- floor(nrow(input.data)*.2)
  X.train <- ts(input.data$total_cases[1:(nrow(input.data)-split)], frequency = 52)
  X.test <- ts(input.data$total_cases[(nrow(input.data)-split):nrow(input.data)], frequency = 52, start=(end(X.train) + c(0,1)))
  xreg.data = input.data %>% select(input.xreg)
  xreg.train <- xreg.data[1:(nrow(input.data)-split),]
  xreg.test <- xreg.data[(nrow(input.data)-split):nrow(input.data),]
  grid = seq(1:grid.len)
  
  
  for (i in grid) {
    X.fit <- nnetar(X.train, p=3, repeats = 2, scale.inputs = TRUE, xreg=xreg.train, MaxNWts=nMaxNWts)
    X.fc <- forecast(X.fit, h=split, xreg=xreg.test)
    model.mae <- accuracy(X.fc)[3]
    
    true.mae <- mean(abs(as.vector(tail(X.fc$fitted, n=(split+1)))-as.vector(tail(input.data$total_cases, n=(split+1)))))
    if (true.mae < best_mae) {
      best_fit <- X.fit
      best_mae <- true.mae
    }
    cat('Iteration:', i, 'model.mae:', model.mae, ' validation mae:', true.mae, '\n')
  }
  X.fit <- nnetar(X.train, xreg=xreg.train, model = best_fit, MaxNWts=nMaxNWts)
  X.fc <- forecast(X.fit, h=split, xreg=xreg.test)
  plot(X.fc, type='l')
  lines(X.test,col='red')
  #Note: the next line retrains the model on the full data set, so that future forecasts will start
  #   from where the "end" of the training data cuts off. If you want to play with "valication forecasts",
  #   Then comment-out the following line so that future forecasts will start from the 20%-end of the
  #   training set, so you can review results with the hold out data.
  X.fit <- nnetar(ts(input.data$total_cases, frequency=52), xreg=xreg.data, model = best_fit, MaxNWts=nMaxNWts)
  return(X.fit)
}


validate_perf <- function(input.model, input.xreg, input.train.X) {
  #best_model <- get_best_model(train.X.sj, cust.xreg)
  #xreg = train.X.sj %>% select(cust.xreg)
  model.predicted <- input.train.X$total_cases
  
  split <- floor(nrow(input.xreg)*.2)
  
  for (i in (split:nrow(input.xreg))) {
    point.fc <- forecast(input.model, h=1, xreg=input.xreg[i,])
    point.fc <- as.numeric(point.fc$mean)
    if(point.fc < 0) {point.fc <- 0}
    model.predicted[i] <- point.fc
    if (i<(nrow(input.xreg)-1)& 'case_lag1' %in% colnames(input.xreg)) {
      input.xreg$case_lag1[i+1] <- point.fc
      if (i<(nrow(input.xreg)-2) & 'case_lag2' %in% colnames(input.xreg)) {
        input.xreg$case_lag2[i+2] <- point.fc
        if (i<(nrow(input.xreg)-3) & 'case_lag3' %in% colnames(input.xreg)) {
          input.xreg$case_lag3[i+3] <- point.fc
        }
      }
    }
  }
  print(mean(abs(model.predicted - input.train.X$total_cases)))
  plot(input.train.X$total_cases, type='l', xlim=c(split,nrow(input.xreg)))
  lines(model.predicted, col='blue')
}

validate_perf_nocaselag <- function(input.model, input.xreg, input.train.X) {
  #best_model <- get_best_model(train.X.sj, cust.xreg)
  #xreg = train.X.sj %>% select(cust.xreg)
  #model.predicted <- input.train.X$total_cases
  
  split <- floor(nrow(input.xreg)*.2)
  #print(paste('split: ', split))
  X.train <- ts(input.train.X$total_cases[1:(nrow(input.train.X)-split)], frequency = 52)
  X.test <- ts(input.train.X$total_cases[(nrow(input.train.X)-split):nrow(input.train.X)], frequency = 52, start=(end(X.train) + c(0,1)))
  xreg.train <- input.xreg[1:(nrow(input.train.X)-split),]
  xreg.test <- input.xreg[(nrow(input.train.X)-split):nrow(input.train.X),]
  
  X.fc <- forecast(input.model, h=split, xreg=xreg.test)
  #print(nrow(xreg.test))
  model.mae <- accuracy(X.fc)[3]
  print(paste('model MAE: ', model.mae))
  plot(X.fc, type='l')
  lines(X.test, col='red')
  return(X.fc)
}


test_forecast <- function(input.model, input.xreg, input.test.X) {
  return.df = input.test.X
  xreg = input.test.X %>% select(input.xreg)
  for (i in (1:(nrow(return.df)))) {
    point.fc <- forecast(input.model, h=1, xreg=xreg[i,])
    point.fc <- as.numeric(point.fc$mean)
    if(point.fc < 0) {point.fc <- 0}
    return.df$predicted[i] <- point.fc
    if (i<(nrow(return.df)-1) & 'case_lag1' %in% colnames(input.xreg)) {
      xreg$case_lag1[i+1] <- point.fc
      if (i<(nrow(return.df)-2) & 'case_lag2' %in% colnames(input.xreg)) {
        xreg$case_lag2[i+2] <- point.fc
        if (i<(nrow(return.df)-3) & 'case_lag3' %in% colnames(input.xreg)) {
          xreg$case_lag3[i+3] <- point.fc
        }
      }
    }
  }
  plot(return.df$predicted, type='l', main=c(input.test.X$city[1],' predicted cases'))
  return(return.df)
}

test_forecast_nocaselag <- function(input.model, input.xreg, input.test.X) {
  return.df = input.test.X
  xreg = input.test.X %>% select(input.xreg)
  
  city.fc <- forecast(input.model, h=(nrow(xreg)), xreg=xreg)
  plot(city.fc, type='l', main=paste(input.test.X$city[1],' predicted cases'))
  return.df$predicted <- tail(city.fc$fitted, n=nrow(xreg))
  return.df %<>% select(city, year, weekofyear, predicted)
  return.df$predicted <- as.integer(round(return.df$predicted))
  #return.df[return.df$predicted < 0] <- 2
  return(return.df)
}


#full_model <- get_best_model(train.X.sj, c(xreg.full,'case_lag1'),1)
#validate_perf(full_model, train.X.sj %>% select(c(xreg.full,'case_lag1')), train.X.sj)
full_model <- get_best_model(train.X.sj, c(xreg.full),1,6000)
validate_perf(full_model, train.X.sj %>% select(c(xreg.full)), train.X.sj)
sj.fc <- test_forecast(full_model, xreg.full, test.X.sj)




#full_model <- get_best_model(train.X.iq, c(xreg.full,'case_lag1'),1)
#validate_perf(full_model, train.X.iq %>% select(c(xreg.full,'case_lag1')), train.X.iq)
full_model <- get_best_model(train.X.iq, c(xreg.full),1,6000)
validate_perf(full_model, train.X.iq %>% select(c(xreg.full)), train.X.iq)
iq.fc <- test_forecast(full_model, xreg.full, test.X.iq)


submission = read.csv('submission_format.csv')
inner_join(submission, rbind(sj.fc,iq.fc), by=c('city', 'year', 'weekofyear')) %>%
  select(city, year, weekofyear, total_cases=predicted) ->
  prediction.df
prediction.df$total_cases <- round(prediction.df$total_cases)

write.csv(prediction.df, file='AdamsStewartSubmission.csv', row.names = FALSE)

#================================================================
#    Now lets tweak to pick out important vars
#================================================================
rf.select <- train.X.iq %>% select(-c(week_start_date, city, year, case_lag1, case_lag2, case_lag3))
                                      # seas_min_temp, seas_max_temp, seas_avg_temp,
                                      # seas_min_temp1, seas_min_temp2, seas_min_temp3, seas_min_temp4,
                                      # seas_max_temp1,seas_max_temp2,seas_max_temp3,seas_max_temp4,
                                      # seas_avg_temp1,seas_avg_temp2,seas_avg_temp3,seas_avg_temp4))
rf.select <- ungroup(rf.select)
split <- floor(nrow(rf.select)*.2)
y.train <- ts(rf.select$total_cases[1:(nrow(rf.select)-split)], frequency = 52)
y.val <- ts(rf.select$total_cases[(nrow(rf.select)-split):nrow(rf.select)], frequency = 52, start=(end(y.train) + c(0,1)))
X.all <- rf.select %>% select(-c(total_cases))
X.train <- X.all[1:(nrow(X.all)-split),]
X.val <- X.all[(nrow(X.all)-split):nrow(X.all),]

library(party)
cf1 <- cforest(y.train ~ ., data=X.train, control=cforest_unbiased(mtry=2, ntree=50))
#varimp(cf1, conditional = TRUE)
v <- varimpAUC(cf1)
v <- sort(abs(v), decreasing = TRUE)
barplot(v)
cat(paste(names(v), round(v,2), sep=':', collapse='\n'))
cust.xreg10 <- names(v)[1:10]
cust.xreg15 <- names(v)[1:15]
cust.xreg20 <- names(v)[1:20]

#=============================================================
#      Does it perform better?
#=============================================================
cust10_model <- get_best_model(train.X.sj, c(cust.xreg10),1)
mean(abs(as.vector(tail(cust10_model$fitted, n=(split+1)))-as.vector(tail(train.X.sj$total_cases, n=(split+1)))))

validate_perf(cust10_model, train.X.sj %>% select(c(cust.xreg10)), train.X.sj)
sj.val.fc <- validate_perf_nocaselag(cust10_model, train.X.sj %>% select(c(cust.xreg10)), train.X.sj)
sj.fc <- test_forecast(cust10_model, cust.xreg10, test.X.sj)
sj.fc <- test_forecast_nocaselag(cust10_model, cust.xreg10, test.X.sj)
sj.fc$predicted[sj.fc$predicted<0] <- 2

cust15_model <- get_best_model(train.X.sj, c(cust.xreg15),1)
validate_perf(cust15_model, train.X.sj %>% select(c(cust.xreg15)), train.X.sj)
validate_perf_nocaselag(cust15_model, train.X.sj %>% select(c(cust.xreg15)), train.X.sj)
sj.fc <- test_forecast(cust15_model, cust.xreg15, test.X.sj)
sj.fc <- test_forecast_nocaselag(cust15_model, cust.xreg15, test.X.sj)

cust20_model <- get_best_model(train.X.sj, c(cust.xreg20),1)
validate_perf(cust20_model, train.X.sj %>% select(c(cust.xreg20)), train.X.sj)
validate_perf_nocaselag(cust20_model, train.X.sj %>% select(c(cust.xreg20)), train.X.sj)
sj.fc <- test_forecast(cust20_model, cust.xreg20, test.X.sj)

cust10_model_iq <- get_best_model(train.X.iq, c(cust.xreg10),1)
mean(abs(as.vector(tail(cust10_model_iq$fitted, n=(split+1)))-as.vector(tail(train.X.iq$total_cases, n=(split+1)))))


validate_perf(cust10_model_iq, train.X.iq %>% select(c(cust.xreg10)), train.X.iq)
iq.val.fc <- validate_perf_nocaselag(cust10_model_iq, train.X.iq %>% select(c(cust.xreg10)), train.X.iq)
iq.fc <- test_forecast(cust10_model_iq, cust.xreg10, test.X.iq)
iq.fc <- test_forecast_nocaselag(cust10_model_iq, cust.xreg10, test.X.iq)
iq.fc$predicted[iq.fc$predicted<0] <- 2
#=============================================================
#      Combine to one...
#=============================================================
submission = read.csv('submission_format.csv')
inner_join(submission, rbind(sj.fc,iq.fc), by=c('city', 'year', 'weekofyear')) %>%
  select(city, year, weekofyear, total_cases=predicted) ->
  prediction.df

write.csv(prediction.df, file='AdamsStewartSubmission.csv', row.names = FALSE)

#=============================================================
#      Arima better?
#=============================================================
#xreg=stats::lag(stocks.ts, 1)

fit.aa <- auto.arima(y.train, xreg=(X.train %>% select(cust.xreg10)))
fit.aa.fc <- forecast(fit.aa, h=split, xreg=(X.val %>% select(cust.xreg10)))
plot(fit.aa.fc)
lines(y.val, col='red')
fit.aa.acc <- accuracy(fit.aa.fc, as.vector(y.val))
noquote(c("aa w/ reg10:",format(fit.aa.acc[2,], digits=2, nsmall=2)))


fit.aa <- auto.arima(y.train, xreg=(X.train %>% select(cust.xreg15)))
fit.aa.fc <- forecast(fit.aa, h=split, xreg=(X.val %>% select(cust.xreg15)))
plot(fit.aa.fc)
lines(y.val, col='red')
fit.aa.acc <- accuracy(fit.aa.fc, as.vector(y.val))
noquote(c("aa w/ reg15:",format(fit.aa.acc[2,], digits=2, nsmall=2)))


fit.aa <- auto.arima(y.train, xreg=(X.train %>% select(cust.xreg20)))
fit.aa.fc <- forecast(fit.aa, h=split, xreg=(X.val %>% select(cust.xreg20)))
plot(fit.aa.fc)
lines(y.val, col='red')
fit.aa.acc <- accuracy(fit.aa.fc, as.vector(y.val))
noquote(c("aa w/ reg20:",format(fit.aa.acc[2,], digits=2, nsmall=2)))





