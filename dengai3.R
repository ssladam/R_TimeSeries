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
#todo: numerical indicator of season, Spring / Summer / Fall / Winter

full.X.comp %<>% add_group_summaries(c('city', 'year', 'weekofyear'),
                               seas_min_temp = mean(reanalysis_min_air_temp_k),
                               seas_max_temp = mean(reanalysis_max_air_temp_k))
#boxplot(full.X.comp$seas_max_temp)



full.sj <- full.X.comp %>% filter(city=='sj')
full.iq <- full.X.comp %>% filter(city=='iq')
nrow.sj.train <- train.X %>% filter(city=='sj') %>% nrow()
nrow.sj.test <- test.X %>% filter(city=='sj') %>% nrow()
nrow.iq.train <- train.X %>% filter(city=='iq') %>% nrow()
nrow.iq.test <- test.X %>% filter(city=='iq') %>% nrow()
sj.train <- head(full.sj, n=nrow.sj.train)
sj.test <- tail(full.sj, n=nrow.sj.test)
sj.cases <- train.y$total_cases

ccf(sj.cases, sj.train$reanalysis_min_air_temp_k, main='CCF min_air_temp', 52) #8 and 9
ccf(sj.cases, sj.train$reanalysis_max_air_temp_k, main='CCF max_air_temp') #7
ccf(sj.cases, sj.train$precipitation_amt_mm, main='CCF precipitation') #2, 5, both very weak
ccf(sj.cases, sj.train$reanalysis_relative_humidity_percent, main='CCF rel humidity') #3
ccf(sj.cases, sj.train$reanalysis_avg_temp_k, main='CCF avg temp') #8
acf(sj.cases)


check_ccf <- function(x, y, input.title) {
  zz <- ccf(x, y, main=input.title)
  cat('lag:',zz$lag[which.max(zz$acf)], 'max corr:', max(zz$acf),'\n')
  cat('lag:',zz$lag[which.min(zz$acf)], 'min corr:', min(zz$acf))
}

check_ccf(sj.cases, sj.train$weekofyear, 'weekofyear') #include at 0
check_ccf(sj.cases, sj.train$ndvi_nw, 'ndvi_nw') #include at 10
check_ccf(sj.cases, sj.train$ndvi_se, 'ndvi_se') #include at 21
check_ccf(sj.cases, sj.train$precipitation_amt_mm, 'precipitation_amt_mm') #include at 2
check_ccf(sj.cases, sj.train$reanalysis_air_temp_k, 'reanalysis_air_temp_k') #include at 8
check_ccf(sj.cases, sj.train$reanalysis_avg_temp_k, 'reanalysis_avg_temp_k') #include at 8
check_ccf(sj.cases, sj.train$reanalysis_dew_point_temp_k, 'reanalysis_dew_point_temp_k') #include at 8
check_ccf(sj.cases, sj.train$reanalysis_max_air_temp_k, 'reanalysis_max_air_temp_k') #include at 7
check_ccf(sj.cases, sj.train$reanalysis_min_air_temp_k, 'reanalysis_min_air_temp_k') #include at 8
check_ccf(sj.cases, sj.train$reanalysis_precip_amt_kg_per_m2, 'reanalysis_precip_amt_kg_per_m2') #include at 2
check_ccf(sj.cases, sj.train$reanalysis_relative_humidity_percent, 'reanalysis_relative_humidity_percent') #include at 3
check_ccf(sj.cases, sj.train$reanalysis_specific_humidity_g_per_kg, 'reanalysis_specific_humidity_g_per_kg') #include at 8
check_ccf(sj.cases, sj.train$reanalysis_tdtr_k, 'reanalysis_tdtr_k') #include at 0
check_ccf(sj.cases, sj.train$station_avg_temp_c, 'station_avg_temp_c') #include at 10
check_ccf(sj.cases, sj.train$station_diur_temp_rng_c, 'station_diur_temp_rng_c') #include at 25
check_ccf(sj.cases, sj.train$station_max_temp_c, 'station_max_temp_c') #include at 11
check_ccf(sj.cases, sj.train$station_min_temp_c, 'station_min_temp_c') #include at 10
check_ccf(sj.cases, sj.train$seas_min_temp, 'seas_min_temp') #include at 8

full.sj %<>% mutate(lag10_ndvi_nw = lag(ndvi_nw, 10))
full.sj %<>% mutate(lag21_ndvi_se = lag(ndvi_se, 21))
full.sj %<>% mutate(lag2_precipitation_amt_mm = lag(precipitation_amt_mm, 2))
full.sj %<>% mutate(lag8_reanalysis_air_temp_k = lag(reanalysis_air_temp_k, 8))
full.sj %<>% mutate(lag8_reanalysis_avg_temp_k = lag(reanalysis_avg_temp_k, 8))
full.sj %<>% mutate(lag8_reanalysis_dew_point_temp_k = lag(reanalysis_dew_point_temp_k, 8))
full.sj %<>% mutate(lag7_reanalysis_max_air_temp_k = lag(reanalysis_max_air_temp_k, 7))
full.sj %<>% mutate(lag8_reanalysis_min_air_temp_k = lag(reanalysis_min_air_temp_k, 8))
full.sj %<>% mutate(lag2_reanalysis_precip_amt_kg_per_m2 = lag(reanalysis_precip_amt_kg_per_m2, 2))
full.sj %<>% mutate(lag3_reanalysis_relative_humidity_percent = lag(reanalysis_relative_humidity_percent, 3))
full.sj %<>% mutate(lag8_reanalysis_specific_humidity_g_per_kg = lag(reanalysis_specific_humidity_g_per_kg, 8))
full.sj %<>% mutate(lag10_station_avg_temp_c = lag(station_avg_temp_c, 10))
full.sj %<>% mutate(lag25_station_diur_temp_rng_c = lag(station_diur_temp_rng_c, 25))
full.sj %<>% mutate(lag11_station_max_temp_c = lag(station_max_temp_c, 11))
full.sj %<>% mutate(lag10_station_min_temp_c = lag(station_min_temp_c, 10))
full.sj %<>% mutate(lag8_seas_min_temp = lag(seas_min_temp, 8))


full.iq %<>% mutate(lag10_ndvi_nw = lag(ndvi_nw, 10))
full.iq %<>% mutate(lag21_ndvi_se = lag(ndvi_se, 21))
full.iq %<>% mutate(lag2_precipitation_amt_mm = lag(precipitation_amt_mm, 2))
full.iq %<>% mutate(lag8_reanalysis_air_temp_k = lag(reanalysis_air_temp_k, 8))
full.iq %<>% mutate(lag8_reanalysis_avg_temp_k = lag(reanalysis_avg_temp_k, 8))
full.iq %<>% mutate(lag8_reanalysis_dew_point_temp_k = lag(reanalysis_dew_point_temp_k, 8))
full.iq %<>% mutate(lag7_reanalysis_max_air_temp_k = lag(reanalysis_max_air_temp_k, 7))
full.iq %<>% mutate(lag8_reanalysis_min_air_temp_k = lag(reanalysis_min_air_temp_k, 8))
full.iq %<>% mutate(lag2_reanalysis_precip_amt_kg_per_m2 = lag(reanalysis_precip_amt_kg_per_m2, 2))
full.iq %<>% mutate(lag3_reanalysis_relative_humidity_percent = lag(reanalysis_relative_humidity_percent, 3))
full.iq %<>% mutate(lag8_reanalysis_specific_humidity_g_per_kg = lag(reanalysis_specific_humidity_g_per_kg, 8))
full.iq %<>% mutate(lag10_station_avg_temp_c = lag(station_avg_temp_c, 10))
full.iq %<>% mutate(lag25_station_diur_temp_rng_c = lag(station_diur_temp_rng_c, 25))
full.iq %<>% mutate(lag11_station_max_temp_c = lag(station_max_temp_c, 11))
full.iq %<>% mutate(lag10_station_min_temp_c = lag(station_min_temp_c, 10))
full.iq %<>% mutate(lag8_seas_min_temp = lag(seas_min_temp, 8))


main.xreg <- c('weekofyear','lag10_ndvi_nw','reanalysis_tdtr_k','lag10_ndvi_nw',
               'lag21_ndvi_se', 'lag2_precipitation_amt_mm', 'lag8_reanalysis_air_temp_k',
               'lag8_reanalysis_avg_temp_k','lag8_reanalysis_dew_point_temp_k',
               'lag7_reanalysis_max_air_temp_k','lag8_reanalysis_min_air_temp_k',
               'lag2_reanalysis_precip_amt_kg_per_m2', 'lag3_reanalysis_relative_humidity_percent',
               'lag8_reanalysis_specific_humidity_g_per_kg','lag10_station_avg_temp_c',
               'lag25_station_diur_temp_rng_c', 'lag11_station_max_temp_c',
               'lag10_station_min_temp_c','lag8_seas_min_temp',
               'case_lag1')


# full.sj %<>% mutate(lag8_min_air_temp = lag(reanalysis_min_air_temp_k, 8))
# full.sj %<>% mutate(lag7_max_air_temp = lag(reanalysis_max_air_temp_k, 7))
# full.sj %<>% mutate(lag5_precip = lag(precipitation_amt_mm, 5))
# full.sj %<>% mutate(lag3_humid = lag(reanalysis_relative_humidity_percent, 3))
# full.sj %<>% mutate(lag8_avg_temp = lag(reanalysis_avg_temp_k, 8))

# full.iq %<>% mutate(lag8_min_air_temp = lag(reanalysis_min_air_temp_k, 8))
# full.iq %<>% mutate(lag7_max_air_temp = lag(reanalysis_max_air_temp_k, 7))
# full.iq %<>% mutate(lag2_precip = lag(precipitation_amt_mm, 2))
# full.iq %<>% mutate(lag5_precip = lag(precipitation_amt_mm, 5))
# full.iq %<>% mutate(lag3_humid = lag(reanalysis_relative_humidity_percent, 3))
# full.iq %<>% mutate(lag8_avg_temp = lag(reanalysis_avg_temp_k, 8))

sj.train <- head(full.sj, n=nrow.sj.train)
sj.test <- tail(full.sj, n=nrow.sj.test)
iq.train <- head(full.iq, n=nrow.iq.train)
iq.test <- tail(full.iq, n=nrow.iq.test)

full.train <- rbind(sj.train, iq.train)
full.test <- rbind(sj.test, iq.test)
full.X <- rbind(full.train, full.test)


#Now break them back into train / test set, now that feature creation has finished
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


get_best_model <- function(input.data, input.xreg, grid.len=3) {
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
    X.fit <- nnetar(X.train, p=1, repeats = 2, scale.inputs = TRUE, xreg=xreg.train)
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
  X.fit <- nnetar(X.train, xreg=xreg.train, maxit=256, model = best_fit)
  X.fc <- forecast(X.fit, h=split, xreg=xreg.test)
  plot(X.fc, type='l')
  lines(X.test,col='red')
  X.fit <- nnetar(ts(input.data$total_cases, frequency=52), xreg=xreg.data, model = best_fit)
  return(X.fit)
}

validate_perf <- function(input.model, input.xreg, input.train.X) {
  #best_model <- get_best_model(train.X.sj, cust.xreg)
  #xreg = train.X.sj %>% select(cust.xreg)
  model.predicted <- input.train.X$total_cases
  
  split <- floor(nrow(input.xreg)*.8)
  
  for (i in (split:nrow(input.xreg))) {
    point.fc <- forecast(input.model, h=1, xreg=input.xreg[i,])
    point.fc <- as.numeric(point.fc$mean)
    if(point.fc < 0) {point.fc <- 0}
    model.predicted[i] <- point.fc
    if (i<(nrow(input.xreg)-1)) {
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


main.xreg <- c('weekofyear','reanalysis_tdtr_k','lag10_ndvi_nw',
               'lag21_ndvi_se', 'lag2_precipitation_amt_mm', 'lag8_reanalysis_air_temp_k',
               'lag8_reanalysis_avg_temp_k','lag8_reanalysis_dew_point_temp_k',
               'lag7_reanalysis_max_air_temp_k','lag8_reanalysis_min_air_temp_k',
               'lag2_reanalysis_precip_amt_kg_per_m2', 'lag3_reanalysis_relative_humidity_percent',
               'lag8_reanalysis_specific_humidity_g_per_kg','lag10_station_avg_temp_c',
               'lag25_station_diur_temp_rng_c', 'lag11_station_max_temp_c',
               'lag10_station_min_temp_c','lag8_seas_min_temp',
               'case_lag1')

cust.xreg <- c('weekofyear','reanalysis_tdtr_k',
               'lag2_precipitation_amt_mm',
               'lag8_reanalysis_avg_temp_k','lag8_reanalysis_dew_point_temp_k',
               'lag7_reanalysis_max_air_temp_k','lag8_reanalysis_min_air_temp_k',
               'lag11_station_max_temp_c',
               'lag10_station_min_temp_c',
               'case_lag1')


best_model <- get_best_model(train.X.sj, main.xreg,1)
best_model <- get_best_model(train.X.sj, cust.xreg,1)
validate_perf(best_model, train.X.sj %>% select(main.xreg), train.X.sj)


# xreg = train.X.sj %>% select(cust.xreg)
# model.predicted <- train.X.sj$total_cases
# 
# split <- floor(nrow(xreg)*.8)
# 
# for (i in (split:nrow(xreg))) {
#   point.fc <- forecast(best_model, h=1, xreg=xreg[i,])
#   point.fc <- as.numeric(point.fc$mean)
#   if(point.fc < 0) {point.fc <- 0}
#   model.predicted[i] <- point.fc
#   if (i<(nrow(xreg)-1)) {
#     xreg$case_lag1[i+1] <- point.fc
#     if (i<(nrow(xreg)-2) & 'case_lag2' %in% colnames(xreg)) {
#       xreg$case_lag2[i+2] <- point.fc
#       if (i<(nrow(xreg)-3) & 'case_lag3' %in% colnames(xreg)) {
#         xreg$case_lag3[i+3] <- point.fc
#       }
#     }
#   }
# }
# mean(abs(model.predicted - train.X.sj$total_cases))
# plot(train.X.sj$total_cases, type='l', xlim=c(750,933))
# lines(model.predicted, col='blue')


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
plot(test.X.sj$predicted, type='l', main='SJ predicted cases')




xreg.iq = test.X.iq %>% select(cust.xreg)
best_model <- get_best_model(train.X.iq, cust.xreg,1)
validate_perf(best_model, train.X.sj %>% select(cust.xreg), train.X.sj)


for (i in (1:(nrow(test.X.iq)))) {
  point.fc <- forecast(best_iq, h=1, xreg=xreg.iq[i,])
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

plot(test.X.iq$predicted, type='l', main='IQ predicted cases')

submission = read.csv('submission_format.csv')
inner_join(submission, rbind(test.X.sj,test.X.iq), by=c('city', 'year', 'weekofyear')) %>%
  dplyr::select(city, year, weekofyear, total_cases=predicted) ->
  prediction.df

prediction.df$total_cases <- round(prediction.df$total_cases)

write.csv(prediction.df, file='AdamsStewartSubmission.csv', row.names = FALSE)


