rm(list=ls())
pkgs <- c('tidyverse','tidyimpute', 'tibbletime', 'corrplot', 'magrittr', 'zoo', 'RColorBrewer', 'gridExtra','MASS')
invisible(lapply(pkgs, require, character.only = T))
rm(pkgs)
select <- dplyr::select

setwd("/temp/NU/413/dengai")
#X <- read.csv("dengue_features_train.csv") #reads in as a list...
X <- read_csv("dengue_features_train.csv")
y <- read_csv("dengue_labels_train.csv")

X %<>% left_join(y, by=c('city', 'year', 'weekofyear')) #brings the response into same frame

X.sj <- X %>% filter(city=='sj')
y.sj <- y %>% filter(city=='sj')
X.iq <- X %>% filter(city=='iq')
y.iq <- y %>% filter(city=='iq')

# impute NAs by the mean value
#X.sj$ndvi_ne <- na.approx(X.sj$ndvi_ne)
#X.iq$ndvi_ne <- na.approx(X.iq$ndvi_ne)
X.sj <- X.sj %>% impute_mean()
X.iq <- X.iq %>% impute_mean()



X.sj.time <- X.sj %>% as_tbl_time(index=week_start_date)
X.sj.monthly <- X.sj.time %>% collapse_by('monthly') %>% group_by(week_start_date) %>%
  summarise(cases = sum(total_cases))
sj.stl <- (stl(ts(X.sj.monthly$cases, frequency=12), s.window='periodic'))
plot(sj.stl)


#Use linear regression to help select most important predictors

X.sj.scale <- X.sj %>% select(-c(city,year,weekofyear,week_start_date)) %>%
  mutate_all(funs(scale(.) %>% as.vector))

fit.lm.sj <- lm(total_cases ~ ., data=X.sj.scale)
plot(fit.lm.sj)
sort.coef <- sort(abs(fit.lm.sj$coefficients[2:21]))
bp <- barplot(sort.coef, horiz = TRUE, axisnames = FALSE, col='lightblue', main='LM most significant coeffs')
text(y=bp, x=0, names(sort.coef),pos=4, cex=0.8)




#auto.arima
#xreg=stats::lag(stocks.ts, 1)
xreg = X.sj %>% select(names(tail(sort.coef,n=7)))
fit.aa.sj<-auto.arima(X.sj[1:(nrow(X.sj)-136),25], xreg=xreg[1:(nrow(X.sj)-136),])
fit.aa.sj.fc<-forecast(fit.aa.sj, h=136,  xreg=xreg[(nrow(X.sj)-136):nrow(X.sj),])
plot(fit.aa.sj.fc)
lines(as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]), x=(nrow(X.sj)-136):nrow(X.sj), col='red')
fit.aa.sj.acc <- accuracy(fit.aa.sj.fc, as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]))
noquote(c("aa w/ xreg:",format(fit.aa.sj.acc[2,], digits=2, nsmall=2)))

#NN with regressors
fit.nn.sj <- nnetar(X.sj$total_cases[1:(nrow(X.sj)-136)],
                    p=1, repeats = 10, maxit=150,
                    scale.inputs = TRUE, xreg=xreg[1:(nrow(X.sj)-136),])
fit.nn.sj.fc <- forecast(fit.nn.sj, h=136, xreg=xreg[(nrow(X.sj)-136):nrow(X.sj),])
plot(fit.nn.sj.fc, type='l')
lines(as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]), x=(nrow(X.sj)-136):nrow(X.sj), col='red')
fit.nn.sj.acc <- accuracy(fit.nn.sj.fc, as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]))
noquote(c("xreg NN:",format(fit.nn.sj.acc[2,], digits=2, nsmall=2)))


#That didn't work very well... let's choose our own regressors
cust.xreg = c('station_max_temp_c', 'station_avg_temp_c', 'reanalysis_sat_precip_amt_mm',
              'reanalysis_specific_humidity_g_per_kg', 'ndvi_se', 'ndvi_sw',
              'ndvi_ne','ndvi_nw')
xreg2 = X.sj %>% select(cust.xreg)

#auto.arima
fit.aa.sj<-auto.arima(X.sj[1:(nrow(X.sj)-136),25], xreg=xreg2[1:(nrow(X.sj)-136),])
fit.aa.sj.fc<-forecast(fit.aa.sj, h=136,  xreg=xreg2[(nrow(X.sj)-136):nrow(X.sj),])
plot(fit.aa.sj.fc)
lines(as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]), x=(nrow(X.sj)-136):nrow(X.sj), col='red')
fit.aa.sj.acc <- accuracy(fit.aa.sj.fc, as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]))
noquote(c("aa w/ xreg:",format(fit.aa.sj.acc[2,], digits=2, nsmall=2)))

#NN with regressors
my.test <- ts(X.sj$total_cases[1:(nrow(X.sj)-136)], frequency = 52)
fit.nn.sj <- nnetar(my.test,
                    p=1, repeats = 10, maxit=150, size=8,
                    scale.inputs = TRUE, xreg=xreg2[1:(nrow(X.sj)-136),])
fit.nn.sj.fc <- forecast(fit.nn.sj, h=136, xreg=xreg2[(nrow(X.sj)-136):nrow(X.sj),])
plot(fit.nn.sj.fc, type='l')
my.train <- ts(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)], frequency = 52, start=c(16,21))
lines(my.train,col='red')
#lines(as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]), x=(nrow(X.sj)-136):nrow(X.sj), col='red')
fit.nn.sj.acc <- accuracy(fit.nn.sj.fc, as.vector(X.sj$total_cases[(nrow(X.sj)-136):nrow(X.sj)]))
fit.nn.sj.acc <- accuracy(fit.nn.sj.fc, my.train)
noquote(c("xreg NN:",format(fit.nn.sj.acc[2,], digits=2, nsmall=2)))







X.iq.time <- X.iq %>% as_tbl_time(index=week_start_date)
X.iq.monthly <- X.iq.time %>% collapse_by('monthly') %>% group_by(week_start_date) %>%
  summarise(cases = sum(total_cases))

iq.stl <- (stl(ts(X.iq.monthly$cases, frequency=12), s.window='periodic'))
plot(iq.stl)



cat('\nSan Juan\n',
    '\t features: ', X.sj %>% ncol, 
    '\t entries: ' , X.sj %>% nrow,
    '\t labels: '  , y.sj %>% nrow)
cat('\nIquitos\n',
    '\t features: ', X.iq %>% ncol, 
    '\t entries: ' , X.iq %>% nrow,
    '\t labels: '  , y.iq %>% nrow)

# Remove `week_start_date` string.
#X.sj %<>% dplyr::select(-week_start_date)
#X.iq %<>% dplyr::select(-week_start_date)
#X.sj %<>% dplyr::select(-c(city, year, weekofyear))
#X.iq %<>% dplyr::select(-c(city, year, weekofyear))



# count missing values (as percent)
X %>% group_by(city) %>%
  summarise_all(funs(round(100*mean(is.na(.)),2))) %>% t

X.sj %>%
  mutate(index = as.numeric(row.names(.))) %>%
  ggplot(aes(index, ndvi_ne)) + 
  geom_line(colour = 'dodgerblue') +
  ggtitle("Vegetation Index over Time")




#Model selection... is Poisson or NB more appropriate?
y %>% group_by(city) %>% dplyr::select(city, total_cases) %>%
  summarise(mean_cases = mean(total_cases), var_cases = var(total_cases), total_cases = sum(total_cases))
#variance and mean are very different. Strong case for NB. Consider ZINB.


# total cases of dengue: histograms
rbind(y.iq, y.sj) %>% 
  ggplot(aes(x = total_cases,fill = ..count..)) + 
  geom_histogram(bins = 12, colour = 'black') + ggtitle('Total Cases of Dengue') +
  scale_y_continuous(breaks = seq(0,700,100)) + facet_wrap(~city)

# corerlations between features
X.sj %<>% mutate('total_cases' = y.sj$total_cases)
X.iq %<>% mutate('total_cases' = y.iq$total_cases)

# plot san juan correlation matrix
X.sj %>% 
  dplyr::select(-city, -year, -weekofyear) %>%
  cor(use = 'pairwise.complete.obs') -> M1
corrplot(M1, type="lower", method="color",
         col=brewer.pal(n=8, name="RdBu"),diag=FALSE)


# plot iq correlation matrix
X.iq %>% 
  dplyr::select(-city, -year, -weekofyear) %>%
  cor(use = 'pairwise.complete.obs') -> M2
corrplot(M2, type="lower", method="color",
         col=brewer.pal(n=8, name="RdBu"),diag=FALSE)


# see the correlations as barplot
sort(M1[21,-21]) %>%  
  as.data.frame %>% 
  `names<-`('correlation') %>%
  ggplot(aes(x = reorder(row.names(.), -correlation), y = correlation, fill = correlation)) + 
  geom_bar(stat='identity', colour = 'black') + scale_fill_continuous(guide = FALSE) + scale_y_continuous(limits =  c(-.15,.25)) +
  labs(title = 'San Jose\n Correlations', x = NULL, y = NULL) + coord_flip() -> cor1

# can use ncol(M1) instead of 21 to generalize the code
sort(M2[21,-21]) %>%  
  as.data.frame %>% 
  `names<-`('correlation') %>%
  ggplot(aes(x = reorder(row.names(.), -correlation), y = correlation, fill = correlation)) + 
  geom_bar(stat='identity', colour = 'black') + scale_fill_continuous(guide = FALSE) + scale_y_continuous(limits =  c(-.15,.25)) +
  labs(title = 'Iquitos\n Correlations', x = NULL, y = NULL) + coord_flip() -> cor2

grid.arrange(cor1, cor2, nrow = 1)

# Now get rid of the label from our training data!
X.sj %<>% dplyr::select(-total_cases)
X.iq %<>% dplyr::select(-total_cases)
X.sj %<>% dplyr::select(-city)
X.iq %<>% dplyr::select(-city)

X.sj.train <- head(X.sj, 800)
X.sj.test <- tail(X.sj, (nrow(X.sj) - 800))
y.sj.train <- head(y.sj, 800)
y.sj.test <- tail(y.sj, (nrow(X.sj) - 800))

mae <- function(error) return(mean(abs(error)))

opt_nb_model <- function(train, test, y.train, y.test) {
  form <- paste("total_cases ~ 1 + ", paste(colnames(train), collapse=" + "))
  grid <- 10^(seq(-8, -1, 1))
  train$total_cases <- y.train$total_cases
  test$total_cases <- y.test$total_cases
  
  best_alpha <- c()
  best_score <- 1000
  
  # Find best model parameters...
  for (i in grid) {
    model <- glm.nb(formula = form, data=train, init.theta = i)
    results <- predict(model, test)
    score <- mae(test$total_cases - results)
    
    if(score < best_score) {
      best_alpha <- i
      best_score <- score
      cat('\nBest score=', best_score)
    }
  }
  
  #refit the best model
  model <- glm.nb(formula = form, data=train, init.theta = best_alpha)
  return (model)
}

sj.model <- opt_nb_model(X.sj.train, X.sj.test, y.sj.train, y.sj.test)

# plot sj

sj.pred <- data.frame(total_cases = y.sj.train$total_cases)
sj.pred$pred_cases <- predict(sj.model, X.sj.train, type = 'response')
sj.pred %>% 
  mutate(index = as.numeric(row.names(.))) %>%
  ggplot(aes(x = index)) + ggtitle("San Jose") +
  geom_line(aes(y = total_cases, colour = "total_cases")) + 
  geom_line(aes(y = pred_cases, colour = "fitted"))

