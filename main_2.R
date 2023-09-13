# Library -----------------------------------------------------------------

source('utils_v2.R')

# Load data ---------------------------------------------------------------

recordings <- get_audio_objects("processed-recs-paper")

# Extract Features --------------------------------------------------------

mel_features_list <- mel_filter_features(recordings, n = 200)
mel_features_list <- sample(mel_features_list)

# CV for bandwidth ----------------------------------------------------------------------

tries <- rep(NA, 50)
grid <- seq(0,10,length.out=50)

for (i in 1:50) {
  LOOCV <- loo_cv(mel_features_list, h=grid[i])
  tries[i] <- LOOCV$L2
}
plot(grid,tries,xlab='Bandwidth',ylab='RMSE',main='LOOCV for bandwidth')
h_hat <- grid[which.min(tries)]
abline(v = h_hat,lty=2,col='red')


# CV for k in kNN ---------------------------------------------------------

kNN_tries[3]
kNN_grid <- seq(2,15,length.out=14)
kNN_tries <- rep(NA,14)

for (i in 1:14) {
  LOOCV <- kNN_loo_cv(mel_features_list, K=kNN_grid[i])
  kNN_tries[i] <- LOOCV$L2
}


plot(kNN_grid,kNN_tries,xlab='K',ylab='RMSE',main='LOOCV for K in kNN')
k_hat <- kNN_grid[which.min(kNN_tries)]
abline(v = k_hat,lty=2,col='red')

#bag <- bagging_cv(mel_features_list, N = 100)
#mean(bag)

cv_result <- loo_cv(mel_features_list, h = grid[which.min(tries)])

cv_ranges_report <- cv_ranges_scores(mel_features_list, 
                                     metric_fun = Accuracy)


# Best models -------------------------------------------------------------


KR <- c('Model'='Kernel regression','Optimized hyperparameter'='Bandwidth',
        'Validated hyperparameter'=2.0408, 'RMSE'=5.8674)

kNN <- c('Model'='kNN regression','Optimized hyperparameter'='K',
         'Validated hyperparameter'=2, 'RMSE'=6.3303)

df <- data.frame(rbind(KR,kNN))

# Predictions: kernel regression ------------------------------------------

best_KR <- min(tries)
c(train,test) := train_test_split(mel_features_list, test_ratio=0.2)
KR_preds <- KR_predict_set(train,test,h=grid[which.min(tries)])
RMSE(KR_preds)
Accuracy(KR_preds, degrees=5)


# Predictions: kNN regression ---------------------------------------------

best_kNN <- min(kNN_tries)
c(train,test) := train_test_split(mel_features_list, test_ratio=0.2)
KR_preds <- kNN_predict_set(train,test,K=grid[which.min(kNN_tries)])
RMSE(kNN_preds)
Accuracy(kNN_preds, degrees=5)

# Plots -------------------------------------------------------------------

plot(1:200, rep(0,200),col = "white", ylim  = c(-20,-5), main = "Exctracted Features hot water", xlab = "x", ylab = expression(phi(x)))
for (i in 1:length(mel_features_list)){
  audio <- mel_features_list[[i]]
  if(audio$Label > 90 && audio$Label < 100){
    points(1:200, audio$Mel,"l",col = rgb( 1, runif(1), 0, 1),lwd = 3)
  }
}

plot(1:200, rep(0,200),col = "white", ylim  = c(-20,-5), main = "Exctracted Features cold water", xlab = "x", ylab = expression(phi(x)))

for (audio in mel_features_list){
  if(audio$Label > 0 && audio$Label <10){
    points(1:200, audio$Mel,"l",col = rgb( 0, runif(1), 1, 0.2),lwd = 3)
  }
}

legend(170, -5, legend=c("Hot water", "Cold water"),
       col=c("red", "blue"), lty=1:1, cex=0.5)


  