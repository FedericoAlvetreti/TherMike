# Libraries and basics -----------------------------------------------------

# Libraries
library(readr)
library(tuneR)
library(fda)
library(seewave)


# Operator used to assign multiple variables in the same line of code
':=' <- function(lhs, rhs) {
  frame <- parent.frame()
  lhs <- as.list(substitute(lhs))
  if (length(lhs) > 1)
    lhs <- lhs[-1]
  if (length(lhs) == 1) {
    do.call(`=`, list(lhs[[1]], rhs), envir=frame)
    return(invisible(NULL)) 
  }
  if (is.function(rhs) || is(rhs, 'formula'))
    rhs <- list(rhs)
  if (length(lhs) > length(rhs))
    rhs <- c(rhs, rep(list(NULL), length(lhs) - length(rhs)))
  for (i in 1:length(lhs))
    do.call(`=`, list(lhs[[i]], rhs[[i]]), envir=frame)
  return(invisible(NULL)) 
}


# Chunk x in n equal sized vectors
chunk_in_n <- function(x, n) split(x, cut(seq_along(x), n, labels = FALSE))


# Split dataset into train and test
train_test_split <- function(dataset, test_ratio = 0.2){
  
  n <- length(dataset)
  train_ratio <- 1 - test_ratio 
  
  train_idx <- sample(c(TRUE, FALSE), n, replace = TRUE, 
                      prob = c(train_ratio, test_ratio))
  
  train  <- dataset[train_idx]
  test   <- dataset[!train_idx]
  
  return(list("Train" = train, "Test"= test))
}


# Data manipulation -----------------------------------------------------


# Define a custom audio object given a .wav file path
get_audio_object <- function(audio_file_path){
  
  # Load audio file
  audio <- readWave(audio_file_path)
  
  # Get sample size
  n <- length(audio@left)
  
  # Get sampling rate
  sampling_rate <- audio@samp.rate
  
  # Get frequencies
  frequencies <- c(0:(n-1)) * (sampling_rate / n)
  frequencies <- frequencies[1:(n %/% 2)] # Halve since it is real valued
  
  # Get power spectrum
  power_spectrum <- Mod(fft(c(audio@left))^2) / n
  power_spectrum <- power_spectrum[1:(n %/% 2)] # Halve since it is real valued
  
  # Get label
  num_part <- unlist(strsplit(audio_file_path, "_"))[4]
  label <- as.numeric(substring(num_part, 1, nchar(num_part) - 4))
  
  # Get author
  author <- unlist(strsplit(audio_file_path, "_"))[3]
  
  return(list("Author" = author,
              "Recording" = audio@left,
              "Label" = label,
              "Frequencies" = frequencies,
              "Power_spectrum" = power_spectrum,
              "Sample_size" = n,
              "Sampling_rate" = sampling_rate))
}


# Get all audio objects given a folder path
get_audio_objects <- function(directory_path){
  
  files <- list.files(path = directory_path) # Get audio file names
  n_iter <- length(files) # Get number of iterations
  
  # Set up a progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style [1, 2, 3]
                       width = 50,   # Progress bar width
                       char = "=")   # Character used to create the bar
  
  i = 1 # Set iterator to update the progress bar
  
  audio_objects <- list()
  
  for(file_path in files){
    audio_objects <- append(audio_objects, 
                            list(get_audio_object(paste(directory_path,
                                                        file_path,
                                                        sep="\\"))))
    setTxtProgressBar(pb, i) # Update progress bar
    i = i + 1 
  }
  
  close(pb) 
  return(audio_objects)
}


# Compute the spectrogram of a given audio object
spectro <- function(data, nfft = 1024, window = 256, overlap = 128, t0=0, 
                   plot_spec = T, normalize = T, return_data = F){
  
  library(signal)
  library(oce)
  
  # Get audio signal
  snd = data$Recording
  
  # Demean to remove DC offset
  snd = snd - mean(snd)
  
  # Calculate duration
  dur = length(snd) / data$Sampling_rate
  
  # Create spectrogram
  spec = specgram(x = snd,
                  n = nfft,
                  Fs = data$Sampling_rate,
                  window = window,
                  overlap = overlap)
  
  
  P = abs(spec$S) # Discard phase info
  
  if(normalize) P = P/max(P) # Normalize
  
  P = 10 * log10(P) # Convert to dB
  
  # Config time axis
  if(t0 == 0) t = as.numeric(spec$t) else t = as.POSIXct(spec$t, origin = t0)
  
  f = spec$f # Rename freq
  
  if(plot_spec){
    
    # Change plot color defaults
    par(bg = "black")
    par(col.lab="white")
    par(col.axis="white")
    par(col.main="white")
    
    # Plot spectrogram
    imagep(t, f, t(P), col = oce.colorsViridis, drawPalette = T,
           ylab = 'Frequency [Hz]', axes = F)
    
    box(col = 'white')
    axis(2, labels = T, col = 'white')
    title(main = data$Label)
    
    # Add x axis
    if(t0 == 0){
      
      axis(1, labels = T, col = 'white')
      
    } else {
      
      axis.POSIXct(seq.POSIXt(t0, t0+dur, 10), side = 1, 
                   format = '%H:%M:%S', col = 'white', las = 1)
      
      mtext(paste0(format(t0, '%B %d, %Y')), side = 1, adj = 0,
            line = 2, col = 'white')
    }
  }
  
  if(return_data){
    
    spec = list(
      t = t,
      f = f,
      p = t(P)
    )
    
    return(spec)  
  }
}


# Distances ---------------------------------------------------------------


# Euclidean distance between f1 and f2
L2 <- function(f1, f2){return(sqrt(sum(((f1 - f2)^2))))}

# Root Mean Squared Error 
RMSE <- function(preds){
  f1 <- preds$Label
  f2 <- preds$Prediction
  return(sqrt(sum(((f1 - f2)^2))/nrow(preds)))
}

# Mean absolute difference between two vectors
Mean_abs_diff <- function(v1,v2){return(mean(abs(v1-v2)))}


# Gaussian kernel
Kernel <- function(x) {return(0.5 * exp(-0.5 * x^2))}


# Metrics -----------------------------------------------------------------


# Euclidean distance between labels and predictions
L2_score <- function(preds) {
  f1 <- preds$Label
  f2 <- preds$Prediction
  return(sqrt(sum(((f1 - f2)^2))))
}



# Mean absolute difference between labels and predictions
Mean_abs_diff_score <- function(preds){
  v1 <- preds$Label
  v2 <- preds$Prediction
  
  return(mean(abs(v1-v2)))
}


# Set up a sensitivity range to evaluate accuracy of each prediction
Accuracy <- function(preds, degrees = 5){
  
  k <- 0
  n <- dim(preds)[1]
  
  if(n == 0) return(NA)
  
  for(i in (1:n)){
    
    prediction <- preds$Prediction[i]
    label <- preds$Label[i]
    
    if(label < prediction + degrees && label > prediction - degrees) k <- k + 1
  }
  return(k / n)
}


# Plots -------------------------------------------------------------------


# Plot audio object in time domain
plot_audio <- function(audio_object){
  
  time <- (1 : audio_object$Sample_size)
  audio <- audio_object$Recording
  
  plot(time,
       audio,
       type = "l",
       ylab = "Amplitude", 
       xlab = "Time",
       main = paste(audio_object$Label, "°", sep = ""))
}


# Plot audio object's power spectrum
plot_power_spectrum <- function(audio_object){
  
  plot(audio_object$Frequencies,
       audio_object$Power_spectrum,
       type = "l",
       ylab = "Power",
       xlab = "Frequency",
       main = paste(audio_object$Label, "°", sep = ""))
}


# Plot labels' distribution given a set of recordings and k number of bins
plot_labels <- function(recordings, k){
  
  labels <- list()
  
  for(recording in recordings){
    labels <- append(labels, recording$Label)
  }
  
  labels <- unlist(labels)
  hist(labels, breaks = k, probability = F, 
       ylab = "Abs Frequency",
       xlab = "Label",
       main = "Recordings distribution")
}


# Functional Feature Extraction -------------------------------------------


# Get a smooth approximation given a power spectrum and new basis parameters
get_smooth_from_ps <- function(basis_func, n_basis, audio_object){
  
  # Set up new basis
  my_basis <- basis_func(range(0, (audio_object$Sampling_rate/2)),
                         nbasis = n_basis)
  
  freq <- audio_object$Frequencies
  power <- audio_object$Power_spectrum
  
  # Get parameters of the power spectrum approximation
  params <- smooth.basis(argvals = freq,
                         y = power ,
                         fdParobj = my_basis)$fd$params
  
  
  return(list("Params" = params, "Label" = audio_object$Label))
}


# Apply smooth approximations to a list of audio objects 
get_smooth_list_from_ps <- function(basis_func, n_basis, list){
  
  n_iter <- length(list) # Get number of iterations
  
  # Set up a progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style [1, 2, 3]
                       width = 50,   # Progress bar width
                       char = "=")   # Character used to create the bar
  
  i = 1 # Set iterator to update the progress bar
  
  smooth_list <- list()
  
  for(audio_object in list){
    smooth_list <- append(smooth_list, list(get_smooth_from_ps(basis_func, 
                                                     n_basis, 
                                                     audio_object)))
    
    setTxtProgressBar(pb, i) # Update progress bar
    i = i + 1 
  }
  
  close(pb) 
  return(smooth_list)
}


# MEL Feature Extraction --------------------------------------------------


# Compute MFCCs of an audio's power spectrum using mel-filterbank analysis
mel_filter_feature <- function(audio_obj, n = 200){
  
  # Retrieve the power spectrum
  power_spectrum <- audio_obj$Power_spectrum
  
  # Calculate the mel filterbank
  mel_filter_bank <- melfilterbank(f = audio_obj$Sampling_rate,
                                   wl = 2 * length(power_spectrum),
                                   m = n)
  
  mel_filter <- mel_filter_bank$amp # Extract mel filter amplitudes 
  mel_freq <- 1000 * mel_filter_bank$central.freq # Extract central frequencies
  
  # Apply mel filterbank to the power spectrum and normalize
  mel_ps <- (c(power_spectrum %*% mel_filter)) / colSums(mel_filter)
  mel_ps <- log(mel_ps / sum(mel_ps * mel_freq) ) 
  
  return(list("Mel" = mel_ps, 
              "Label" = audio_obj$Label,
              "Author"= audio_obj$Author))
}


# Apply mel_filter_feature to a list of audio objects
mel_filter_features <- function(audio_obj_list, n = 200){
  
  n_iter <- length(audio_obj_list) # Get number of iterations
  
  # Set up a progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style [1, 2, 3]
                       width = 50,   # Progress bar width
                       char = "=")   # Character used to create the bar
  
  i = 1 # Set iterator to update the progress bar
  
  mel_features <- list()
  
  for(audio_object in audio_obj_list){
    mel_features <- append(mel_features,
                           list(mel_filter_feature(audio_object, n)))
    
    setTxtProgressBar(pb, i) # Update progress bar
    i = i + 1 
  }
  
  close(pb) 
  return(mel_features)
}


# kNN regression ----------------------------------------------------------

kNN_predict_audio <- function(train_set, new_data, K) {
  
  # Initialize data structure to store distances and labels
  distances_matrix <- matrix(NA,nrow=length(train_set),ncol=2)
  
  # Main loop for pair-wise distances
  for (i in 1:length(train_set)) {
    audio <- train_set[[i]]
    distances_matrix[i,1] <- L2(audio$Mel, new_data$Mel)
    distances_matrix[i,2] <- audio$Label
  }
  
  # Sort by increasing distance and retrieve the kNN
  idxs <- order(distances_matrix[,1])[1:K]
  
  # Retrieve the mean of the labels of the kNN
  return(mean(distances_matrix[idxs,2]))
}

kNN_predict_set <- function(train_set, test_set, K) {
  
  # Set up a dataframe to store predictions
  preds <- data.frame()
  
  # Predict for each element in the test_set
  for(audio in test_set){
    preds <- rbind(preds, 
                   c(audio$Label, kNN_predict_audio(train_set, audio, K)))
  }
  
  names(preds) <- c("Label", "Prediction") # Add column names
  preds <- preds[order(preds$Label), ] # Sort by Label
  
  return(preds)
}

kNN_loo_cv <- function(features_list, K) {
  
  n <- length(features_list)  # Get size
  
  #LOOCV is k-fold with k = n
  return(kNN_k_folds_cv(features_list, f = n, K))
}

# K-folds cross validation
kNN_k_folds_cv <- function(features_list, f = 10, K){
  
  n <- length(features_list)  # Get size
  folds <- chunk_in_n(1:n, f) # Split data in folds
  
  # Initialize metrics
  L2_err <- 0
  acc <- 0
  abs_err <- 0
  
  for(fold in folds){
    
    test <- features_list[fold] # Get test features
    train <- features_list[-fold] # Get train features
    
    # Evaluate predictions
    preds <- kNN_predict_set(train, test, K)
    
    # Increment metrics
    L2_err <- L2_err + L2(preds$Label, preds$Prediction)
    acc <- acc + Accuracy(preds)
    abs_err <- abs_err + Mean_abs_diff(preds$Label, preds$Prediction)
  }
  
  # Compute final metrics
  L2_err <- L2_err / f
  acc <- acc / f
  abs_err <- abs_err / f
  
  return(list("L2" = L2_err, "Accuracy" = acc, "Mean_abs" = abs_err))
}

# KR prediction ----------------------------------------------------------


# Single prediction
KR_predict_audio <- function(train_set, new_data, h){
  
  # Set up a dataframe to save distances
  weights_df <- data.frame()
  
  # Calculate L2 distance between new_data and each element in train
  for(train_audio in train_set){
    weight <-  Kernel(L2(train_audio$Mel, new_data$Mel)/h) 
    weights_df <- rbind(weights_df, c(weight, train_audio$Label))
  }
  
  names(weights_df) <- c("Weight", "Label")
  
  # Evaluate prediction
  tot_weight <- sum(weights_df$Weight)
  prediction <- sum(weights_df$Weight * weights_df$Label) / tot_weight
  
  if(is.na(prediction)) prediction <- 0
  
  return(prediction)
}

# Multiple prediction
KR_predict_set <- function(train_set, test_set, h){
  
  # Set up a dataframe to store predictions
  preds <- data.frame()
  
  # Predict for each element in the test_set
  for(audio in test_set){
    preds <- rbind(preds, 
                   c(audio$Label, KR_predict_audio(train_set, audio, h)))
  }
  
  names(preds) <- c("Label", "Prediction") # Add column names
  preds <- preds[order(preds$Label), ] # Sort by Label
  
  return(preds)
}


# CV ----------------------------------------------------------------------

# Leave-one-out-cross-validation

loo_cv <- function(features_list, h = 0.5, k) {
  
  n <- length(features_list)  # Get size
  
  #LOOCV is k-fold with k = n
  return(k_folds_cv(features_list, k = n, h = 0.5))
}

# K-folds cross validation
k_folds_cv <- function(features_list, k = 10, h = 0.5){
  
  n <- length(features_list)  # Get size
  folds <- chunk_in_n(1:n, k) # Split data in folds
  
  # Initialize metrics
  L2_err <- 0
  acc <- 0
  abs_err <- 0
  
  for(fold in folds){
    
    test <- features_list[fold] # Get test features
    train <- features_list[-fold] # Get train features
    
    # Evaluate predictions
    preds <- KR_predict_set(train, test, h)
    
    # Increment metrics
    L2_err <- L2_err + L2(preds$Label, preds$Prediction)
    acc <- acc + Accuracy(preds)
    abs_err <- abs_err + Mean_abs_diff(preds$Label, preds$Prediction)
  }
  
  # Compute final metrics
  L2_err <- L2_err / k 
  acc <- acc / k
  abs_err <- abs_err / k
  
  return(list("L2" = L2_err, "Accuracy" = acc, "Mean_abs" = abs_err))
}

# Leave-one-out-cross-validation
loo_cv <- function(features_list, h = 0.5) {
  
  n <- length(features_list)  # Get sample size
  
  # LOOCV is k-fold with k = n
  return(k_folds_cv(features_list, h, k = n))
}

# Evaluate a metric function over a set of ranges
evaluate_ranges <- function(preds, metric_fun, degrees = 5) {
  
  preds <- preds[order(preds$Label), ] # Sort by label
  
  # Calculate lower bound for each degree range
  lower_bounds <- seq(0, 100, by = degrees)
  lower_bounds <- lower_bounds[1 : (length(lower_bounds)-1)]
  
  report_df <- data.frame() # Set up a dataframe to store results
  
  for(lower_bound in lower_bounds) {
    
    upper_bound = lower_bound  + degrees
    
    # Get current range predictions subset
    preds_subset <- preds[preds$Label > lower_bound,]
    preds_subset <- preds_subset[preds_subset$Label < upper_bound,]
    
    # Write lower-upper bound interval
    interval <- paste(as.character(lower_bound),
                      as.character(upper_bound),
                      sep = "-")
    
    # Add row to report
    row <- data.frame(interval = interval, 
                      score = round(metric_fun(preds_subset),3))
    
    report_df <- rbind(report_df, row)
  }
  
  return(report_df)
}


# Apply k-fold cross validation and custom metric function to ranges of y
cv_ranges_scores <- function(features_list,
                            metric_fun = Accuracy,
                            k = 10,
                            h = 0.5,
                            degrees = 20){
  
  n <- length(features_list) # Get size
  folds <- chunk_in_n(1:n, k) # Split data in folds
  
  n_int <- floor(100 / degrees) # Get number of intervals
  int_scores <- rep(0, n_int)  # Set up a list of scores
  
  for(fold in folds){
    
    test <- features_list[fold] # Get train features
    train <- features_list[-fold] # Get test features
    
    # Evaluate predictions
    preds <- KR_predict_set(train, test, h)
    
    # Apply metric function to 
    int_scores <- int_scores + evaluate_ranges(preds, metric_fun, degrees)$score
  }
  
  # Set up a report based on last fold
  report <-  evaluate_ranges(preds, metric_fun, degrees)
  
  int_scores <- int_scores / k # Adjust scores
  report$score <- int_scores   # Update report
  
  return(report)
}


bagging_cv <- function(features_list, 
                       grid_vals = seq(0, 10, length.out = 50)[-1],
                       N = 20) {
  
  # Optimal sub-sampling size 
  n <- length(features_list)  # Get the sample size
  r <- 20                     # Naive suboptimal subsampling size
  
  # Initialize a vector to store the sub-optimal values
  baggable <- rep(NA, N)
  
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = N,      # Maximum value of the progress bar
                       style = 3,    # Progress bar style [1, 2, 3]
                       width = 50,   # Progress bar width
                       char = "=")   # Character used to create the bar
  # Main loop 
  for (i in 1:N) {
    
    # Random sub-sampling
    idxs <- sample(1:n, r)
    
    # Initialize a vector to store grid-search scores
    scores <- rep(NA, length(grid_vals))

    # Inner LOOCV loop
    for (j in 1:length(grid_vals)) {
      
      LOOCV <- loo_cv(features_list[idxs], h=grid_vals[j])
      scores[j] <- LOOCV$L2
    }
    
    # Store the best bandwidth
    baggable[i] <- grid_vals[which.min(scores)]
    
    
    setTxtProgressBar(pb, i) # Update progress bar
  }

  # Montecarlize
  return(baggable)
}

# Conformal prediction  ---------------------------------------------------

train_test_conformal_split <- function(dataset, ratio = c(0.8,0.1,0.1)) {
  
  n <- length(dataset)
  idxs <- 1:n
  train_idx <- sample(idxs, ratio[1]*n)
  
  train  <- dataset[train_idx]
  
  test_conf_idxs <- idxs[-train_idx]
  test_idx <- sample(test_conf_idxs, ratio[2]*n)
  
  test <- dataset[test_idx]
  conf <- dataset[-c(train_idx,test_idx)]
  return(list("Train" = train, "Test"= test, "Conf"=conf))
}

Mean_abs_diff <- function(v1,v2){return(mean(abs(v1-v2)))}

calibration_quantile <- function(dataset, ratio=c(0.8,0.1,0.1),h,alpha=0.05) {
  
  # Splitting the data
  c(train, test, conf) := train_test_conformal_split(dataset,ratio)
  
  # Gathering prediction
  preds <- KR_predict_set(train,conf,h)
  scores <- Mean_abs_diff(preds[,1], preds[,2])
  
  # Retrieving an empirical quantile of the scores
  p <- ceiling(length(conf)+1)*(1-alpha)/length(conf)
  return(quantile(scores,p))
}

conformal_prediction <- function(dataset, h, quantile) {
  datamatrix <- matrix(NA, nrow=length(dataset), ncol=3)
  for (i in 1:length(dataset)) {
    pred <- KR_predict_audio(dataset[-i],dataset[i], h)
    score <- Mean_abs_diff(pred[,1], pred[,2])
    datamatrix[i,1] <- pred[,2] - (score+quantile)
    datamatrix[i,2] <- pred[,2]
    datamatrix[i,3] <- pred[,2] + (score+quantile)
  }
  return(datamatrix)
}
