# Libraries and tools -----------------------------------------------------

# Libraries
library(readr)
library(tuneR)
library(fda)

# Operator to assign multiple variables in the same line of code
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

# Data manipulation and feature extraction --------------------------------

# Define audio object from a .wav file
get_audio_object <- function(file_audio_path){
  
  # actual recording
  audio <- readWave(file_audio_path)
  
  # sample size
  n <- length(audio@left)

  # label
  num_part <- unlist(strsplit(file_audio_path, "_"))[4]
  label <- as.numeric(substring(num_part, 1, nchar(num_part) - 4))
  
  # author
  author <- unlist(strsplit(file_audio_path, "_"))[3]
  
  # get frequencies and halve it (since it is real valued) (MUST CHECK IF RIGHT)
  sampling_rate <- audio@samp.rate
  frequencies <- c(0:(n-1)) * (sampling_rate / n)
  frequencies <- frequencies[1:(n %/% 2)]
  
  # get power spectrum and halve it (since it is real valued)
  my_ps <- Mod(fft(c(audio@left))^2) / n
  my_ps <- my_ps[1:(n %/% 2)]
  
  
  return(list("Author" = author,
              "Recording" = audio@left,
              "Label" = label,
              "Frequencies" = frequencies,
              "Power_spectrum" = my_ps,
              "Sample_size" = n,
              "Sampling_rate" = sampling_rate))
}

# Get audio objects from a folder
get_audio_objects <- function(directory_path){
  
  files <- list.files(path = directory_path)
  n_iter <- length(files)
  
  # setup progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  i = 1 
  
  result <- list()
  for(file_path in files){
    result <- append(result, list(get_audio_object(paste(directory_path,file_path,sep="\\"))))
    setTxtProgressBar(pb, i)
    i = i + 1 
  }
  
  close(pb) 
  return(result)
}

# Thresholding algorithm, detects peaks and compue moving mean given a signal
ThresholdingAlgo <- function(y,lag,threshold,influence) {
  signals <- rep(0,length(y))
  filteredY <- y[0:lag]
  avgFilter <- NULL
  stdFilter <- NULL
  avgFilter[lag] <- mean(y[0:lag], na.rm=TRUE)
  stdFilter[lag] <- sd(y[0:lag], na.rm=TRUE)
  for (i in (lag+1):length(y)){
    if (abs(y[i]-avgFilter[i-1]) > threshold*stdFilter[i-1]) {
      if (y[i] > avgFilter[i-1]) {
        signals[i] <- 1;
      } else {
        signals[i] <- -1;
      }
      filteredY[i] <- influence*y[i]+(1-influence)*filteredY[i-1]
    } else {
      signals[i] <- 0
      filteredY[i] <- y[i]
    }
    avgFilter[i] <- mean(filteredY[(i-lag):i], na.rm=TRUE)
    stdFilter[i] <- sd(filteredY[(i-lag):i], na.rm=TRUE)
  }
  return(list("signals"=signals,"avgFilter"=avgFilter,"stdFilter"=stdFilter))
}

# Divide into test and train sets
train_test_split <- function(tot_set,test_ratio){
  n <- length(tot_set)
  train_ratio <- 1 - test_ratio 
  sample <- sample(c(TRUE, FALSE), n, replace = TRUE, prob = c(train_ratio,test_ratio))
  train  <- tot_set[sample]
  test   <- tot_set[!sample]
  return(list("Train" = train, "Test"= test))
}

# Computes the spectogram of a audio object
spectro = function(data, nfft=1024, window=256, overlap=128, t0=0, 
                   plot_spec = T, normalize = T, return_data = F){
  
  library(signal)
  library(oce)
  
  # extract signal
  snd = data$Recording
  
  # demean to remove DC offset
  snd = snd-mean(snd)
  
  # determine duration
  dur = length(snd)/data$Sampling_rate
  
  # create spectrogram
  spec = specgram(x = snd,
                  n = nfft,
                  Fs = data$Sampling_rate,
                  window = window,
                  overlap = overlap
  )
  
  # discard phase info
  P = abs(spec$S)
  
  # normalize
  if(normalize){
    P = P/max(P)  
  }
  
  # convert to dB
  P = 10*log10(P)
  
  # config time axis
  if(t0==0){
    t = as.numeric(spec$t)
  } else {
    t = as.POSIXct(spec$t, origin = t0)
  }
  
  # rename freq
  f = spec$f
  
  if(plot_spec){
    
    # change plot colour defaults
    par(bg = "black")
    par(col.lab="white")
    par(col.axis="white")
    par(col.main="white")
    
    # plot spectrogram
    imagep(t,f, t(P), col = oce.colorsViridis, drawPalette = T,
           ylab = 'Frequency [Hz]', axes = F)
    
    box(col = 'white')
    axis(2, labels = T, col = 'white')
    title(main = data$Label)
    # add x axis
    if(t0==0){
      
      axis(1, labels = T, col = 'white')
      
    }else{
      
      axis.POSIXct(seq.POSIXt(t0, t0+dur, 10), side = 1, format = '%H:%M:%S', col = 'white', las = 1)
      mtext(paste0(format(t0, '%B %d, %Y')), side = 1, adj = 0, line = 2, col = 'white')
      
    }
  }
  
  if(return_data){
    
    # prep output
    spec = list(
      t = t,
      f = f,
      p = t(P)
    )
    
    return(spec)  
  }
}


# Distances ---------------------------------------------------------------

L2 <- function(f1, f2){return(sqrt(sum((f1 - f2)^2)))}
Kernel <- function(x) {
  return(0.5*exp(-0.5*x^2))
}

# Plots -------------------------------------------------------------------

# Plot recording from audio object
simple_plot <- function(audio_object){
  
  time <- (1 : audio_object$Sample_size)
  audio <- audio_object$Recording
  
  # plot 
  plot(time,audio,"l",ylab = "Audio",xlab = "Time",main = paste(audio_object$Label, "°", sep = ""))
}

# Plot power spectrum from audio object 
power_plot <- function(audio_object){
  plot(audio_object$Frequencies,
       audio_object$Power_spectrum,
       "l",ylab = "Power",xlab = "Frequency",
       main = paste(audio_object$Label, "°", sep = ""))
}

# Plot histogram of labels given a list of recordings objects
labels_hist <- function(recordings){
  labels <- list()
  for(recording in recordings){
    labels <- append(labels,recording$Label)
  }
  labels <- unlist(labels)
  hist(labels)
  
}


# Functional features extraction and prediction ---------------------------

# Get the smooth approximation of a the power spectrum of an audio object defined a basis
get_smooth_from_ps <- function(basis_func, n_basis, audio_object){
  
  my_basis <- basis_func(range(0,(audio_object$Sampling_rate/2)),nbasis = n_basis)
  freq <- audio_object$Frequencies
  power <- audio_object$Power_spectrum
  appr <- smooth.basis(argvals = freq, y = power ,fdParobj = my_basis)$fd
  appr <- eval.fd(0:(audio_object$Sampling_rate/2),appr)
  appr <- appr / sum(appr)
  
  return(list("PS_smooth_func" = appr,
              "Label" = audio_object$Label))
}

# Get smooth approximations from a list of audio objects 
get_smooth_list_from_ps <- function(basis_func, n_basis, list){
  
  n_iter <- length(list)
  
  # setup progress bar
  pb <- txtProgressBar(min = 0,      # Minimum value of the progress bar
                       max = n_iter, # Maximum value of the progress bar
                       style = 3,    # Progress bar style (also available style = 1 and style = 2)
                       width = 50,   # Progress bar width. Defaults to getOption("width")
                       char = "=")   # Character used to create the bar
  i = 1 
  
  result <- list()
  for(audio_object in list){
    result <- append(result, list(get_smooth_from_ps(basis_func, n_basis, audio_object)))
    setTxtProgressBar(pb, i)
    i = i + 1 
  }
  close(pb) 
  return(result)
}


# KNN prediction ----------------------------------------------------------

# Predict a single audio 
KNN_predict_audio <- function(train_set, new_data, k){
  
  weight_label_dataframe <- data.frame()
  
  for(train_audio in train_set){
    weight <-  Kernel( L2(train_audio$PS_smooth_func,new_data$PS_smooth_func) ) 
    weight_label_dataframe <- rbind(weight_label_dataframe,c(weight,train_audio$Label))
  }
  
  names(weight_label_dataframe) <- c("Weight","Label")
  weight_label_dataframe <- weight_label_dataframe[order(weight_label_dataframe$Weight,decreasing = T),]
  
  weight_label_dataframe <- weight_label_dataframe[1:k,]
  
  tot_weight <- sum(weight_label_dataframe$Weight)
  
  prediction <- sum(weight_label_dataframe$Weight * weight_label_dataframe$Label) / tot_weight
  
  return(prediction)
}

# Predict a set of audios
KNN_predict_set <- function(train_set,test_set,k){
  preds <- data.frame()
  for(audio in test_set){
    preds <- rbind(preds,
                   c(audio$Label,KNN_predict_audio(train_set, audio, k)))
  }
  names(preds) <- c("Label", "Prediction")
  preds <- preds[order(preds$Label),]
}

# All NN prediction  ------------------------------------------------------

# Predict a single audio 
All_NN_predict_audio <- function(train_set, new_data){
  
  weight_label_dataframe <- data.frame()
  
  for(train_audio in train_set){
    weight <-  Kernel( L2(train_audio$PS_smooth_func,new_data$PS_smooth_func) )
    weight_label_dataframe <- rbind(weight_label_dataframe,c(weight,train_audio$Label))
  }
  
  names(weight_label_dataframe) <- c("Weight","Label")

  tot_weight <- sum(weight_label_dataframe$Weight)
  
  prediction <- sum(weight_label_dataframe$Weight * weight_label_dataframe$Label) / tot_weight
  
  return(prediction)
}

# Predict a set of audios
All_NN_predict_set <- function(train_set,test_set){
  preds <- data.frame()
  for(audio in test_set){
    preds <- rbind(preds,
                   c(audio$Label,All_NN_predict_audio(train_set, audio)))
  }
  names(preds) <- c("Label", "Prediction")
  preds <- preds[order(preds$Label),]
}

