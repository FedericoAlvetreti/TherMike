
# Libraries
library(readr)
library(tuneR)


# operator to assign multiple variables in the same line of code
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


# Extract file audios with label 
get_audio_object <- function(file_audio_path){
  
  # actual recording
  audio <- readWave(file_audio_path)
  
  # sample size
  n <- length(audio@left)
  
  # label
  num_part <- unlist(strsplit(file_audio_path, "_"))[4]
  label <- as.numeric(substring(num_part, 1, nchar(num_part) - 4))
  
  # frequencies (MUST CHECK IF RIGHT)
  sampling_rate <- audio@samp.rate
  frequencies <- c(0:(n-1)) * (sampling_rate / n)
  
  return(list("Recording" = audio@left,
              "Label" = label,
              "Frequencies" = frequencies,
              "Sample_size" = n))
}
get_audio_objects <- function(directory_path){
  files <- list.files(path = directory_path)
  result <- list()
  for(file_path in files){
    result <- append(result, list(get_audio_object(paste(directory_path,file_path,sep="\\"))))
  }
  return(result)
}



# Plot power spectrum from audio object defined above
# ( the way I got power and the frequencies must be absolutely checked ! ! ! )
pw_plot <- function(audio_object){
  
  # get power spectrum and halve it (since it is real valued)
  my_ps <- Mod(fft(c(audio_object$Recording))^2) / audio_object$Sample_size
  my_ps <- my_ps[1:(audio_object$Sample_size %/% 2)]
  
  # get frequencies and halve it (since it is real valued)
  f <- audio_object$Frequencies
  f <- f[1:(audio_object$Sample_size %/% 2)]
  
  # plot 
  plot(f,my_ps,"l",ylab = "Power",xlab = "Frequency",main = paste(audio_object$Label, "°", sep = ""))
}
