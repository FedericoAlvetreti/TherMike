# Library -----------------------------------------------------------------
source('utils.R')

# Loading data ------------------------------------------------------------
my_recordings <- get_audio_objects("SL-splash-recordings")

# An example
pw_plot(my_recordings[[84]])
