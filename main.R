# Library -----------------------------------------------------------------
source('utils.R')


# Loading and merging data ------------------------------------------------

our_recordings <- get_audio_objects("processed-recs")
paper_recordings <- get_audio_objects("processed-paper")
tot_recordings <- c(our_recordings,paper_recordings)


# Kind of a feature extraction --------------------------------------------

tot_smooth_funcs <- get_smooth_list_from_ps(create.bspline.basis, n_basis = 100, tot_recordings)

# Split in test and train -------------------------------------------------

c(train,test) := train_test_split(tot_smooth_funcs,0.1)


# Predictions -------------------------------------------------------------

preds <- predict_set(train,test)