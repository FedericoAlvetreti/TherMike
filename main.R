# Library -----------------------------------------------------------------
source('utils.R')


# Loading and merging data ------------------------------------------------

our_recordings <- get_audio_objects("SL-splash-recordings")
paper_recordings <- get_audio_objects("paper-processed-recordings")
tot_recordings <- c(our_recordings,paper_recordings)

smooth_funcs <- get_smooth_list_from_ps(create.bspline.basis, n_basis = 200, our_recordings)

# Predictions -------------------------------------------------------------

c(train,test) := train_test_split(smooth_funcs,0.1)

KNN_preds <- KNN_predict_set(train,test,5)
All_NN_preds <- All_NN_predict_set(train,test)

All_NN_predict_audio(train,test[[1]])
