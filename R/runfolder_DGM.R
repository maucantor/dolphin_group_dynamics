# dolphin_group_metrics.R
# title: "Dolphin-group-dynamics-metrics"
# author: "Mauricio Cantor and João Valle-Pereira"
# date: "2025-12-04, 05"
# 
# compute metrics for all csv files in a folder and prepare database

# for Julia, the selected variables that seem to make sense are these 5:
# MPD and MST for spatial cohesion
# MRL and angular velocity for heading coordination
# Breath synchrony

# 
# ---------------------------
# Setting up
# ---------------------------


# 1. Loading all functions from the other file
source(paste0(getwd(), '/R/functions_DGM_julia.R'))
 

# 2. Read all csv files in a folder
# 
# Pick a folder a working directory in the cetacean server
# you can run the code below for all videos in a site and given flight altitude:
# e.g.: "/Volumes/LABIRINTO/Data_processing/JuliaPierry/Cananeia/output/20m"
# or all videos in a site:
# e.g.,: "/Volumes/LABIRINTO/Data_processing/JuliaPierry/Cananeia/output/"

folder_path <- "/Volumes/LABIRINTO/Data_processing/JuliaPierry/Cananeia/output/20m"

# pick all the .csv files in that folder (including subfolders)
file_list <- list.files(path = folder_path,
                        pattern = "\\.csv$",
                        recursive = TRUE,     # Search within subdirectories
                        full.names = TRUE)

# these are the files (e.g., should be 10 for Cananeia/20m)
file_list


# ---------------------------
# Running
# ---------------------------

# 1. Run master compute_group_metrics() through all files
# This function in R/functions_DGM_julia.R only calculates those 5 metrics above

result_list = list()

# Loop across files, saving only the dataframe with the aggregated metrics per entire videos
for(i in 1:length(file_list)){
 aux1 = compute_group_metrics(csv_path = file_list[[i]],
                                           whole_video = TRUE,
                                           fps = 30, 
                                           sliding = TRUE, 
                                           window_sec = NA,  
                                           slide_step_frames = NULL,
                                           interp_gap = 0, 
                                           Nmin = 1,
                                           angle_col = "MovingAvgAngle_deg",
                                           breath_col = NULL)
 # add video name, site name, flight altitude and save the video-summary row
 aux2 = cbind(Video = gsub("_output.csv", "", basename(file_list[i])),
              Site = extract_path_info_regex(file_list[i])[1],
              FlightAltitude = extract_path_info_regex(file_list[i])[2],
              aux1$window_metrics)
 # save video summary in a list
 result_list[[i]] = aux2
}

# Combine all elements of the list into a dataframe
result_full_videos = as.data.frame(dplyr::bind_rows(result_list))

result_full_videos
view(result_full_videos)


