# dolphin_group_metrics.R
# title: "Dolphin-group-dynamics-metrics"
# author: "Mauricio Cantor and João Valle-Pereira"
# date: "2025-12-04"
# 
# groundtruthing group metrics: closer look at 1 video and the csv outputs
# 

# loading all functions
source(paste0(getwd(), '/R/functions_DGM.R'))


# Example: the first 50 seconds of this video from Cananeia at low altitude (20m), few individuals (3), with good detection rate and 2 individuals in close proximity most of the time, and another that branches off; there are some good turning angles too"
# Video file: /Volumes/LABIRINTO/Data_processing/JuliaPierry/Cananeia/output/20m/E01_CN01_03ago23_(2)_corte_126/E01_CN01_03ago23_(2)_corte_126_predictions.mp4
# frame range for 50 secs: 000001 to 001495


csv_file = "/Volumes/LABIRINTO/Data_processing/JuliaPierry/Cananeia/output/20m/E01_CN01_03ago23_(2)_corte_126/E01_CN01_03ago23_(2)_corte_126_output.csv"


test = compute_group_metrics(csv_path = csv_file,
                             # should output summary for the whole video? (that is, 1 single sliding window)
                             whole_video = FALSE,
                             # If not, then define these parameters for sliding window analyses:
                             fps = 30, 
                             # whether to summarize metrics per sliding window
                             sliding = TRUE, 
                             # length, in secs, of sliding window of analyses; ==NA, if whole_video=TRUE or
                             window_sec = 10,  
                             # how many frames to jump between windows
                             slide_step_frames = NULL,
                             # interpolation of detections per frame to smooth out misdetections; ==0 to use csv as is
                             interp_gap = 0, 
                             # minimum number of detections to calculate metrics (Nmin=1 to include all individuals)
                             Nmin = 1,
                             # for angle analyses, give the column names if present
                             angle_col = "MovingAvgAngle_deg",
                             # for breath sync analyses, give column name it present 
                             # (normally no pre-defined breath variable, so this will be computed from detections)
                             breath_col = NULL, 
                             # should summary be plotted? (doesnt work all the time)
                             return_plots = FALSE)

str(test)

#data input
view(test[[1]])

# metrics per frame
str(test[[2]])
as.data.frame(test$frame_metrics)
view(test$frame_metrics)

# summary metrics per X-sec window
str(test[[3]])
as.data.frame(test$window_metrics)
view(test$window_metrics)

