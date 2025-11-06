# Surface proportion of each habitat within buffer --> we need to apply transformation to avoid compositional data bias (see here for details : https://docs.google.com/document/d/1cN9vJ6I4fHzPXZfXjOm77Klk5hCSFBBkOj6Mhgum_S8/edit?tab=t.0 and https://medium.com/@nextgendatascientist/a-guide-for-data-scientists-log-ratio-transformations-in-machine-learning-a2db44e2a455) 

###### Distance between seabed and depth sampling ############

# Check negative dist_seabed_depthsampling
dsamp %>%
  filter(dist_seabed_depthsampling < 0) %>%
  pull (method) # dive_transect, dive_motionless and motionless_descended_from_surface
  
  # Set negative dist_seabed_depthsampling to 0