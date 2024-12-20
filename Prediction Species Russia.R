####################Predict Tree Species########################
#####################by Jacob Schladebach#######################
##########################18.10.2024############################



library(lidR) 
library(dplyr)
library(sf) 
library(terra) 
library(alphashape3d)
require(geosphere)
require(tidyr)
require(purrr)
require(caret)

#########functions
#######the following function is according to Liam Irwin (https://liamirwin.github.io/CompTreeR/)
get_crown_attributes <- function(X, Y, Z, export_ashape = FALSE) {
  if (length(X) <= 3 || length(Y) <= 3 || length(Z) <= 3) {
    print('Cannot compute a 3D hull from 3 or fewer points')
    return(0)
  }
  
  # Create a matrix of 3D points
  a3d <- cbind(X, Y, Z)
  
  # Get treetop location
  top_x <- a3d[which.max(a3d[,3]),1]
  top_y <- a3d[which.max(a3d[,3]),2]
  
  # Normalize X and Y
  a3d[,1] <- a3d[,1] - mean(a3d[,1])
  a3d[,2] <- a3d[,2] - mean(a3d[,2])
  
  alpha <- 1
  
  ashape <- alphashape3d::ashape3d(x = a3d, alpha = alpha, pert = TRUE)
  
  df <- data.frame(
    Zmax = max(Z),
    Zq999 = as.numeric(quantile(Z, 0.999)),
    Zq99 = as.numeric(quantile(Z, 0.990)),
    Z_mean = mean(Z),
    n_points = length(Z),
    vol_concave = alphashape3d::volume_ashape3d(ashape, indexAlpha = 1),
    CV_Z = sd(Z) / mean(Z),
    CRR = (mean(Z) - min(Z)) / (max(Z) - min(Z)),
    X = top_x,
    Y = top_y
  )
  
  if (export_ashape) {
    return(list(metrics = df, ashape = ashape, a3d = a3d))
  } else {
    return(df)
  }
}
triangle_area <- function(x1, z1, x2, z2, x3, z3) {
  return(abs((x1 * (z2 - z3) + x2 * (z3 - z1) + x3 * (z1 - z2)) / 2))
}




#####loading data
rf_model<-readRDS("")
out<-""

test<-list.files(pattern=".laz")

for(n in 143:length(test)) {
  result <- tryCatch({  
    las<-readLAS(test[n])
    if (is.null(las$R)){next}
    predict_data<-data.frame(Tree = NULL, LAT= NULL, LONG= NULL, Zq999 = NULL, Zq99=NULL, vol_concave= NULL, CV_Z=NULL, CRR=NULL, density_concave =NULL,  vertical_variability=NULL, pointedness =NULL, top_angle =NULL,Z_widest_distance =NULL, relative_height_widest=NULL, widest_distance = NULL, Probability=NULL)    
    points <- las@data
    unique_values_count <- unique(points$Tree)
    for (i in 1:length(unique_values_count)) {
      tree<-subset(las, las$Tree == unique_values_count[i] )
      tree@data$Z<-tree@data$Z-min(tree@data$Z)
      if(unique_values_count[i]!=1 & unique_values_count[i]!=2 & tree@header$`Number of point records`>50 & max(tree@data$Z)>1 ){ 
        ashape <- get_crown_attributes(tree@data$X,
                                       tree@data$Y,
                                       tree@data$Z,
                                       export_ashape = TRUE)
        
        ############
        ashape_df<-ashape[[2]]
        # plot(ashape_df, indexAlpha = 1, transparency = 0.3, axes = TRUE)
        
        df<-as.data.frame(ashape_df$x)

        ####locate widest point of crown
        slizes<-length(tree$X)/40
        segment_size <- (max(df$Z)-min(df$Z))/slizes
        breaks <- seq(min(df$Z), max(df$Z), by = segment_size)
        df$z_segment <- cut(df$Z, breaks = breaks, include.lowest = TRUE)
        
        widest_segments <- df %>%
          group_by(Y, z_segment) %>%
          summarise(
            min_x = min(X),
            max_x = max(X),
            distance = max_x - min_x,
            min_height = min(Z)  # Assuming 'z' is the height variable
          ) %>%
          filter(distance > 0)
        
        widest_segments2 <- df %>%
          group_by(X, z_segment) %>%
          summarise(
            min_y = min(Z),
            max_y = max(Z),
            distance = max_y - min_y,
            min_height = min(Z)  # Assuming 'z' is the height variable
          ) %>%
          filter(distance > 0)
        if (nrow(widest_segments)>0){        
          widest<-max(widest_segments$distance)
          widest2<-max(widest_segments2$distance)
          
          if (widest>widest2) {
            widest_coords<-subset(widest_segments, widest_segments$distance==widest)
          } else {
            widest_coords<-subset(widest_segments2, widest_segments2$distance==widest2)
          }
          top<-subset(df, df$Z == max(df$Z))
          B1<-data.frame(X=widest_coords[,3], Z = widest_coords$min_height)
          B2<-data.frame(X=widest_coords[,4], Z = widest_coords$min_height)
          vertices <- data.frame(
            x = c(B1$min_x, B2$max_x, top$X),  
            z = c(widest_coords$min_height, widest_coords$min_height, top$Z)   
          )

          #########construct triangle in crown
          A <- vertices[1, ]  
          B <- vertices[2, ]  
          C <- vertices[3, ]
          a <- sqrt((B$x - C$x)^2 + (B$z - C$z)^2)  
          b <- sqrt((A$x - C$x)^2 + (A$z - C$z)^2)  
          c <- sqrt((A$x - B$x)^2 + (A$z - B$z)^2) 
          right_dist <- sqrt((B$x - C$x)^2 + (B$z - C$z)^2)  
          left_dist <- sqrt((A$x - C$x)^2 + (A$z - C$z)^2)  
          widest_distance <- sqrt((A$x - B$x)^2 + (A$z - B$z)^2) 
          
          ####get angle at top of crown triangle
          top_angle<- acos((a^2 + b^2 - c^2) / (2 * a * b)) * (180 / pi)
          
          # P<-data.frame(x=df$X,z=df$Z)
          #  P<-subset(P, P$z>widest_coords$min_height)
          #  PL<-subset(P,P$x<top$X)
          # PR<-subset(P,P$x>top$X)
          
          # P$inside <- apply(P[, c("x", "z")], 1,
          # function(p) is_point_in_triangle(A, B, C, p))
          # inside_count <- sum(P$inside)
          # outside_count <- nrow(P) - inside_count
          
          #####get pointedness coefficient 
          height<-max(tree@data$Z)-min(tree@data$Z)
          average_height <- mean(tree@data$Z)
          pointedness<-(height-average_height)/height

          ####get ashape metrics and put df together
          result<-as.data.frame(ashape$metrics)
          crs<-st_crs(tree)
          sf<- st_as_sf(df, coords = c("X", "Y"), crs = crs)
          sf<-st_transform(sf,4326)
          coords<-st_coordinates(sf)
          result$LONG<-mean(coords[,1])
          result$LAT<-mean(coords[,2])
          result$pointedness<-pointedness
          result$top_angle<-top_angle
          result$widest_distance<-widest_distance
          widest_coords<-widest_coords[1,]
          result$Z_widest_distance<-widest_coords$min_height
          result$relative_height_widest<-result$Z_widest_distance/max(tree$Z)
          result$density_concave <- result$n_points/result$vol_concave
          result$vertical_variability <- sd(tree$Z)       
          
          predicted_labels <- predict(rf_model, newdata = result,type = "prob")
          predicted_classes <- apply(predicted_labels, 1, function(x) colnames(predicted_labels)[which.max(x)])
          code <- ifelse(predicted_classes == "LARIX", 3,
                         ifelse(predicted_classes == "PISY", 4,
                                ifelse(predicted_classes =="BETU",5)))
                                       
                                                                                        
          
          result$Tree<-unique_values_count[i]
          result$Species<- code
          result$Probability<-as.numeric(max(predicted_labels))
        }} else {
          if (unique_values_count[i]!=2){
            p<-1
          }else{
            p<-2
          }
          result <- unclassified <- data.frame(Tree= unique_values_count[i], Species = p, LAT= 0, LONG=0,   Zmax = 0,X=0, Y=0, Z_mean=0, n_points=0, Zq999 = 0 , Zq99=0,  vol_concave= 0, CV_Z=0, CRR=0, density_concave =0,  vertical_variability=0, pointedness =0, top_angle =0,Z_widest_distance =0, relative_height_widest=0, widest_distance = 0, Probability=0)
        }
      predict_data<-rbind(predict_data,result)
    }
    predict_data<-predict_data[order(predict_data$Tree),]
    labels<-c(predict_data$Species)
    prob<-c(predict_data$Probability)
    las_df<-as.data.frame(las@data)
    unique_trees <- sort(unique(las_df$Tree))
    aggregated_data <- las_df %>%
      group_by(Tree) %>%
      summarise(mean_height = mean(Z, na.rm = TRUE), .groups = 'drop')
    aggregated_data$Species<-labels
    aggregated_data$Probability<-prob
    
    
    las@data <- merge(las@data, aggregated_data[, c("Tree", "Species")], by = "Tree", all.x = TRUE)
    las@data <- merge(las@data, aggregated_data[, c("Tree", "Probability")], by = "Tree", all.x = TRUE)
    las <- add_lasattribute(las, las@data$Species, "Species", "Predicted Species")
    las <- add_lasattribute(las, las@data$Probability, "Probability", "Prediction Probability")
    name<-basename(test[n])
    name<-gsub("_segmented","_predicted",name)
    name2<-basename(test[n])
    name2<-gsub("_segmented.laz","_predicted_meta.csv",name2)
    write.csv2(predict_data,paste0(out,name2))
    writeLAS(las,paste0(out,name))
  }, error = function(e) {
    cat("Error on iteration", n, ": ", conditionMessage(e), "/n")
    next 
  })}

