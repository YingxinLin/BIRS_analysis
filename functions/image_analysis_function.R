# image function
library(raster)
library(ggplot2)
library(ggthemes)

plotExpression <- function(raster_obj, sce_cellLabel, exprsMat,
                           gene, scale = TRUE) {
  
  
  ddf <- rasterToPoints(raster_obj)
  ddf <- data.frame(ddf)
  colnames(ddf) <- c("X", "Y")
  cell_label <- raster::values(r)
  
  notInLabel <- unique(cell_label[!cell_label %in% sce_cellLabel])
  
  # exprsMat <- assay(sce, exprs_value)
  
  new_values <- mapvalues((cell_label), 
                          from = notInLabel, 
                          to = rep(min(exprsMat[gene, ]),
                                   length(notInLabel)))
  
  
  new_values <- mapvalues((new_values), 
                          from = sce_cellLabel, 
                          to = exprsMat[gene, ])
  
  if (scale) {
    new_values <- zero_one_scale(new_values)
  }
  
  
  
  g <- ggplot(NULL) + 
    geom_raster(data = ddf, aes(X, Y, fill = new_values)) +
    theme_minimal() +
    scale_fill_viridis_c() +
    coord_quickmap() +
    theme(aspect.ratio = 1) +
    labs(fill = gene)
  g
  #return(g)
}


zero_one_scale <- function(mat) {
  (mat - min(mat))/(max(mat) - min(mat))
}


# A function to map the sce value to raster object
mapValueToCoord <- function(raster_obj, 
                            sce_cellLabelInImage, 
                            sce_newValue,
                            cont = FALSE) {
  
  cell_label <- raster::values(raster_obj)
  
  
  notInLabel <- unique(cell_label[!cell_label %in% sce_cellLabelInImage])
  
  if (cont) {
    new_values <- mapvalues((cell_label), 
                            from = notInLabel, 
                            to = rep(min(sce_newValue),length(notInLabel)))
  } else {
    new_values <- mapvalues((cell_label), 
                            from = notInLabel, 
                            to = rep("Background",length(notInLabel)))
  }
  
  
  new_values <- mapvalues(new_values, 
                          from = sce_cellLabelInImage, 
                          to = sce_newValue)
  
  return(new_values)
}


rearrange_string <- function(str) {
  unlist(lapply(strsplit(str, "_"), function(x) paste(sort(x), collapse = "_")))
}

L_stats <- function(ppp_obj = NULL, from = NULL, to = NULL, L_dist = NULL) {
  L <- spatstat::Lcross(ppp_obj, 
                        from = from,
                        to = to,
                        verbose = FALSE,
                        correction = "best")
  
  
  L_theo <- L$theo[L$r <= L_dist]
  L_iso <- L$iso[L$r <= L_dist]
  L_res <- mean(L_iso - L_theo)
  
  return(L_res)
}
