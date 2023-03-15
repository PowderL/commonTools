gen_raster_define <- function(value_cells, grid.extent = c(102.5, 105, 30, 31.5), grid.res=0.01) {
  # degree; west, east, south, north; longitude, latitude
  value_raster <- raster(extent(grid.extent), resolution=grid.res, crs="+proj=longlat +datum=WGS84")
  
  row_min <- ceiling((90 - grid.extent[4] + grid.res/2) / grid.res)  #  5908
  col_min <- ceiling((grid.extent[1] + 180 + grid.res/2) / grid.res)
  
  n.rows <- nrow(value_raster)
  n.cols <- ncol(value_raster)
  r.values <- rep(NA, n.rows*n.cols)
  i.loc <- (value_cells[, 1]-row_min)*n.cols + (value_cells[,2]-col_min) + 1
  r.values[i.loc] <- value_cells[, 3]
  values(value_raster) <- r.values
  return (value_raster)
}

gen_grid_data <- function(ext, resolution){
  lon <- seq(ext[1] + resolution/2, ext[2] - resolution/2, by = resolution)
  lat <- seq(ext[3] + resolution/2, ext[4] - resolution/2, by = resolution)
  grid <- expand.grid(lon, lat)
  names(grid) <- c("lon", "lat")
  grid$col001 <- ceiling((180 + grid$lon)/resolution)
  grid$row001 <- ceiling((90 - grid$lat)/resolution)
  return (as.data.table(grid))
}

gen_grids_in_bounds <- function(map, resolution) {
  ext <- raster::extent(map)
  ext_outbound <- c(floor(ext[1]), ceiling(ext[2]), floor(ext[3]), ceiling(ext[4]))
  r <- raster(extent(ext_outbound), resolution=resolution, crs="+proj=longlat +datum=WGS84")
  values(r) <- 1
  mask_r <- mask(map, r)
  grids.all <- gen_grid_data(ext_outbound, resolution)
  setorder(grids.all, row001, col001)
  grids_in_bouds <- grids.all[!is.na(values(mask_r)), ]
  return(grids_in_bouds)
}

gen_maps <- function(r.list, map, unit = "a", magnitude = 10){
  if (is.list(r.list)){
    r.season <- stack(r.list)
    names.attr <- names(r.list)
  }else{
    r.season <- r.list
  }
  val <- values(r.season)
  values_quantile <- quantile(val, probs = c(0.0005, 0.9995), na.rm = T)
  low_limit <- values_quantile[1]
  high_limit <- values_quantile[2]
  ## to remove white grid
  min_value <- (low_limit %/% magnitude + 1) * magnitude
  
  max_value <- high_limit %/% magnitude * magnitude
  
  breaks <- seq(low_limit, high_limit, 0.01*magnitude)
  labels <- seq(min_value, max_value, magnitude)
  if (length(labels) > 10) labels <- seq(min_value, max_value, 5*magnitude)
  if (length(labels) > 10) labels <- seq(min_value, max_value, 10*magnitude)
  breaks_label <- labels
  labels_label <- labels
  labels_1 <- as.character(labels_label[1])
  if (unit == 'a'){
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mu*g*"/"*m^3)))
  }else if (unit == 'b') {
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mg*"/"*m^3)))
  }else if(unit == 'c'){
    labels_1_quote <- paste0(labels_1, " (ppm)")
  }else if (unit == "d"){
    labels_1_quote <- paste0(labels_1, " (ppb)")
  }else{
    labels_1_quote <- paste0(labels_1)
  }
  
  labels_label <- c(labels_1_quote, labels_label[2:length(labels_label)])
  cols <- matlab.like(length(breaks))
  low_value <- min(breaks)
  high_value <- max(breaks)
  r.season[r.season < low_value] <- low_value
  r.season[r.season > high_value] <- high_value
  tmp <- spplot(r.season, names.attr= names.attr, 
                col.regions=cols, at=breaks, maxpixels=500000,
                colorkey=list(labels=list(labels= labels_label, at=breaks_label), space="right", width=1.2),
                panel=function(...) {
                  panel.gridplot(...) 
                  sp.polygons(map, lwd=0.4, col="black")
                }
  )
  return (tmp)
}
gen_maps_bars <- function(r.list, map, unit = "a", breaks = seq(-1.5-0.01, 1+0.01, 0.001), labels = seq(-1.5, 1, 0.5)){
  if (is.list(r.list)){
    r.season <- stack(r.list)
    names.attr <- names(r.list)
  }else{
    r.season <- r.list
  }
  breaks_label <- labels
  labels_label <- labels
  labels_1 <- as.character(labels_label[1])
  if (unit == 'a'){
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mu*g*"/"*m^3)))
  }else if (unit == 'b') {
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mg*"/"*m^3)))
  }else if(unit == 'c'){
    labels_1_quote <- paste0(labels_1, " (ppm)")
  }else if (unit == "d"){
    labels_1_quote <- as.expression(bquote(.(labels_1)(molec*"/"*m^2)))
  }else if (unit == "e"){
    labels_1_quote <- as.expression(bquote(atop(NA, atop(textstyle(paste(.(labels_1), "      ")), textstyle((mu*"mol/"*m^2))))))
  }else if (unit == 'f') {
    labels_1_quote <- as.expression(bquote(.(labels_1), (mu*"mol/"*m^2)))
  }else if (unit == "p") {
    labels_1_quote <- paste0(labels_1, " (%)")
  }
    
  cut.val <- 0
  theme.novpadding <-
    list(layout.heights =
           list(top.padding = cut.val,
                main.key.padding = 0.5,
                key.axis.padding = cut.val,
                axis.xlab.padding = cut.val,
                xlab.key.padding = cut.val,
                key.sub.padding = cut.val,
                bottom.padding = 1),
         layout.widths =
           list(left.padding = 0.5,
                key.ylab.padding = cut.val,
                ylab.axis.padding = cut.val,
                axis.key.padding = cut.val,
                right.padding = 0.5))  
  labels_label <- c(labels_1_quote, labels_label[2:length(labels_label)])
  cols <- matlab.like(length(breaks))
  low_value <- min(breaks)
  high_value <- max(breaks)
  r.season[r.season < low_value] <- low_value
  r.season[r.season > high_value] <- high_value
  tmp <- spplot(r.season, names.attr= names.attr, 
                col.regions=cols, at=breaks, maxpixels=500000,
                colorkey=list(labels=list(labels= labels_label, at=breaks_label), space="right", width=1.2),
                panel=function(...) {
                  panel.gridplot(...) 
                  sp.polygons(map, lwd=0.4, col="black")
                },
                par.settings=theme.novpadding
  )
  return (tmp)
}

gen_maps_bars_list <- function(r.list, map, unit = "a", breaks = seq(-1.5-0.01, 1+0.01, 0.001), labels = seq(-1.5, 1, 0.5)){
  breaks_label <- labels
  labels_label <- labels
  labels_1 <- as.character(labels_label[1])
  if (unit == 'a'){
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mu*g*"/"*m^3)))
  }else if (unit == 'b') {
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (mg*"/"*m^3)))
  }else if(unit == 'c'){
    labels_1_quote <- paste0(labels_1, "\n(ppm)")
  }else if (unit == "d"){
    labels_1_quote <- as.expression(bquote(.(labels_1) ~ (molec*"/"*m^2)))
  }else{
    labels_1_quote <- paste0(labels_1)
  }
  cut.val <- 0
  theme.novpadding <-
    list(layout.heights =
           list(top.padding = cut.val,
                main.key.padding = 0.5,
                key.axis.padding = cut.val,
                axis.xlab.padding = cut.val,
                xlab.key.padding = cut.val,
                key.sub.padding = cut.val,
                bottom.padding = cut.val),
         layout.widths =
           list(left.padding = 0.5,
                key.ylab.padding = cut.val,
                ylab.axis.padding = cut.val,
                axis.key.padding = cut.val,
                right.padding = 0.5))
  
  labels_label <- c(labels_1_quote, labels_label[2:length(labels_label)])
  cols <- matlab.like(length(breaks))
  low_value <- min(breaks)
  high_value <- max(breaks)
  name.attrs <- names(r.list)
  map.list <- foreach (i = 1:length(r.list))%do%{
    r.season <- r.list[[i]]
    r.season[r.season < low_value] <- low_value
    r.season[r.season > high_value] <- high_value
    tmp <- spplot(r.season, main=name.attrs[i], 
                  col.regions=cols, at=breaks, maxpixels=500000,
                  colorkey=list(labels=list(labels= labels_label, at=breaks_label), space="right", width=1.2),
                  panel=function(...) {
                    panel.gridplot(...) 
                    sp.polygons(map, lwd=0.4, col="black")},
                  par.settings=theme.novpadding
    )
    return (tmp)
  }
  return (map.list)
}