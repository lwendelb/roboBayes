##########################################################
## Purpose: Do cloudmasking for each point for each sensor,
##          combine to single list of datatables
## Input: Raw band values output from GEE
## Variable created: sat_list
## Files created: "data/myanmar/sat_data.Rdata" + "data/peru/sat_data.Rdata"
## Creator: Ian McGregor, imcgreg@ncsu.edu
## System: R Version 3.6.3, May 2020
##########################################################
library(data.table)

# cloudmasking functions ####
parseMOD09 <- function(x) {
  ## REMINDER - this function is based on the state_1km band from mod09ga. This is used
  ## for both mod09ga and mod09gq, as in GEE it says gq is meant to be used in conjunction
  ## with ga "where important quality information is stored".
  
  # cloud state
  cloud_state <- bitwAnd(x, 3) # 0=clear, 1=cloudy, 2=mixed, 3=not set, assumed clear
  # is it cloud shadow
  is_cloud_shadow <- ifelse(bitwAnd(bitwShiftR(x, 2), 1), TRUE, FALSE)
  # land/water flag
  lw_flag <- bitwAnd(bitwShiftR(x, 3), 7) # 0=shallow ocean, 1=land, 2=ocean coastlines and lake shorelines, 3=shallow inland water, 4=ephemeral water, 5=deep inland water, 6=continental/moderate ocean, 7=deep ocean
  # aerosol quantity
  aero_quan <- bitwAnd(bitwShiftR(x, 6), 3) # 0=climatology, 1=low, 2=average, 3=high
  # cirrus detected
  cirrus_detected <- bitwAnd(bitwShiftR(x, 8), 3) # 0=None, 1=small, 2=average, 3=high
  # internal cloud algorithm flag
  is_cloud <- ifelse(bitwAnd(bitwShiftR(x, 10), 1), TRUE, FALSE)
  # internal fire algorithm flag
  is_fire <- ifelse(bitwAnd(bitwShiftR(x, 11), 1), TRUE, FALSE)
  # MOD35 snow/ice flag
  is_mod35_snow <- ifelse(bitwAnd(bitwShiftR(x, 12), 1), TRUE, FALSE)
  # is it adjacent to cloud
  is_adj_cloud <- ifelse(bitwAnd(bitwShiftR(x, 13), 1), TRUE, FALSE)
  # is it brdf corrected
  is_brdf_corr <- ifelse(bitwAnd(bitwShiftR(x, 14), 1), TRUE, FALSE)
  # is it internal snow mask
  is_inter_snow <- ifelse(bitwAnd(bitwShiftR(x, 15), 1), TRUE, FALSE)
  
  return(list(
    cloud_state = cloud_state, cloud_shadow = is_cloud_shadow, lw_flag = lw_flag, aero_quan = aero_quan,
    cirrus_detected = cirrus_detected, cloud = is_cloud, fire = is_fire, mod35_snow = is_mod35_snow,
    adj_cloud = is_adj_cloud, brdf_corr = is_brdf_corr, inter_snow = is_inter_snow
  ))
}

# Functions for QA/QC filtering and calculating ndvi of matrix outputs
parseL8SR_pixel <- function(x){
  # QA_PIXEL: We only want to keep the best quality images, so only those that are
  ## clear and have low confidences
  
  # Binary
  ## Bit 0 - if pixel is fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  ## Bit 1 - if dilated cloud, then true
  dilatedCloud <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  ## Bit 2 - if cirrus, then true
  cirrus <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  ## Bit 3 - if cloud, then true
  cloud <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  ## Bit 4 - if cloud shadow, then true
  cloudShadow <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  ## Bit 5 - if snow, then true
  snow <- ifelse(bitwAnd(bitwShiftR(x,5),1), TRUE, FALSE)
  ## Bit 6 - if clear, then true
  clear <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  ## Bit 7 - if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x,7),1), TRUE, FALSE)
  
  # Confidences
  ## Confidences should be interpreted as the answer to the question: "What are the chances
  ### I will see X outside?", with X being cloud, cloud shadow, etc. 
  
  ## Bits 8-9 - if cloud conf low or no level set, then false
  ### 0=no level set, 1=low, 2=medium, 3=high
  cloudConf <- ifelse(bitwAnd(bitwShiftR(x,8), 3) == 1, FALSE, TRUE)
  ## Bits 10-11 - if cloud shadow confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cloudShadowConf <- ifelse(bitwAnd(bitwShiftR(x,10), 3) == 1, FALSE, TRUE)
  ## Bits 12-13 - if snow/ice confidence low, then false
  ### 0=no level set, 1=low, 2=reserved, 3=high
  snowConf <- ifelse(bitwAnd(bitwShiftR(x, 12), 3) == 1, FALSE, TRUE)
  ## Bits 14-15 - if low cirrus confidence, then false; 
  ### 0=no level set, 1=low, 2=reserved, 3=high
  cirrusConf <- ifelse(bitwAnd(bitwShiftR(x,14), 3) == 1, FALSE, TRUE)
  
  return(list(fill=fill, dilatedCloud=dilatedCloud, cirrus=cirrus, cloud=cloud, 
              cloudShadow=cloudShadow, snow=snow, clear=clear, water=water,
              cloudConf=cloudConf, cloudShadowConf=cloudShadowConf, 
              snowConf=snowConf, cirrusConf=cirrusConf)
  )
}
parseL8SR_radsat <- function(x){
  # QA_RADSAT: We only want best images, so no saturation and no occlusion
  
  #is saturated?
  b1 <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  b2 <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  b3 <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  b4 <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  b5 <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  b6 <- ifelse(bitwAnd(bitwShiftR(x,5), 1), TRUE, FALSE)
  b7 <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  #band 8 is not used
  b9 <- ifelse(bitwAnd(bitwShiftR(x,8), 1), TRUE, FALSE)
  terrainOcclusion <- ifelse(bitwAnd(bitwShiftR(x,11), 1), TRUE, FALSE)
  
  return(list(
    b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, b6=b6, b7=b7, b9=b9,
    terrainOcclusion=terrainOcclusion
  ))
}
parseL8SR_aerosol <- function(x){
  # SR_QA_AEROSOL: We want best images, so no fill, no water, and low aerosol
  ## difference with climatology if correction applied (see user guide link in
  ## L8 section below)
  
  # Bit 0: if fill, then true
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  # Bit 2: if water, then true
  water <- ifelse(bitwAnd(bitwShiftR(x, 2), 1), TRUE, FALSE)
  #aerosol level; 0=climatology (no correction), 1=low, 2=med, 3=high)
  aerosolLow <- ifelse(bitwAnd(bitwShiftR(x, 6), 3) < 2, TRUE, FALSE)
  return(list(fill=fill, water=water, aerosolLow=aerosolLow))
}

qa_remap <- function(G, qa=""){
  if(qa=="scl"){
    # https://sentinels.copernicus.eu/web/sentinel/technical-guides/sentinel-2-msi/level-2a/algorithm
    return(ifelse(G %in% c(4,5), "good", "bad"))
  }
  
  if(qa=="pixel"){
    bit_output <- parseL8SR_pixel(G)
    cond_false <- c("fill", "dilatedCloud", "cirrus", "cloud", "cloudShadow", 
                    "snow", "water", "cloudConf", "cloudShadowConf", 
                    "snowConf", "cirrusConf")
    cond_true <- c("clear")
    
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
  
  if(qa=="radsat"){
    bit_output <- parseL8SR_radsat(G)
    cond_false <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7",
                    "b9", "terrainOcclusion")
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0, "good", "bad"))
  }
  
  if(qa=="aerosol"){
    ### pixel = 0 && radsat = 0. We are not filtering out based on aerosol due to low numbers of obs
    bit_output <- parseL8SR_aerosol(G)
    cond_false <- c("fill", "water")
    cond_true <- c("aerosolLow")
    #in other words, the data is valid (TRUE), if aerosol is low.
    return(ifelse(sum(sapply(bit_output[cond_false], sum))==0 &
                    sum(bit_output[[cond_true]]) == 1, "good", "bad"))
  }
}

# assuming the bands are in order (e.g., Band 4 then Band 5)
calcIndex <- function(dataMat, s=sat, b=bandNames){
  
  # no scaling for S2 bc already done in L2A conversion
  bandA <- dataMat[[b[1]]]
  bandB <- dataMat[[b[2]]]
  bandC <- dataMat[[b[3]]]
  
  # scaling comes from Table 6-1 https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1619_Landsat-8-Collection2_Level-2_Science-Product-Guide-v3.pdf
  if(grepl("landsat", s)){
    bandA <- bandA*0.0000275 + (-0.2)
    bandB <- bandB*0.0000275 + (-0.2)
    bandC <- bandC*0.0000275 + (-0.2)
  }
  
  return(cbind((bandB - bandA) / (bandB + bandA),bandC))
}


# apply QA/QC to the bands individually, then calculate ndvi. This returns a matrix of ndvi
filterThenIndex <- function(input, sat=names(vrts)[X], bandNames=bands){
  qaBands <- bandNames[!grepl("B", bandNames)]
  indexBands <- bandNames[grepl("B", bandNames)]
  
  valid <- rep(T,nrow(input))
  
  if(grepl("landsat", sat)){
    qaBands <- qaBands[!grepl("AEROSOL", qaBands)]
    qaNames <- gsub("QA_", "", qaBands)
    
    #filter out values that are not valid (outside of valid range based on user guide - see above in calcIndex)
    for(j in 1:length(indexBands)){
      valid[input[[bandNames[j]]] < 7273 | input[[bandNames[j]]] > 43636] <- F
    }
  } else {
    qaNames <- qaBands
  }
  
  for(i in 1:length(qaBands)){
    look <- data.table(v=sort(unique(as.numeric(input[[qaBands[i]]]))))
    look[, qa := sapply(v, qa_remap, qa=tolower(qaNames[i]))] # apply qa_map function
    
    valid[input[[qaBands[i]]] %in% look[qa == "bad", v]] <- F
    #for(j in 1:length(indexBands)){
    #  input[[bandNames[j]]][input[[qaBands[i]]] %in% look[qa == "bad", v]] <- NA #screen out qa values that are invalid (1)
    #}
  }
  
  return(data.frame(valid,calcIndex(input, s=sat, b=bandNames)))
}

parseL8SR <- function(x){
  #is the pixel fill?
  fill <- ifelse(bitwAnd(x, 1), TRUE, FALSE)
  #is the pixel clear?
  clear <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  #is the pixel water?
  water <- ifelse(bitwAnd(bitwShiftR(x,2),1), TRUE, FALSE)
  #is cloud shadow?
  shadow <- ifelse(bitwAnd(bitwShiftR(x,3),1), TRUE, FALSE)
  #is snow?
  snow <- ifelse(bitwAnd(bitwShiftR(x,4),1), TRUE, FALSE)
  #is cloud?
  cloud <- ifelse(bitwAnd(bitwShiftR(x,5),1), TRUE, FALSE)
  #cloud confidence; 0=none, 1=low, 2=medium, 3=high
  conf_cloud <- bitwAnd(bitwShiftR(x, 6), 3)
  #cirrus confidence; 0=not set, 1=low, 2=medium, 3=high
  conf_cirrus <- bitwAnd(bitwShiftR(x, 8), 3)
  #is terrain occluded? (hidden terrain?)
  occl <- ifelse(bitwAnd(bitwShiftR(x,10),1), TRUE, FALSE)

  return(list(
    fill=fill, clear=clear, water=water, shadow=shadow, snow=snow,
    cloud=cloud, confidence_cloud=conf_cloud, confidence_cirrus=conf_cirrus,
    terrain_occlusion=occl
  ))
}
parseL8SR_radsat <- function(x){
  #is fill?
  fill <- ifelse(bitwAnd(x, 1), "invalid", "valid")
  #is saturated?
  b1 <- ifelse(bitwAnd(bitwShiftR(x,1), 1), TRUE, FALSE)
  b2 <- ifelse(bitwAnd(bitwShiftR(x,2), 1), TRUE, FALSE)
  b3 <- ifelse(bitwAnd(bitwShiftR(x,3), 1), TRUE, FALSE)
  b4 <- ifelse(bitwAnd(bitwShiftR(x,4), 1), TRUE, FALSE)
  b5 <- ifelse(bitwAnd(bitwShiftR(x,5), 1), TRUE, FALSE)
  b6 <- ifelse(bitwAnd(bitwShiftR(x,6), 1), TRUE, FALSE)
  b7 <- ifelse(bitwAnd(bitwShiftR(x,7), 1), TRUE, FALSE)
  #band 8 is not used
  b9 <- ifelse(bitwAnd(bitwShiftR(x,9), 1), TRUE, FALSE)
  b10 <- ifelse(bitwAnd(bitwShiftR(x,10), 1), TRUE, FALSE)
  b11 <- ifelse(bitwAnd(bitwShiftR(x,11), 1), TRUE, FALSE)

  return(list(
    is_good_data = fill, b1=b1, b2=b2, b3=b3, b4=b4, b5=b5, b6=b6,
    b7=b7, b9=b9, b10=b10, b11=b11
  ))
}
parseL8SR_aerosol <- function(x){
  #aerosol level; 0=climatology (data not obtained, 1=low, 2=med, 3=high)
  aerosol <- bitwAnd(bitwShiftR(x, 6), 3)
  return(list(aerosol = aerosol))
}
parseL8TOA <- function(x){
  # fill
  fill <- ifelse(bitwAnd(x,1), TRUE, FALSE)
  #occlusion
  occl <- ifelse(bitwAnd(bitwShiftR(x,1),1), TRUE, FALSE)
  #rad_sat; 0=none, 1=1-2 bands, 2=3-4 bands, 3=5+
  rad_sat <- bitwAnd(bitwShiftR(x,2), 3)
  #cloud
  cloud <- ifelse(bitwAnd(bitwShiftR(x,4),1), TRUE, FALSE)
  #cloud confidence; 0=na or no, 1=low, 2=med, 3=high
  conf_cloud <- bitwAnd(bitwShiftR(x,5), 3)
  #cloud shadow confidence; 0=na/no, 1=low, 2=med, 3=high
  conf_shad <- bitwAnd(bitwShiftR(x,7), 3)
  #snow/ice confidence; 0=na/no, 1=low, 2=med, 3=high
  conf_snow <- bitwAnd(bitwShiftR(x,9), 3)
  #cirrus confidence; 0=na/no, 1=low, 2=med, 3=high
  conf_cirrus <- bitwAnd(bitwShiftR(x,11), 3)

  return(list(
    fill=fill, terrain_occlusion = occl, rad_sat=rad_sat, cloud=cloud,
    conf_cloud=conf_cloud, conf_shadow=conf_shad, conf_snow=conf_snow,
    conf_cirrus=conf_cirrus
  ))
}
parseS2TOA <- function(x){
  #clouds
  clouds <- ifelse(bitwAnd(bitwShiftR(x,10),1), TRUE, FALSE)
  #cirrus
  cirrus <- ifelse(bitwAnd(bitwShiftR(x,11),1), TRUE, FALSE)
  
  return(list(clouds=clouds, cirrus=cirrus))
}
cloudmask <- function(dt = dt, sat_list=sat_list, var="aerosol", sat="mod09ga", nm="state", qcnm="state_1km",...){
  
  if(grepl("s2", sat)){
    # classSCL <- c(1,3,8:11)
    # dt <- dt[, valid_true :=  ifelse(SCL %in% classSCL, FALSE, valid)
    #          ][, valid := NULL]
    # setnames(dt, old="valid_true", new="valid")
    
    ## s2cloudless
    ### we only need the valid column, so this is enough
    dt <- dt[, valid := ifelse(cloudmask == 1, FALSE, TRUE)]
    return(dt)
  } else {
    look <- 
      data.table(
        bitval= 
          unique(dt[,get(names(dt)[grepl(nm, names(dt))])]))
    names(look) <- "bitval"
    
    applyfilter <- function(x, var){
      result <- #return TRUE if data is valid (QA means clear conditions)
        if(is.na(x)){FALSE} else {
          if(sat=="l8sr"){ #L8SR =======================================
            if(qcnm=="pixel_qa"){
              bit_output <- parseL8SR(x)
              cond_false <- 
                c("fill", "water", "shadow", "snow", "cloud", 
                  "terrain_occlusion")
              cond_true <- c("clear")
              
              ifelse(
                all(lapply(bit_output[cond_false], isFALSE)) &
                  bit_output$clear == TRUE &
                  bit_output$confidence_cloud %in% c(0,1) &
                  bit_output$confidence_cirrus %in% c(0,1),
                TRUE, FALSE)
            } else if(qcnm=="radsat_qa"){
              bit_output <- parseL8SR_radsat(x)
              cond_false <- c("b1", "b2", "b3", "b4", "b5", "b6", "b7",
                              "b9", "b10")
              ifelse(all(lapply(bit_output[cond_false], isFALSE)) &
                       bit_output$is_good_data == "valid",
                     TRUE, FALSE)
            } else if(qcnm=="sr_aerosol"){
              bit_output <- parseL8SR_aerosol(x)
              
              #in other words, the data is valid (TRUE), if aerosol is low.
              ifelse(bit_output$aerosol == 1, TRUE, FALSE)
            }
          } else if(grepl("l8toa", sat)){ #L8TOA =================================
            bit_output <- parseL8TOA(x)
            cond_false <- c("fill", "terrain_occlusion", "cloud")
            cond_zero <- c("rad_sat", "conf_cloud", "conf_shadow", "conf_snow",
                           "conf_cirrus")
            ifelse(all(lapply(bit_output[cond_false], isFALSE)) &
                     all(lapply(bit_output[cond_zero], function(x) x<=1)),
                   TRUE, FALSE)
            
          } else if(grepl("mod", sat)){ #MOD09 ==================================
            bit_output <- parseMOD09(x)
            
            if(var == "aerosol"){
              ifelse(bit_output$aero_quan==1, TRUE, FALSE)
            } else {
              cond_false <- c("cloud_shadow", "cloud")
              ifelse(bit_output$cloud_state %in% c(0,3) &
                       bit_output$cirrus_detected == 0 &
                       bit_output$lw_flag == 1 &
                       all(lapply(bit_output[cond_false], isFALSE)),
                     TRUE, FALSE)
            }
          } 
          # else if(grepl("s2", sat)){ #S2 =======================================
          #   bit_output <- parseS2TOA(x)
          #   ifelse(bit_output$clouds == FALSE &
          #            bit_output$cirrus == FALSE,
          #          TRUE, FALSE)
          # }
        }
      return(result)
    }
    
    look[,("valid") := sapply(look$bitval, applyfilter, var)]
    
    setkeyv(dt, qcnm)
    setkeyv(look, "bitval")
    dt <- dt[look]
    setkey(dt, NULL)
    
    return(dt[,.SD,.SDcols=c("date","pointid", qcnm, "valid")])
    }
  }
  #
  ### original s2 - modis bi-cloudmasking ####
  # else if(sat=="s2toa"){ #for time being (July 2020), filter out points if modis shows cloud
  #   nm <- "state"
  #   sat <- "mod"
  #   qcnm <- "state_1km"
  #   var <- "nope"
  #   mod <- sat_list[["mod09ga"]][,c("date", "pointid", "state_1km")]
  #   look <-
  #     data.table(
  #       bitval=
  #         unique(mod[,.SD, .SDcols=grep(nm, names(mod))]))
  #   names(look) <- "bitval"
  #   
  #   look[,("valid") := sapply(look$bitval, applyfilter, var=var)]
  #   
  #   setkeyv(mod, qcnm)
  #   setkeyv(look, "bitval")
  #   mod <- mod[look]
  #   setkey(mod, NULL)
  #   
  #   dt <- dt[,`:=` (date = as.Date(date), matchid=paste0(date, pointid))]
  #   mod <- mod[, date := as.Date(date)
  #   ][date %in% dt[,date], 
  #   ][, matchid := paste0(date, pointid)]
  #   
  #   dt <- merge(dt, mod[,.(valid, matchid)], by="matchid", all.x=TRUE)
  #   setnames(dt, old=c("valid.x", "valid.y"), new=c("valid", "modval"))
  #   
  #   # tp <- mod[,valid][match(dt[,matchid], mod[,matchid])]
  #   
  #   dt <- dt[, valid := ifelse(modval==FALSE, FALSE, valid)
  #   ][, modval := NULL]
  #   qcnm <- "QA60"
  #   nm <- "QA"
  # }
  # 
  # if(grepl("s2", sat)){ 
  #   #after using the mod09ga cloudmask to also mask s2toa, there are duplicate days
  #   #(byproduct of overlapping graticules) that have contrasting TRUE/FALSE 
  #   #designations for being valid. 
  #   #thus I consolidate those and say, if there are multiple readings for one day 
  #   #and any one of those has a reading of FALSE, then assign FALSE.
  #   dt <- do.call(rbind, lapply(unique(dt[,pointid]), function(x){
  #     d <- dt[pointid==x,]
  #     d <- do.call(rbind, lapply(unique(d[,date]), function(y){
  #       f <- d[date==as.Date(y),]
  #       
  #       if(nrow(f)==1){
  #         return(f)
  #       } else {
  #         f <- f[,validno :=
  #                  lapply(c(1:nrow(f)),
  #                         function(z){
  #                           out <- ifelse(is.na(f[,valid][z]) |
  #                                           is.na(f[,get(qcnm)][z]), NA,
  #                                         ifelse(isTRUE(f[,valid][z]), 1, 0))
  #                           return(out)
  #                         })]
  #         
  #         f <- f[, valid := ifelse(all(f[,validno]==1) |
  #                                    all(f[validno == !is.na(validno),validno]==1) &
  #                                    !(0 %in% f[,validno]), TRUE, FALSE)
  #         ][, validno := NULL]
  #         return(f)}
  #     }))
  #     return(d)
  #   })
  #   )
  # }

## onboard data and use the cloudmask functions above ####
onboard_cloud <- function(p ="myanmar", inputFolder="gee", 
                          pointData = "training_locations.csv",
                          saveRdata=TRUE){
  files <- list.files(paste0(dataFolder, inputFolder), pattern=".csv")
  files <- files[!grepl("_", files)]
  fname <- c("l8sr")#, "l8toa", "mod09ga", "s1", "s2sr", "s2toa") #order of files
  # fname <- c("l8sr", "l8toa", "mod09ga", "mod09gq", "s1", 
  #            "s2sr", "s2toa")
  
  sat_list <- list()
  for(i in 1:length(fname)){
    data <- fread(paste0(dataFolder, inputFolder, "/", files[i]))
    data$date <- as.Date(data$date,format="%m/%d/%Y")
    sat_list[[i]] <- data
  }
  names(sat_list) <- fname
  
  ## run the cloudmask function
  var <- "nope"
  totals <- NULL
  nonS1 <- which(!grepl("s1", fname))
  for(j in nonS1){
    sat <- names(sat_list)[[j]]
    dt <- sat_list[[j]]
    #
    # MODIS ===========================================================
    if(grepl("mod", sat)){
      #see table 9 (gq) and 10 (ga) 
      #https://lpdaac.usgs.gov/documents/306/MOD09_User_Guide_V6.pdf
      
      ##reminder that we are using state_1km for QC for both ga and gq because 
      ##a., gq doesn't have its own equivalent and
      ##b., ga is meant to be used in conjunction with gq according to GEE
      
      nm <- "state"
      # sq <- c(3,9:10,11,17:32) #cirrus (bits 9-10 (real 8-9) need to be 00, aka none)
      qcnm <- "state_1km"
      dt <- sat_list[["mod09ga"]][,c("date", "pointid", "state_1km")]
      vec <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
      
      ##aerosol mask
      var <- "aerosol" 
      # sq <- c(8) ##only doing with low aerosol (bits 6-7 = 01)
      vec1 <- cloudmask(dt, sat_list, var, sat, nm, qcnm) ## combine
      vec <- vec[,valid_aerosol := vec1$valid
      ][, date := as.Date(date)]
      
      sat_list[[j]] <- sat_list[[j]][, date := as.Date(date)]
      
      setkeyv(vec, c("pointid", "date"))
      setkeyv(sat_list[[j]], c("pointid", "date"))
      sat_list[[j]] <- sat_list[[j]][vec]
      setkey(sat_list[[j]], NULL)
      
      var <- "nope"
      
      #mod09gq (250m): QC band should be 0 for first bit 
      ##(indicating a produced product),
      ## then Bits 4-7 and 8-11 are a 4-bit series for 
      ##each band, indicating quality.
      # Highest quality is having 0000, so we say sum(bits)==0. TABLE 9
      
      #mod09ga (500m): QC band same as gq above, and 
      ##the 4-bit series is also the same,
      # but this time there are 7 bands. Highest quality condition still applies, so
      # we want sum(bits)==0 (Bits 2-29).At this point, not worries about atmospheric 
      # or adjaency correction. TABLE 10
    }
    
    # Landsat 8 =======================================================
    if(grepl("l8", sat)){
      ## https://prd-wret.s3.us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1574_L8_Data_Users_Handbook-v5.0.pdf
      ## https://prd-wret.s3-us-west-2.amazonaws.com/assets/palladium/production/atoms/files/LSDS-1368_L8_SurfaceReflectanceCode-LASRC_ProductGuide-v2.pdf
      
      setkeyDup <- function(DT){
        setkeyv(DT, c("pointid", "date")) #order by pointid + date
        DT <- DT[, `:=` (uniqueid = 1:nrow(DT))] #add a unique id
        setkeyv(DT, c("uniqueid", "pointid", "date")) #reorder using the unique id
        return(DT)
      }
      
      if(sat=="l8toa"){
        nm <- qcnm <- "BQA"
        # sq <- c(1,2,4,5,7,9,11,13)
        vec <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
        vec <- vec[, date := as.Date(date)][, "BQA" := NULL]
      } else if(sat=="l8sr"){
        for(k in 1:3){
          if(k==1){
            nm <- qcnm <- "pixel_qa"
            # sq <- c(1,3:6,8,10:11)
            vec_pix <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
          } else if(k==2){
            nm <- qcnm <- "radsat_qa"
            # sq <- c(1:32)
            vec_rad <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
          } else if(k==3){
            nm <- qcnm <- "sr_aerosol"
            # sq <- c(8)
            vec_aer <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
          }
        }
        
        #order by a unique id in order to combine across all three at once
        setnames(vec_pix, old="valid", new="valid_qa")
        setnames(vec_rad, old="valid", new="valid_radsatqa")
        setnames(vec_aer, old="valid", new="valid_aerosol")
        
        vec_pix <- setkeyDup(vec_pix)
        vec_rad <- setkeyDup(vec_rad)
        vec_aer <- setkeyDup(vec_aer)
        vec_pix <- vec_pix[vec_rad[vec_aer]]
        
        #remove unnecessary cols
        vec_pix <- vec_pix[,radsat_qa := NULL]
        vec_pix <- vec_pix[,sr_aerosol := NULL]
        setkey(vec_pix, NULL)
        
        vec_pix[, `:=` (valid_cloudAndRadqa = 
                          ifelse(valid_qa == TRUE & valid_radsatqa == TRUE,
                                 TRUE, FALSE),
                        date = as.Date(date))]
        vec <- vec_pix[,pixel_qa := NULL]
      }
      
      #make sure the tables are ordered the same, then combine
      vec <- setkeyDup(vec)
      dt <- setkeyDup(dt)
      dt <- dt[vec]
      setkey(dt, NULL)
      
      #assign back to list
      sat_list[[j]] <- dt
      
      
      
      #l8toa: we say Bits 0,1,3,4 (fill, terrain occlusion, saturation, and cloud) should
      # sum==0. Plus, have the confidence for cloud/shadow/snow/cirrus be low
      # (sum(bits 6,8,10,12)) == 0. (remember, these bits are read left to right, so
      # in other words should be read rev(as.integer(intToBits(2720))[6:7]).
      
      #l8sr: for pixel_qa, we want bits 0,2,3,4,5 (fill,water,shadow,snow,cloud) to 
      # be 0 and bit 1 (clear) to be 1. Plus, we want bit 10 (occlusion) to be 0.
      #l8sr: for radsat_qa, we want everything to be 0 (valid data)
      #l8sr: for sr_aerosol, set bits 0,2,3,4,6 to be 0, allowing for aerosol level 
      # of 00 or 01 (bits 6-7).
    }
    
    # Sentinel 2 ======================================================
    if(grepl("s2", sat)){ 
      #focus on QA60 flag
      ## the internal S2 QA60 is pretty terrible, so initially we used the
      ## SCL from s2sr and then used mod09 to mask out s2toa.
      
      ### there is now s2cloudless algorithm in GEE, which is much better,
      ### so we are using that for both.
      
      # nm <- "QA"
      # # sq <- 1:32
      # qcnm <- "QA60"
      # classSCL <- c(1,3,8:11) #scene classification for surf reflec
      
      #with s2 cloudless, we have a simple job
      dt <- cloudmask(dt, sat_list, var, sat, nm, qcnm)
      dt <- setkeyv(dt, c("pointid", "date"))
      sat_list[[j]] <- dt
      
      ### PRE-s2cloudless algorithm ###
      #s2toa: sum(bits)=0, there are only two conditions, Bit 10 (dense cloud)
      # and bit 11 (cirrus).
      
      #s2sr: same as s2toa, but it has additional classification 
      # table (SCL), from which we exclude saturated (1), 
      # cloud shadow (3), cloud medium prob (8),
      # cloud high prob (9), cirrus (10), and snow/ice (11) classes.
      ## this removes a further 99 data points than if we didn't exclude.
    }
    
    # calculate % filtered data ====================
    if(sat=="l8sr"){
      tot <- data.table(satellite = sat, 
                        orig = nrow(sat_list[[j]]), 
                        valid_qa = nrow(sat_list[[j]][valid_cloudAndRadqa==TRUE, ]),
                        valid_aerosol = nrow(sat_list[[j]][valid_aerosol==TRUE, ]),
                        valid_total = nrow(sat_list[[j]][
                          valid_cloudAndRadqa==TRUE &
                            valid_aerosol==TRUE, ]))
    } else if (grepl("mod", sat)){
      tot <- data.table(satellite = sat, 
                        orig = nrow(sat_list[[j]]), 
                        valid_qa = nrow(sat_list[[j]][valid==TRUE,]),
                        valid_aerosol = nrow(sat_list[[j]][valid_aerosol==TRUE, ]),
                        valid_total = nrow(sat_list[[j]][valid==TRUE &
                                                           valid_aerosol==TRUE,]))
    } else {
      tot <- data.table(satellite = sat, 
                        orig = nrow(sat_list[[j]]), 
                        valid_qa = nrow(sat_list[[j]][valid==TRUE,]),
                        valid_aerosol = NA,
                        valid_total = NA)
    }
    tot[,("perc_bad_qa") := 1-valid_qa/orig] #only focus on one main for now
    totals <- rbind(totals, tot)
  }
  
  #bring in dates to be used for model runs
  base <- fread(paste0(dataFolder, pointData))
  base <- base[,.(pointid, datePre, dateDist)]
  setnames(base, old=c("datePre", "dateDist"),
           new=c("date_model", "date_dist"))
  
  ids <- sort(unique(base[,pointid]))
  
  #update dates
  date_fn <- function(dae){
    date_true <- ifelse(is.na(as.Date(dae)) | (Sys.Date() - as.Date(dae)) > 10000, 
                        as.Date(dae, format="%d/%m/%Y"),
                        as.Date(dae))
    date_true <- as.Date(date_true, origin="1970-01-01")
    return(date_true)
  }
  
  base <- base[, `:=` (date_model = date_fn(date_model),
                       date_dist = date_fn(date_dist))]
  
  get_moddate <- function(sat_list=sat_list){
    dat <- sat_list
    setkey(dat, pointid)
    setkey(base, pointid)
    dat <- dat[base]
  }
  sat_list <- lapply(sat_list, get_moddate)
  
  if(saveRdata){
    save(sat_list, file=paste0(dataFolder, "/sat_data_", p, ".Rdata"))
  }
  return(list(totals=as.data.table(totals), sat_list = sat_list, ids=ids))
}
