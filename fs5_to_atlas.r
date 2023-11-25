## To extract atlas ROI values from fsaverage5 vertex-wise data
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################

fs5_to_atlas=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360
{
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
  nregions=max(ROImap[[1]][,atlas])
  if(length(data)%%20484!=0)
  {
    stop("Length of data is not a multiple of 20484")
  }
  data[is.na(data)]=0
  if(length(data)==20484)
  {
  data=matrix(data,ncol=20484,nrow=1)  
  ROI=rep(NA,nregions)
    for (region in 1:nregions)
    {
      ROI[region]=mean(data[which(ROImap[[1]][,atlas]==region)])
    }
  } else 
  {
    ROI=matrix(NA, nrow=NROW(data), ncol=nregions)
    for (region in 1:nregions)
    {
      ROI[,region]=rowMeans(data[,which(ROImap[[1]][,atlas]==region)])
    }
  }
  return(ROI)
}

atlas_to_fs5=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360
  {
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
  nregions=max(ROImap[[1]][,atlas])
  fs5_dat=rep(NA,20484)
  for (region in 1:nregions)
    {
    fs5_dat[which(ROImap[[1]][,atlas]==region]]=data[region]
    }
  return(as.numeric(fs5_dat))
  }
