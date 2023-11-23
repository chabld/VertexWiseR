## WRAPPERS FOR VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################

fs5_to_sch100=function(data)
{
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/ROImap.rdata?raw=TRUE"))
  if(length(data)%%20484!=0)
  {
    stop("Length of data is not a multiple of 20484")
  }
  data[is.na(data)]=0
  if(length(data)==20484)
  {
  data=matrix(dat,ncol=20484,nrow=1)  
  ROI=rep(NA,100)
    for (region in 1:100)
    {
      ROI[region]=mean(data[which(ROImap[[1]]==region)])
    }
  } else 
  {
    ROI=matrix(NA, nrow=NROW(data), ncol=100)
    for (region in 1:100)
    {
      ROI[,region]=rowMeans(data[,which(ROImap[[1]]==region)])
    }
  }
  return(ROI)
}
