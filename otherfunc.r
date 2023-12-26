## OTHER VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
## Efficient way to extract t statistics from linear regression models
extract.t=function(mod,row)
{
  p = mod$rank
  rdf = mod$df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr$qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se                          
}
############################################################################################################################
############################################################################################################################
##find clusters using edgelist
getClusters=function(data)
{ 
  if(!exists("fs5_edgelist", inherit=F))  {load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs5edgelist.rdata?raw=TRUE"))} 
  
  vert=which(data!=0)
  
  fs5_edgelist0=fs5_edgelist[which(!is.na(match(fs5_edgelist[,1],vert))),]
  edgelist=fs5_edgelist0[which(!is.na(match(fs5_edgelist0[,2],vert))),]
  if(length(edgelist)>2)
  {
    com=igraph::components(igraph::graph.data.frame(edgelist, directed = F))
    clust.size=com$csize
    clust.map=rep(NA,20484)
    #cluster mappings
    for(clust.no in 1:com$no)
    {
      clust.map[as.numeric(names(which(com$membership==clust.no)))]=clust.no
    } 
  } else if(length(edgelist)==2)
  {
    clust.size=2
    clust.map=rep(NA,20484)
    clust.map[edgelist]=1
  } else 
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size))
}

############################################################################################################################
############################################################################################################################
## To extract atlas ROI values from fsaverage5 vertex-wise data
fs5_to_atlas=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
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

atlas_to_fs5=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
  {
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
  nregions=max(ROImap[[1]][,atlas])
  fs5_dat=rep(NA,20484)
  for (region in 1:nregions)
    {
    fs5_dat[which(ROImap[[1]][,atlas]==region)]=data[region]
    }
  return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################
##smoothing fsaverage5 and fsaverage6 data
smooth=function(data,FWHM=10)
{
  brainstat.mesh.data=reticulate::import("brainstat.mesh.data")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  col0=which(colSums(CT_dat==0) == nrow(CT_dat))
  if(ncol(CT_dat)==20484)
    {
    surftemp=brainstat.datasets$fetch_template_surface("fsaverage5", join=T)
    vert_mm=3.5
    }
  if(ncol(CT_dat)==81924)
    {
    surftemp=brainstat.datasets$fetch_template_surface("fsaverage6", join=T)
    vert_mm=2
    }
  smooth=brainstat.mesh.data$mesh_smooth(Y=CT_dat,surf=surftemp, FWHM = FWHM/vert_mm)
  smooth[,col0]=0
  return(smooth)  
}
############################################################################################################################
############################################################################################################################
##smoothing fsaverage5 and fsaverage6 data
##CT surface plots
plotCT=function(data, filename,title="",surface="inflated",cmap,fs_path, range=NULL)
{
  if(length(data) != 20484)
  {
    stop("Data has to be a numeric vector with 20484 values")
  } 
  if(!missing("fs_path")){cat("The fs_path parameter and the fsaverage5 files are no longer needed in the updated plotCT function\n")}
  
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.plotting=reticulate::import("brainspace.plotting")  

  left=brainstat.datasets$fetch_template_surface("fsaverage5", join=F, layer=surface)[1]
  right=brainstat.datasets$fetch_template_surface("fsaverage5", join=F, layer=surface)[2]
  
  if(missing("cmap"))
  {
    if(range(data,na.rm = T)[1]>=0)
      {
      cmap="Reds"
      range=NULL
      }
    else if (range(data,na.rm = T)[2]<=0)
      {cmap="Blues_r"
      range=NULL
      }
    else
      {
      cmap="RdBu_r"
      range="sym"
      }  
  }
  
  CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(data),cmap=cmap, 
                                              size=reticulate::tuple(as.integer(c(1920,400))),nan_color=reticulate::tuple(c(0.7, 0.7, 0.7, 1)),
                                              return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=range,
                                              label_text=list('left'=list(title)),interactive=F, color_bar=T,  transparent_bg=FALSE)
  CTplot$screenshot(filename=filename,transparent_bg = F)

}
############################################################################################################################
############################################################################################################################
##CT image decoding
decode_img=function(img,contrast="positive")
{
  ##input checks
  if(length(img) != 20484)
  {
    stop("img object has to be a numeric vector with 20484 values")
  } 
  if(contrast != "positive" & contrast != "negative" )
  {
    stop("contrast has to be either positive or negative")
  } 
  tflow=reticulate::import("templateflow.api")
  np=reticulate::import("numpy")
  interpolate=reticulate::import("brainstat.mesh.interpolate")
  discrete=reticulate::import("nimare.decode")
  nimare.dataset=reticulate::import("nimare.dataset")
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    img[is.na(img)]=0
    img[img<0]=0
    img[img>0]=1
    stat_labels=r_to_py(img)
    stat_nii = interpolate$`_surf2vol`("fsaverage5", stat_labels$flatten())
  } else if (contrast=="negative")
  {
    img[is.na(img)]=0
    img[img>0]=0
    img[img<0]=1
    stat_labels=r_to_py(img)
    stat_nii = interpolate$`_surf2vol`("fsaverage5", stat_labels$flatten())
  }
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/neurosynth_dataset.pkl.gz?raw=TRUE",destfile = "neurosynth_dataset.pkl.gz")
  } 
  ##running the decoding process
  neurosynth_dset = nimare.dataset$Dataset$load("neurosynth_dataset.pkl.gz")
  cat("Correlating input image with images in the neurosynth database. This may take a while\n")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)
  decoder_df = data.matrix(decoder$transform())
  
  ##compiling the results
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  
  return(result)
}  
############################################################################################################################
############################################################################################################################
