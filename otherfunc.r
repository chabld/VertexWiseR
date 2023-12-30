## OTHER VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY

############################################################################################################################
############################################################################################################################
## Efficient way to extract t statistics from linear regression models
## adapted from https://stackoverflow.com/questions/15820623/obtain-t-statistic-for-regression-coefficients-of-an-mlm-object-returned-by-l
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
  n_vert=length(data)
  
  #listing out non-zero vertices
  vert=which(data!=0)

  #matching non-zero vertices with adjacency matrices to obtain list of edges connecting between the non-zero vertices
  edgelist0=edgelist[!is.na(match(edgelist[,1],vert)),]
  edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]
  
  if(length(edgelist1)>2) #if at least 2 edges are identified
  {
    #extracting cluster-related info from list of non-zero edges
    com=igraph::components(igraph::graph.data.frame(edgelist1, directed = F))
    clust.size=com$csize
    
    #cluster mappings
    clust.map=rep(NA,n_vert)
    clust.map[as.numeric(names(com$membership))]=com$membership
  
  } else if(length(edgelist1)==2) #bypass cluster extraction procedure if only 1 edge is identified
  {
    clust.size=2
    clust.map=rep(NA,n_vert)
    clust.map[edgelist1]=1
  } else #bypass cluster extraction procedure if no edges are identified
  {
    clust.map="noclusters"
    clust.size="noclusters"
  }
  return(list(clust.map,clust.size))
}
############################################################################################################################
############################################################################################################################
## To extract atlas ROI values from fsaverage5 vertex-wise data and vice-versa
fs5_to_atlas=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
{  
  #check length of vector
  if(length(data)%%20484!=0) {stop("Length of data is not a multiple of 20484")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
  
  #init variables
  nregions=max(ROImap[[1]][,atlas])
  data[is.na(data)]=0

  #mapping fsaverage5 space vertice to atlas regions if data is a 1x20484 vector
  if(length(data)==20484) 
  {
    data=matrix(data,ncol=20484,nrow=1)  
    ROI=rep(NA,nregions)
    for (region in 1:nregions)  {ROI[region]=mean(data[which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  } else 
  {
  #mapping fsaverage5 space vertice to atlas regions if data is a Nx20484 matrix
    ROI=matrix(NA, nrow=NROW(data), ncol=nregions)
    for (region in 1:nregions)  {ROI[,region]=rowMeans(data[,which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  }
  return(ROI)
}

atlas_to_fs5=function(data,atlas) ## atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
  {
    #load atlas mapping data
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap.rdata?raw=TRUE"))
  
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    fs5_dat=rep(NA,20484)
  
    #mapping atlas label to fsaverag5 space
    for (region in 1:nregions)  {fs5_dat[which(ROImap[[1]][,atlas]==region)]=data[region]}
    return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################
##smoothing fsaverage5 and fsaverage6 data
smooth=function(data,FWHM=10)
{
  ##import python libraries
  brainstat.mesh.data=reticulate::import("brainstat.mesh.data")
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  
  ##identify NA vertices
  col0=which(colSums(data==0) == nrow(data))
  
  if(ncol(data)==20484) ##fsaverage5 parameters
    {
      surftemp=brainstat.datasets$fetch_template_surface("fsaverage5", join=T)
      vert_mm=3.5
    } else if (ncol(data)==81924) ##fsaverage6 parameters
    {
      surftemp=brainstat.datasets$fetch_template_surface("fsaverage6", join=T)
      vert_mm=2
    } 
    else {stop("data vector should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}

  ##smoothing
  smooth=brainstat.mesh.data$mesh_smooth(Y=CT_dat,surf=surftemp, FWHM = FWHM/vert_mm)
  
  ##the smoothing process might eat into previously identified NA vertices, hence we need to recode the previously identified NA vertices to NA; 
  smooth[,col0]=NA
  return(smooth)  
}
############################################################################################################################
############################################################################################################################
##CT surface plots
plotCT=function(data, filename,title="",surface="inflated",cmap,fs_path, range=NULL , colorbar=T)
{
  #check length of vector
  n_vert=length(data)
  if(n_vert==20484) {template="fsaverage5"}
  else if (n_vert==81924) {template="fsaverage6"} 
  else {stop("data vector should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}

  #legacy input message
  if(!missing("fs_path")){cat("The fs_path parameter and the fsaverage5 files are no longer needed in the updated plotCT function\n")}

  #import python libraries
  brainstat.datasets=reticulate::import("brainstat.datasets")  
  brainspace.plotting=reticulate::import("brainspace.plotting")  

  #loading fsaverage surface
  left=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[1]
  right=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[2]

  #setting color maps
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

  #plot object
  CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(data),cmap=cmap, 
                                              size=reticulate::tuple(as.integer(c(1920,400))),nan_color=reticulate::tuple(c(0.7, 0.7, 0.7, 1)),
                                              return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=range,
                                              label_text=list('left'=list(title)),interactive=F, color_bar=colorbar,  transparent_bg=FALSE)
  #output plot as a .png image
  CTplot$screenshot(filename=filename,transparent_bg = F)
}
############################################################################################################################
############################################################################################################################
##CT image decoding
decode_img=function(img,contrast="positive")
{
  ##checks
    #check length of vector
    n_vert=length(img)
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {template="fsaverage6"} 
    else {stop("data vector should only contain 20484 (fsaverage5) or 81924 (fsaverage6) columns")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  ##import python libraries
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
    stat_nii = interpolate$`_surf2vol`(template, stat_labels$flatten())
  } else if (contrast=="negative")
  {
    img[is.na(img)]=0
    img[img>0]=0
    img[img<0]=1
    stat_labels=r_to_py(img)
    stat_nii = interpolate$`_surf2vol`(template, stat_labels$flatten())
  }
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/neurosynth_dataset.pkl.gz?raw=TRUE",destfile = "neurosynth_dataset.pkl.gz")
  } 
  ##running the decoding procedure
  neurosynth_dset = nimare.dataset$Dataset$load("neurosynth_dataset.pkl.gz")
  cat("Correlating input image with images in the neurosynth database. This may take a while\n")
  decoder = discrete$ROIAssociationDecoder(stat_nii)
  decoder$fit(neurosynth_dset)

  ##compiling the results
  decoder_df = data.matrix(decoder$transform())
  row.names(decoder_df)=gsub(pattern = "terms_abstract_tfidf__",x=row.names(decoder_df), replacement = "")
  result=data.frame(row.names(decoder_df),round(as.numeric(decoder_df),3))
  colnames(result)=c("keyword","r")
  result=result[order(-result$r),]
  
  return(result)
}  
############################################################################################################################
############################################################################################################################
