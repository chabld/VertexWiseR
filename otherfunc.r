## OTHER VERTEX-WISE FUNCTIONS
## FOR USE IN THE COGNITIVE AND BRAIN HEALTH LABORATORY
############################################################################################################################
############################################################################################################################
## permutation function for random subject effects
## Paired/grouped data points are first shuffled within subjects, then these pairs/groups are shuffled between subjects
perm_within_between=function(random)
{
  ##for groups of 2 or more (subjects with 2 or more measurements)
  perm.idx=rep(NA, length(random))
  for(count in 2:max(table(random)))
  {
    if(length(which(table(random)==count)>0))
    {
      sub.id=as.numeric(which(table(random)==count))
        if(length(sub.id)>1)
        {
          ##between group shuffling
          recode.vec=sample(sub.id)
          vec.idx=1
          for(sub in sub.id)
          {
            perm.idx[which(random==sub)]=sample(which(random==recode.vec[vec.idx])) ##sample— within subject shuffling
            vec.idx=vec.idx+1
          }   
          remove(vec.idx,recode.vec)  
        } else 
        {
          ##if only one subject has a certain count, between subject shuffling will not be possible, only within-subject shuffling will be carried out
          perm.idx[which(random==sub.id)]=sample(which(random==sub.id)) ##sample— within subject shuffling
        }
    }
  }
  ##for subjects with a single measurement
  sub.idx=which(is.na(perm.idx))
  if(length(sub.idx)>1)
  {
    perm.idx[sub.idx]=sample(sub.idx)  
  } else 
  {
    perm.idx[sub.idx]=sub.idx
  }
  return(perm.idx)
}

############################################################################################################################
############################################################################################################################
## smooth fsaverage5/6 images
smooth=function(data, FWHM)
  {
  reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/python/smooth.py?raw=TRUE")
  
  ##setting default FWHM values for fsaverage5/6 space
  if(missing("FWHM")) 
    {
    if(NCOL(data)==20484) {FWHM=10}
    else if(NCOL(data)==81924) {FWHM=5}
    }
  smoothed=mesh_smooth(data, FWHM)
  smoothed[is.na(smoothed)]=0
  return(smoothed)
  }
############################################################################################################################
############################################################################################################################
## Efficient way to extract t statistics from linear regression models to speed up the permutation process
## adapted from https://stackoverflow.com/questions/15820623/obtain-t-statistic-for-regression-coefficients-of-an-mlm-object-returned-by-l
extract.t=function(mod,row)
{
  p = mod$rank
  df.residual=NROW(mod$residuals)-NROW(mod$coefficients)
  rdf = df.residual
  Qr = mod$qr
  p1 = 1L:p
  r = mod$residuals
  rss = colSums(r^2)
  resvar = rss/rdf
  R = chol2inv(Qr[p1, p1, drop = FALSE])  
  se = (sqrt(diag(R) %*% t(resvar)))[row,]
  est = mod$coefficients[row,]
  tval = est/se 
  return(tval)
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
  if(length(edgelist0)>2)  {edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]} 
  else if (length(edgelist0)==2)  ##if only a single edge was identified, edgelist will no longer be a Nx2 matrix, hence need to reshape it into a matrix
    { 
    edgelist0=matrix(edgelist0,ncol=2,nrow=1)
    edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]
    } else {edgelist1=0}
  remove(data,vert,edgelist0)
  
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
## Atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
## ROI to vertex mapping data for 1 to 4 obtained from enigmatoolbox https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations
## ROI to vertex mapping data for 5 obtained from nilearn.datasets.fetch_atlas_surf_destrieux https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py
fs5_to_atlas=function(data,atlas) 
{  
  #check length of vector
  if(length(data)%%20484!=0) {stop("Length of data is not a multiple of 20484")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
  
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

atlas_to_fs5=function(data,atlas) 
  {
    #load atlas mapping data
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
  
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    fs5_dat=rep(NA,20484)
  
    #mapping atlas label to fsaverag5 space
    for (region in 1:nregions)  {fs5_dat[which(ROImap[[1]][,atlas]==region)]=data[region]}
    return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################
#convert between fsaverage5 and fsaverage6 spacing
fs5_to_fs6=function(data)
{
  #check length of vector
  if(length(data)%%20484!=0) {stop("Length of data is not a multiple of 20484")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  #mapping fsaverage5 to fsaverage6 space if data is a Nx20484 matrix
  if(length(data)==20484) {data.fs6=data[fs6_to_fs5]} 
  #mapping fsaverage5 to fsaverage6 space if data is a Nx20484 matrix
  else {data.fs6=data[,fs6_to_fs5]}
  return(data.fs6)
}

fs6_to_fs5=function(data)
{
  #check length of vector
  if(length(data)%%81924!=0) {stop("Length of data is not a multiple of 81924")}
  
  #load atlas mapping data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  
  if(length(data)==81924) #mapping fsaverage6 to fsaverage5 space if data is a Nx20484 matrix
  {
    data=matrix(data,ncol=81924,nrow=1)  
    data.fs5=matrix(NA,ncol=20484,nrow=1)
    
    for (vert in 1:20484)  {data.fs5[vert]=mean(data[fs6_to_fs5==vert],na.rm = T)} 
  } else #mapping fsaverage6 to fsaverage5 space if data is a Nx20484 matrix
  {
    data.fs5=matrix(NA,ncol=20484,nrow=NROW(data))
    for (vert in 1:20484)  {data.fs5[,vert]=rowMeans(data[,fs6_to_fs5==vert],na.rm = T)} 
  }
  return(data.fs5)
}
############################################################################################################################
############################################################################################################################
##Cortical surface/hippocampal plots
##input data can be a matrix with multiple rows, for multiple plots in a single .png file
plotCT=function(data, filename,title="",surface="inflated",cmap,fs_path, limits, colorbar=T)
{
  #format title for single row
  if(is.null(nrow(data)))  {title=list('left'=list(title))}
  
  #check length of vector
  n_vert=length(data)
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524!=0) {stop("data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
  #legacy input message
  if(!missing("fs_path"))  {cat("The fs_path parameter and the fsaverage5 files are no longer needed in the updated plotCT function\n")}
  
  #setting color maps
  if(missing("cmap"))
  {
    if(range(data,na.rm = T)[1]>=0)  {cmap="Reds"}
    else if (range(data,na.rm = T)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
  
  #setting limits
  if(missing("limits")) {limits=range(data,na.rm = T)}
  limits=reticulate::tuple(limits[1],limits[2])
  
  if(n_vert%%14524!=0)
  {
    ##cortical surface fplots
      #import python libraries
      brainstat.datasets=reticulate::import("brainstat.datasets")  
      brainspace.plotting=reticulate::import("brainspace.plotting")  
      
      #loading fsaverage surface
      left=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[1]
      right=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[2]
      
      CTplot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(data),cmap=cmap, 
                                                  size=reticulate::tuple(as.integer(c(1920,400))),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                  return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=limits,
                                                  label_text=title,interactive=F, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
  } else
  {
    ##hippocampal plots
      #import python libraries
      reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/python/hipp_plot.py?raw=TRUE")
      
      #reshaping data into a 7262 x 2 x N array
      if(is.null(nrow(data)))  {data=cbind(data[1:7262],data[7263:14524])} 
      else  {data=array(cbind(data[,1:7262],data[,7263:14524]),c(7262,2,nrow(data)))}
      
      CTplot=surfplot_canonical_foldunfold(data,color_bar=colorbar,share="row",nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                           cmap=cmap,color_range=limits,label_text=title, return_plotter=T,interactive=F) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  CTplot$screenshot(filename=filename,transparent_bg = F)
}
############################################################################################################################
############################################################################################################################
##converting surface to volumetric data and exporting it as a .nii file

surf_to_vol=function(surf_data, filename="output.nii")
  {
  #check length of vector
    n_vert=length(surf_data)
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {template="fsaverage6"} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) or 81924 (fsaverage6) is accepted")}
  
  #load python libraries
    interpolate=reticulate::import("brainstat.mesh.interpolate")
    nibabel=reticulate::import("nibabel")

  #convert and export .nii file
    stat_nii = interpolate$`_surf2vol`(template, surf_data)
    nibabel$save(stat_nii,filename)
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
    else if (n_vert==81924) {stop("decoding of fsaverage6-space image is current not implemented, please resample the image to fsaverage5 space")} 
    else {stop("Only an img vector with a length of 20484 (fsaverage5) is accepted")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  ##import python libraries
  interpolate=reticulate::import("brainstat.mesh.interpolate")
  discrete=reticulate::import("nimare.decode")
  nimare.dataset=reticulate::import("nimare.dataset")
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    img[is.na(img)]=0
    img[img<0]=0
    img[img>0]=1
  } else if (contrast=="negative")
  {
    img[is.na(img)]=0
    img[img>0]=0
    img[img<0]=1
  }

  ##convert img vector to nii image
  stat_labels=reticulate::r_to_py(img)
  stat_nii = interpolate$`_surf2vol`(template, stat_labels)
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/neurosynth_dataset.pkl.gz",destfile = "neurosynth_dataset.pkl.gz")
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
