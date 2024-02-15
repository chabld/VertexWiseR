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

#' @title Smooth surface
#'
#' @description Smooths surface data at defined full width at half maximum (FWHM) as per the corresponding template of surface data
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format
#' @param FWHM A numeric vector object containing the desired smoothing width in mm 
#'
#' @returns A matrix object with smoothed vertex-wise values
#' @examples
#' smooth_surf(CT_data, 10)
#' @importFrom reticulate source_python
#' @export

## smooth surface data 
## FWHM input is measured in mm, which is subsequently converted into mesh units
smooth_surf=function(surf_data, FWHM)
{
  ##source python function
  reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/python/smooth.py?raw=TRUE")
  
  n_vert=ncol(surf_data)
  ##select template, set its FWHM parameter and load its edgelist file
  if(n_vert==20484) 
  {
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs5.rdata?raw=TRUE"))
    FWHM=FWHM/3.5 #converting mm to mesh units
  } else if(n_vert==81924) 
  {
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistfs6.rdata?raw=TRUE"))
    FWHM=FWHM/2 #converting mm to mesh units
  } else if(n_vert==14524) 
  {
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/edgelistHIP.rdata?raw=TRUE"))
    FWHM=FWHM/0.5 #converting m to mesh units
  } else {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}
  
  smoothed=mesh_smooth(surf_data,edgelist, FWHM)
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
getClusters=function(surf_data)
{ 
  n_vert=length(surf_data)
  
  #listing out non-zero vertices
  vert=which(surf_data!=0)

  #matching non-zero vertices with adjacency matrices to obtain list of edges connecting between the non-zero vertices
  edgelist0=edgelist[!is.na(match(edgelist[,1],vert)),]
  if(length(edgelist0)>2)  {edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]} 
  else if (length(edgelist0)==2)  ##if only a single edge was identified, edgelist will no longer be a Nx2 matrix, hence need to reshape it into a matrix
    { 
    edgelist0=matrix(edgelist0,ncol=2,nrow=1)
    edgelist1=edgelist0[!is.na(match(edgelist0[,2],vert)),]
    } else {edgelist1=0}
  remove(surf_data,vert,edgelist0)
  
  if(length(edgelist1)>2) #if at least 2 edges are identified
  {
    #extracting cluster-related info from list of non-zero edges
    com=igraph::components(igraph::graph.surf_data.frame(edgelist1, directed = F))
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

#' @title Fs5 to atlas
#'
#' @description Returns the mean vertex-wise surface data in fsaverage5 space for each ROI of a selected atlas
#' @details The function currently works with the Desikan-Killiany, Schaefer-100, Schaefer-200, Glasser-360, or Destrieux-148 atlases. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{enigmatoolbox} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{nilearn.datasets.fetch_atlas_surf_destrieux}
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148. 
#'
#' @returns A matrix object with ROI as column and corresponding average vertex-wise values as row
#' @seealso \code{\link{atlas_to_fs5}}
#' @examples
#' fs5_to_atlas(CTv, 1)
#' @export

## To extract atlas ROI values from fsaverage5 vertex-wise data and vice-versa
## Atlas: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148
## ROI to vertex mapping data for 1 to 4 obtained from enigmatoolbox https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations
## ROI to vertex mapping data for 5 obtained from nilearn.datasets.fetch_atlas_surf_destrieux https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py
fs5_to_atlas=function(surf_data,atlas) 
{  
  #check length of vector
  if(length(surf_data)%%20484!=0) {stop("Length of surf_data is not a multiple of 20484")}
  
  #load atlas mapping surf_data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
  
  #init variables
  nregions=max(ROImap[[1]][,atlas])
  surf_data[is.na(surf_data)]=0

  #mapping fsaverage5 space vertice to atlas regions if surf_data is a 1x20484 vector
  if(length(surf_data)==20484) 
  {
    surf_data=matrix(surf_data,ncol=20484,nrow=1)  
    ROI=rep(NA,nregions)
    for (region in 1:nregions)  {ROI[region]=mean(surf_data[which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  } else 
  {
  #mapping fsaverage5 space vertice to atlas regions if surf_data is a Nx20484 matrix
    ROI=matrix(NA, nrow=NROW(surf_data), ncol=nregions)
    for (region in 1:nregions)  {ROI[,region]=rowMeans(surf_data[,which(ROImap[[1]][,atlas]==region)])} #vertices are averaged within the atlas ROI
  }
  return(ROI)
}


#' @title Atlas to fsaverage5
#'
#' @description Returns the vertex-wise surface data mapped in fsaverage5 space from data parcellated with a selected atlas
#' @details The function currently works with the Desikan-Killiany, Schaefer-100, Schaefer-200, Glasser-360, or Destrieux-148 atlases. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{enigmatoolbox} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{nilearn.datasets.fetch_atlas_surf_destrieux}
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @seealso \code{\link{fs5_to_atlas}}
#' @examples
#' atlas_to_fs5(CTv, 1)
#' @export


atlas_to_fs5=function(surf_data,atlas) 
  {
    #load atlas mapping surf_data
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/ROImap_fs5.rdata?raw=TRUE"))
  
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    fs5_dat=rep(NA,20484)
  
    #mapping atlas label to fsaverage5 space
    for (region in 1:nregions)  {fs5_dat[which(ROImap[[1]][,atlas]==region)]=surf_data[region]}
    return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################


#' @title fsaverage5 to fsaverage6
#'
#' @description Remaps vertex-wise surface data in fsaverage5 space to fsaverage6 space 
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage6 space
#' @seealso \code{\link{fs6_to_fs5}}
#' @examples
#' fs5_to_fs6(CTv)
#' @export

#convert between fsaverage5 and fsaverage6 spacing
fs5_to_fs6=function(surf_data)
{
  #check length of vector
  if(length(surf_data)%%20484!=0) {stop("Length of surf_data is not a multiple of 20484")}
  
  #load atlas mapping surf_data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  #mapping fsaverage5 to fsaverage6 space if surf_data is a vector length of 20484
  if(length(surf_data)==20484) {surf_data.fs6=surf_data[fs6_to_fs5]} 
  #mapping fsaverage5 to fsaverage6 space if surf_data is a Nx20484 matrix
  else {surf_data.fs6=surf_data[,fs6_to_fs5]}
  return(surf_data.fs6)
}

#' @title fsaverage6 to fsaverage5
#'
#' @description Remaps vertex-wise surface data in fsaverage6 space to fsaverage5 space 
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @seealso \code{\link{fs5_to_fs6}}
#' @examples
#' fs6_to_fs5(CTv)
#' @export

fs6_to_fs5=function(surf_data)
{
  #check length of vector
  if(length(surf_data)%%81924!=0) {stop("Length of surf_data is not a multiple of 81924")}
  
  #load atlas mapping surf_data
  load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/fs6_to_fs5.rdata?raw=TRUE"))
  
  if(length(surf_data)==81924) #mapping fsaverage6 to fsaverage5 space if surf_data is a Nx81924 matrix
  {
    surf_data=matrix(surf_data,ncol=81924,nrow=1)  
    surf_data.fs5=matrix(NA,ncol=20484,nrow=1)
    
    for (vert in 1:20484)  {surf_data.fs5[vert]=mean(surf_data[fs6_to_fs5==vert],na.rm = T)} 
  } else #mapping fsaverage6 to fsaverage5 space if surf_data is a Nx20484 matrix
  {
    surf_data.fs5=matrix(NA,ncol=20484,nrow=NROW(surf_data))
    for (vert in 1:20484)  {surf_data.fs5[,vert]=rowMeans(surf_data[,fs6_to_fs5==vert],na.rm = T)} 
  }
  return(surf_data.fs5)
}
############################################################################################################################
############################################################################################################################

#' @title Surface plotter
#'
#' @description Plots surface data in a grid with one or multiple rows, for multiple plots in a .png file
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param filename A string object containing the desired name of the output .png file.
#' @param title A string object containing the title wanted in the plot. Default is none. 
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param cmap A string object containing the colormap for the plot. Options are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{Matplotlib plotting library}. 
#' @param limits Combined pair of numeric vectors composed of the lower limit and the upper limit value of the surface data. Default is min and max values across all the vertices. 
#' @param colorbar A logical object stating whether to include a color bar in the plot or not (default is TRUE).
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @examples
#' plot_surf(CTv, filename = 'output', title = 'Cortical thickness', surface = 'white', cmap = 'RdBu_r', colorbar = FALSE)
#' @importFrom reticulate tuple import np_array source_python
#' @export

##Cortical surface/hippocampal plots
##input surf_data can be a matrix with multiple rows, for multiple plots in a single .png file
plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=T)
{
  #format title for single row
  if(nrow(surf_data)==1) 
  {
    title=list('left'=list(title))
    rows=1
    surf_data=as.numeric(surf_data)
  } else {rows=nrow(surf_data)}
  
  #check length of vector
  n_vert=length(surf_data)
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524!=0) {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}

  #setting color maps
  if(missing("cmap"))
  {
    if(range(surf_data,na.rm = T)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = T)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
  
  #setting limits
  if(missing("limits")) {limits=range(surf_data,na.rm = T)}
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
    
    surf_plot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(surf_data),cmap=cmap, 
                                                size=reticulate::tuple(as.integer(c(1920,rows*400))),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=1.25,color_range=limits,
                                                label_text=title,interactive=F, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
  } else
  {
  ##hippocampal plots
    #import python libraries and hippocampal template data
    reticulate::source_python("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/python/hipp_plot.py?raw=TRUE")
    load(file = url("https://github.com/CogBrainHealthLab/VertexWiseR/blob/main/data/hip_points_cells.rdata?raw=TRUE"))
    
    #reshaping surf_data into a 7262 x 2 x N array
    if(is.null(nrow(surf_data)))  {surf_data=cbind(surf_data[1:7262],surf_data[7263:14524])} #if N=1
    else  
      {
        surf_data.3d=array(NA,c(7262,2,nrow(surf_data))) #if N>1
        for (row in 1:nrow(surf_data))  {surf_data.3d[,,row]=cbind(surf_data[row,1:7262],surf_data[row,7263:14524])}
        surf_data=surf_data.3d
      }
    
    surf_plot=surfplot_canonical_foldunfold(surf_data,hipdat =hip_points_cells,color_bar=colorbar,share="row",nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),size=as.integer(c(350,300)),
                                         cmap=cmap,color_range=limits,label_text=title, return_plotter=T,interactive=F) ##disabling interactive mode because this causes RStudio to hang
  }
  #output plot as a .png image
  surf_plot$screenshot(filename=filename,transparent_bg = F)
}
############################################################################################################################
############################################################################################################################

#' @title Surface to volume
#'
#' @description Converts surface data to volumetric data (.nii file)
#'
#' @param surf_data A matrix object containing the surface data, either in fsaverage5 or fsaverage6 space. See SURFvextract() output format. 
#' @param filename A string object containing the desired name of the output .nii file (default is 'output.nii').
#'
#' @returns A .nii volume file
#' @examples
#' surf_to_vol(CTv, filename = 'volume.nii')
#' @importFrom reticulate import
#' @export

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
#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric studies and identifies their neuropsychological correlates
#'
#' @details The ]href{https://nimare.readthedocs.io/en/stable/index.html}{NiMARE} python module is used for the imaging decoding and is imported via the reticulate package. It also downloads the \href{https://neurosynth.org/}{neurosynth} database (~9 Mb) for correlation.
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param contrast A string object indicating whether to decode positive or negative clusters ('positive' or 'negative')
#'
#' @returns A data.frame object listing the images that correlate the most with the clusters, indicating the pearson r and names their neuropsychological correlate
#' @examples
#' decode_surf_data(CTv, 'positive')
#' @importFrom reticulate import r_to_py
#' @export

##CT image decoding
decode_surf_data=function(surf_data,contrast="positive")
{
  ##checks
    #check length of vector
    n_vert=length(surf_data)
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {stop("decoding of fsaverage6-space image is current not implemented, please resample the image to fsaverage5 space")} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) is accepted")}

    #check contrast
    if(contrast != "positive" & contrast != "negative")  {stop("contrast has to be either positive or negative")} 
  
  ##import python libraries
  interpolate=reticulate::import("brainstat.mesh.interpolate")
  discrete=reticulate::import("nimare.decode")
  nimare.dataset=reticulate::import("nimare.dataset")
  
  ##selecting contrasts
  if(contrast=="positive")
  {
    surf_data[is.na(surf_data)]=0
    surf_data[surf_data<0]=0
    surf_data[surf_data>0]=1
  } else if (contrast=="negative")
  {
    surf_data[is.na(surf_data)]=0
    surf_data[surf_data>0]=0
    surf_data[surf_data<0]=1
  }

  ##convert surf_data vector to nii image
  stat_labels=reticulate::r_to_py(surf_data)
  stat_nii = interpolate$`_surf2vol`(template, stat_labels)
  
  ##download neurosynth database if necessary 
  if(file.exists("neurosynth_dataset.pkl.gz")==F)
  {
    cat("\neurosynth_surf_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/data/neurosynth_dataset.pkl.gz",destfile = "neurosynth_dataset.pkl.gz")
  } 
  ##running the decoding procedure
  neurosynth_dset = nimare.dataset$dataset$load("neurosynth_dataset.pkl.gz")
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
