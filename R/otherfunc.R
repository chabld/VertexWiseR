## OTHER VERTEX-WISE FUNCTIONS
############################################################################################################################
############################################################################################################################
	
## permutation functions for random subject effects
  ## Paired/grouped data points are first shuffled within subjects, then these pairs/groups are shuffled between subjects
  perm_within_between=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count))>0)
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

  ## Paired/grouped data points are shuffled within subjects, order of subjects in the dataset remains unchanged
  perm_within=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
  
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count)>0))
      {
        sub.id=as.numeric(which(table(random)==count))
        for(sub in sub.id)
        {
          perm.idx[which(random==sub)]=sample(which(random==sub))
        }  
      }
    }
    return(perm.idx)
  }

  ## Paired/grouped data points are shuffled between subjects, order of data points within subjects remains unchanged.
  perm_between=function(random)
  {
    ##for groups of 2 or more (subjects with 2 or more measurements)
    perm.idx=rep(NA, length(random))
    for(count in 2:max(table(random)))
    {
      if(length(which(table(random)==count))>0)
      {
        sub.id=as.numeric(which(table(random)==count))
        if(length(sub.id)>1)
        {
          ##between group shuffling
          recode.vec=sample(sub.id)
          vec.idx=1
          for(sub in sub.id)
          {
            perm.idx[which(random==sub)]=which(random==recode.vec[vec.idx])
            vec.idx=vec.idx+1
          }   
          remove(vec.idx,recode.vec)  
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
#' @param surf_data A matrix object containing the surface data, see SURFvextract() or HIPvextract() output format
#' @param FWHM A numeric vector object containing the desired smoothing width in mm 
#'
#' @returns A matrix object with smoothed vertex-wise values
#' @examples
#' surf_data = readRDS(file = url(paste0("https://github.com",
#'"/CogBrainHealthLab/VertexWiseR/blob/main/inst/demo_data/",
#'"SPRENG_CTv.rds?raw=TRUE")))[1:3,]
#' smooth_surf(surf_data, 10)
#' @importFrom reticulate source_python
#' @export

## smooth surface data 
## FWHM input is measured in mm, which is subsequently converted into mesh units
smooth_surf=function(surf_data, FWHM)
{
  #Check if required python dependencies and libraries are  imported
  VWRrequirements()
  
  #Solves the "no visible binding for global variable" issue
  . <- mesh_smooth <- NULL 
  internalenv <- new.env()
  assign("mesh_smooth", mesh_smooth, envir = internalenv)
  
  #mesh_smooth() fails if surf_data is not a matrix object
  if (class(surf_data)[1] != 'matrix') {
  surf_data = as.matrix(surf_data) }
  
  ##source python function
  reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/smooth.py'))
  
  n_vert=ncol(surf_data)
  ##select template, set its FWHM parameter and load its edgelist file
  
  if(n_vert==20484) 
  {
    edgelist<- get('edgelistfs5') 
    FWHM=FWHM/3.5 #converting mm to mesh units
  } else if(n_vert==81924) 
  {
    edgelist<- get('edgelistfs6') 
    FWHM=FWHM/2 #converting mm to mesh units
  } else if(n_vert==14524) 
  {
    edgelist<- get('edgelistHIP') 
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
#' @importFrom igraph graph_from_data_frame
#' @importFrom igraph components

##find clusters using edgelist
getClusters=function(surf_data)
{ 
  n_vert=length(surf_data)
  
  #listing out non-zero vertices
  vert=which(surf_data!=0)
  
  #visible binding for edgelist object
  edgelist <- get("edgelist")
  
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
    com=igraph::components(igraph::graph_from_data_frame(edgelist1, directed = F))
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
#' @details The function currently works with the Desikan-Killiany-70, Schaefer-100, Schaefer-200, Glasser-360, or Destrieux-148 atlases. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{enigmatoolbox} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{nilearn.datasets.fetch_atlas_surf_destrieux}
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#' @param atlas A numeric integer object corresponding to the atlas of interest. 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148. 
#'
#' @returns A matrix object with ROI as column and corresponding average vertex-wise values as row
#' @seealso \code{\link{atlas_to_fs5}}
#' @examples
#' CTv = runif(20484,min=0, max=100)
#' fs5_to_atlas(CTv, 1)
#' @export

fs5_to_atlas=function(surf_data,atlas) 
{  
  #check length of vector
  if(length(surf_data)%%20484!=0) {stop("Length of surf_data is not a multiple of 20484")}
  
  if(missing("atlas")) {stop("Please specify an atlas number: 1=Desikan, 2=Schaefer-100, 3=Schaefer-200, 4=Glasser-360, 5=Destrieux-148")}
  
  #load atlas mapping surf_data
  ROImap <- get('ROImap_fs5')
  
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
#' @description Maps average parcellation surface values (e.g. produced with the fs5_to_atlas() function) to the fsaverage5 space
#' @details The function currently works with the Desikan-Killiany-70, Schaefer-100, Schaefer-200, Glasser-360, or Destrieux-148 atlases. ROI to vertex mapping data for 1 to 4 were obtained from the \href{https://github.com/MICA-MNI/ENIGMA/tree/master/enigmatoolbox/datasets/parcellations}{enigma toolbox} ; and data for 5 from \href{https://github.com/nilearn/nilearn/blob/a366d22e426b07166e6f8ce1b7ac6eb732c88155/nilearn/datasets/atlas.py}{nilearn.datasets.fetch_atlas_surf_destrieux} . atlas_to_fs5() will automatically detect the atlas based on the number of columns.
#'
#' @param parcel_data A matrix object containing average surface measures for each region of interest, see the fs5_to_atlas() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @seealso \code{\link{fs5_to_atlas}}
#' @examples
#' atlas_data = runif(100,min=0, max=100)
#' atlas_to_fs5(atlas_data)
#' @export


atlas_to_fs5=function(parcel_data) 
  {
  
    #load atlas mapping surface data
    ROImap <- get('ROImap_fs5')
    
    if (ncol(parcel_data) == 70) {atlas=1} 
    else if (ncol(parcel_data) == 100) {atlas=2} 
    else if (ncol(parcel_data) == 200) {atlas=3} 
    else if (ncol(parcel_data) == 360) {atlas=4} 
    else if (ncol(parcel_data) == 148) {atlas=5} 
    else { stop('The function could not identify what atlas your data was parcellated with, based on the number of columns (parcels). The function currently works with the Desikan-Killiany-70, Schaefer-100, Schaefer-200, Glasser-360, or Destrieux-148 atlases.')}
    
    #init variables
    nregions=max(ROImap[[1]][,atlas])
    fs5_dat=rep(NA,20484)
  
    #mapping atlas label to fsaverage5 space
    for (region in 1:nregions)  {fs5_dat[which(ROImap[[1]][,atlas]==region)]=parcel_data[region]}
    return(as.numeric(fs5_dat))
  }
############################################################################################################################
############################################################################################################################


#' @title fsaverage5 to fsaverage6
#'
#' @description Remaps vertex-wise surface data in fsaverage5 space to fsaverage6 space using the nearest neighbor approach 
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage6 space
#' @seealso \code{\link{fs6_to_fs5}}
#' @examples
#' CTv = runif(20484,min=0, max=100)
#' fs5_to_fs6(CTv)
#' @export

#convert between fsaverage5 and fsaverage6 spacing
fs5_to_fs6=function(surf_data)
{
  #check length of vector
  if(length(surf_data)%%20484!=0) {stop("Length of surf_data is not a multiple of 20484")}
  
  #load atlas mapping surf_data
  fs6_to_fs5_map <- get('fs6_to_fs5_map')
  
  #mapping fsaverage5 to fsaverage6 space if surf_data is a vector length of 20484
  if(length(surf_data)==20484) {surf_data.fs6=surf_data[fs6_to_fs5_map]} 
  #mapping fsaverage5 to fsaverage6 space if surf_data is a Nx20484 matrix
  else {surf_data.fs6=surf_data[,fs6_to_fs5_map]}
  return(surf_data.fs6)
}

#' @title fsaverage6 to fsaverage5
#'
#' @description Remaps vertex-wise surface data in fsaverage6 space to fsaverage5 space using the nearest neighbor approach
#'
#' @param surf_data A matrix object containing the surface data, see SURFvextract() output format. 
#'
#' @returns A matrix object containing vertex-wise surface data mapped in fsaverage5 space
#' @seealso \code{\link{fs5_to_fs6}}
#' @examples
#' CTv = runif(81924,min=0, max=100);
#' fs6_to_fs5(CTv)
#' 
#' @export

fs6_to_fs5=function(surf_data)
{
  #check length of vector
  if(length(surf_data)%%81924!=0) {stop("Length of surf_data is not a multiple of 81924")}
  
  #load atlas mapping surf_data
  fs6_to_fs5_map <- get('fs6_to_fs5_map')
  
  if(length(surf_data)==81924) #mapping fsaverage6 to fsaverage5 space if surf_data is a Nx81924 matrix
  {
    surf_data=matrix(surf_data,ncol=81924,nrow=1)  
    surf_data.fs5=matrix(NA,ncol=20484,nrow=1)
    
    for (vert in 1:20484)  {surf_data.fs5[vert]=mean(surf_data[fs6_to_fs5_map==vert],na.rm = T)} 
  } else #mapping fsaverage6 to fsaverage5 space if surf_data is a Nx20484 matrix
  {
    surf_data.fs5=matrix(NA,ncol=20484,nrow=NROW(surf_data))
    for (vert in 1:20484)  {surf_data.fs5[,vert]=rowMeans(surf_data[,fs6_to_fs5_map==vert],na.rm = T)} 
  }
  return(surf_data.fs5)
}
############################################################################################################################
############################################################################################################################

#' @title Surface plotter
#'
#' @description Plots surface data in a grid with one or multiple rows in a .png file
#'
#' @param surf_data  A numeric vector (length of V) or a matrix (N rows x V columns), where N is the number of subplots, and V is the number of vertices. It can be the output from SURFvextract() as well as masks or vertex-wise results outputted by analyses functions.
#' @param filename A string object containing the desired name of the output .png file.
#' @param title A string object for setting the title in the plot. Default is none. For titles that too long to be fully displayed within the plot, we recommend splitting them into multiple lines by inserting "\\n".
#' @param surface A string object containing the name of the type of cortical surface background rendered. Possible options include "white", "smoothwm","pial" and "inflated" (default). The surface parameter is ignored for hippocampal surface data.
#' @param cmap A string object containing the colormap for the plot. Options are listed in the \href{https://matplotlib.org/stable/gallery/color/colormap_reference.html}{Matplotlib plotting library}. 
#' 
#' Default cmap is set to `"Reds"` for positive values, `"Blues_r"` for negative values and `"RdBu"` when both positive and negative values exist. 
#' @param limits A combined pair of numeric vector composed of the lower and upper color scale limits of the plot. If the limits are specified, the same limits will be applied to all subplots. When left unspecified, the same symmetrical limits c(-max(abs(surf_dat),max(abs(surf_dat))) will be used for all subplots. If set to NULL, each subplot will have its own limits corresponding to their min and max values
#' @param colorbar A logical object stating whether to include a color bar in the plot or not (default is TRUE).
#' @param size A combined pair of numeric vector indicating the image dimensions (width and height in pixels). Default is c(1920,400) for whole-brain surface and c(400,200) for hippocampal surface.
#' @param zoom A numeric value for adjusting the level of zoom on the figures. Default is 1.25 for whole-brain surface and 1.20 for hippocampal surface.
#'
#' @returns Outputs the plot as a .png image
#' @examples
#'if(interactive()){
#'results = runif(20484,min=0, max=1)
#'plot_surf(results, filename='output.png',title = 
#' 'Cortical thickness', surface = 'inflated', cmap = 'Blues')
#'}
#' @importFrom reticulate tuple import np_array source_python
#' @export

plot_surf=function(surf_data, filename, title="",surface="inflated",cmap,limits, colorbar=T, size, zoom)
{
  #Check if required python dependencies and libraries are  imported
  VWRrequirements()
  
  if (missing("filename")) {
    cat('No filename argument was given. The plot will be saved as "plot.png" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/plot.png')
  }
  
  #format title for single row
  if(is.null(nrow(surf_data)))
  {
    title=list('left'=list(title))
    rows=1
    surf_data=as.numeric(surf_data)
  } else {rows=nrow(surf_data)}

  #in a multi-row data scenario: insert a dummy title if title is missing  or repeat the title nrow times
  if(rows>1) 
    {
       if(missing("title")) {title=rep(NULL,rows)}
       else if (missing("title")) {title=rep(title,rows)}
    }
	  
  #check length of vector
  n_vert=length(surf_data)
  if(n_vert%%20484==0) {template="fsaverage5"}
  else if (n_vert%%81924==0) {template="fsaverage6"} 
  else if (n_vert%%14524!=0) {stop("surf_data vector should only contain 20484 (fsaverage5), 81924 (fsaverage6) or 14524 (hippocampal vertices) columns")}

  #if cmap is missing, select cmaps depending on whether the image contains positive only or negative only values
  if(missing("cmap"))
  {
    if(range(surf_data,na.rm = T)[1]>=0)  {cmap="Reds"}
    else if (range(surf_data,na.rm = T)[2]<=0)  {cmap="Blues_r"}
    else  {cmap="RdBu_r"}  
  }
  
  #setting color scale limits
   maxlimit=max(abs(range(surf_data,na.rm = T)))
    if(missing("limits")) 
    {
      if(range(surf_data,na.rm = T)[1]>=0) {limits=reticulate::tuple(0,range(surf_data,na.rm = T)[2])} ##if image contains all positive values
      else if(range(surf_data,na.rm = T)[2]<=0) {limits=reticulate::tuple(range(surf_data,na.rm = T)[1],0)} ##if image contains all negative values
      else {limits=reticulate::tuple(-maxlimit,maxlimit)} ##symmetrical limits will be used if image contains both positive and negative values
    } else {
      ##user specified limits
      if(!is.null(limits))
      {
        limits=reticulate::tuple(limits[1],limits[2])  
      }   
    }
  
  if(n_vert%%14524!=0)
  {
  ##cortical surface fplots
    #import python libraries
    brainstat.datasets=reticulate::import("brainstat.datasets")  
    brainspace.plotting=reticulate::import("brainspace.plotting")  
    
    #loading fsaverage surface
    left=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[1]
    right=brainstat.datasets$fetch_template_surface(template, join=F, layer=surface)[2]
    
    #default cortical size and zoom parametes
    if(missing("size")) { size=c(1920,rows*400)}
    if(missing("zoom")) { zoom=1.25 }
    
    surf_plot=brainspace.plotting$plot_hemispheres(left[[1]], right[[1]],  array_name=reticulate::np_array(surf_data),cmap=cmap, 
                                                size=reticulate::tuple(as.integer(size)),nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),
                                                return_plotter=T,background=reticulate::tuple(as.integer(c(1,1,1))),zoom=zoom,color_range=limits,
                                                label_text=title,interactive=F, color_bar=colorbar,  transparent_bg=FALSE)  ##disabling interactive mode because this causes RStudio to hang
  } else
  {
    #Solves the "no visible binding for global variable" issue
    . <- surfplot_canonical_foldunfold  <- NULL 
    internalenv <- new.env()
    assign("surfplot_canonical_foldunfold", surfplot_canonical_foldunfold, envir = internalenv)
    
  ##hippocampal plots
    #import python libraries
    reticulate::source_python(paste0(system.file(package='VertexWiseR'),'/python/hipp_plot.py'))
    
    #default hippocampal size and zoom parametes
    if(missing("size")) { size=c(400,200)}
    if(missing("zoom")) { zoom=1.2 }

    #reshaping surf_data into a 7262 x 2 x N array
    if(is.null(nrow(surf_data)))  {surf_data=cbind(surf_data[1:7262],surf_data[7263:14524])} #if N=1
    else  
      {
        surf_data.3d=array(NA,c(7262,2,nrow(surf_data))) #if N>1
        for (row in 1:nrow(surf_data))  {surf_data.3d[,,row]=cbind(surf_data[row,1:7262],surf_data[row,7263:14524])}
        surf_data=surf_data.3d
      }
    
    surf_plot=surfplot_canonical_foldunfold(surf_data,hipdat =get('hip_points_cells'),color_bar=colorbar,share="row",nan_color=reticulate::tuple(0.7, 0.7, 0.7, 1),size=as.integer(size), zoom=zoom,
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
#' if(interactive()){
#' CTv = runif(20484,min=0, max=100);
#' surf_to_vol(CTv, filename = 'volume.nii')
#' }
#' @importFrom reticulate import
#' @export

##converting surface to volumetric data and exporting it as a .nii file

surf_to_vol=function(surf_data, filename)
  {
  #Check if required python dependencies and libraries are  imported
  VWRrequirements()
  
  if (missing("filename")) {
    cat('No filename argument was given. The volume will be saved as "vol.nii" in R temporary directory (tempdir()).\n')
    filename=paste0(tempdir(),'/vol.nii')
  }
  
  #check length of vector
    n_vert=ncol(surf_data)
    if(n_vert==20484) {template="fsaverage5"}
    else if (n_vert==81924) {template="fsaverage6"} 
    else {stop("Only an surf_data vector with a length of 20484 (fsaverage5) or 81924 (fsaverage6) is accepted")}
  
  #load python libraries
    interpolate=reticulate::import("brainstat.mesh.interpolate")
    nibabel=reticulate::import("nibabel")

  #convert and export .nii file
    stat_nii = interpolate$`_surf2vol`(template, surf_data)
    nibabel$save(stat_nii,filename)
    cat(filename)
  }

############################################################################################################################
############################################################################################################################
#' @title Decode surface data
#'
#' @description Correlates the significant clusters of an earlier vertex-wise analysis with a database of task-based fMRI and voxel-based morphometric statistical maps and associate them with relevant key words
#'
#' @details The \href{https://nimare.readthedocs.io/en/stable/index.html}{NiMARE} python module is used for the imaging decoding and is imported via the reticulate package. The function also downloads the \href{https://github.com/neurosynth/neurosynth-data}{neurosynth} database in the package's inst/extdata direcotry (~8 Mb) for the analysis.
#'
#' @param surf_data a numeric vector with a length of 20484
#' @param contrast A string object indicating whether to decode the positive or negative mask ('positive' or 'negative')
#'
#' @returns A data.frame object listing the keywords and their Pearson's R values
#' @examples
#' CTv = rbinom(20484, 1, 0.001) 
#' decode_surf_data(CTv, 'positive')
#' @importFrom reticulate import r_to_py
#' @importFrom RCurl url.exists
#' @export

##CT image decoding
decode_surf_data=function(surf_data,contrast="positive")
{
  #Check if required python dependencies and libraries are  imported
  VWRrequirements()
  
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
  if(file.exists(system.file('extdata','neurosynth_dataset.pkl.gz', package='VertexWiseR'))==F)
  {
    cat("\nneurosynth_dataset.pkl.gz is not detected in the current working directory. The neurosynth database will be downloaded\n")
    
    #Check if URL works and avoid returning error but only print message as requested by CRAN:
    url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz"
    if(RCurl::url.exists(url)) {
        download.file(url="https://raw.githubusercontent.com/CogBrainHealthLab/VertexWiseR/main/inst/extdata/neurosynth_dataset.pkl.gz",destfile = paste0(system.file(package='VertexWiseR'),'/extdata/neurosynth_dataset.pkl.gz'))
    } else { 
      return("The neurosynth database (neurosynth_dataset.pkl.gz) could not be downloaded from the github VertexWiseR directory. Please check your internet connection or visit https://github.com/CogBrainHealthLab/VertexWiseR/tree/main/inst/extdata to download the object.") #ends function
    } 
  }
  
  ##running the decoding procedure
  neurosynth_dset = nimare.dataset$Dataset$load(system.file("extdata/neurosynth_dataset.pkl.gz", package='VertexWiseR'))
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
#' @title VertexWiseR system requirements installation
#'
#' @description Helps the user install all system requirements for VertexWiseR functions to work (miniconda, brainstat toolbox and libraries). If they are installed already nothing will be overwritten. 
#'
#' @details VertexWiseR imports and makes use of the R package reticulate. reticulate is a package that allows R to borrow or translate python functions into R. Using reticulate, the package calls functions from the brainstat python module. For reticulate to work properly with VertexWiseR, the latest version of miniconda needs to be installed with it — miniconda is a lightweight version of python, specifically for use within RStudio. Likewise, analyses of cortical surface require fsaverage templates as imported by brainstat.
#'
#' @examples
#' VWRfirstrun()
#' @importFrom reticulate conda_binary py_module_available
#' @importFrom fs path_home
#' @export

VWRfirstrun=function() 
{
  cat('Checking for VertexWiseR system requirements ... \n')
  if (is(tryCatch(reticulate::conda_binary(), error=function(e) e))[1] == 'simpleError')
  {
    cat('Miniconda could not be found in the environment. \n')
    prompt = utils::menu(c("Yes", "No"), title=" Do you want miniconda to be installed now?")
    if (prompt==1){reticulate::install_miniconda()} else {stop('VertexWiseR will not work properly without miniconda. reticulate::conda_list() should detect it on your system.')}

    
  } else if(!reticulate::py_module_available("brainstat")) 
  {
    cat('Brainstat could not be found in the environment. \n')
    prompt = utils::menu(c("Yes", "No"), title=" Do you want brainstat to be installed now?")
    if (prompt==1){reticulate::py_install("brainstat",pip=TRUE)} else {stop('VertexWiseR will not work properly without brainstat.')}
    
    
  } else if (!file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5'))) 
  {     
    cat('VertexWiseR could not find brainstat fsaverage5 templates in $home/brainstat_data/.')  
    prompt = utils::menu(c("Yes", "No"), title=" Do you want the templates to be downloaded now?")
    if (prompt==1){    
      brainstat.datasets.base=reticulate::import("brainstat.datasets.base")
      brainstat.datasets.base$fetch_template_surface("fsaverage5")
    } else {stop('VertexWiseR will not be able to analyse fsaverage5 data without the brainstat templates.')}
    
  } else if (!file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6')))
  {
    cat('VertexWiseR could not find brainstat fsaverage6 templates in $home/brainstat_data/..')
    prompt = utils::menu(c("Yes", "No"), title=" Do you want the templates to be downloaded now?")
    if (prompt==1){    
      brainstat.datasets.base=reticulate::import("brainstat.datasets.base")
      brainstat.datasets.base$fetch_template_surface("fsaverage6")
    } else {stop('VertexWiseR will not be able to analyse fsaverage6 data without the brainstat templates.')}
  } else
  {
    cat('All system requirements are installed.')
  }
}

############################################################################################################################
############################################################################################################################
#This function checks that any external system requirement is fulfilled before using the functions that need python libraries
#' @importFrom utils menu

VWRrequirements=function() 
{
  if (is(tryCatch(reticulate::conda_binary(), error=function(e) e))[1] == 'simpleError')
  {
    cat('VertexWiseR requires miniconda to be installed on your system for this function to run.')
    prompt = utils::menu(c("Yes", "No"), title=" VWRfirstrun() let you install all VertexWiseR requirements automatically. Do you want to run it?")
    if (prompt==1){VWRfirstrun()} else {stop('VertexWiseR cannot run this function without miniconda. reticulate::conda_list() should detect it on your system.')}
    
  } else if(!reticulate::py_module_available("brainstat")) 
  {
    cat('VertexWiseR requires brainstat to be installed on your system for this function to run.')
    prompt = utils::menu(c("Yes", "No"), title=" VWRfirstrun() let you install all VertexWiseR requirements automatically. Do you want to run it?")
    if (prompt==1){VWRfirstrun()} else {stop('VertexWiseR cannot run this function without brainstat.')}
    
  } else if (!file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage5')))
  {
    cat('VertexWiseR requires fsaverage5 templates to be imported by brainstat in $home/brainstat_data/ for this function to run.')
    prompt = utils::menu(c("Yes", "No"), title=" VWRfirstrun() let you install all VertexWiseR requirements automatically. Do you want to run it?")
    if (prompt==1){VWRfirstrun()} else {stop('VertexWiseR cannot run this function without brainstat fsaverage templates.')}
      
} else if  (!file.exists(paste0(fs::path_home(),'/brainstat_data/surface_data/tpl-fsaverage/fsaverage6')))
{
  cat('VertexWiseR requires fsaverage6 templates to be imported by brainstat in $home/brainstat_data/ for this function to run.')
  prompt = utils::menu(c("Yes", "No"), title=" VWRfirstrun() let you install all VertexWiseR requirements automatically. Do you want to run it?")
  if (prompt==1){VWRfirstrun()} else {stop('VertexWiseR cannot run this function without brainstat fsaverage6 templates.')} 
}
}