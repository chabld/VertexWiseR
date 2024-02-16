#' @title HIPvextract
#'
#' @description Extracts hippocampal vertex-wise surface-based measures for each subject in the hippunfold subjects directory, and stores it into a single matrix.
#' @details The function searches for the hippocampal surface data by listing out files with certain suffixes, and then read these files and save both the left and right hippocampal vertex data for each subject as rows in a N x 14524 data matrix. 
#'
#' @param sdirpath A string object containing the path to the hipunfold subjects directory.
#' @param filename A string object containing the desired name of the output RDS file.
#' @param measure A string object containing the name of the measure of interest. Options are thickness and area. Default is thickness.
#' @param subj_ID A logical object stating whether to return a list object containing both subject ID and data matrix.
#'
#' @returns A .rds file with a list containing the list of subject IDs (first element) and a surface data matrix object (second element), or a data matrix object. The matrix can be used independently by VertexWiseR statistical analysis functions. Each row corresponds to a subject (in the same order as 1) and contains the left to right hemispheres' vertex-wise values.
#' @examples
#' HIPvextract(sdirpath = "./", filename = "hip_data.RDS", measure = "thickness") 
#' @importFrom gifti readgii
#' @export

HIPvextract=function(sdirpath, filename, measure="thickness", subj_ID = T)
{
  setwd(sdirpath)
  
  ## get filelists and subject lists
  lh.filelist=list.files(pattern=paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""), recursive=T)
  rh.filelist=gsub(paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""),
                   paste("_hemi-R_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""),
                   lh.filelist)
  sublist=gsub(paste("_hemi-L_space-T1w_den-0p5mm_label-hipp_",measure,".shape.gii",sep=""), 
               "",
               basename(lh.filelist))

  ##read data and save data for each subject as rows in a data matrix
  hip_dat=matrix(NA, nrow=NROW(sublist), ncol=14524)
  
  for (sub in 1:NROW(sublist))
  {
    lh=gifti::readgii(lh.filelist[sub])
    rh=gifti::readgii(rh.filelist[sub])
    hip_dat[sub,]=c(lh$data$normal,rh$data$normal)
  }

  ##output file depending on subj_ID==T
  hip_dat=hip_dat[order(sublist),]
    if(subj_ID==T)
    {
    saveRDS(list(sublist,hip_dat), file=filename)
    } else
    {
    saveRDS(hip_dat, file=filename)
    }
}

