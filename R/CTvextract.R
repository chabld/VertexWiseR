#' @title CTvextract
#'
#' @description Extracts whole-brain vertex-wise cortical thickness values for each subject in a Freesurfer output subjects directory, resamples the data to a common surface template, and stores it into a single matrix. This function requires the Freesurfer environment to be preset; and the 'freesurferformats' R package.
#'
#' @param sdirpath A string object containing the path to the Freesurfer output subjects directory.
#' @param filename A string object containing the desired name of the output RDS file.
#' @param template A string object containing the name of surface template (available: 'fsaverage5', 'fsaverage6'). Default is fsaverage5.
#' @param measure A string object containing the name of the measure of interest. Options are thickness, curv, sulc, area. Default is thickness.
#'
#' @return A matrix object that can be used independently by VertexWiseR to compute statistical analyses. Each row corresponds to a subject's left and right hemisphere vertex-wise cortical thickness values
#' @examples
#' CTvextract('myfreesurfer_output_path/subjects_directory/', template='fsaverage5') 
#' @importFrom freesurferformats read.fs.mgh
#' @export

CTvextract=function(sdirpath, filename, template='fsaverage5', measure = 'thickness') 
{ 
#finds specifically subject folders in the directory (checks if a surf folder is present) and stores their ordered IDs in a list  
system(paste0("export SUBJECTS_DIR=", sdirpath))
system("find $SUBJECTS_DIR -maxdepth 1 -type d -exec test -e '{}/surf' \\; -exec basename {} > ./sublist.txt \\; \n sort -n ./sublist.txt -o ./sublist.txt")
  
#Calls Freesurfer to extract vertex-wise thickness data from the sample and resample it to the fsaverage5 common-space surface; and concatenate it into mgh files
system(paste0("ln -s $FREESURFER_HOME/subjects/", template, " -t $SUBJECTS_DIR \n
       mris_preproc --f ./sublist.txt --target fsaverage5 --hemi lh --meas", measure, " --out lh.mgh \n 
       mris_preproc --f ./sublist.txt --target fsaverage5 --hemi rh --meas", measure, " --out rh.mgh"))

#Reads mgh files to stores and assign the thickness values to each subject in a matrix object usable by VertexWiseR
CT=t(rbind(drop(freesurferformats::read.fs.mgh("lh.mgh")),drop(freesurferformats::read.fs.mgh("rh.mgh"))))
saveRDS(CT, file=paste0(filename,".rds"))
}