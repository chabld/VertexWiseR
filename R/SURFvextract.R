#' @title SURFvextract
#'
#' @description Extracts whole-brain vertex-wise surface-based measures for each subject in a Freesurfer output subjects directory, resamples the data to a common surface template, and stores it into a single matrix. This function requires the Freesurfer environment to be preset.
#' @details The function runs system shell commands that will produce in the set subjects directory: 1) a sorted list of subjects "sublist.txt"; 2) a link file to the selected surface fsaverage template. 3) left and right hemisphere .mgh maps outputted by FreeSurfer's mris_preproc.
#'
#' @param sdirpath A string object containing the path to the Freesurfer subjects directory.
#' @param filename A string object containing the desired name of the output RDS file.
#' @param template A string object containing the name of surface template (available: 'fsaverage5', 'fsaverage6'). Default is fsaverage5.
#' @param measure A string object containing the name of the measure of interest. Options are thickness, curv, sulc, area. Default is thickness.
#' @param subj_ID A logical object stating whether to include subject IDs (folder names in the subjects directory) as a first column to the output matrix. Default is TRUE.
#'
#' @returns A .rds file containing a matrix object. The matrix can be used independently by VertexWiseR to compute statistical analyses. Each row corresponds to a subject's left and right hemisphere vertex-wise values
#' @examples
#' SURFvextract(sdirpath = "myfreesurfer_output_path/subjects_directory/", filename = "CTv", template="fsaverage5", measure = "curv") 
#' @importFrom freesurferformats read.fs.mgh
#' @export

SURFvextract=function(sdirpath, filename, template='fsaverage5', measure = 'thickness', subj_ID = T) 
{ 
#finds specifically subject folders in the directory (checks if a surf folder is present) and stores their ordered IDs in a list  
Sys.setenv(SUBJECTS_DIR=sdirpath)
system("find $SUBJECTS_DIR -maxdepth 1 -type d -exec test -e '{}/surf' \\; -exec basename {} > $SUBJECTS_DIR/sublist.txt \\;");
system("sort -n $SUBJECTS_DIR/sublist.txt -o $SUBJECTS_DIR/sublist.txt");
  
#Calls Freesurfer to extract vertex-wise thickness data from the sample and resample it to the fsaverage5 common-space surface; and concatenate it into mgh files
system(paste0("ln -s $FREESURFER_HOME/subjects/", template, " -t $SUBJECTS_DIR \n
       mris_preproc --f $SUBJECTS_DIR/sublist.txt --target fsaverage5 --hemi lh --meas ", measure, " --out $SUBJECTS_DIR/lh.mgh \n 
       mris_preproc --f $SUBJECTS_DIR/sublist.txt --target fsaverage5 --hemi rh --meas ", measure, " --out $SUBJECTS_DIR/rh.mgh"));

#Reads mgh files to stores and assign the thickness values to each subject in a matrix object usable by VertexWiseR
if (subj_ID == T) 
{
sublist = read.delim(paste0(sdirpath,"/sublist.txt"));
SURFdata= t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
SURFdata=cbind(sublist,SURFdata); 
} else {
SURFdata=t(rbind(drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"lh.mgh"))),drop(freesurferformats::read.fs.mgh(paste0(sdirpath,"rh.mgh")))));
}

saveRDS(SURFdata, file=filename)
}