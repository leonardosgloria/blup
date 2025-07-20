#' download_BLUPF90
#'
#' Description:
#' This function downloads all the BLUPF90 software files from their official repository
#' (http://nce.ads.uga.edu/html/projects/programs/) and saves them to the specified
#' destination folder. The BLUPF90 software is used for statistical analysis
#' and genetic evaluation of animal and plant breeding data.
#' By default, the function will save the software files in the R user folder,
#' but you can provide a different destination folder using the dest_folder parameter.
#' If the update parameter is set to TRUE, the function will replace any existing
#' files in the destination folder with the latest versions from the repository.
#' If update is set to FALSE, the function will skip the download for files that already exist in the destination folder.
#'
#'
#'
#' Usage:
#' \code{download_BLUPF90(dest_folder = NULL)}
#'
#' Arguments:
#' @param dest_folder (optional) A character string specifying the destination folder where the BLUPF90 files will be saved.
#' @param update (optional): Specifies whether to update the existing BLUPF90 software files if they already exist in the destination folder. If set to TRUE, the function will download and replace any existing files. If set to FALSE, the function will skip the download if the files already exist. The default value is FALSE.
#'
#' @return
#' Details:
#' The \code{download_BLUPF90} function downloads BLUPF90 software files from the appropriate URL based on the operating system. The function identifies the operating system (Linux, Mac_OSX, or Windows) and constructs the URL accordingly. It retrieves the list of available BLUPF90 files from the URL and compares them to the files in the local destination folder. It then downloads the missing files from the URL to the specified destination folder.
#'
#' @import RCurl httr
#'
#' @examples
#' \code{
#' # Download BLUPF90 files to the default destination folder
#' download_BLUPF90()
#'
#' # Download BLUPF90 files to a specific destination folder
#' download_BLUPF90(dest_folder = "~/blupf90_files")
#' }
#'
#' Note:
#' The \code{download_BLUPF90} function requires the \code{RCurl} package to be installed.
#'
#' Please ensure that you have proper permissions to access the destination folder and
#' that you have an active internet connection during the execution of this function.
#'
download_BLUPF90 <- function(dest_folder=NULL,update=F){
 
  if(is.null(dest_folder)==T){
    dest_folder=paste0(.libPaths()[1],"/blupf90")
  }
  if(gsub("([0-9]|\\.)","",version$os)=="linux-gnu"){
    S_OP <- "Linux"
  } else if(gsub("([0-9]|\\.)","",version$os)=="darwin"){
    S_OP <- "Mac_OSX"}else{
      S_OP <- "Windows"
    }
  
  url_blupf90 <- paste0("https://nce.ads.uga.edu/html/projects/programs/",S_OP,"/64bit/")
  alltext <- RCurl::getURL(url_blupf90,ssl.verifypeer = FALSE)
  pattern <- "([A-z]+)(f90|f90[+])"
  m <- base::gregexec(pattern, alltext)
  list_blupf90 <- base::regmatches(alltext, m)[[1]][1,]|>unique()|>data.frame()
  colnames(list_blupf90) <- c("stw")
  if(S_OP=="Windows"){
    list_blupf90$stw <- paste0(list_blupf90$stw,".exe")
  }
  
  d_f <- function(stw_BLUPF90,dest_folder=dest_folder){
       
    if(is.null(dest_folder)==T){
      destfile <- paste0(getwd(),"/",stw_BLUPF90)
      dest_folder <- getwd()
    }else{
      destfile <- paste0(dest_folder,"/",stw_BLUPF90)
    }
    dir.create(dirname(dest_folder), recursive = TRUE, showWarnings = FALSE)
    
    # Apply download.file function in R
    if(S_OP=="Windows"){
      download.file(paste0(url_blupf90,stw_BLUPF90),destfile,mode="wb",ssl_verifypeer = FALSE)
    }else{
      # Download using httr and bypass SSL check
      res <- httr::GET(paste0(url_blupf90,stw_BLUPF90), httr::config(ssl_verifypeer = FALSE))
      
      # Write to disk
      writeBin(httr::content(res, "raw"), destfile)
    }
    if(S_OP!="Windows"){
      Sys.chmod(paste0(dest_folder,"/",stw_BLUPF90),  # Apply Sys.chmod function
                mode = "0755")}
  }
  
  stw_local <- data.frame(stw=list.files(getwd()))
  if(is.null(dest_folder)==T){
    stw_local <- data.frame(stw=list.files(getwd()))
  }else{
    stw_local <- data.frame(stw=list.files(dest_folder))
  }
  
  diff_stw <- data.frame(stw=setdiff(list_blupf90$stw,stw_local$stw))
  if(nrow(diff_stw)==0 & update==F){print("You have downloaded all BLUPF90 softwares")
  }else{
    apply(list_blupf90, 1, d_f,dest_folder=dest_folder)
    }
}

