#' download_BLUPF90
#'
#' Description:
#' Download all BLUPF90 software from their official repository (http://nce.ads.uga.edu/html/projects/programs/),
#' and saves them to the destination folder.
#'
#' Usage:
#' \code{download_BLUPF90(dest_folder = NULL)}
#'
#' Arguments:
#' @param dest_folder (optional) A character string specifying the destination folder where the BLUPF90 files will be saved.
#'
#' @return
#' Details:
#' The \code{download_BLUPF90} function downloads BLUPF90 software files from the appropriate URL based on the operating system. The function identifies the operating system (Linux, Mac_OSX, or Windows) and constructs the URL accordingly. It retrieves the list of available BLUPF90 files from the URL and compares them to the files in the local destination folder. It then downloads the missing files from the URL to the specified destination folder.
#'
#' @import RCurl
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
#' Please ensure that you have proper permissions to access the destination folder and that you have an active internet connection during the execution of this function.
#'
download_BLUPF90 <- function(dest_folder=paste0(.libPaths()[1],"/blupf90")){
  #download_BLUPF90(dest_folder = paste0(.libPaths()[1],"/blupf90"))

  if(gsub("([0-9]|\\.)","",version$os)=="linux-gnu"){
    S_OP <- "Linux"
  } else if(gsub("([0-9]|\\.)","",version$os)=="darwin"){
    S_OP <- "Mac_OSX"}else{
      S_OP <- "Windows"
    }

  url_blupf90 <- paste0("http://nce.ads.uga.edu/html/projects/programs/",S_OP,"/64bit/")
  alltext <- RCurl::getURL(url_blupf90)
  pattern <- "([A-z]+)(f90|f90[+])"
  m <- base::gregexec(pattern, alltext)
  list_blupf90 <- base::regmatches(alltext, m)[[1]][1,]|>unique()|>data.frame()
  colnames(list_blupf90) <- c("stw")
  if(S_OP=="Windows"){
    list_blupf90$stw <- paste0(list_blupf90$stw,".exe")
  }

  d_f <- function(stw_BLUPF90,dest_folder=dest_folder){
    # Specify destination where file should be saved
    if(is.null(dest_folder)==T){
      destfile <- paste0(getwd(),"/",stw_BLUPF90)
      dest_folder <- getwd()
    }else{
      destfile <- paste0(dest_folder,"/",stw_BLUPF90)
    }
    dir.create(dest_folder)
    #dir.create(paste0(.libPaths()[1]),"/blupf90")
    # Apply download.file function in R
    if(S_OP=="Windows"){
      download.file(paste0(url_blupf90,stw_BLUPF90),destfile,mode="wb")
    }else{
      download.file(paste0(url_blupf90,stw_BLUPF90),destfile)
    }
    if(S_OP!="Windows"){
      Sys.chmod(paste0(dest_folder,"/",stw_BLUPF90),  # Apply Sys.chmod function
                mode = "0777")}
  }

  stw_local <- data.frame(stw=list.files(getwd()))
  if(is.null(dest_folder)==T){
    stw_local <- data.frame(stw=list.files(getwd()))
  }else{
    stw_local <- data.frame(stw=list.files(dest_folder))
  }

  diff_stw <- data.frame(stw=setdiff(list_blupf90$stw,stw_local$stw))
  if(nrow(diff_stw)==0){print("You have downloaded all BLUPF90 softwares")
  }else{
    apply(diff_stw, 1, d_f,dest_folder=dest_folder)}
}




