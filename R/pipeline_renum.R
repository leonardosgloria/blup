#' @title blup
#' @author Leonardo S. Gloria, Luiz F. Brito, Alencar Xavier, Keith A. Cherkauer, Victor B. Pedrosa, Katy Martin Rainey
#' @description
#'  blup is a standard analytical function in R for single and multi-trait methods and single-step analysis in the breeding pipeline using
#' BLUPF90 software
#'
#' @param datarenum The dataset used for the analysis.
#' @param formula The formula specifying the model for the analysis.
#' @param fields_output (optional) A character vector specifying the fields to be included in the output.
#' @param weights_object (optional) The weights object to be used in the analysis.
#' @param residual_start The starting value for residual matrix.
#' @param VCA_RRM (optional) The genetic variance matrix for analysis using random regression model.
#' @param VCD_RRM (optional) The random effect (non-genetic effect) variance matrix for analysis using random regression model.
#' @param VCA The genetic variance matrix for analysis.
#' @param VCD The random effect (non-genetic effect) variance matrix for analysis.
#' @param ped_name The name of the pedigree file.
#' @param PED_DEPTH The depth of the pedigree.
#' @param genotype_file The genotype file used in the analysis.
#' @param missing_values The value representing missing data in the dataset. Default value is -99.
#' @param RRM_option The random regression option option for the analysis.
#' @param het_res_variance (optional) The heterogeneity of residual variance parameters.
#' @param fit_option The fit options for the analysis.
#' @param run_model A list specifying the method and options for running the model.
#' @param cross_validation A list specifying the cross-validation options.
#' @param keep_files (optional) Logical value indicating whether to keep intermediate files generated during the analysis.
#' @param extra.option.blup (optional) Additional options specific to the BLUP analysis. See more options in: http://nce.ads.uga.edu/wiki/doku.php?id=readme.blupf90plus
#'
#' @details
#' The `blup` function performs Best Linear Unbiased Prediction (BLUP) analysis on the given dataset using the specified model formula.
#' It estimates genetic parameters, such as additive genetic variance and non-genetic random effects, using the provided input parameters.
#'  Various options are available to customize the analysis, including random regression option using Legendre polynomials,
#'  residual heritability, pedigree information, genotype data, single-step GBLUP.
#'
#' @examples
#' \dontrun{
#' # Example usage of blup function
#' blup(datarenum, formula, fields_output = c("field1", "field2"), weights_object = weights,
#'      residual_start = 1, VCA_RRM = VCA_matrix, VCD_RRM = VCD_matrix, VCA = 0.1, VCD = 0.1,
#'      ped_name = "pedigree.txt", PED_DEPTH = 3, genotype_file = "genotypes.txt", missing_values = -999,
#'      RRM_option = "option1", het_res_variance = 0.2, fit_option = "option2",
#'      run_model = list(method = "blupf90+", dense = FALSE, Gibbs_option = NULL, slurm_option = NULL,
#'                       n_threads = 1),
#'      cross_validation = list(group_name_cv, LOO_cv, prop, nrep, eff_name, RRM_option,
#'                             keep_folder = FALSE, parameter_file = "renf90.par", logname = "blupf90.log"),
#'      keep_files = TRUE, extra.option.blup = NULL)
#' }
#'
#' @references
#' \itemize{
#'   \item Lourenco, D., S. Tsuruta, I. Aguilar, Y. Masuda, M. Bermann, A. Legarra, and I. Misztal. "Recent updates in the BLUPF90 software suite." In Proceedings of 12th World Congress on Genetics Applied to Livestock Production (WCGALP) Technical and species orientated innovations in animal breeding, and contribution of genetics to solving societal challenges, pp. 1530-1533. Wageningen Academic Publishers, 2022.
#'
#'   \item Lourenco, Daniela, Andres Legarra, Shogo Tsuruta, Yutaka Masuda, Ignacio Aguilar, and Ignacy Misztal. "Single-step genomic evaluations from theory to practice: using SNP chips and sequence data in BLUPF90." Genes 11, no. 7 (2020): 790.
#'
#' }
#'
#'
#'
#' @export
#'
blup <- function(datarenum,formula,fields_output=NULL,weights_object=NULL,residual_start=1,VCA_RRM=NULL,VCD_RRM=NULL,VCA=0.1,VCD=0.1,
                           ped_name,PED_DEPTH=3,genotype_file,missing_values=0,RRM_option,het_res_variance=NULL,
                           fit_option,
                           run_model=list(method="blupf90+",dense=F,
                                          Gibbs_option= NULL,
                                          slurm_option=NULL,
                                          n_threads=1),
                           cross_validation=list(group_name_cv,LOO_cv,prop,nrep,eff_name,RRM_option,
                                                 keep_folder=F,parameter_file="renf90.par",logname = "blupf90.log"),
                           keep_files=T,
                           extra.option.blup=NULL
){
  "%ni%" <- Negate("%in%")
  #CHECK and download blupf90 family software
  #download_BLUPF90(dest_folder = paste0(.libPaths()[1],"/blupf90"))
  blupf90_folder <- paste0(.libPaths()[1],"/blupf90")
  ###########
  files_init <- list.files()

  if(inherits(model, "formula")==TRUE){ #check if is single or multi trait
    pipeline_renum <- pipeline_renum_single(datarenum,formula,fields_output,weights_object,residual_start,VCA_RRM,VCD_RRM,VCA,VCD,
                                            ped_name,PED_DEPTH,genotype_file,missing_values,RRM_option=RRM_option,
                                            het_res_variance,fit_option,extra.option.blup) #call function for single trait
    #####Call function to fit the model
    if(is.null(run_model)==F){

      fit_default <- list(method="blupf90+",dense=F,
                          Gibbs_option= NULL,
                          slurm_option=NULL,n_threads=2)

      run_model <- c(fit_default[names(fit_default)[!(names(fit_default) %in% names(run_model))]],
                     run_model)


      pipeline_fit_model(method=run_model$method,Gibbs_option=run_model$Gibbs_option,dense=run_model$dense,
                         slurm_option=run_model$slurm_option,
                         parameter_file="renf90.par",n_threads=run_model$n_threads)

    ############################
    #Call function to extract VC
    ############################
    out_variance_c <- pipeline_extract_vc(datarenum=datarenum,formula,logname="blupf90.log",parms_card="renf90.par",
                                          method=run_model$method,RRM_option=RRM_option)


    #    Creating renf90.par with using output
    fields_output <- c(fields_output, unique(c(cross_validation$group_name_cv,cross_validation$LOO_cv)))

    pipeline_renum_single(datarenum,formula,fields_output,
                          weights_object,VCA=out_variance_c$GENETIC,VCD=out_variance_c$RANDOM,
                          residual_start=out_variance_c$RESIDUAL,
                          VCA_RRM=out_variance_c$RRM_GENETIC,
                          VCD_RRM=out_variance_c$RRM_RANDOM,
                          ped_name,
                          PED_DEPTH,genotype_file,missing_values,
                          RRM_option,het_res_variance,
                          fit_option=list(VCE=F,maxrounds=1),extra.option.blup)
    }
    ###################################
    #####CROSS VALIDATION
    ##################################
    if(is.null(cross_validation)==F){

      CROSS_default <- list(LOO_cv=NULL,
                            prop=0.2,nrep=5,
                            keep_folder=F,parameter_file="renf90.par",
                            logname="blupf90.log",
                            method="blupf90+")

      cross_validation <- c(CROSS_default[names(CROSS_default)[!(names(CROSS_default) %in% names(cross_validation))]],
                      cross_validation)


      pipeline_cross_validation (datarenum=datarenum,formula=formula,
                                 group_name_cv=cross_validation$group_name_cv,
                                 LOO_cv=cross_validation$LOO_cv,
                                 prop=cross_validation$prop,nrep=cross_validation$nrep,
                                 eff_name=cross_validation$eff_name,
                                 RRM_option=RRM_option,keep_folder=cross_validation$keep_folder,
                                 parameter_file=cross_validation$parameter_file,
                                 logname = cross_validation$logname,genotype_file=genotype_file,
                                 n_threads=run_model$n_threads,
                                 method=run_model$method)
    }

###################################
#####SUMMARY
##################################
if(is.null(RRM_option)==F){
  h2_vc_time_var <- read.table("VC_Time_var.txt",h=T)
  solutions <- read.table("solutions.orig",h=T)
  sol_cross_val <- read.table("cross_validation_result.txt",h=T)
  total_EBV <- read.table("EBV_t.txt",h=T)
  continuous_EBV <- read.table("EBV_Time_var.txt",h=T) %>%
    mutate(Time_var=rep(h2_vc_time_var$Time_var,length(total_EBV$ID))) %>%
    select(ID,Time_var,EBV)

##ADD th acc
  write.table(continuous_EBV,"EBV_Time_var.txt",row.names = F,quote = F)

  summary_blup <- list(formula,method=run_model$method,h2_vc_time_var,
                       out_variance_c,solutions,sol_cross_val,total_EBV,continuous_EBV)

  names(summary_blup) <- c("formula","Method","Heritability","Variance components",
                           "Solutions","Validation","EBV_Total","EBV_time")

  summary_blup$`Variance components`$RRM_RANDOM <-
  matrix(summary_blup$`Variance components`$RRM_RANDOM,ncol=RRM_option$poly+1)

summary_blup$`Variance components`$RRM_GENETIC <-
  matrix(summary_blup$`Variance components`$RRM_GENETIC,ncol=RRM_option$poly+1)

summary_blup$polynomial <- data.table::fread("fi.txt")
}else{
  #ADD heritability single trait nonRRM
  solutions <- read.table("solutions.orig",h=T)

  if(is.null(cross_validation)==T){
    sol_cross_val <- NULL}else{
    sol_cross_val <-read.table("cross_validation_result.txt",h=T)
    }
  summary_blup <- list(formula,method=run_model$method,
                       out_variance_c,solutions,sol_cross_val)

  names(summary_blup) <- c("formula","Method","Variance components",
                           "Solutions","Validation")

}
###################################################################
#Check output files and delete what was not before the analysis
files_end <- list.files()
nonINIT_file <- files_end[files_end%ni%files_init]

if(keep_files==F){
  unlink(nonINIT_file, recursive = TRUE)
}

return(summary_blup)
  }else{
    ###################################
    #####MULTI-TRAIT MODEL
    ##################################
    pipeline_renum <- pipeline_renum_multi(datarenum,formula,fields_output=fields_output,weights_object,residual_start,VCA_RRM,VCD_RRM,VCA,VCD,
                                           ped_name,PED_DEPTH,genotype_file,missing_values,RRM_option,het_res_variance,
                                           fit_option,extra.option.blup) #call function for multi trait

    #####Call function to fit the model
    if(is.null(run_model)==F){

      fit_default <- list(method="blupf90+",dense=F,
                          Gibbs_option= NULL,
                          slurm_option=NULL,n_threads=2)

      run_model <- c(fit_default[names(fit_default)[!(names(fit_default) %in% names(run_model))]],
                     run_model)

      pipeline_fit_model(method=run_model$method,Gibbs_option=run_model$Gibbs_option,dense=run_model$dense,
                         slurm_option=run_model$slurm_option,
                         parameter_file="renf90.par",n_threads=run_model$n_threads)


    ############################
    #Call function to extract VC
    ############################
    out_variance_c <- pipeline_extract_vc_multi(datarenum=datarenum,formula,logname="blupf90.log",parms_card="renf90.par",
                                                method=run_model$method,RRM_option=RRM_option)


    #    Creating renf90.par with using output
    fields_output <- c(fields_output, unique(c(cross_validation$group_name_cv,cross_validation$LOO_cv)))

    pipeline_renum_multi(datarenum,formula,fields_output=fields_output,
                         weights_object,VCA=out_variance_c$GENETIC,VCD=out_variance_c$RANDOM,
                         residual_start=out_variance_c$RESIDUAL,
                         VCA_RRM=out_variance_c$RRM_GENETIC,
                         VCD_RRM=out_variance_c$RRM_RANDOM,
                         ped_name,
                         PED_DEPTH,genotype_file,missing_values,
                         RRM_option,het_res_variance,
                         fit_option=list(VCE=F,maxrounds=1),extra.option.blup)

    output_multi<- list(formula,method=run_model$method)
    output_multi$Var_Cor <- out_variance_c
    ###################################
    #####Genetic Correlation
    ##################################
    if(is.null(RRM_option)==F){

    out_variance_c$RESIDUAL=NULL

    gpl <- list()
    i=1
    ngpl <- NULL
    for (gpi in 1:length(out_variance_c)) {

      if(length(out_variance_c[[gpi]]) == 0) {
        gpl[[gpi]]=NULL
        i=i
      }else{
        gpl[[i]]=out_variance_c[[gpi]]
        ngpl[i] <- names(out_variance_c)[gpi]
        i=i+1
      }
    }

    names(gpl) <- ngpl

    gen_cor <- RRM_cor(log_name="blupf90.log",vc_names=ngpl,parms_card="renf90.par",DAPmin=RRM_option$Pmin,RRM_option$Pmax)
    output_multi$gen_cor <- gen_cor
    output_multi$polynomial <- data.table::fread("fi.txt")

    }

    ###################################
    #####Heritability
    ##################################
    if(is.null(RRM_option)==F){
      h2multi <- h2_multi(datarenum=datarenum,formula=model,logname="blupf90.log",parms_card="renf90.par",
             method="blupf90+",RRM_option=RRM_option)
      output_multi$h2_multi <- h2multi

    }
    ###################################
    #####CROSS VALIDATION
    ##################################
    if(is.null(cross_validation)==F){

      CROSS_default <- list(LOO_cv=NULL,
                            prop=0.2,nrep=5,
                            keep_folder=F,parameter_file="renf90.par",
                            logname="blupf90.log")

      cross_validation <- c(CROSS_default[names(CROSS_default)[!(names(CROSS_default) %in% names(cross_validation))]],
                            cross_validation)
    pipeline_cross_validation_multi(datarenum=datarenum,formula=formula,group_name_cv=cross_validation$group_name_cv,LOO_cv=cross_validation$LOO_cv,
                                    prop=cross_validation$prop,nrep=cross_validation$nrep,eff_name=cross_validation$eff_name,keep_folder=cross_validation$keep_folder,
                                    genotype_file=genotype_file,
                                    parameter_file="renf90.par",logname = "blupf90.log",n_threads=run_model$n_threads,RRM_option=RRM_option)

      }
     }
    return(output_multi)
  }

#Check output files and delete what was not before the analysis
files_end <- list.files()
nonINIT_file <- files_end[files_end%ni%files_init]


if(keep_files==F){
  unlink(nonINIT_file, recursive = TRUE)
    }
}
