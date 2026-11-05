# Internal deterministic module blocks. Run the project entry point for a full rebuild.
if(.Platform$OS.type=='windows')Sys.setlocale('LC_CTYPE','English_United States.utf8')
a<-commandArgs(trailingOnly=FALSE);this<-sub('^--file=','',a[grep('^--file=',a)][1]);root<-normalizePath(file.path(dirname(this),'..'),winslash='/');Sys.setenv(ALPS_REPRO_ROOT=root)
block<-commandArgs(trailingOnly=TRUE)[1]
if(is.na(block)||!block%in%c('clinical','prediction','expression'))stop('Use run_all.R or ALPS_RUN("models"); run_block.R is an internal CLI entry.')
blocks<-list(
 clinical=c('01_longitudinal_cognition_prognosis.R','02_mediation_genetics_BBB.R','03_network_regularisation_and_amyloid.R','04_additional_reporting_models.R','08_additional_reporting.R','12_matched_matrices.R','13_additional_models.R','14_target_longitudinal.R','17_anatomical_terms.R','18_regional_models.R','19_structural_network.R','20_review_sensitivities.R'),
 prediction=c('09_linear_prediction.R','10_nested_prediction.R'),
 expression=c('05_transcriptome_reproduction_and_diagnosis.R','06_transcriptome_GSEA_from_refitted_model.R','11_GSEA_bootstrap.R','16_CAMERA.R','15_fixed_mapping_GSEA.R'))
stopifnot(block%in%names(blocks))
for(n in blocks[[block]]){cat('START',n,format(Sys.time()),'\n');flush.console();source(file.path(root,'R/models',n),local=new.env(parent=globalenv()),encoding='UTF-8');cat('FINISHED',n,format(Sys.time()),'\n');flush.console()}
cat('BLOCK_COMPLETED',block,format(Sys.time()),'\n')
