# Central path setup for console, RStudio Source, and Rscript.
if(.Platform$OS.type=='windows')Sys.setlocale('LC_CTYPE','English_United States.utf8')
.alps_boot_files<-unlist(lapply(sys.frames(),function(z){x<-get0('ofile',envir=z,inherits=FALSE);if(is.character(x))x else character()}),use.names=FALSE)
.alps_boot_files<-as.character(.alps_boot_files)
.alps_boot_files<-.alps_boot_files[basename(.alps_boot_files)=='bootstrap_local.R']
.alps_boot_root<-if(length(.alps_boot_files))dirname(dirname(normalizePath(tail(.alps_boot_files,1),winslash='/',mustWork=TRUE)))else Sys.getenv('ALPS_LOCAL_PROJECT')
if(!nzchar(.alps_boot_root)).alps_boot_root<-getwd()
ALPS_SETUP<-local({default_root<-.alps_boot_root;function(project=default_root,quiet=FALSE){
 project<-normalizePath(project,winslash='/',mustWork=TRUE)
 cfg_file<-file.path(project,'project_paths.json')
 if(!file.exists(cfg_file))stop('Missing configuration: ',cfg_file,call.=FALSE)
 if(!requireNamespace('jsonlite',quietly=TRUE))stop('R package jsonlite is missing. Install/restore the documented R environment first.',call.=FALSE)
 cfg<-jsonlite::fromJSON(cfg_file)
 path<-function(x){if(is.null(x)||!nzchar(x))return('');if(!grepl('^[A-Za-z]:[/\\\\]|^/|^\\\\\\\\',x))x<-file.path(project,x);normalizePath(x,winslash='/',mustWork=FALSE)}
 program<-function(x,kind){if(!is.null(x)&&nzchar(x)&&!grepl('[/\\]',x)){v<-Sys.which(x);if(nzchar(v))return(normalizePath(v,winslash='/',mustWork=TRUE));if(kind=='Rscript')return(file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe'else'Rscript'))};path(x)}
 paths<-list(project=project,data=file.path(project,'data'),results=file.path(project,'results'),figures=file.path(project,'figures'),logs=file.path(project,'logs'),raw=path(cfg$numeric_source),python=program(cfg$python_executable,'python'),python_library=path(cfg$python_library),rscript=program(cfg$rscript_executable,'Rscript'))
 for(n in c('results','figures','logs'))dir.create(paths[[n]],recursive=TRUE,showWarnings=FALSE)
 Sys.setenv(ALPS_LOCAL_PROJECT=project,ALPS_REPRO_ROOT=project,ALPS_RAW_DATA=paths$raw,PYTHONIOENCODING='utf-8')
 if(nzchar(paths$python))Sys.setenv(ALPS_PYTHON=paths$python)
 if(nzchar(paths$rscript))Sys.setenv(ALPS_RSCRIPT=paths$rscript)
 if(dir.exists(paths$python_library)){
  old<-Sys.getenv('PYTHONPATH');entries<-unique(c(paths$python_library,if(nzchar(old))strsplit(old,.Platform$path.sep,fixed=TRUE)[[1]]));Sys.setenv(PYTHONPATH=paste(entries,collapse=.Platform$path.sep))
 }
 if(!quiet){cat('PROJECT :',paths$project,'\nDATA    :',paths$data,'\nRESULTS :',paths$results,'\nFIGURES :',paths$figures,'\nNUMERIC_SOURCE:',paths$raw,'\n')}
 invisible(paths)
}})
ALPS_REQUIRE_FILES<-function(files,folder='data'){
 paths<-ALPS_SETUP(quiet=TRUE);base<-paths[[folder]];if(is.null(base))stop('Unknown path group: ',folder)
 expected<-file.path(base,files);bad<-expected[!file.exists(expected)]
 if(length(bad))stop('Required files are missing:\n',paste(bad,collapse='\n'),'\nRun prepare mode with the configured corrected numerical source, or restore the complete local folder.',call.=FALSE)
 invisible(expected)
}
ALPS_READ<-function(name,...){p<-ALPS_REQUIRE_FILES(name);utils::read.csv(p,check.names=FALSE,fileEncoding='UTF-8-BOM',...)}
ALPS_CHECK<-function(check_packages=TRUE,stop_on_error=TRUE){
 paths<-ALPS_SETUP();files<-c('baseline_master.csv','longitudinal_with_site_protocol_proxy.csv','structural_PET216.csv','prognostic_covariates_preindex365.csv')
 tab<-data.frame(File=file.path(paths$data,files),Exists=file.exists(file.path(paths$data,files)));print(tab,row.names=FALSE)
 pkg<-NULL
 if(check_packages){expected<-utils::read.csv(file.path(paths$project,'R_PACKAGE_VERSIONS.csv'),stringsAsFactors=FALSE,fileEncoding='UTF-8-BOM');pkg<-data.frame(Package=expected$Package,Expected=expected$Version,Installed=vapply(expected$Package,function(p)if(requireNamespace(p,quietly=TRUE))as.character(packageVersion(p))else'MISSING',character(1)));print(pkg,row.names=FALSE);utils::write.csv(pkg,file.path(paths$logs,'R_package_check.csv'),row.names=FALSE)}
 ok<-all(tab$Exists)&&(is.null(pkg)||!any(pkg$Installed=='MISSING'))
 utils::write.csv(tab,file.path(paths$logs,'Input_file_check.csv'),row.names=FALSE)
 if(!ok&&stop_on_error)stop('Preflight failed. The missing paths/packages are listed above.',call.=FALSE)
 if(ok)cat('PREFLIGHT_OK. R version:',as.character(getRversion()),'\n')
 invisible(list(ok=ok,paths=paths,files=tab,packages=pkg))
}
ALPS_MODEL<-function(module){
 paths<-ALPS_SETUP();file<-file.path(paths$project,'R/models',module)
 if(!file.exists(file))stop('Unknown model module: ',file,call.=FALSE)
 ALPS_REQUIRE_FILES(c('baseline_master.csv','longitudinal_with_site_protocol_proxy.csv'))
 cat('Running:',file,'\n');env<-new.env(parent=globalenv());env$read.csv<-utils::read.csv;source(file,local=env,encoding='UTF-8');cat('MODEL_COMPLETED:',module,'\n');invisible(env)
}
ALPS_REDRAW<-function(){
 paths<-ALPS_SETUP();ALPS_REQUIRE_FILES(c('baseline_master.csv','longitudinal_with_site_protocol_proxy.csv'))
 source(file.path(paths$project,'R/plots/draw_figures.R'),local=new.env(parent=globalenv()),encoding='UTF-8');cat('PDF_FIGURES_COMPLETED:',paths$figures,'\n');invisible(paths$figures)
}
ALPS_RUN<-function(mode=c('check','all','prepare','models','figures','validate','paired')){
 mode<-match.arg(mode);paths<-ALPS_SETUP()
 if(mode=='check')return(invisible(ALPS_CHECK()))
 if(mode=='figures')return(invisible(ALPS_REDRAW()))
 if(!nzchar(paths$python))stop('Configure python_executable in project_paths.json.',call.=FALSE)
 status<-system2(shQuote(paths$python),c(shQuote(file.path(paths$project,'run_pipeline.py')),mode))
 if(status!=0)stop('Pipeline failed. Inspect logs/pipeline_*.log.',call.=FALSE)
 invisible(status)
}
ALPS_PATHS<-ALPS_SETUP(quiet=TRUE)
