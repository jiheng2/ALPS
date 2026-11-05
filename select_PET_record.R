# Prospective merger. Exact existing PET date is required unless nearest selection is explicitly requested.
# Returns synchronized ROI values, source record ID, date and both interval definitions.
select_PET_record<-function(master,source,allow_nearest=FALSE,window_days=183L){
 stopifnot(all(c('RID','DTI_DATE','FDG_Date')%in%names(master)),'Subject_ID'%in%names(source))
 source$RID<-as.integer(sub('.*_S_([0-9]+).*','\\1',source$Subject_ID))
 source$PET_SOURCE_DATE<-as.Date(sub('.*_([0-9]{4}-[0-9]{2}-[0-9]{2})$','\\1',source$Subject_ID))
 if(anyNA(source$RID))stop('Source ID lacks an unambiguous participant RID')
 if(anyNA(source$PET_SOURCE_DATE))stop('Source ID lacks an unambiguous acquisition date')
 if(anyDuplicated(paste(source$RID,source$PET_SOURCE_DATE)))stop('Resolve multiple source records on the same participant/date using documented image identifiers')
 roi<-grep('^AAL3_',names(source),value=TRUE);out<-master
 out$PET_SOURCE_RECORD_ID<-NA_character_;out$PET_SOURCE_DATE<-as.Date(NA);out$PET_minus_DTI_days<-NA_real_;out$PET_absolute_gap_days<-NA_real_
 for(i in seq_len(nrow(master))){
  if(is.na(master$FDG_Date[i])&&!allow_nearest)next
  z<-source[source$RID==master$RID[i],];gap<-as.numeric(z$PET_SOURCE_DATE-as.Date(master$DTI_DATE[i]));z<-z[is.finite(gap)&abs(gap)<=window_days,,drop=FALSE]
  if(!nrow(z))next
  exact<-which(z$PET_SOURCE_DATE==as.Date(master$FDG_Date[i]));j<-exact
  if(!length(j)&&allow_nearest){g<-abs(as.numeric(z$PET_SOURCE_DATE-as.Date(master$DTI_DATE[i])));j<-which(g==min(g));if(length(j)!=1)stop('Equidistant PET records: resolve explicitly')}
  if(length(j)!=1)stop('No unique exact PET record; refusing to replace values while keeping an old date')
  out[i,roi]<-z[j,roi];if('WholeBrain_SUVR_AAL3'%in%names(source))out$WholeBrain_SUVR_AAL3[i]<-z$WholeBrain_SUVR_AAL3[j]
  out$PET_SOURCE_RECORD_ID[i]<-z$Subject_ID[j];out$PET_SOURCE_DATE[i]<-z$PET_SOURCE_DATE[j];out$FDG_Date[i]<-as.character(z$PET_SOURCE_DATE[j]);out$PET_minus_DTI_days[i]<-as.numeric(z$PET_SOURCE_DATE[j]-as.Date(master$DTI_DATE[i]));out$PET_absolute_gap_days[i]<-abs(out$PET_minus_DTI_days[i]);out$DTI_FDG_gap_days[i]<-out$PET_absolute_gap_days[i]
  if('FDG_DATE'%in%names(out))out$FDG_DATE[i]<-as.character(z$PET_SOURCE_DATE[j])
  if('PET_visit_id'%in%names(out))out$PET_visit_id[i]<-z$Subject_ID[j]
  if('DTI_PET_gap_days'%in%names(out))out$DTI_PET_gap_days[i]<-out$PET_absolute_gap_days[i]
 }
 out
}
