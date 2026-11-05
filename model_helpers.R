options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(sandwich))
hc3<-function(d,formula,term,min_n=30){
 vars<-all.vars(formula);z<-droplevels(d[complete.cases(d[,vars,drop=FALSE]),,drop=FALSE]);n<-nrow(z)
 ans<-data.frame(N=n,Estimate=NA_real_,SE=NA_real_,CI_low=NA_real_,CI_high=NA_real_,P=NA_real_,Status='Nonestimable')
 if(n<min_n||!is.finite(sd(z[[vars[1]]]))||sd(z[[vars[1]]])<1e-10)return(ans)
 f<-lm(formula,z);X<-model.matrix(f)
 if(f$rank<ncol(X)||!term%in%names(coef(f)))return(ans)
 b<-coef(f)[term];se<-sqrt(diag(vcovHC(f,type='HC3')))[term]
 if(!is.finite(se)||se<=0)return(ans)
 ans$Estimate<-b;ans$SE<-se;ans$CI_low<-b-qnorm(.975)*se;ans$CI_high<-b+qnorm(.975)*se;ans$P<-2*pnorm(-abs(b/se));ans$Status<-'OK';ans
}
bh_rows<-function(d,by,p='P',q='q',family_n=NULL){
 d[[q]]<-NA_real_;keys<-interaction(d[,by,drop=FALSE],drop=TRUE,lex.order=TRUE)
 for(k in levels(keys)){i<-which(keys==k&is.finite(d[[p]]));if(length(i))d[[q]][i]<-p.adjust(d[[p]][i],'BH',n=if(is.null(family_n))length(i)else family_n)};d
}
partial_rank<-function(d,x,y,covars,min_n=15){
 z<-droplevels(d[complete.cases(d[,unique(c(x,y,covars)),drop=FALSE]),,drop=FALSE]);n<-nrow(z)
 constant<-covars[vapply(z[,covars,drop=FALSE],function(a)length(unique(a))<2,logical(1))]
 out<-data.frame(N=n,rho=NA_real_,P=NA_real_,df=NA_real_,CN=sum(z$Group=='CN'),MCI=sum(z$Group=='MCI'),AD=sum(z$Group=='AD'),Constant_covariates_omitted=paste(constant,collapse=';'))
 if(n<min_n||length(unique(z[[x]]))<3||length(unique(z[[y]]))<3)return(out)
 covars<-setdiff(covars,constant)
 zz<-z;for(k in covars)if(is.numeric(zz[[k]]))zz[[k]]<-rank(zz[[k]],ties.method='average')
 X<-model.matrix(reformulate(covars),zz)
 if(qr(X)$rank<ncol(X))return(out)
 rx<-qr.resid(qr(X),rank(z[[x]],ties.method='average'));ry<-qr.resid(qr(X),rank(z[[y]],ties.method='average'));r<-cor(rx,ry);df<-n-ncol(X)-1
 out$rho<-r;out$df<-df;out$P<-2*pt(-abs(r*sqrt(df/(1-r*r))),df);out
}
