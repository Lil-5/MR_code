mr_scatter_plot_modified <- function(mr_results, dat)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  requireNamespace("ggplot2", quietly=TRUE)
  requireNamespace("plyr", quietly=TRUE)
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point(size = 2 ) +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE,linewidth=0.8) +
      ggplot2::scale_colour_manual(values=c("#ff0000", "#16982b", "#4798b3", "#ffd700",  "#e11ae3", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      
      ggplot2::labs(colour="MR Test", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(legend.position="top", legend.direction="vertical") +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2))
  })
  mrres
}
library(data.table)
library(R.utils)
setwd("E:\\micro")
library(TwoSampleMR)
index <- read.table("1.txt",as.is = T)
index
index <- as.vector(t(index))
sensitive=data.frame()

for (ebi in index) {
  #exposure_dat <-extract_instruments(ebi,
  #                                    p1=1e-05,clump=TRUE,
  #                                    access_token = NULL)  
  exposure_dat <-fread(file = paste0(ebi,".csv"),header = T)
  if(length(exposure_dat[,1])>0){
    #exposure_dat<-snp_add_exposure_eaf(exposure_dat)
    outcome_dat <-extract_outcome_data(snps=exposure_dat$SNP, 
                                       outcomes="NAFLD.csv",access_token = NULL)
    dat <- harmonise_data(exposure_dat,outcome_dat)
    res <- generate_odds_ratios(mr_res = mr(dat,method_list = c("mr_ivw",
                                                                "mr_egger_regression",
                                                                "mr_weighted_median",
                                                                "mr_weighted_mode",
                                                                "mr_simple_mode",
                                                                "mr_wald_ratio")))
    write.csv(res,file = paste0("dres",ebi,".csv"), row.names = FALSE)
    write.csv(dat,file = paste0("ddat",ebi,".csv"), row.names = FALSE)
    if(res[1,9]<0.05){q=mr_heterogeneity(dat)
    egger=mr_pleiotropy_test(dat)
    dat$samplesize.outcome=463010
    direct=directionality_test(dat)
    #library(MRPRESSO)
    presso=run_mr_presso(dat = dat,NbDistribution = 1000,SignifThreshold = 0.05)
    sensitive <- rbind(sensitive,data.frame(row.names = ebi,
                                            Q=q[2,8],
                                            MRegger_interpreter=egger$egger_intercept,
                                            MRegger_interpreter_pval=egger$pval,
                                            steiger_pval=direct$steiger_pval,
                                            MRPRESSO_GLOBAL=presso[[1]][["MR-PRESSO results"]]$`Global Test`$Pvalue))
    
    
    ##
    pdf(file = paste0("scatter",ebi,".pdf"),width = 6,height = 6)
    p1=mr_scatter_plot_modified(res,dat)
    print(p1)
    dev.off()
    pdf(file = paste0("plot",ebi,".pdf"),width = 6,height = 6)
    p2=mr_funnel_plot(singlesnp_results = mr_singlesnp(dat))
    print(p2)
    dev.off()
    pdf(file = paste0("leave one out",ebi,".pdf"),width = 6,height = 6)
    p3=mr_leaveoneout_plot(leaveoneout_results = mr_leaveoneout(dat))
    print(p3)
    dev.off()}}}
write.csv(sensitive,"sensitive.csv")
a = list.files(pattern = "dres")
dir = paste("./",a,sep="")
n = length(dir) 
merge_data = read.csv(file = dir[1],sep=",",header=T) 
for (i in 2:n){
  new_data = read.csv(file = dir[i],sep=",",header=T)
  merge_data = merge(merge_data,new_data,all=TRUE)
}
write.csv(merge_data,file="1.csv")
a = list.files(pattern = "ddat")
dir = paste("./",a,sep="")
n = length(dir) 
merge_data = read.csv(file = dir[1],sep=",",header=T) 
for (i in 2:n){
  new_data = read.csv(file = dir[i],sep=",",header=T)
  merge_data = merge(merge_data,new_data,all=TRUE)
}
write.csv(merge_data,file="2.csv")
#totalsnp<-merge_data[,c(26,1,2,3,8,6,26,27,7,16,18,29)]
totalsnp<-merge_data[,c("id.exposure","SNP","effect_allele.exposure",
                        "other_allele.exposure","eaf.exposure","beta.exposure",
                        "se.exposure","pval.exposure","beta.outcome",
                        "se.outcome","pval.outcome","mr_keep")]
colnames(totalsnp)<-c("id.exposure","SNP","effect_allele",
                      "other_allele",
                      "eaf",
                      "beta.exposure",
                      "se.exposure",
                      "pval.exposure",
                      "beta.outcome",	"se.outcome",	"pval.outcome",
                      "mr_keep"
)
totalsnp<-totalsnp[totalsnp$mr_keep=="TRUE",]
totalsnp<-totalsnp[,-12]
totalsnp=dplyr::arrange(totalsnp,id.exposure)
totalsnp<-transform(totalsnp,R2=2*eaf*(1-eaf)*(beta.exposure)^2)
totalsnp<-transform(totalsnp,F=18338*R2/(1-R2))
which(totalsnp$F<10)
write.csv(totalsnp,"ivs.csv",row.names = FALSE)

#library(openxlsx)
s=read.csv("ivs.csv")
q=read.csv("1.csv")
#index=as.vector(t(read.table("11.txt",as.is = T)))
index=rownames(sensitive)
df=s[s$id%in%index,]
dq=q[q$id.exposure%in%index,]
write.csv(df,"positive IVs.csv",row.names = FALSE)
write.csv(dq,"positive MR.csv",row.names = FALSE)

library(phenoscanner)
n<-length(df[,"SNP"])
if(n>100){
  mulsnp1<-df[1:50,"SNP"]
  mulsnp2<-df[51:n,"SNP"]
  pheno1<-phenoscanner(snpquery = mulsnp1)
  pheno2<-phenoscanner(snpquery = mulsnp2)
  p<-rbind(pheno1[["results"]],pheno2[["results"]])
  p1<-p[,c("snp","trait")]
  write.csv(p1,"phenoscanner.csv",row.names = FALSE)
}else{mulsnp<-df[,"SNP"]
pheno<-phenoscanner(snpquery = mulsnp)
p<-pheno[["results"]]
p1<-p[,c("snp","trait")]
write.csv(p1,"phenoscanner.csv",row.names = FALSE)}
df<-read.csv("API.csv")
df<-df[,2:3]
d<-read.csv("positive MR.csv")
d<-d[,-1]
result <- merge(d, df, by.x="id.exposure",by.y = "id",
                all.x=TRUE, all.y=FALSE)
result$id.outcome<-"NAFLD"
write.csv(result,"positive MR.csv",row.names = FALSE)

d<-read.csv("sensitive.csv")
result <- merge(d, df, by.x="X",by.y = "id",
                all.x=TRUE, all.y=FALSE)
write.csv(result,"sensitive.csv",row.names = FALSE)
d<-read.csv("1.csv")
result <- merge(d, df, by.x="id.exposure",by.y = "id",
                all.x=TRUE, all.y=FALSE)
write.csv(result,"1.csv",row.names = FALSE)