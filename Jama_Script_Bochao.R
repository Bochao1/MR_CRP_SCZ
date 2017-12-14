##-----------------------------------
## Script name: replicate table 1 from Hartwig et al with 18 liberal threshold CRP SNPs. 
##      author: Bochao Lin
##        date: 7th Dec 2017
##-------------------------------------
rm(list = ls(all.names = TRUE))
library(MendelianRandomization)
library(MRInstruments)
library(TwoSampleMR)
#-----------------
#-----------------Part 1 extract intrument variables from GWAScatalog for both CRP and SCZ GWAS by {MRInstruments} packages and do MR analyses.
#-----------------
##--get CRP summary statistics from Deghan et al.
snp18=c('rs2794520','rs4420638','rs1183910','rs4420065','rs4129267','rs1260326',
        'rs12239046','rs6734238','rs9987289','rs10745954','rs1800961','rs340029',
        'rs10521222','rs12037222','rs13233571','rs2847281','rs6901250','rs4705952')
gwas=data(gwas_catalog)
crpgwas=subset(gwas_catalog,gwas_catalog$PubmedID==21300955)
crp_18 <- TwoSampleMR::format_data(crpgwas, type="exposure")
crp_18=subset(crp_18,crp_18$SNP %in% snp18)
##--get SCZ summary statistics from PGC
ao <- available_outcomes()
psy=subset(ao,ao$subcategory=='Psychiatric / neurological')
outcome_dat_scz <- extract_outcome_data(snps=snp18, outcomes=22, proxies = TRUE, rsq = 0.8,
                                        align_alleles = 1, palindromes = 1, maf_threshold = 0.3,
                                        access_token = get_mrbase_access_token())
pubdata <- harmonise_data(crp_18, outcome_dat_scz,action=1)
pubdata$beta.outcome=round(pubdata$beta.outcome,4)
subset(pubdata,select=c(SNP,effect_allele.exposure,other_allele.exposure,effect_allele.outcome,other_allele.outcome,
                        beta.exposure,se.exposure,beta.outcome,se.outcome))
## after data harmonization, the effect allele of rs9987289 is different from Hartwig.
### Issues
crp_18[crp_18$SNP=="rs9987289",]
outcome_dat_scz[outcome_dat_scz$SNP=="rs9987289",]
##--do MR-base analysis in {TwosampleMR}.
result1=TwoSampleMR::mr(pubdata,method = c('mr_ivw','mr_egger_regression','mr_weighted_median'))
res_TwosampleMR=data.frame(package=rep('TwosampleMR',3),
  model=c("ivw_random","egger_random",'weighted_mdeilan'),
  OR=round(exp(result1$b),2),
  SE.low=round(exp(result1$b-1.96*result1$se),2),
  SE.up=round(exp(result1$b+1.96*result1$se),2))
##--do MR-base analysis in {MendelianRandomization}
MRInputObject <- mr_input(bx = pubdata$beta.exposure,
                          bxse = pubdata$se.exposure,
                          by = pubdata$beta.outcome,
                          byse = pubdata$se.outcome)
x=MendelianRandomization::mr_ivw(MRInputObject,model = 'fixed')
mr_ivm_fixed_MR2=round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
x=MendelianRandomization::mr_ivw(MRInputObject,model = 'random')
mr_ivm_random_MR2=round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
x=MendelianRandomization::mr_egger(MRInputObject)
mr_egger_random_MR2=round(exp(c(x@Estimate,x@CILower.Est,x@CIUpper.Est)),2)
x=MendelianRandomization::mr_median(MRInputObject,weighting = "weighted", distribution = "normal")
mr_median_MR2=round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
resMR2=data.frame(package=rep('MendelianRandomization',4),
                  model=c("ivw_fixed","ivw_random","egger_random",'weighted_mdeilan'),
                  rbind(mr_ivm_fixed_MR2,mr_ivm_random_MR2,mr_egger_random_MR2,mr_median_MR2))
colnames(resMR2)=c('package','model','OR','SE.low','SE.up')
##--do MR-base analysis using our own R scripts
#data entry
betaXG=pubdata$beta.exposure
betaYG=pubdata$beta.outcome
sebetaXG=pubdata$se.exposure
sebetaYG=pubdata$se.outcome
# Inverse-variance weighted method (FE and then RE)
betaIVW = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1]
sebetaIVW.fixed = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,2]/
  summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$sigma
sebetaIVW.random = summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$coef[1,2]/
  min(summary(lm(betaYG~betaXG-1, weights=sebetaYG^-2))$sigma,1)
mr_ivw_random_manual=round(exp(c(betaIVW,betaIVW-1.96*sebetaIVW.random,betaIVW+1.96*sebetaIVW.random)),2)
mr_ivw_fixed_manual=round(exp(c(betaIVW,betaIVW-1.96*sebetaIVW.fixed,betaIVW+1.96*sebetaIVW.fixed)),2)
# Egger regression method (RE)
betaYG = betaYG*sign(betaXG); betaXG = abs(betaXG)
betaEGGER = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,1]
sebetaEGGER.random = summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$coef[2,2]/
min(summary(lm(betaYG~betaXG, weights=sebetaYG^-2))$sigma, 1)
mr_egger_random_manual=round(exp(c(betaEGGER,betaEGGER-1.96*sebetaEGGER.random,betaEGGER+1.96*sebetaEGGER.random)),2)
resManual=data.frame(package=rep('none',3),model=c('ivw_random','ivw_fixed','egger_random'),
  rbind(mr_ivw_random_manual,mr_ivw_fixed_manual,mr_egger_random_manual))
colnames(resManual)=c('package','model','OR','SE.low','SE.up')
##Compare estimates from three different methods:
com.res=rbind(res_TwosampleMR,resMR2,resManual)
View(com.res)
# Conclusions: 
# 1. We used a different  effect allele for rs9987289 on CRP than Hartwig et al.
# 2. The estimations from two packages and our own script are the same.

#-----------------
#-----------------Part 2:  Use Hartwig et al supplement etable 2 as input data (no rounding), and do MR analyses.
#-----------------
bx=c(0.1600,0.2360,0.1490,0.0900,0.0790,0.0720,
     0.0470,0.0500,0.0690,0.0390,0.0880,0.0320,
     0.1040,0.0450,0.0540,0.0310,0.0350,0.0420)
bxse=c(0.0060,0.0090,0.0060,0.0050,0.0050,0.0050,
       0.0060,0.0060,0.0110,0.0060,0.0150,0.0060,
       0.0150,0.0070,0.0090,0.0060,0.0060,0.0070)
by=c(-0.0220,-0.0086,-0.0278,-0.0224,-0.0257,-0.0061,
     -0.0034,-0.0027,-0.0832,-0.0251,0.0307,0.0043,
     -0.0091,-0.0095,0.0131,0.0113,-0.0133,-0.0033)
byse=c(0.0110,0.0146,0.0112,0.0110,0.0108,0.0108,
       0.0110,0.0110,0.0188,0.0106,0.0300,0.0112,
       0.0325,0.0128,0.0167,0.0108,0.0113,0.0121)
hartiwigetable2=cbind(snp18,bx,bxse,by,byse)
hartiwigetable2

MRInputObject2 <- mr_input(bx = bx,
                          bxse = bxse,
                          by = by,
                          byse = byse)
x=MendelianRandomization::mr_ivw(MRInputObject2,model = 'fixed')
round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
x=MendelianRandomization::mr_ivw(MRInputObject2,model = 'random')
round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
x=MendelianRandomization::mr_egger(MRInputObject2)
round(exp(c(x@Estimate,x@CILower.Est,x@CIUpper.Est)),2)
round(exp(c(x@Intercept,x@CILower.Int,x@CIUpper.Int)),2)
x=MendelianRandomization::mr_median(MRInputObject2)
round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)

# Even when using the same input data as Hartwig et al, the estimates we get are different from Hartwig table 1.

#-----------------
#-----------------Part 3 replicate table 3 in Prins et al.
#-----------------
snp15=c('rs2794520','rs4420065','rs4129267','rs1260326',
        'rs12239046','rs6734238','rs10745954','rs1800961','rs340029',
        'rs10521222','rs12037222','rs13233571','rs2847281','rs6901250','rs4705952')
gwas=data(gwas_catalog)
crpgwas=subset(gwas_catalog,gwas_catalog$PubmedID==21300955)
crp_15 <- TwoSampleMR::format_data(crpgwas, type="exposure")
crp_15=subset(crp_15,crp_15$SNP %in% snp15)

##--get SCZ summary statistics from PGC
ao <- available_outcomes()
psy=subset(ao,ao$subcategory=='Psychiatric / neurological')
outcome_dat_scz <- extract_outcome_data(snps=snp15, outcomes=22, proxies = TRUE, rsq = 0.8,
                                        align_alleles = 1, palindromes = 1, maf_threshold = 0.3,
                                        access_token = get_mrbase_access_token())
pubdata <- harmonise_data(crp_15, outcome_dat_scz,action=1)
pubdata$beta.outcome=round(pubdata$beta.outcome,4)
subset(pubdata,select=c(SNP,effect_allele.exposure,other_allele.exposure,effect_allele.outcome,other_allele.outcome,
                        beta.exposure,target_a1.outcome,beta.outcome,se.outcome))
##--do MR-base analysis in {TwosampleMR}.
head(pubdata)
result1=TwoSampleMR::mr(pubdata,method = c('mr_ivw','mr_egger_regression','mr_weighted_median'))
exp(result1$b)
exp(result1$b-1.96*result1$se)
exp(result1$b+1.96*result1$se)
##--do MR-base analysis in {MendelianRandomization}
MRInputObject <- mr_input(bx = pubdata$beta.exposure,
                          bxse = pubdata$se.exposure,
                          by = pubdata$beta.outcome,
                          byse = pubdata$se.outcome)
x=MendelianRandomization::mr_ivw(MRInputObject,model = 'fixed')
round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)
x=MendelianRandomization::mr_ivw(MRInputObject,model = 'random')
round(exp(c(x@Estimate,x@CILower,x@CIUpper)),2)

#IVW model results are consistent with Prins et al results.

