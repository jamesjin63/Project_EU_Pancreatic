# R1.1: pancreatic cancer YLL/YLD decomposition and YLL-based VLW sensitivity.
# Run from EuropeanUnion project root. No original inputs/outputs overwritten.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
root <- getwd()
out <- file.path(root,'analysis','R1_1_YLL_YLD')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
newfile <- file.path(root,'data/GBD2023_Pancreatic_cancer_32locations_YLLs_YLDs_1990_2023.csv')
eu <- fread(file.path(root,'EU28_location_list.csv'))[group_type=='EU-28 member',.(location_id,location_name)]
stopifnot(nrow(eu)==28,uniqueN(eu$location_id)==28)
ids <- eu$location_id; ages <- c(1L,6:20,30:32,235L)
rawcols <- c('location_id','location_name','sex_id','sex_name','age_id','age_name','year','measure_id','metric_id','cause_id','val','lower','upper')
old <- fread(file.path(root,'data/merged.csv'),select=rawcols)
new <- fread(newfile,select=rawcols)
stopifnot(all(new$cause_id==456),setequal(new$year,1990:2023))
keys <- c('location_id','sex_id','age_id','year')
old <- old[cause_id==456 & measure_id==2 & metric_id==1]
# Keep official aggregate locations only for separate diagnostics, never sum with countries.
d <- old[age_id %in% ages & location_id %in% new$location_id]
setnames(d,c('val','lower','upper'),c('DALY','DALY_lo','DALY_hi'))
for(meas in c(3,4)) {
 label <- if(meas==3)'YLD' else 'YLL'
 x <- new[measure_id==meas & metric_id==1,c(keys,'val','lower','upper'),with=FALSE]
 stopifnot(!anyDuplicated(x[,..keys]))
 setnames(x,c('val','lower','upper'),paste0(label,c('','_lo','_hi')))
 d <- merge(d,x,by=keys,all.x=TRUE,sort=FALSE)
}
stopifnot(!anyDuplicated(d[,..keys]))
missing <- d[is.na(YLL)|is.na(YLD),.(location_id,location_name,sex_id,age_id,age_name,year,DALY,YLD,YLL)]
fwrite(missing,file.path(out,'QA_missing_component_rows.csv'))
# No blanket imputation: only infer a missing nonnegative YLL from exact zero DALY and YLD.
d[,YLL_zero_inferred:=is.na(YLL) & age_id %in% c(1,6,7) & DALY==0 & DALY_lo==0 & DALY_hi==0 & !is.na(YLD) & YLD==0 & YLD_lo==0 & YLD_hi==0]
stopifnot(all(!is.na(d$YLD)),all(!is.na(d$YLL)|d$YLL_zero_inferred))
d[YLL_zero_inferred==TRUE,c('YLL','YLL_lo','YLL_hi'):=list(0,0,0)]
d[,`:=`(components=YLL+YLD,residual=DALY-YLL-YLD)]
d[,relative_error:=abs(residual)/pmax(abs(DALY),1)]
fwrite(d[,c(keys,'location_name','age_name','DALY','YLL','YLD','residual','relative_error','YLL_zero_inferred'),with=FALSE],file.path(out,'QA_DALY_component_identity.csv'))
cat('Maximum DALY identity absolute error:',max(abs(d$residual)),'relative:',max(d$relative_error),'\n')
stopifnot(max(d$relative_error)<1e-6)
# Select only countries for the substantive EU-28 analysis.
d <- d[location_id %in% ids]
stopifnot(nrow(d)==28*20*3*34)
allcheck <- d[,.(age_sum=sum(DALY)),by=.(location_id,sex_id,year)]
a <- old[location_id %in% ids & age_id==22,.(location_id,sex_id,year,all_age_DALY=val)]
allcheck <- merge(allcheck,a,by=c('location_id','sex_id','year'))
allcheck[,relative_error:=abs(age_sum-all_age_DALY)/pmax(all_age_DALY,1)]
fwrite(allcheck,file.path(out,'QA_age_sum_vs_original_All_ages.csv'))
stopifnot(nrow(allcheck)==2856,max(allcheck$relative_error)<1e-6)
g <- fread('data/gdp.csv')[year==2023,.(location_name=country,GDPpc=NY.GDP.PCAP.PP.CD,GDP=NY.GDP.MKTP.PP.CD)]
g[location_name=='Slovak Republic',location_name:='Slovakia']
h <- fread('data/HALE.csv')[year==2023 & age_name=='All ages' & metric_name=='Years' & location_id %in% ids,.(location_id,sex_id,HALE=val)]
stopifnot(nrow(h)==84,!anyDuplicated(h[,.(location_id,sex_id)]))
meta <- merge(eu,g,by='location_name');stopifnot(nrow(meta)==28,all(is.finite(meta$GDP)),all(meta$GDP>0))
# Preserve R0 subregions solely for like-for-like sensitivity comparison; not an official regrouping.
reg <- fread(file.path(root,'results_VLW','pancreatic_EU_IE1','T2_country_VLW_2023.csv'))[,.(location_name=Country,subregion=Subregion)]
meta <- merge(meta,reg,by='location_name')
d <- merge(d,meta[,.(location_id,GDPpc,GDP,subregion)],by='location_id');d <- merge(d,h,by=c('location_id','sex_id'))
stopifnot(nrow(d)==57120,all(is.finite(d$HALE)),all(d$HALE>0))
d[,VSLY:=13.2e6*(GDPpc/82304.62)/(HALE/2)]
d[,`:=`(VLW_DALY=VSLY*DALY/1e9,VLW_YLL=VSLY*YLL/1e9,VLW_YLD=VSLY*YLD/1e9)]
metrics <- c('DALY','YLL','YLD','VLW_DALY','VLW_YLL','VLW_YLD')
country <- d[,lapply(.SD,sum),by=.(location_id,location_name,subregion,sex_id,sex_name,year),.SDcols=metrics]
country <- merge(country,meta[,.(location_id,GDP)],by='location_id')
# Aggregate country totals once per stratum; never count GDP once per age.
agg <- function(x,bycols) {
 z <- x[,lapply(.SD,sum),by=bycols,.SDcols=c(metrics,'GDP')]
 z[,`:=`(YLD_share_pct=100*YLD/DALY,YLL_share_pct=100*YLL/DALY,
          VLW_reduction_pct=100*(VLW_DALY-VLW_YLL)/VLW_DALY,
          VLW_DALY_GDP_pct=100*VLW_DALY*1e9/GDP,VLW_YLL_GDP_pct=100*VLW_YLL*1e9/GDP)]
 z
}
country <- agg(country,c('location_id','location_name','subregion','sex_id','sex_name','year'))
eu_trend <- agg(country,c('sex_id','sex_name','year'))
sub_trend <- agg(country,c('subregion','sex_id','sex_name','year'))
age <- d[,lapply(.SD,sum),by=.(sex_id,sex_name,age_id,age_name,year),.SDcols=metrics]
age[,`:=`(YLD_share_pct=fifelse(DALY>0,100*YLD/DALY,NA_real_),VLW_reduction_pct=fifelse(VLW_DALY>0,100*VLW_YLD/VLW_DALY,NA_real_))]
# Independent reference check: original all-age EU VLW archive.
reference <- fread(file.path(root,'results_VLW','forecast_EU_IE1.0','Forecast_EU_raw.csv'))[period=='Historical',.(year,original_VLW=VLW)]
check <- merge(eu_trend[sex_id==3,.(year,VLW_DALY)],reference,by='year')
check[,error:=VLW_DALY-original_VLW];fwrite(check,file.path(out,'QA_VLW_vs_original.csv'))
stopifnot(max(abs(check$error))<1e-6)
# Rank stability for absolute and GDP-normalized results.
c2023 <- country[year==2023 & sex_id==3]
c2023[,`:=`(rank_DALY=frank(-VLW_DALY,ties.method='min'),rank_YLL=frank(-VLW_YLL,ties.method='min'),rank_GDP_DALY=frank(-VLW_DALY_GDP_pct,ties.method='min'),rank_GDP_YLL=frank(-VLW_YLL_GDP_pct,ties.method='min'))]
c2023[,`:=`(rank_change=rank_YLL-rank_DALY,rank_GDP_change=rank_GDP_YLL-rank_GDP_DALY)]
setorder(c2023,rank_DALY)
fwrite(eu_trend[order(sex_id,year)],file.path(out,'Table_S6A_EU28_annual_decomposition.csv'))
fwrite(c2023,file.path(out,'Table_S6B_country_2023_VLW_comparison.csv'))
fwrite(sub_trend[order(subregion,sex_id,year)],file.path(out,'Table_S6C_subregion_annual_decomposition.csv'))
fwrite(age[year==2023][order(sex_id,age_id)],file.path(out,'Table_S6D_age_sex_2023.csv'))
fwrite(country[order(location_id,sex_id,year)],file.path(out,'Table_S6E_country_annual_decomposition.csv'))
fwrite(d[,c('location_id','location_name','sex_id','sex_name','age_id','age_name','year','DALY','YLL','YLD','HALE','GDPpc','VSLY','VLW_DALY','VLW_YLL','VLW_YLD','YLL_zero_inferred'),with=FALSE],file.path(out,'Analysis_age_country_year_inputs_and_VLW.csv'))
# IE sensitivity of metric substitution, maintaining the existing transfer framework.
ies <- rbindlist(lapply(c(.5,1,1.5),function(ie){
 x <- copy(d[sex_id==3 & year==2023]);x[,v:=13.2e6*(GDPpc/82304.62)^ie/(HALE/2)]
 x[,.(IE=ie,VLW_DALY=sum(v*DALY/1e9),VLW_YLL=sum(v*YLL/1e9),VLW_YLD=sum(v*YLD/1e9))]
}))
ies[,reduction_pct:=100*VLW_YLD/VLW_DALY];fwrite(ies,file.path(out,'Table_S6F_IE_metric_sensitivity_2023.csv'))
# Compact report evidence; no fabricated 95% interval for shares or aggregate components.
summary <- data.table(item=c('new_file_rows','matched_country_age_rows','inferred_zero_YLL_country_rows','max_component_relative_error_country','max_age_sum_relative_error','max_VLW_absolute_error','rank_spearman_VLW','rank_spearman_GDP','countries_changed_absolute_rank','countries_changed_GDP_rank','max_YLD_pct_country_year_sex','min_YLD_pct_country_year_sex'),value=c(nrow(new),nrow(d),sum(d$YLL_zero_inferred),max(d$relative_error),max(allcheck$relative_error),max(abs(check$error)),cor(c2023$VLW_DALY,c2023$VLW_YLL,method='spearman'),cor(c2023$VLW_DALY_GDP_pct,c2023$VLW_YLL_GDP_pct,method='spearman'),sum(c2023$rank_change!=0),sum(c2023$rank_GDP_change!=0),max(country$YLD_share_pct),min(country$YLD_share_pct)))
fwrite(summary,file.path(out,'QA_summary.csv'))
cat('EU Both key years:\n');print(eu_trend[sex_id==3 & year %in% c(1990,2023)])
print(summary)
saveRDS(list(eu=eu_trend,country=c2023,sub=sub_trend,age=age[year==2023],ie=ies),file.path(out,'plot_data.rds'))
capture.output(sessionInfo(),file=file.path(out,'sessionInfo.txt'))
