# Completes remaining requested sensitivity analyses; retains DOT/PPP/fixed-2023 GDP.
# All PSA distributions/dependence choices are assumptions, not recovered GBD draws.
suppressPackageStartupMessages(library(data.table))
root<-getwd();base<-file.path(root,'analysis','Review_completion')
for(z in c('R1_3_R2_4_5_age','R2_7_discount','R2_8_PSA','R1_4_ratio','Manuscript_alignment','qa'))dir.create(file.path(base,z),recursive=TRUE,showWarnings=FALSE)
savecsv<-function(x,dir,name)fwrite(x,file.path(base,dir,name))
eu<-fread('EU28_location_list.csv')[group_type=='EU-28 member',.(location_id,location_name)]
ids<-eu$location_id;old<-fread('data/merged.csv');g<-fread('data/gdp.csv')[year==2023,.(location_name=country,GDPpc=NY.GDP.PCAP.PP.CD,GDP=NY.GDP.MKTP.PP.CD)]
g[location_name=='Slovak Republic',location_name:='Slovakia'];g<-merge(eu,g,by='location_name')
hall<-fread('data/HALE.csv')[year==2023 & metric_name=='Years' & location_id %in% ids]
h<-hall[age_id==22,.(location_id,sex_id,H=val,Hlo=lower,Hhi=upper)]
reg<-fread('results_VLW/pancreatic_EU_IE1/T2_country_VLW_2023.csv')[,.(location_name=Country,subregion=Subregion)]
g<-merge(g,reg,by='location_name');stopifnot(nrow(g)==28,nrow(h)==84)
d<-old[cause_id==456 & location_id %in% ids & measure_id==2 & metric_id==1 & age_id==22,.(location_id,sex_id,sex_name,year,D=val,Dlo=lower,Dhi=upper)]
d<-merge(merge(d,g,by='location_id'),h,by=c('location_id','sex_id'))
stopifnot(nrow(d)==2856,!anyDuplicated(d[,.(location_id,sex_id,year)]),all(d$Dlo>0),all(d$Dhi>=d$D),all(d$Dlo<=d$D),all(d$Hlo>0))
d[,`:=`(VSL=13.2e6*GDPpc/82304.62,VSLY=13.2e6*GDPpc/82304.62/(H/2))];d[,VLW:=VSLY*D/1e9]
savecsv(d,'R2_8_PSA','Input_country_year_sex_summary_bounds.csv')
# Age valuation: apply all four explicit annualization assumptions to the same counts.
a<-fread(file.path(root,'analysis/R1_1_YLL_YLD/Analysis_age_country_year_inputs_and_VLW.csv'))[year==2023]
ageh<-hall[,.(location_id,sex_id,age_id,H_age=val)]
a<-merge(a,ageh,by=c('location_id','sex_id','age_id'),all.x=TRUE)
stopifnot(nrow(a)==1680,all(is.finite(a$H_age)|a$DALY==0))
a[,VSL:=13.2e6*GDPpc/82304.62]
a[,`:=`(base=VLW_DALY,full_H=VSL/HALE*DALY/1e9,
 remaining_H=fifelse(DALY==0,0,VSL/H_age*DALY/1e9),
 annuity_Hhalf=VSL/((1-1.03^(-HALE/2))/.03)*DALY/1e9)]
stopifnot(all(is.finite(a$remaining_H)),all(abs(a$full_H-a$base/2)<1e-9))
ms<-c('base','full_H','remaining_H','annuity_Hhalf')
savecsv(a,'R1_3_R2_4_5_age','Table_S8A_country_age_valuation.csv')
aa<-a[,lapply(.SD,sum),by=.(age_id,age_name,sex_id,sex_name),.SDcols=ms]
aa[,remaining_to_base:=fifelse(base>0,remaining_H/base,NA_real_)]
savecsv(aa,'R1_3_R2_4_5_age','Table_S8B_EU_age_valuation.csv')
ac<-a[,lapply(.SD,sum),by=.(location_id,location_name,sex_id,sex_name),.SDcols=ms]
ac[,remaining_to_base:=remaining_H/base]
savecsv(ac,'R1_3_R2_4_5_age','Table_S8C_country_total_valuation.csv')
at<-a[,lapply(.SD,sum),by=.(sex_id,sex_name),.SDcols=ms];at[,remaining_to_base:=remaining_H/base]
savecsv(at,'R1_3_R2_4_5_age','Table_S8D_EU_total_valuation.csv')
peaks<-rbindlist(lapply(ms,function(m)aa[,.(specification=m,age_peak=age_name[which.max(get(m))],VLW_peak=max(get(m))),by=.(sex_id,sex_name)]))
savecsv(peaks,'R1_3_R2_4_5_age','Table_S8E_age_peak_comparison.csv')
# No clinical remaining-life claim: 0.4*HALE is a chosen horizon coefficient.
d23<-d[year==2023]
disc<-rbindlist(lapply(c(.2,.4,.5,.6,1),function(k)rbindlist(lapply(c(0,.015,.03,.05),function(r){
 x<-copy(d23);x[,factor:=if(r==0)1 else (1-(1+r)^(-k*H))/(r*k*H)]
 x[,.(horizon_fraction=k,discount_rate=r,VLW=sum(VLW),VLW_discounted=sum(VLW*factor)),by=.(sex_id,sex_name)]
}))))
disc[,reduction_pct:=100*(1-VLW_discounted/VLW)]
savecsv(disc,'R2_7_discount','Table_S9_discount_rate_horizon_grid.csv')
stopifnot(abs(disc[sex_id==3 & horizon_fraction==.4 & discount_rate==.03]$VLW_discounted-426.490774199837)<1e-8)
# Numerical validation of income cancellation, 3 elasticities and all 28 countries.
ratio<-rbindlist(lapply(c(.5,1,1.5),function(e){
 x<-copy(d23[sex_id==3]);x[,population_implied:=GDP/GDPpc]
 x[,`:=`(IE=e,direct_pct=100*13.2e6*(GDPpc/82304.62)^e/(H/2)*D/GDP,
 algebra_pct=200*13.2e6/82304.62^e*(D/(population_implied*H))*GDPpc^(e-1))]
 x[,.(location_id,location_name,IE,GDPpc,population_implied,D,H,direct_pct,algebra_pct)]
}));ratio[,error:=direct_pct-algebra_pct];stopifnot(max(abs(ratio$error))<1e-12)
savecsv(ratio,'R1_4_ratio','Table_R1_4_income_cancellation.csv')
# Correct Figure1 input geography without changing monetary estimates.
main<-old[location_id %in% ids & cause_id==456 & age_id==22 & sex_id==3 & metric_id==1 & measure_id %in% c(1,2),.(val=sum(val)),by=.(year,measure_id)]
main<-dcast(main,year~measure_id,value.var='val');setnames(main,c('1','2'),c('Deaths','DALYs'))
v<-d[sex_id==3,.(VLW=sum(VLW),GDP=sum(GDP)),by=year];main<-merge(main,v,by='year');main[,VLW_GDP:=100*VLW*1e9/GDP]
savecsv(main,'Manuscript_alignment','Table_Figure1_corrected_EU28.csv')
rm(old);gc()
# PSA: fixed GDP; uncertain health summaries plus explicitly assumed VSL/IE distributions.
# Primary reconstruction: piecewise log-normal quantile map. Point estimate is a median assumption.
# The alternate log-normal matches endpoints and uses geometric midpoint, illustrating shape/center sensitivity.
M<-50000L;seed<-20260914L;set.seed(seed)
qmap<-function(z,v,lo,hi,shape='split'){
 if(length(z)==1L && length(v)>1L) z<-rep(z,length(v))
 if(shape=='split') exp(log(v)+ifelse(z<0,log(v/lo),log(hi/v))*z/qnorm(.975)) else exp((log(lo)+log(hi))/2+(log(hi)-log(lo))*z/(2*qnorm(.975)))
}
tri<-function(u,lo,mode,hi){c<-(mode-lo)/(hi-lo);ifelse(u<c,lo+sqrt(u*(hi-lo)*(mode-lo)),hi-sqrt((1-u)*(hi-lo)*(hi-mode)))}
# Shared random numbers across scenarios and years enable paired comparisons; not empirical dependence.
zD0<-rnorm(M);zH0<-rnorm(M);eD<-matrix(rnorm(M*28),M,28);eH<-matrix(rnorm(M*28),M,28)
uA<-runif(M);uE<-runif(M)
qstats<-function(v)c(mean=mean(v),median=median(v),lower=unname(quantile(v,.025)),upper=unname(quantile(v,.975)))
configs<-data.table(scenario=c('health_only','joint_triangular','joint_rho0','joint_rho1','joint_uniform','joint_lognormal','joint_kappa_negative','joint_kappa_positive'),rho=c(.5,.5,0,1,.5,.5,.5,.5),shape=c('split','split','split','split','split','lognormal','split','split'),prior=c('fixed','triangular','triangular','triangular','uniform','triangular','triangular','triangular'),kappa=c(0,0,0,0,0,0,-.5,.5))
all_res<-list();country_res<-list();batch_res<-list();stage_res<-list();eu_draws<-list();rank_res<-list();conv<-list()
for(sc in seq_len(nrow(configs))){
 cfg<-configs[sc];cat('PSA scenario',cfg$scenario,'\n')
 zd<-sqrt(cfg$rho)*zD0+sqrt(1-cfg$rho)*eD
 zh0<-sqrt(cfg$rho)*zH0+sqrt(1-cfg$rho)*eH
 zh<-cfg$kappa*zd+sqrt(1-cfg$kappa^2)*zh0
 A<-if(cfg$prior=='fixed')rep(13.2e6,M) else if(cfg$prior=='uniform') (7+7*uA)*1e6 else tri(uA,7,13.2,14)*1e6
 IE<-if(cfg$prior=='fixed')rep(1,M) else if(cfg$prior=='uniform') .5+uE else tri(uE,.5,1,1.5)
 years<-if(cfg$scenario %in% c('health_only','joint_triangular'))1990:2023 else 2023
 for(yr in years){
  xx<-d[year==yr & sex_id==3][order(location_id)];stopifnot(nrow(xx)==28)
  Dm<-qmap(zd,rep(xx$D,each=M),rep(xx$Dlo,each=M),rep(xx$Dhi,each=M),cfg$shape)
  Hm<-qmap(zh,rep(xx$H,each=M),rep(xx$Hlo,each=M),rep(xx$Hhi,each=M),cfg$shape)
  inc<-exp(outer(IE,log(xx$GDPpc/82304.62)))
  W<-Dm/Hm*2*inc*A/1e9
  # Matrix recycling repeats draw-specific common parameters down each country column.
  stopifnot(all(is.finite(W)),all(W>=0))
  sums<-list('EU-28'=rowSums(W))
  for(regn in unique(xx$subregion))sums[[regn]]<-rowSums(W[,xx$subregion==regn,drop=FALSE])
  for(regn in names(sums)){
   gd<-sum(xx[subregion==regn | regn=='EU-28']$GDP)
   basev<-sum(xx[subregion==regn | regn=='EU-28']$VLW)
   st<-qstats(sums[[regn]])
   all_res[[length(all_res)+1]]<-data.table(scenario=cfg$scenario,year=yr,region=regn,base_VLW=basev,t(as.matrix(st)),GDP=gd)
  }
  if(yr==2023){
   eu_draws[[cfg$scenario]]<-sums[['EU-28']]
   cr<-rbindlist(lapply(1:28,function(i)data.table(location_id=xx$location_id[i],location_name=xx$location_name[i],GDP=xx$GDP[i],base_VLW=xx$VLW[i],t(as.matrix(qstats(W[,i]))))))
   cr[,scenario:=cfg$scenario];country_res[[sc]]<-cr
   for(n in c(5000,10000,20000,50000))conv[[length(conv)+1]]<-data.table(scenario=cfg$scenario,n=n,t(as.matrix(qstats(sums[['EU-28']][1:n]))))
   for(bi in 1:10)batch_res[[length(batch_res)+1]]<-data.table(scenario=cfg$scenario,batch=bi,t(as.matrix(qstats(sums[['EU-28']][((bi-1)*5000+1):(bi*5000)]))))
   if(cfg$scenario=='joint_triangular'){
    # Report rank uncertainty; common VSL factor cancels in within-draw ranking.
    rankabs<-max.col(W,ties.method='first');rankgdp<-max.col(sweep(W,2,xx$GDP,'/'),ties.method='first')
    rank_res[[1]]<-data.table(location_id=xx$location_id,location_name=xx$location_name,P_highest_VLW=tabulate(rankabs,28)/M,P_highest_GDP_ratio=tabulate(rankgdp,28)/M)
    stages<-list('DALY only'=rowSums(sweep(Dm,2,xx$VSLY/1e9,'*')),
                 'DALY + HALE'=rowSums(sweep(Dm/Hm*2,2,13.2e6*xx$GDPpc/82304.62/1e9,'*')),
                 'DALY + HALE + VSL'=rowSums(sweep(Dm/Hm*2,2,xx$GDPpc/82304.62/1e9,'*')*A),
                 'DALY + HALE + VSL + IE'=rowSums(W))
    for(sn in names(stages))stage_res[[length(stage_res)+1]]<-data.table(stage=sn,t(as.matrix(qstats(stages[[sn]]))))
   }
  }
 }
}
res<-rbindlist(all_res);res[,`:=`(median_GDP_pct=median*1e11/GDP,lower_GDP_pct=lower*1e11/GDP,upper_GDP_pct=upper*1e11/GDP)]
cr<-rbindlist(country_res);cr[,`:=`(median_GDP_pct=median*1e11/GDP,lower_GDP_pct=lower*1e11/GDP,upper_GDP_pct=upper*1e11/GDP)]
cv<-rbindlist(conv);bt<-rbindlist(batch_res)
# Between-batch quantile variability estimates Monte Carlo error, not model uncertainty.
mcse<-bt[,.(median_MCSE=sd(median)/sqrt(.N),lower_MCSE=sd(lower)/sqrt(.N),upper_MCSE=sd(upper)/sqrt(.N)),by=scenario]
mcse<-merge(mcse,res[year==2023 & region=='EU-28',.(scenario,median,lower,upper)],by='scenario');mcse[,max_relative_MCSE:=pmax(median_MCSE/median,lower_MCSE/lower,upper_MCSE/upper)]
savecsv(configs,'R2_8_PSA','Table_S11A_PSA_scenarios.csv')
savecsv(res,'R2_8_PSA','Table_S12A_EU_subregion_PSA.csv')
savecsv(cr,'R2_8_PSA','Table_S12B_country_2023_PSA.csv')
savecsv(cv,'R2_8_PSA','Table_S12C_convergence.csv')
savecsv(mcse,'R2_8_PSA','Table_S12D_Monte_Carlo_error.csv')
savecsv(rbindlist(stage_res),'R2_8_PSA','Table_S12E_sequential_uncertainty.csv')
savecsv(rank_res[[1]],'R2_8_PSA','Table_S12F_rank_probabilities.csv')
fwrite(as.data.table(eu_draws),file.path(base,'R2_8_PSA','EU2023_draws_by_scenario.csv'))
# Quantile-map deterministic boundary test using exact normal quantiles.
xx<-d[year==2023 & sex_id==3]
stopifnot(max(abs(qmap(qnorm(.025),xx$D,xx$Dlo,xx$Dhi)-xx$Dlo))<1e-7,max(abs(qmap(qnorm(.975),xx$D,xx$Dlo,xx$Dhi)-xx$Dhi))<1e-7)
savecsv(data.table(check=c('country_year_sex_rows','age_2023_rows','zero_DALY_missing_age_HALE','ratio_identity_max_abs','primary_PSA_draws','seed','max_quantile_relative_MCSE'),value=c(nrow(d),nrow(a),sum(is.na(a$H_age)&a$DALY==0),max(abs(ratio$error)),M,seed,max(mcse$max_relative_MCSE))),'qa','checks.csv')
capture.output(sessionInfo(),file=file.path(base,'sessionInfo.txt'))
cat('AGE TOTALS\n');print(at);cat('PSA 2023 EU\n');print(res[year==2023 & region=='EU-28']);cat('MC error\n');print(mcse)
