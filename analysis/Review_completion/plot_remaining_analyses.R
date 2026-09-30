suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
b<-file.path(getwd(),'analysis','Review_completion')
read<-function(d,f)fread(file.path(b,d,f))
blue<-'#1B4F72';green<-'#117A65';red<-'#C0392B';purple<-'#6C3483'
th<-function()theme_minimal(base_size=8,base_family='Helvetica')+theme(axis.title=element_text(size=8,face='bold',colour='grey20'),axis.text=element_text(size=7,colour='grey30'),legend.text=element_text(size=6.5),legend.title=element_text(size=7),panel.grid.major=element_line(colour='grey92',linewidth=.3),panel.grid.minor=element_blank(),axis.line=element_line(colour='grey30',linewidth=.3),axis.ticks=element_line(colour='grey30',linewidth=.3),legend.position='bottom',plot.tag=element_text(size=11),plot.margin=margin(6,7,5,5))
save<-function(p,d,n,w=7,h=5.5){ggsave(file.path(b,d,paste0(n,'.pdf')),p,width=w,height=h,device=grDevices::cairo_pdf);ggsave(file.path(b,d,paste0(n,'.png')),p,width=w,height=h,dpi=300,device=ragg::agg_png,bg='white')}
# Age annualization diagnostic, no alternative treated as the correct willingness-to-pay model.
a<-read('R1_3_R2_4_5_age','Table_S8B_EU_age_valuation.csv')[sex_id==3];ord<-c(1,6:20,30:32,235);a[,index:=match(age_id,ord)];setorder(a,index)
ac<-read('R1_3_R2_4_5_age','Table_S8A_country_age_valuation.csv')[sex_id==3]
h<-ac[,.(remaining=median(H_age,na.rm=TRUE),half=median(HALE/2),full=median(HALE)),by=.(age_id,age_name)];h[,index:=match(age_id,ord)]
h<-melt(h,id.vars=c('age_id','age_name','index'),variable.name='Basis');h<-h[is.finite(value)]
ax<-scale_x_continuous(limits=c(1,20),breaks=1:20,labels=sub(' years','',a$age_name))
ath<-th()+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
p1<-ggplot(h,aes(index,value,colour=Basis))+geom_line(linewidth=.7)+ax+scale_colour_manual(values=c(remaining=red,half=blue,full=green),labels=c(full='All-age H',half='All-age H/2',remaining='Age-specific h'),name=NULL)+labs(x='Age (years)',y='Median denominator (years)',tag='a')+ath
l<-melt(a,id.vars=c('age_id','age_name','index'),measure.vars=c('base','full_H','remaining_H','annuity_Hhalf'),variable.name='Basis')
cols<-c(base=blue,full_H=green,remaining_H=red,annuity_Hhalf=purple);labels<-c(base='H/2 (base)',full_H='H',remaining_H='Age-specific h',annuity_Hhalf='Annuity over H/2')
p2<-ggplot(l,aes(index,value,colour=Basis))+geom_line(linewidth=.7)+ax+scale_colour_manual(values=cols,labels=labels,name=NULL)+guides(colour=guide_legend(nrow=2))+labs(x='Age (years)',y='VLW (billion USD, 2023 PPP)',tag='b')+ath
p3<-ggplot(a[is.finite(remaining_to_base)],aes(index,remaining_to_base))+geom_line(colour=red,linewidth=.7)+geom_hline(yintercept=1,linetype=2,colour='grey50')+ax+labs(x='Age (years)',y='Age-specific / base VLW',tag='c')+ath
at<-read('R1_3_R2_4_5_age','Table_S8D_EU_total_valuation.csv')[sex_id==3];at<-melt(at,measure.vars=names(cols),variable.name='Basis')
at[,Basis:=factor(Basis,levels=rev(names(cols)))]
p4<-ggplot(at,aes(value,Basis,fill=Basis))+geom_col(width=.6)+geom_text(aes(label=sprintf('%.1f',value)),hjust=-.1,size=2.5)+scale_fill_manual(values=cols,guide='none')+scale_y_discrete(labels=labels)+scale_x_continuous(expand=expansion(mult=c(0,.2)))+labs(x='Total VLW (billion USD)',y=NULL,tag='d')+th()
save((p1|p2)/(p3|p4),'R1_3_R2_4_5_age','Figure_R1_3_age_annualization',7.5,6)
# Discount sensitivity: optional multiplier, not clinical life expectancy.
d<-read('R2_7_discount','Table_S9_discount_rate_horizon_grid.csv')[sex_id==3]
d[,fraction:=factor(horizon_fraction)]
p1<-ggplot(d,aes(discount_rate*100,VLW_discounted,colour=fraction,group=fraction))+geom_line(linewidth=.7)+geom_point(size=1)+scale_colour_manual(values=c(blue,green,red,purple,'#5D6D7E'),name='Horizon / HALE')+labs(x='Annual discount rate (%)',y='Discounted VLW (billion USD)',tag='a')+th()
p2<-ggplot(d,aes(factor(discount_rate*100),factor(horizon_fraction),fill=reduction_pct))+geom_tile(colour='white')+geom_text(aes(label=sprintf('%.1f%%',reduction_pct)),size=2.4)+scale_fill_gradient(low='white',high='#A6BDDB',guide='none')+labs(x='Annual discount rate (%)',y='Horizon / HALE',tag='b')+th()
save(p1|p2,'R2_7_discount','Figure_R2_7_discount_sensitivity',7,3.3)
# Probabilistic sensitivity, all intervals explicitly simulation intervals.
r<-read('R2_8_PSA','Table_S12A_EU_subregion_PSA.csv')[year==2023 & region=='EU-28']
scenario_labels<-c(health_only='DALY + HALE; fixed valuation',joint_triangular='Joint: triangular assumptions',joint_rho0='Joint: between-country rho = 0',joint_rho1='Joint: between-country rho = 1',joint_uniform='Joint: uniform assumptions',joint_lognormal='Joint: lognormal health inputs',joint_kappa_negative='Joint: DALY-HALE kappa = -0.5',joint_kappa_positive='Joint: DALY-HALE kappa = +0.5')
rr<-read('R2_8_PSA','EU2023_draws_by_scenario.csv');rr<-melt(rr,measure.vars=c('health_only','joint_triangular','joint_uniform'),variable.name='Scenario')
p1<-ggplot(rr,aes(value,colour=Scenario))+geom_density(linewidth=.7)+geom_vline(xintercept=634.135197,linetype=2,colour='grey40')+scale_colour_manual(values=c(health_only=blue,joint_triangular=green,joint_uniform=red),labels=c('Health only','Joint triangular','Joint uniform'),name=NULL)+labs(x='2023 VLW (billion USD)',y='Simulation density',tag='a')+th()+guides(colour=guide_legend(nrow=2))
r[,scenario:=factor(scenario,levels=rev(names(scenario_labels)))]
p2<-ggplot(r,aes(median,scenario))+geom_segment(aes(x=lower,xend=upper,yend=scenario),colour=blue,linewidth=.65)+geom_point(colour=green,size=1.5)+geom_vline(xintercept=634.135197,linetype=2,colour='grey50')+scale_y_discrete(labels=scenario_labels)+labs(x='Median and 95% simulation interval',y=NULL,tag='b')+th()+theme(axis.text.y=element_text(size=6))
st<-read('R2_8_PSA','Table_S12E_sequential_uncertainty.csv');st[,width:=upper-lower];st[,stage:=factor(stage,levels=rev(stage))]
p3<-ggplot(st,aes(width,stage))+geom_col(fill=purple,width=.6)+labs(x='95% simulation interval width (billion USD)',y=NULL,tag='c')+th()+theme(axis.text.y=element_text(size=6))
cv<-read('R2_8_PSA','Table_S12C_convergence.csv')[scenario=='joint_triangular'];cv<-melt(cv,id.vars='n',measure.vars=c('lower','median','upper'),variable.name='Quantile')
p4<-ggplot(cv,aes(n/1000,value,colour=Quantile))+geom_line(linewidth=.7)+geom_point(size=1)+scale_colour_manual(values=c(lower=blue,median=green,upper=red),labels=c('2.5%','50%','97.5%'),name=NULL)+labs(x='Monte Carlo draws (thousands)',y='VLW quantile (billion USD)',tag='d')+th()
save((p1|p2)/(p3|p4),'R2_8_PSA','Figure_S15_PSA_overview',9.5,6.3)
cr<-read('R2_8_PSA','Table_S12B_country_2023_PSA.csv')[scenario %in% c('health_only','joint_triangular')]
lev<-cr[scenario=='health_only'][order(base_VLW)]$location_name;cr[,Country:=factor(location_name,levels=lev)]
cr[,y:=as.numeric(Country)+ifelse(scenario=='health_only',-.13,.13)]
p1<-ggplot(cr,aes(median,y,colour=scenario))+geom_segment(aes(x=lower,xend=upper,yend=y),linewidth=.55)+geom_point(size=1)+scale_y_continuous(breaks=1:28,labels=lev)+scale_colour_manual(values=c(health_only=blue,joint_triangular=red),labels=c('Health only','Joint triangular'),name=NULL)+labs(x='VLW (billion USD, 2023 PPP)',y=NULL,tag='a')+th()
p2<-ggplot(cr,aes(median_GDP_pct,y,colour=scenario))+geom_segment(aes(x=lower_GDP_pct,xend=upper_GDP_pct,yend=y),linewidth=.55)+geom_point(size=1)+scale_y_continuous(breaks=1:28,labels=rep('',28))+scale_colour_manual(values=c(health_only=blue,joint_triangular=red),labels=c('Health only','Joint triangular'),name=NULL)+labs(x='VLW relative to 2023 GDP (%)',y=NULL,tag='b')+th()
save(p1|p2,'R2_8_PSA','Figure_S16_country_PSA',7,7.3)
# Corrected original Figure1: all four panels share 28-country geography; point estimates.
f<-read('Manuscript_alignment','Table_Figure1_corrected_EU28.csv');plots<-list();ys<-c('DALYs','Deaths','VLW','VLW_GDP');yl<-c('DALYs (thousands)','Deaths (thousands)','VLW (billion USD, 2023 PPP)','VLW relative to 2023 GDP (%)')
for(i in 1:4){ff<-copy(f);ff[,y:=get(ys[i])/if(i<=2)1000 else 1];plots[[i]]<-ggplot(ff,aes(year,y))+geom_line(colour=c(blue,red,green,purple)[i],linewidth=.8)+scale_x_continuous(breaks=c(1990,2000,2010,2023))+labs(x='Year',y=yl[i],tag=letters[i])+th()}
save((plots[[1]]|plots[[2]])/(plots[[3]]|plots[[4]]),'Manuscript_alignment','Figure_1_R1_corrected_EU28',7,5.5)
# Existing US-reference sensitivity visual; no new European anchors.
v<-fread(file.path(getwd(),'analysis/Response_support/Table_R1_reference_VSL_sensitivity.csv'));v[,IE:=factor(IE)]
p1<-ggplot(v,aes(VSL_reference_million_2023_int_dollars,VLW_billion,colour=IE))+geom_line(linewidth=.7)+geom_point(size=1.3)+scale_colour_manual(values=c(blue,green,red),name='Income elasticity')+labs(x='US reference VSL (million USD, 2023)',y='VLW (billion USD, 2023 PPP)',tag='a')+th()
p2<-ggplot(v,aes(VSL_reference_million_2023_int_dollars,VLW_GDP_percent,colour=IE))+geom_line(linewidth=.7)+geom_point(size=1.3)+scale_colour_manual(values=c(blue,green,red),name='Income elasticity')+labs(x='US reference VSL (million USD, 2023)',y='VLW relative to 2023 GDP (%)',tag='b')+th()
save(p1|p2,'Manuscript_alignment','Figure_R1_2_US_reference_sensitivity',7,3)
cat('Six PDF/PNG figure pairs saved.\n')
