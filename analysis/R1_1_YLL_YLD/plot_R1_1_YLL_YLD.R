# Publication figures using the R0 theme, typography and palette.
suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork);library(scales)})
out <- file.path(getwd(),'analysis','R1_1_YLL_YLD')
x <- readRDS(file.path(out,'plot_data.rds'))
blue <- '#1B4F72'; green <- '#117A65';red <- '#C0392B';purple <- '#6C3483'
sexcol <- c(Male='#2E86C1',Female=red,Both='#5D6D7E')
regcol <- c('Western Europe'=blue,'Central Europe'=green,'Eastern Europe'=red)
th <- function(){theme_minimal(base_size=8,base_family='Helvetica')+theme(axis.title=element_text(size=8,face='bold',colour='grey20'),axis.text=element_text(size=7,colour='grey30'),legend.title=element_text(size=7,face='bold'),legend.text=element_text(size=7),legend.key.size=unit(.35,'cm'),panel.grid.major=element_line(colour='grey92',linewidth=.3),panel.grid.minor=element_blank(),strip.text=element_text(size=8,face='bold'),plot.margin=margin(6,7,5,5,'pt'),plot.tag=element_text(size=11),axis.line=element_line(colour='grey30',linewidth=.3),axis.ticks=element_line(colour='grey30',linewidth=.3),legend.position='bottom')}
xt <- scale_x_continuous(breaks=c(1990,2000,2010,2023))
ys <- scale_y_continuous(labels=label_number(accuracy=.1),expand=expansion(mult=c(0,.07)))
savefig <- function(p,name,w,h){
 ggsave(file.path(out,paste0(name,'.pdf')),p,width=w,height=h,device=grDevices::cairo_pdf)
 ggsave(file.path(out,paste0(name,'.png')),p,width=w,height=h,dpi=300,device=ragg::agg_png,bg='white')
}
e <- x$eu[sex_id==3]
long <- melt(e,id.vars='year',measure.vars=c('DALY','YLL'),variable.name='Metric',value.name='value')
long[,Metric:=factor(Metric,levels=c('DALY','YLL'),labels=c('DALYs','YLLs'))]
p1 <- ggplot(long,aes(year,value/1e3,colour=Metric,linetype=Metric))+geom_line(linewidth=.8)+scale_colour_manual(values=c(DALYs=blue,YLLs=red),name=NULL)+scale_linetype_manual(values=c(DALYs='solid',YLLs='dashed'),name=NULL)+xt+scale_y_continuous(labels=label_number(accuracy=1,big.mark=','))+labs(x='Year',y='Health loss (thousand years)',tag='a')+th()
p2 <- ggplot(x$eu,aes(year,YLD_share_pct,colour=sex_name))+geom_line(linewidth=.8)+scale_colour_manual(values=sexcol,name='Sex')+xt+scale_y_continuous(limits=c(0,1.3),breaks=c(0,.4,.8,1.2))+labs(x='Year',y='YLDs / DALYs (%)',tag='b')+th()
lv <- melt(e,id.vars='year',measure.vars=c('VLW_DALY','VLW_YLL'),variable.name='Metric',value.name='VLW')
lv[,Metric:=factor(Metric,levels=c('VLW_DALY','VLW_YLL'),labels=c('DALY-based','YLL-based'))]
p3 <- ggplot(lv,aes(year,VLW,colour=Metric,linetype=Metric))+geom_line(linewidth=.8)+scale_colour_manual(values=c('DALY-based'=green,'YLL-based'=red),name=NULL)+scale_linetype_manual(values=c('DALY-based'='solid','YLL-based'='dashed'),name=NULL)+xt+scale_y_continuous(labels=label_number(accuracy=1))+labs(x='Year',y='VLW (billion USD, 2023 PPP)',tag='c')+th()
p4 <- ggplot(e,aes(year,VLW_YLD))+geom_line(colour=purple,linewidth=.8)+xt+scale_y_continuous(limits=c(0,7),breaks=c(0,2,4,6))+labs(x='Year',y='YLD-based VLW (billion USD)',tag='d')+th()
savefig((p1|p2)/(p3|p4),'Figure_R1_1_1_EU_decomposition_4panel',7,5.6)
# Country comparison: dollar estimates side by side plus the small proportional reduction.
c <- copy(x$country);c[,Country:=factor(location_name,levels=rev(location_name))]
p1 <- ggplot(c,aes(y=Country))+geom_segment(aes(x=VLW_YLL,xend=VLW_DALY,yend=Country),colour='grey50',linewidth=.65)+geom_point(aes(x=VLW_DALY,colour='DALY-based'),shape=1,size=2,stroke=.55)+geom_point(aes(x=VLW_YLL,colour='YLL-based'),shape=16,size=1.15)+scale_colour_manual(values=c('DALY-based'=green,'YLL-based'=red),name=NULL)+scale_x_continuous(expand=expansion(mult=c(.02,.08)))+labs(x='VLW (billion USD, 2023 PPP)',y=NULL,tag='a')+th()
p2 <- ggplot(c,aes(VLW_reduction_pct,Country,fill=subregion))+geom_col(width=.65)+geom_text(aes(label=sprintf('%.2f',VLW_reduction_pct)),hjust=-.15,size=2.2,colour='grey25')+scale_fill_manual(values=regcol,name=NULL)+scale_x_continuous(limits=c(0,1.85),breaks=c(0,.5,1,1.5),expand=expansion(mult=c(0,0)))+labs(x='Reduction with YLL-based VLW (%)',y=NULL,tag='b')+th()+theme(axis.text.y=element_blank(),legend.text=element_text(size=6),legend.key.size=unit(.22,'cm'))+guides(fill=guide_legend(nrow=2))
savefig(p1|p2,'Figure_R1_1_2_country_VLW_comparison',7,7.1)
# Age analyses: exclude zero-denominator ages from ratio panel explicitly.
a <- copy(x$age);ageids<-c(1,6:20,30:32,235)
a[,age_index:=match(age_id,ageids)];a[,age_label:=sub(' years','',age_name)]
ax <- scale_x_continuous(limits=c(1,20),breaks=1:20,labels=a[sex_id==3][order(age_index)]$age_label)
atheme <- th()+theme(axis.text.x=element_text(angle=60,hjust=1,size=6))
p1<-ggplot(a,aes(age_index,DALY/1e3,colour=sex_name))+geom_line(linewidth=.7)+geom_point(size=.7)+ax+scale_colour_manual(values=sexcol,name='Sex')+labs(x='Age group (years)',y='DALYs (thousands)',tag='a')+atheme
p2<-ggplot(a,aes(age_index,YLD,colour=sex_name))+geom_line(linewidth=.7)+geom_point(size=.7)+ax+scale_colour_manual(values=sexcol,name='Sex')+scale_y_continuous(labels=label_number(accuracy=1,big.mark=','))+labs(x='Age group (years)',y='YLDs (years)',tag='b')+atheme
p3<-ggplot(a[!is.na(YLD_share_pct)],aes(age_index,YLD_share_pct,colour=sex_name))+geom_line(linewidth=.7)+geom_point(size=.7)+ax+scale_colour_manual(values=sexcol,name='Sex')+scale_y_continuous(limits=c(0,NA))+labs(x='Age group (years)',y='YLDs / DALYs (%)',tag='c')+atheme
p4<-ggplot(a,aes(age_index,VLW_YLD,colour=sex_name))+geom_line(linewidth=.7)+geom_point(size=.7)+ax+scale_colour_manual(values=sexcol,name='Sex')+labs(x='Age group (years)',y='YLD-based VLW (billion USD)',tag='d')+atheme
savefig(((p1|p2)/(p3|p4)),'Figure_R1_1_3_age_sex_decomposition',7.5,6)
# R0 study-defined subregions; matched geography for all components.
s <- x$sub[sex_id==3]
p1<-ggplot(s,aes(year,YLD_share_pct,colour=subregion))+geom_line(linewidth=.8)+xt+scale_y_continuous(limits=c(0,1.4),breaks=c(0,.4,.8,1.2))+scale_colour_manual(values=regcol,name=NULL)+labs(x='Year',y='YLDs / DALYs (%)',tag='a')+th()
p2<-ggplot(s,aes(year,VLW_reduction_pct,colour=subregion))+geom_line(linewidth=.8)+xt+scale_y_continuous(limits=c(0,1.4),breaks=c(0,.4,.8,1.2))+scale_colour_manual(values=regcol,name=NULL)+labs(x='Year',y='Reduction with YLL-based VLW (%)',tag='b')+th()
savefig((p1|p2),'Figure_R1_1_4_subregional_sensitivity',7,3)
cat('Saved four PDF and PNG figure pairs.\n')
