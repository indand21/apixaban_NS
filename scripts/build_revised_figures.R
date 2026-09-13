source('scripts/load_project.R')
stopifnot(file.exists('output/revised/completed_at.txt') || Sys.getenv('NS_FIGURE_PREVIEW')=='1')
figdir <- file.path('output', 'revised', 'figures')
dir.create(figdir,recursive=TRUE,showWarnings=FALSE)
read_result<-function(n)read.csv(file.path('output/revised',paste0(n,'.csv')))
stages<-c('Normal','NS_Mild','NS_Moderate','NS_Severe')
labels<-c('Reference','NS mild','NS moderate','NS severe')
stage_factor<-function(x)factor(x,levels=stages,labels=labels)
pal<-setNames(c('#2964A1','#B38A22','#B64B80','#454545'),labels)
binary<-c('Feedforward'='#2964A1','Feedback'='#CF762F')
theme_set(theme_bw(base_size=11,base_family='Arial')+
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_line(colour='#E9E9E9',linewidth=.25),
    plot.title=element_text(face='bold',size=13),plot.subtitle=element_text(size=10),
    strip.background=element_rect(fill='#F3F3F3'),legend.position='top',legend.title=element_blank(),
    plot.caption=element_text(hjust=0,size=9),plot.margin=margin(8,10,8,8)))
save_fig<-function(p,id,w=7.5,h=5.2){
  ggsave(file.path(figdir,paste0(id,'.png')),p,width=w,height=h,dpi=300,device=ragg::agg_png)
  ggsave(file.path(figdir,paste0(id,'.pdf')),p,width=w,height=h,device=cairo_pdf)
}
pk<-read_result('pk_profiles');pk$stage<-stage_factor(pk$stage)
long<-rbind(data.frame(pk,metric='Total concentration',value=pk$total_ng_mL),data.frame(pk,metric='Unbound concentration',value=pk$free_ng_mL))
p<-ggplot(long,aes(time_h,value,colour=stage,linetype=stage))+geom_line(linewidth=.75)+
  facet_wrap(~metric,scales='free_y')+scale_colour_manual(values=pal)+scale_linetype_manual(values=c('solid','dashed','dotdash','longdash'))+
  scale_x_continuous(breaks=seq(0,12,3))+labs(title='Steady state apixaban profiles',
  subtitle='5 mg every 12 hours; reference age 40 years, weight 70 kg, eGFR 100',x='Time after dose (hours)',y='Concentration (ng/mL)',
  caption='Illustrative NS scenarios. Panels have different vertical scales; curves are simulations.')
save_fig(p,'Figure_1_PK',h=4)
q<-subset(read_result('qsp_profiles'),mode=='feedforward');q$stage<-stage_factor(q$stage)
q$dose_label<-factor(q$dose,levels=c(0,5),labels=c('No drug','5 mg BID trough'))
p<-ggplot(q,aes(time_s,IIa,colour=dose_label,linetype=dose_label))+geom_line(linewidth=.75)+
  facet_wrap(~stage,ncol=2)+scale_colour_manual(values=c('#2964A1','#CF762F'))+
  scale_linetype_manual(values=c('solid','dashed'))+coord_cartesian(xlim=c(0,800))+
  labs(title='Simulated thrombin generation',subtitle='Feedforward model; TF 5 pM; implicit protein C activation disabled',
    x='Assay time (seconds)',y='Thrombin (nM)',caption='Integration extends to 1200 seconds. Curves use matched no-drug baselines.')
save_fig(p,'Figure_2_Thrombin')
d<-read_result('score_curves');d$stage<-stage_factor(d$stage)
d$mode<-factor(d$mode,levels=c('feedforward','coupled'),labels=c('Feedforward','Feedback'))
p<-ggplot(subset(d,score>0),aes(dose,score,colour=mode,linetype=mode))+geom_line(linewidth=.75)+
  facet_wrap(~stage,ncol=2)+scale_colour_manual(values=binary)+scale_linetype_manual(values=c('solid','dashed'))+
  scale_y_log10()+labs(title='Exploratory composite score',subtitle='Peak-thrombin suppression multiplied by matched platelet aggregate retention',
    x='Simulated dose per administration (mg every 12 hours)',y='Composite score (log scale)',
    caption='Zero scores are omitted on the logarithmic axis. This is not a clinical therapeutic index.')
save_fig(p,'Figure_3_Score')
s<-subset(read_result('structural_sensitivity'),TF_pM==5 & exposure=='trough')
s$stage<-stage_factor(s$stage);s$activation<-ifelse(s$k39_nM_s==0,'Disabled','Legacy activation')
p<-ggplot(s,aes(factor(duration_s),stage,fill=100*ETP_suppression))+geom_tile(colour='white',linewidth=.8)+
  geom_text(aes(label=sprintf('%.2f',100*ETP_suppression)),size=3.5,colour='black')+
  facet_wrap(~activation)+scale_fill_gradient(low='#FFF7FA',high='#C56992',name='Suppression (%)')+
  labs(title='Dependence on assay duration and protein C activation',subtitle='ETP suppression at 5 mg BID trough, fixed TF 5 pM',
    x='Integration duration (seconds)',y=NULL,caption='Legacy activation assumes thrombomodulin action without an explicit thrombomodulin state.')+
  theme(legend.position='bottom',legend.title=element_text(size=10),panel.grid=element_blank())
save_fig(p,'Figure_4_Structure',h=4)
if(file.exists('output/revised/synthetic_population.csv')) {
v<-read_result('synthetic_population');v$stage<-stage_factor(v$stage)
v$mode<-factor(v$mode,levels=c('feedforward','coupled'),labels=c('Feedforward','Feedback'))
p<-ggplot(v,aes(stage,score,fill=mode))+geom_boxplot(outlier.size=.7,linewidth=.4)+scale_fill_manual(values=binary)+
  labs(title='Composite score in synthetic individuals',subtitle=paste(length(unique(v$id)),'matched draws per scenario and mode; maximum over the 0 to 20 mg search interval'),
    x=NULL,y='Maximum composite score',caption='Design-distribution spread, not clinical population uncertainty or a confidence interval.')
save_fig(p,'Figure_S1_Population',h=4)
}
cat('FIGURES COMPLETE\n')
