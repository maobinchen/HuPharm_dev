#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
#library(tidyverse)
library(ggplot2)
library(ggpubr)
library(PMCMRplus)
library(RCurl)
library(DT)
library(rmarkdown)
library(zip)

mytheme <- theme(plot.title=element_text(face="bold.italic",size=16, color="brown",hjust=0.5),
                 axis.title.x=element_text(face="bold.italic",size=14, color="brown"),
                 axis.text=element_text(face="bold", size=10,color="darkblue"),
                 axis.title.y=element_text(face="bold",size=14, color="brown"),
                 strip.text = element_text(face="bold", size=12,color="brown"),
                 panel.background=element_rect(fill="white",color="darkblue"),
                 panel.grid.major.y=element_line(color="grey",linetype=1),
                 panel.grid.minor.y=element_line(color="grey",linetype=2),
                 panel.grid.minor.x=element_blank(),
                 legend.position="top")

grps_cols = c("gray0", "red","dark green","yellow","purple","blue","grey",
             "chartreuse1","light blue","palevioletred1","salmon4","peachpuff","gray48",
              "greenyellow","orange","salmon","purple2","lightgoldenrod1","brown4","darkorchid4")

vecSummary=function(v){
    mv=sprintf('%.2f',mean(v))
    n=length(v)
    sem=sprintf('%.2f',sd(v)/sqrt(n))
    paste0(mv,' +/- ',sem,'(',n,')')
}

#convert a matrix a vector
matrix2vec=function(m){
    vec=as.vector(m)
    names(vec)=as.vector(sapply(colnames(m),paste,rownames(m),sep=" - "))
    vec
}

twoSampleTest=function(v1,v2){
    v1=v1[!is.na(v1)]
    v2=v2[!is.na(v2)]
    if(n_distinct(v1)<2 | n_distinct(v2)<2) return(NA)
    pb=bartlett.test(list(v1,v2))$p.val #bartlett test 
    pval=NA
    if(pb<0.05){
        #pval=round(wilcox.test(v1,v2)$p.val,3)
        pval=sprintf("%.3g",wilcox.test(v1,v2)$p.val)
    }else{
        pval=sprintf("%.3g",t.test(v1,v2)$p.val)
    }
    return(as.numeric(pval))
}
#pairwise comparison of tumor volume between groups for some specific day
cmpDay=function(df,var_sel){
    groups=as.character(levels(df$Group))
    gcs=combn(groups,2)
    pvals=apply(gcs,2,function(x) twoSampleTest(df[df$Group==x[1],var_sel],df[df$Group==x[2],var_sel]))
    names(pvals)=apply(gcs,2,function(x) paste0(x[1],' vs ',x[2]))
    pvals
}

#calculate TGI based on tumor volume data
getTGI=function(tv_long,start_d,ref_g,def=3){
    tvs=subset(tv_long,day>=start_d)
    mean_tv_gd=tapply(tvs$TV,list(tvs$Group,tvs$day),mean,na.rm=T)
    t_c=apply(mean_tv_gd,1,function(x) x/mean_tv_gd[1,])
    #print(mean_tv_gd)
    tv_diff=mean_tv_gd-mean_tv_gd[,1]
    tv_diff=t(tv_diff[,-1])
    dt_dc=tv_diff/tv_diff[,ref_g]
    tgi=1-dt_dc
    if(def==1) tgi=t_c
    if(def==2) tgi=dt_dc
    tgi=tgi[,colnames(tgi) != ref_g,drop=F]
    tgi
}

Logged = FALSE;
PASSWORD <- data.frame(user = "crownbio", pwd = "unagi")

# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {
    query <- reactive({validate(need(session$clientData$url_search!='','Not logged in'))
        parseQueryString(session$clientData$url_search)})
    output$uiLogin=renderUI({
        if (query()[['user']]==PASSWORD$user & query()[['pwd']]==PASSWORD$pwd) {
            column(4,
                   radioButtons('source','Tumor volume data from:',choices = c('Paste csv'='paste','Upload csv'='upload')),
                   conditionalPanel("input.source=='upload'",
                                    fileInput(inputId='tv_file',
                                              label='Select tumor volume file for statistical analysis',
                                              accept = c(
                                                  "text/csv",
                                                  "text/comma-separated-values,text/plain",
                                                  ".csv")
                                    )),
                   conditionalPanel("input.source=='paste'",
                                    #textInput('tv_csv',label='Paste tumor volume csv data',value=''))
                                    tags$textarea(id='tv_csv',rows=20,cols=80,'')),
                   br(),
                   actionButton("submit",'Start Analysis')
            )
        }else{
            "Login failed!"
        }    
    })

  tvf = reactive({
     req(input$tv_file)
     inFile <- input$tv_file
     if (is.null(inFile)) return(NULL)
     inFile$datapath
  })
  
  observeEvent(input$submit,updateNavbarPage(session, "main",
                                             selected = 'Data'))
  
  tv = eventReactive(input$submit,{
      
      if(input$source=='upload'){
          tv=read.csv(tvf(),header=T,stringsAsFactors = F) #read raw tumor volume data
      }else{
          tv=read.csv(text=input$tv_csv,header=T,stringsAsFactors = F) #read raw tumor volume data
      }        
      #write.csv(tv,'tv.csv')
      names(tv)[1:2]=c('Group','Mouse')
      names(tv)=gsub('^X','Day',names(tv))
      mouseTV=apply(tv[,c(-1,-2)],1,function(x) any(x>0)) #tv volume of a mouse >0 for any given day
      tv=tv[mouseTV,]
      tv$Group = gsub("^.*?(\\d+).*","G\\1",tv$Group)
      tv$Group = as.factor(tv$Group)
      tv$Mouse = as.factor(tv$Mouse)
      tv
  })
  
  #data sanity check, report non-valid TV (should be numeric >=0), only keep days with more at least two groups with at least two meaningful TV
  datap = reactive({
      #to add sanity test;last day should have at least two groups with more than two meaningful observations
      #req(input$tv_file)
      qc=list()
      tv1=melt(tv(),id.vars=c(1,2),value.name="TV")
      tv1=subset(tv1,!is.na(TV))
      tv1$day=as.numeric(substr(tv1$variable,4,10))
      qc$invalid = subset(tv1, !(is.numeric(TV) & TV>=0))
      tv1 = subset(tv1, is.numeric(TV) & TV>=0)
      tvn = with(tv1,tapply(tv1$TV,list(day,Group),n_distinct))
      #print(tvn)
      dmo = apply(tvn,1,function(x) sum(x>=2,na.rm=T))
      valid_days=as.numeric(names(dmo[dmo>=2]))
      del_days=names(dmo[dmo<2])
      del_days_str='Statistical analysis requires at least two groups with two or more distinct observations.</br>'
      if(length(del_days)>0) del_days_str=paste0(del_days_str,'Removing Day <b>',paste(del_days,collapse = ','),'</b> from further analysis')
      qc$msg=del_days_str
      tv1 = subset(tv1,day %in% valid_days)
      #print(del_days_str)
      list('data'=tv1,'qc'=qc)
  })
  
  tv1 = reactive(datap()$data)

  output$invalid_tv = renderTable(datap()$qc$invalid[,c('Group','Mouse','day','TV')])
  
  output$del_days = renderUI(HTML(datap()$qc$msg))
  
  days=reactive(sort(unique(tv1()$day)))
  #days=reactive(1:10)
  
  #qc = reactiveValues('msg'='','invalid'=NULL)
  output$sd1 = renderUI(selectInput(inputId = "start_day", 
                                    label = "Select first day of efficacy study:",
                                    choices = days(), 
                                    selected = days()[1]))
  
  #assume TV0 > 0
  tv_long = reactive({
      tv_long = tv1()
      tv0=tv_long[tv_long$day==input$start_day,c('Mouse','TV')]
      names(tv0)[2]='TV0'
      tv_long=merge(tv_long,tv0,by="Mouse")
      tv_long$RTV=ifelse(tv_long$TV0==0,NA,tv_long$TV/tv_long$TV0)
      tv_long$DeltaTV=tv_long$TV-tv_long$TV0
      #print(tv_long)
      tv_long[,-3]
  })
  
  groups=reactive(levels(tv_long()$Group))
  
  pw=600 #plot width
  asp_ratio=0.618 #aspect ratio
  n_col=2 #number of subgraphs per column
  
  ph=function(){
      pw*asp_ratio*ceiling(length(days())/n_col)
  }

  #A summary of tumor volume statistics,in (mean +/- sd(n) format)
  tv_summary=reactive(with(tv1(),tapply(TV,list(Group,day),vecSummary)))
  
  #compare 
  pw_grp_pvals=reactive({
      days=split(tv_long(),tv_long()$day)
      pval_table=do.call(rbind,lapply(days,cmpDay,input$var_cmp)) 
      pval_table
  })

  output$tv_table = DT::renderDataTable({
      DT::datatable(data = tv(),caption = 'Tumor Volume Data Table',extensions = c('Buttons','FixedColumns'),
                    options = list( dom = 'Blfrtip',autoWidth=T,pageLength = 50,digits=3,
                                    buttons=c('copy','csv','excel','pdf','print'),fixedColumns=list(leftColumns=2)), 
                    rownames = FALSE) %>% formatStyle(3:ncol(tv()),color=styleInterval(c(0,20,3000),c('red','green','black','red'))) %>%
                           formatRound(3:ncol(tv()), digits=2)
  })
  
 output$tv_summary_table = DT::renderDataTable(DT::datatable(tv_summary(),extensions = c('Buttons'),
                                               options = list( dom = 'Blfrtip',autoWidth=T,pageLength = 20,buttons=c('copy','csv','excel','pdf','print'))) %>%
                                               formatStyle(0,fontWeight='bold'))
 
 growth_cruve_plot = reactive({
     p=ggplot(tv1())+aes(day,TV,group=Mouse,color=Group)+scale_color_manual(values=grps_cols)+
         labs(title='Tumor Growth Curve',x='Day',y=expression(paste('Tumor Volume ( ',mm^3,')')))+
         geom_point()+geom_line()+facet_wrap(~Group,ncol=3)+mytheme+theme(legend.position = 'none')
     p
 })
 
 output$growth_curve = renderPlot(growth_cruve_plot())
 
 group_stat_plot=reactive({
     p=ggline(tv1(),x='day',y='TV',add=c('mean_se'),group = 'Group',title='Tumor Volume Summary Plot',
              xlab='Day',ylab=expression(paste('Tumor Volume ( ',mm^3,')')),size=1.5,color='Group',
              add.params =   list(size=1.5,alpha=0.8,position='dodge'))+scale_color_manual(values=grps_cols)+mytheme
     p
 })     
 
 output$group_stat_plot = renderPlot(group_stat_plot())
 
 #grp_cmp_method=reactive(ifelse(input$group_cmp_mtd=='ANOVA','anova','kruskal.test'))
 group_box_cmp=reactive({
     p=ggboxplot(tv_long(),x='Group',y=input$var_cmp,title=paste0(input$var_cmp,'~Group|Day Boxplot'),add=c('jitter'),color='Group',
                 xlab='Day', facet.by = 'day')+scale_color_manual(values=grps_cols)+
         stat_compare_means(method='kruskal.test',aes(label = paste0("p=", ..p.format..)),label.x.npc = 'left',label.y.npc = 'top')+mytheme
     facet(p,facet.by = 'day',scales='free',ncol=n_col)
     p
 })      
 output$group_box_cmp = renderPlot(group_box_cmp(),width=pw*n_col,height=ph)
 
 
  
 output$pw_group_cmp_tbl = DT::renderDataTable({
     DT::datatable(data = t(pw_grp_pvals()), extensions = c('Buttons'),
                   options = list( dom = 'Blfrtip',autoWidth=T,pageLength = 20, buttons=c('copy','csv','excel','pdf','print'))) %>%
                   formatStyle(1:nrow(pw_grp_pvals()),color=styleInterval(c(0.01,0.05),c('darkred','red','black'))) %>%
                   formatStyle(0,fontWeight='bold')
 })
     

  
  output$ad1 = renderUI(selectInput(inputId = "analysis_day", 
                                    label = "Pick which day to compare:",
                                    choices = days(), 
                                    selected = rev(days())[1]))

  #group selection,keep groups with at least two mice
  output$gs1 = renderUI(checkboxGroupInput(inputId = "group_cmps", 
                                    label = "Select which groups should be included:",
                                    choices = groups(),
                                    selected = groups(),
                                    inline = T)) 
  
  #reference group
  output$rg1 = renderUI(selectInput(inputId = "ref_group", 
                                    label = "Select vehicle group:",
                                    choices = input$group_cmps, 
                                    selected = input$group_cmps[1]))
  
  #change group selection depending on day (at least two distinct observations)
 observeEvent(c(input$analysis_day,input$var_cmp),{
      tds = subset(tv_long(),day == input$analysis_day)
      tds = tds[!is.na(tds[,input$var_cmp]),]
      nd_g = table(tds$Group) 
      grps=names(nd_g[nd_g>=2])
      updateCheckboxGroupInput(session,'group_cmps',choices = grps,selected=grps,inline=T)
      updateSelectInput(session,'ref_group',choices=grps,selected=grps[1])
  })
 
 observeEvent(input$var_cmp,{
     updateRadioButtons(session,'var_cmp_day',selected = input$var_cmp)
 })
 
 observeEvent(input$run_ana,{
     updateRadioButtons(session,'var_cmp',selected = input$var_cmp_day)
 })
 
 
  
  #TO DO:remove groups with less than two distinct tv
  tv_day = reactive({
      req(input$analysis_day)
      req(input$var_cmp)
      req(input$group_cmps)
      req(input$ref_group)
      tv_day = subset(tv_long(),day == as.numeric(input$analysis_day) & Group %in% input$group_cmps)
      #tv_day = tv_day[!is.na(tv_day[,input$var_cmp]),]
      tv_day$Group = factor(tv_day$Group)
      tv_day$Group=relevel(tv_day$Group,ref=input$ref_group)
      tv_day
  })
  
  tgi_table = reactive({
      req(input$ref_group)
      req(input$start_day)
      tgi_table=getTGI(tv_long(),as.numeric(input$start_day),input$ref_group,as.numeric(input$tgi_def))
      tgi_table
  })
  
  output$tgi_table = renderTable(tgi_table(),rownames = T,striped = T)

  #output for the by day analysis
  outday=eventReactive(input$run_ana, {
      #cat(input$analysis_day)
      grps=as.character(sort(unique(tv_day()$Group)))
      #run bartlett test to test homogeniety of variance and normality
      pb_day=bartlett.test(tv_day()[,input$var_cmp],tv_day()$Group)$p.val
      all_cmp=as.list(as.data.frame(combn(grps,2)))
      all_cmp=lapply(all_cmp,as.character)
      sel_cmp=rep(T,length(all_cmp))
      if(input$posthoc_cmps == 'treatment_vs_vehicle') sel_cmp = sapply(all_cmp,function(x) input$ref_group %in% x)
      test=ifelse(pb_day<0.05,'wilcox.test','t.test')
      p = ggbarplot(tv_day(), x = "Group", y = input$var_cmp, color='Group',fill=NULL,size=2,
                    add=c('mean_se','jitter'),add.params = list(size=1.5),
                    title=paste0('Day ',input$analysis_day,':', test,' nominal p value'))
      #p = ggboxplot(tv_day, x = "Group", y = "TV", add='jitter')
      p = p +scale_color_manual(values=grps_cols)
      p = p+stat_compare_means(method=test,comparisons = all_cmp[sel_cmp],label='p.signif')+mytheme+theme(legend.position = 'none')

      grp_n=length(grps)
      rps=NULL #p-values of two group comprison
      omnibus_p=NULL #p-value of omni-bus test (ANOVA/Kruskal-wallis test)
      #cmps_names=sapply(all_cmp,function(x) paste0(x[1],"-",x[2]))
      #cmps=c()#user defined cmps, formatted as names of cmps
      posthoc_test=''
      if(grp_n==1){
          print('Only one group,No statistical analysis will be done')      
      }else if(grp_n==2){
          rps=cmpDay(tv_day(),input$var_cmp)
      }else{
          if(pb_day<0.05){
              omnibus_p = kruskal.test(tv_day()[,input$var_cmp]~tv_day()$Group)$p.val
              names(omnibus_p)='Kruskal-Wallis'
              posthoc_test = ifelse(input$posthoc_cmps == 'treatment_vs_vehicle','kwManyOneConoverTest','kwAllPairsConoverTest')
          }else{
              aov_m=aov(tv_day()[,input$var_cmp]~tv_day()$Group)
              omnibus_p = summary(aov_m)[[1]]$Pr[1]
              names(omnibus_p)='ANOVA'
              posthoc_test = ifelse(input$posthoc_cmps == 'treatment_vs_vehicle','dunnettTest','tukeyTest')
          }
          pvals=do.call(posthoc_test,list(tv_day()[,input$var_cmp],tv_day()$Group))$p.val
          rps=matrix2vec(pvals)
          rps=rps[!is.na(rps)]
      }
      all_ps=c(pb_day,omnibus_p,rps)
      test_names=c("Bartlett's test",names(omnibus_p),names(rps))
      print(all_ps)
      #significance category, <0.05 *, <0.01 **, <0.001 ***
      sigs=rep('ns',length(all_ps))
      sigs[all_ps<0.05]='*'
      sigs[all_ps<0.01]='**'
      sigs[all_ps<0.001]='***'
      output_df=data.frame('Test'=test_names,'p values'=sprintf("%.3g",all_ps),'significance level'=sigs,stringsAsFactors = F)
      rownames(output_df)=NULL
      output_df$p.values=as.numeric(output_df$p.values)
      list('table'=output_df,'figure'=p)
  })
  
  observeEvent(input$run_ana, {
      updateTabsetPanel(session,'byday',selected='Plots')
      output$tv_day_table = DT::renderDataTable({
          DT::datatable(data = tv_day(), 
                        options = list(pageLength = 10,digits=3), 
                        rownames = FALSE) %>% formatRound(c('TV','TV0','RTV'), digits=2)
      })
      output$day_cmp_plot=renderPlot(print(outday()$figure))
      
      output$day_cmp_pval=DT::renderDataTable({DT::datatable(data = outday()$table, extensions = c('Buttons'),
                                                             options = list( dom = 'Blfrtip',autoWidth=T,pageLength = 20, buttons=c('copy','csv','excel','pdf','print'))) %>%
              formatStyle('p.values',target='row',color=styleInterval(c(0.01,0.05),c('darkred','red','black'))) %>%
              formatStyle(0,fontWeight='bold')
      })
  })
  
  output$download = downloadHandler(
      filename = 'stat_analysis.html',
      content = function(file){
          tempd=tempdir()
          tempReport <- file.path(tempd,'tv_analysis.Rmd')
          file.copy(c('tv_analysis.Rmd',"crownbio-logo.png"), tempd, overwrite = TRUE)
          params = list(source=input$source,p_adj_method = input$pval_adj,start_day = input$start_day,analysis_day=input$analysis_day,tgi_def=input$tgi_def,
                        sel_group=input$group_cmps,ref_group=input$ref_group, posthoc_cmps=input$posthoc_cmps, analysis_var=input$var_cmp)
          if(input$source=='upload'){
              params$tv_file = tvf()
          }else{      
              params$tv_str=input$tv_csv
          }      
          withProgress(message='Please wait',detail='Generating Report..',value=0,
                       {render(tempReport,output_format = 'html_document', output_file=file, params = params,envir = new.env(parent = globalenv()))})
                       
      }
  )
  
  output$download_pdf = downloadHandler(
      filename = 'stat_analysis.pdf',
      content = function(file){
          tempd=tempdir()
          tempReport <- file.path(tempd,'tv_analysis_pdf.Rmd')
          file.copy(c('tv_analysis_pdf.Rmd',"header.tex","crownbio-logo.png"), tempd, overwrite = TRUE)
          params = list(source=input$source,p_adj_method = input$pval_adj,start_day = input$start_day,analysis_day=input$analysis_day,tgi_def=input$tgi_def,
                        sel_group=input$group_cmps,ref_group=input$ref_group, posthoc_cmps=input$posthoc_cmps, analysis_var=input$var_cmp)
          if(input$source=='upload'){
              params$tv_file = tvf()
          }else{      
              params$tv_str=input$tv_csv
          }      
          withProgress(message='Please wait',detail='Generating Report..',value=0,
                       {render(tempReport,output_format = 'pdf_document', output_file=file, 
                               params = params,envir = new.env(parent = globalenv()))})
          
      }
  )
  
 
  #download zip file of high resolution images and tables 
  output$downloadData <- downloadHandler(
      filename = 'stat_analysis.zip',
      content = function(fname) {
          tables <- list('TV_summary'=tv_summary(),'all_pairwise_pvalue'=pw_grp_pvals(),
                         'TGI'=tgi_table(),'Day'=outday()$table)
          imgs=list('Growth_curve'=growth_cruve_plot(),'Group_mean_sem'=group_stat_plot(),
                    'Daily_pairwise_cmparison'=outday()$figure)
                    #'Overall_group_comparison'=group_box_cmp(),
          names(tables)[4]=paste0('Day',input$analysis_day)
          tmpdir <- tempdir()
          setwd(tempdir())
          fs=c()
          for (kw in names(tables)) {
              path <- paste0(kw, ".csv")
              fs <- c(fs, path)
              write.csv(tables[[kw]], path)
          }
          for (kw in names(imgs)){
              path <- paste0(kw, ".png")
              fs <- c(fs, path)
              ggsave(path,imgs[[kw]],width=8,height=8,dpi=150)
          }
          zip(fname,fs)
      },
      contentType = "application/zip"
  )
  
})
