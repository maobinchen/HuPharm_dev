#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(reshape2)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
  theme=  shinytheme('spacelab'),
  #shinythemes::themeSelector(),
  # Application title
  #titlePanel("Tumor Volume Statistical Analysis",windowTitle = 'TV Analysis'),

  # Sidebar with a slider input for number of bins 
  navbarPage("Tumor Volume Statistical Analysis",id='main',
    tabPanel("Start",
          textOutput('urlText'),     
          uiOutput('uiLogin')
          #uiOutput('pass')
    ),
    
    # Show a plot of the generated distribution
    
   tabPanel(title = "Data", value='Data',
            br(),
            #tableOutput('invalid_tv'),
            htmlOutput('del_days'),
            DT::dataTableOutput(outputId = "tv_table")),
   
   tabPanel(title = "Summary", 
            h3("Tumor Volume Summary Statistics Mean +/- SEM(n))"),
            br(),
            DT::dataTableOutput(outputId = "tv_summary_table"),
            checkboxInput(inputId='show_growth_curve',
                          label="Display tumor growth curves",
                          value=FALSE),
            conditionalPanel('input.show_growth_curve',
                             plotOutput(outputId = "growth_curve")),
            checkboxInput(inputId='show_group_curve',
                         label="Display group mean +/- sem plot",
                         value=FALSE),
            conditionalPanel('input.show_group_curve',
                            plotOutput(outputId = "group_stat_plot"))
   ),
                   
                   
   tabPanel(title="Analysis",
            br(),
            #compare groups
            tabsetPanel(tabPanel(title='Overall',
                                 br(),
            wellPanel(                     
                uiOutput('sd1'),#select first day of efficacy study to calculate RTV & DeltaTV
                radioButtons(inputId = "var_cmp", 
                             label = "Perform statistical analysis on:",
                             choices = c('TV','RTV','DeltaTV'), 
                             selected = 'TV',inline=T)),
            checkboxInput(inputId='overall_cmp',
                          label="Check overall difference among groups",
                          value = FALSE),
            conditionalPanel('input.overall_cmp',
                             #radioButtons(inputId='group_cmp_mtd',
                                          #label='Method',
                                          #choices = c('ANOVA','Kruskal-Wallis')),
                             plotOutput(outputId = "group_box_cmp"))),
            tabPanel(title='Pairwise',
                     br(),
            checkboxInput(inputId='pw_group_cmp',
                          label="Perform pairwise comparisons between individual groups",
                          value=FALSE),
            conditionalPanel('input.pw_group_cmp',
                             DT::dataTableOutput(outputId = "pw_group_cmp_tbl"))),
            tabPanel(title='By Day',
                     br(),
                     tabsetPanel(id='byday',
                         tabPanel(title="Parameters",
                             checkboxInput(inputId='stat_ana_day',
                                           label="Perform statistical analysis on specific day",
                                           value=FALSE),
                             conditionalPanel('input.stat_ana_day',
                                              wellPanel(
                                                  radioButtons(inputId = "tgi_def", 
                                                               label = "TGI defined by:",
                                                               choices = c('T/C'='1','Delta(T)/Delta(C)'='2','1-Delta(T)/Delta(C)'='3'), 
                                                               selected = '3',inline=T),
                                                  radioButtons(inputId = "var_cmp_day", 
                                                               label = "Perform statistical analysis on:",
                                                               choices = c('TV','RTV','DeltaTV'), 
                                                               selected = 'TV',inline=T),
                                                  uiOutput('ad1'),
                                                  uiOutput('gs1'),
                                                  uiOutput('rg1'),
                                                  #selectInput(inputId = "pval_adj", 
                                                   #           label = "Select method for multiple comparison correction:",
                                                   #           choices = p.adjust.methods, 
                                                   #           selected = 'BH'),
                                                  #define group comparisons;all pairwise, one against all or customized
                                                  radioButtons(inputId = "posthoc_cmps", 
                                                              label = "Define which groups to compare:",
                                                              choices = c('all_pairwise','treatment_vs_vehicle'), 
                                                              selected = 'all_pairwise',inline=T),
                                                  DT::dataTableOutput(outputId = "tv_day_table")),
                                               wellPanel(
                                                  actionButton(inputId = 'run_ana',
                                                               label = "Run analysis"))
                                               )),
                          tabPanel(title='Plots',
                                   plotOutput(outputId = "day_cmp_plot")),
                          tabPanel(title='Tables',
                                   wellPanel(conditionalPanel('input.run_ana > 0',
                                                              downloadButton(
                                                                  outputId = 'download',
                                                                  label = 'Download HTML Report'
                                                              ),
                                                              downloadButton(
                                                                  outputId = 'download_pdf',
                                                                  label = 'Download PDF Report'
                                                              ),
                                                              downloadButton(
                                                                  outputId = 'downloadData',
                                                                  label = 'Download Figures and Tables'
                                                              ))),
                                   h3("Tumor Growth Inhibition(TGI) table"),
                                   tableOutput(outputId = 'tgi_table'),
                                   br(),
                                   h3("Statistical analysis at specific day(p.value adjusted for multiple comparison)"),
                                   DT::dataTableOutput(outputId = 'day_cmp_pval'))
                       
                    )
                )
        )
   )
       
   )
))
