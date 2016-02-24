
# This is the user-interface definition of a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

library(shiny)

shinyUI(fluidPage(theme = 'new_theme.css',

  # Application title
  titlePanel("SimFish"),
  img(src = "snap_logo.png", height = 72, width = 72),
  h5('This is a fishery simulater based on  the General MSE code by Cody Szuwalski'),
  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
#       h6('Select number of years for simulation'),
      h6('Select species to use'),
      selectInput(inputId = 'species',label = 'Species',
                  choices = c('Chita', "Lorna Drum",'Sea Cucumber'), selected = 'Lorna Drum'),
      sliderInput("Years",
                  "Number of Years of Simulation:",
                  min = 2,
                  max = 50,
                  value = 30),
      br(),
      sliderInput('sigmaR', 'Recruitment Variability',
                  min = 0, max = 2, value = 0.05, step = 0.05),
      br(),
      sliderInput('fish_select','Fishery Selectivity (% of Max Length)', min = 0, max = 100, value = c(10,10), step = 0.5, pre = '%'),
      br(),
      sliderInput('f_select','Fishing Mortality Rate', min = 0, max = 3, value = c(0.3), step = 0.05),
      br(),
      sliderInput('surv_select','Survey Selectivity (% of Max Length)', min = 0,
                  max = 100, value = c(10,10), step = 0.5, pre = '%')),
    mainPanel(
      tabsetPanel(type = "tabs",
                  tabPanel("Life History", plotOutput('life_plot'), plotOutput('recruitment_plot')),
                  tabPanel("Froese Indicators",plotOutput("froese_trend")),
                  tabPanel("CPUE Trends",plotOutput("cpue_trend")),
                  tabPanel("Catch Trends",plotOutput("catch_trend")),
                  tabPanel("Biomass Trends",plotOutput("biomass_trend")),
                  tabPanel("Length Frequency",plotOutput("length_freq")),
                  tabPanel("Raw Data", tableOutput('bio_table'))
      ) #close tabsetPanel
    ) # close mainPanel
  ) #close sidebarPanel
) #close fluidpage
) #close ShinyUI

