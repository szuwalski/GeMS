
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#



# a = simfish()
set.seed(123)
library(shiny)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggthemes)
library(roxygen2)
sapply(list.files(pattern="[.]R$", path="simfish/R/", full.names=TRUE), source)
# devtools::load_all('simfish')

shinyServer(function(input, output) {

  simmedpop <- reactive({

    if (input$species == 'Lorna Drum')
    {

      MaxAge <- 10

      vbk <- 0.195

      linf <- 550

      t0 <- -0.808

      mat50 <- 2

      mat95 <- 3

      alpha <- 0.00013

      beta <- 2.63

      Rzero <- 3e6

      MaxMove <- 2

      Move50 <- 3

      Move95 <- 4

      steepness <- 0.8

      sigmaR <- 0.1

      NatM <- 0.52


    }

    out <- simfish(SimYear = input$Years, VonKn = vbk, VonKs = vbk,
                   LinfN = linf, LinfS = linf, t0n = t0, t0s = t0,
                   mat50n = mat50, mat50s = mat50, mat95n = mat95,
                   mat95s = mat95, alphaN = alpha, betaN = beta,
                   alphaS = alpha, betaS = beta, MaxMovingN = MaxMove,
                   MaxMovingS = MaxMove, Move50n = Move50,
                   Move50s = Move50, Move95n = Move95, Move95s = Move95,
                   steepnessN = steepness, steepnessS = steepness,
                   sigmaRn = input$sigmaR, sigmaRs = input$sigmaR, RzeroN = Rzero,
                   RzeroS = Rzero, MaxAge = MaxAge, NatMn = NatM, NatMs = NatM,
                   sel50n = input$fish_select[1]/100*linf, sel50s = input$fish_select[1]/100*linf,
                   sel95n = input$fish_select[2]/100*linf,  sel95s = input$fish_select[2]/100*linf,
                   surv50n = input$surv_select[1]/100*linf, surv50s = input$surv_select[1]/100*linf,
                   surv95n = input$surv_select[2]/100*linf,  surv95s = input$surv_select[2]/100*linf,
                   HistoricalF = input$f_select)


#     out <- simfish(SimYear = input$Years)

    return(out)
  })


  output$biomass_trend <- renderPlot({
    #     out <- simfish(SimYear = input$Years)
    simmedpop()$plots$biomass_plot
    #     out$plots$biomass_plot
  })

  output$length_freq <- renderPlot({
    simmedpop()$plots$length_plot
  })

  output$life_plot <- renderPlot({
    simmedpop()$plots$life_plot
  })

  output$recruitment_plot <- renderPlot({
    simmedpop()$plots$recruitment_plot
  })

  output$froese_trend <- renderPlot({
    simmedpop()$plots$cope_punt_plot
  })

  output$cpue_trend <- renderPlot({
    simmedpop()$plots$cpue_plot
  })

  output$catch_trend <- renderPlot({
    simmedpop()$plots$catch_plot
  })

  output$bio_table <- renderTable({
    simmedpop()$biomass_data
  })




})
