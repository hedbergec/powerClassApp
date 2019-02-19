#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(plotly)
#library(shinycssloaders)

# Define UI for application that draws a histogram
ui <- fluidPage(
   # Application title
   titlePanel("Power Estimation Tool"),
   h1("Warning--still in testing--do not use for actual research!"),
   sidebarLayout(
      sidebarPanel(
        selectInput("type",
                    "Sample design:",
                    c("Simple random sample randomized trial" = "srs",
                      "Cluster randomized trial" = "crt",
                      "Multisite randomized trial" = "msrt")),
        selectInput("out",
                    "What to plot:",
                    c("Exact power" = "ep",
                      "Power by sample size" = "yPower",
                      "Minimum detectable effect size (MDES) by sample size" = "yMDES"
                    )),
        radioButtons("cov",
                     "Use covariates?",
                     c("No", 
                       "Yes"),
                     "Yes"),
        numericInput("alpha",
                    "pick two tailed alpha",
                    value = 0.05,
                    min = .01,
                    max = 0.1,
                    step = .01)
      ),
      mainPanel(
        tabsetPanel(
          tabPanel("Plot", 
                   uiOutput("showThePlot")
                   
                   ),
          tabPanel("Table", 
                   downloadLink("downloadTable", "Download Output as CSV"),
                   uiOutput("showTheTable"),
                   uiOutput("showTheLegend")
                   )
        ),
        h4("Parameters"),
        fluidRow(
          column(5,
                 uiOutput("effectSize"),
                 uiOutput("powerSlider"),
                 uiOutput("iccSlider"),
                 uiOutput("covR"),
                 uiOutput("unitSize"),
                 uiOutput("clusterSize")
          ),
          column(5, offset = 1,
                 uiOutput("upsilonSlider"),
                 uiOutput("r2unit"),
                 uiOutput("r2cluster"),
                 uiOutput("r2treat")
          )
        )
      )
   )
)


server <- function(input, output) {
  output$effectSize <- renderUI({
    if (input$out != "yMDES") {
      sliderInput("es",
                  "Effect size",
                  min = 0,
                  max = 1,
                  value = .3,
                  step = .1)
    }
  })
  
  output$powerSlider <- renderUI ({
    if (input$out == "yMDES") {
      sliderInput("power",
                  "Power",
                  min = 0,
                  max = .9,
                  step = .1,
                  value = .8)
    }
  })
  
  output$unitSize <- renderUI({
    if (input$type == "srs" & input$out == "ep") {
      lab <- "Number of units per group"
      sliderInput("n",
                  lab,
                  min = 4,
                  max = 50,
                  value = 30)
    }
    else if (input$type == "crt") {
      lab <- "Number of units per cluster"
      sliderInput("n",
                  lab,
                  min = 4,
                  max = 50,
                  value = 30)
    }
    else if (input$type == "msrt") {
      lab <- "Number of units per group per cluster"
      sliderInput("n",
                  lab,
                  min = 4,
                  max = 50,
                  value = 30)
    }
   })
  
  output$iccSlider <- renderUI ({
    if (input$type != "srs") {
      sliderInput("icc",
                  "Intraclass correlation",
                  min = .00,
                  max = .5,
                  step = .05,
                  value = .2)
    }
  })
  
  output$upsilonSlider <- renderUI ({
    if (input$type == "msrt") {
      sliderInput("upsilon",
                  "Ratio of treatment variance to cluster mean variance",
                  min = .00,
                  max = 3,
                  step = .25,
                  value = .25)
    }
  })
  
  output$clusterSize <- renderUI({
    if (input$type != "srs") {
      if (input$out == "ep") {
        if (input$type == "crt") {
          lab <- "Number of clusters per group"
        }
        if (input$type == "msrt") {
          lab <- "Number of clusters"
        }
        sliderInput("m",
                    lab,
                    min = 2,
                    max = 25,
                    value = 10)
      }
    }
  })

  output$covR <- renderUI({
    if (input$type == "srs" & input$cov == "Yes") {
      sliderInput("ryx",
                  "Population correlation between outcome and covariate",
                  min = -1,
                  max = 1,
                  step = .025,
                  value = 0)
    }
  })
  
   output$r2unit <- renderUI({
     if (input$type != "srs" & input$cov == "Yes") {
       sliderInput("r2u",
                   "Population squared-correlation between outcome and covariate at unit level",
                   min = 0,
                   max = 1,
                   step = .025,
                   value = .25)
     }
   })
   
   output$r2cluster <- renderUI({
     if (input$type == "crt" & input$cov == "Yes") {
       sliderInput("r2c",
                   "Population squared-correlation between outcome and covariate at cluster level",
                   min = 0,
                   max = 1,
                   step = .025,
                   value = .25)
     }
   })
   
   output$r2treat <- renderUI({
     if (input$type == "msrt" & input$cov == "Yes") {
       sliderInput("r2t",
                   "Population squared-correlation between variance of treatment effect and covariate at cluster level",
                   min = 0,
                   max = 1,
                   step = .025,
                   value = .25)
     }
   })
   
   computeNCP <- reactive({
     if (input$type == "srs") {
       ncp <- input$es*sqrt(input$n/2)
       if (input$cov == "Yes") {
         ncp <- input$es*sqrt(input$n/2)*sqrt(1/(1-input$ryx^2))*sqrt((2*input$n-3)/(2*input$n-2))
       }
     }
     else if (input$type == "crt") {
       ncp <- input$es*sqrt(input$n*input$m/2)*sqrt(1/(1+(input$n-1)*input$icc))
       if (input$cov == "Yes") {
         ncp <- input$es*sqrt(input$n*input$m/2)*sqrt(1/(1+(input$n-1)*input$icc-(input$r2u+(input$n*input$r2c-input$r2u)*input$icc)))
       }
     }
     else if (input$type == "msrt") {
       ncp <- input$es*sqrt(input$n*input$m/2)*sqrt(1/(1+(input$upsilon*input$n/2-1)*input$icc))
       if (input$cov == "Yes") {
         ncp <- input$es*sqrt(input$n*input$m/2)*sqrt(1/(1+(input$upsilon*input$n/2-1)*input$icc-(input$r2u+(input$upsilon*input$n/2*input$r2t-input$r2u)*input$icc)))
       }
     }
     return(ncp)
   })
   
   computeDF <- reactive({
     if (input$type == "srs") {
       df <- input$n*2-2
     }
     else if (input$type == "crt") {
       df <- input$m*2-2
     }
     else if (input$type == "msrt") {
       df <- input$m-1
     }
     if (input$cov == "Yes") {
       df <- df-1
     }
     return(df)
   })
   
   computeCurves <- reactive({
     if (input$out == "ep") {
       req(computeDF())
       req(computeNCP())
       retList <- list()
       critR <- qt(input$alpha/2, computeDF(), lower.tail = FALSE)
       critL <- qt(input$alpha/2, computeDF(), lower.tail = TRUE)
       
       retList$t <- seq(-5, 5, length = 100)
       retList$cden <- dt(seq(-5, 5, length = 100), computeDF())
       retList$alphaR.x <- c(critR, seq(critR, 5, 0.01), 5)
       retList$alphaR.y <- c(0, dt(seq(critR, 5, 0.01), computeDF()), 0)
       
       retList$alphaL.x <- c(critL, seq(-5, critL, 0.01), critL)
       retList$alphaL.y <- c(0, dt(seq(-5, critL, 0.01), computeDF()), 0)
       
       retList$beta.x <- c(critL, seq(critL, critR, 0.01), critR)
       retList$beta.y <- c(0, dt(seq(critL, critR, 0.01), computeDF(), computeNCP()), 0)
       
       retList$power.x <- c(critR, seq(critR, 5, 0.01), 5)
       retList$power.y <- c(0, dt(seq(critR, 5, 0.01), computeDF(), computeNCP()), 0)
       
       retList$ncp.x <- computeNCP()
       retList$ncp.y <- dt(computeNCP(), computeDF(), computeNCP())
       
       retList$beta <- pt(qt(input$alpha/2,computeDF(), lower.tail = FALSE),computeDF(),computeNCP()) - 
         pt(-qt(input$alpha/2,computeDF(), lower.tail = FALSE),computeDF(),computeNCP())
       retList$ct <- qt(input$alpha/2,computeDF(), lower.tail = FALSE)
       retList$ncp <- computeNCP()
       retList$power <- 1-(pt(qt(input$alpha/2,computeDF(), lower.tail = FALSE),computeDF(),computeNCP()) - 
                             pt(-qt(input$alpha/2,computeDF(), lower.tail = FALSE),computeDF(),computeNCP()))
       return(retList)
     }
     else {
       return(list())
     }
   })
   
   computePower <- reactive({
     req(input$type)
     req(input$cov)
     req(input$es)
     if (input$type == "srs") {
       whatisX <- "Units per group"
       x <- seq(2,50, by=2)
       N <- 2*x
       df <- 2*x-2
       ncp <- input$es*sqrt(x/2)
       if (input$cov == "Yes") {
         df <- df -1
         ncp <- input$es*sqrt(x/2)*sqrt(1/(1-input$ryx^2))*sqrt((2*x-3)/(2*x-2))
       } 
     }
     else if (input$type == "crt") {
       whatisX <- "Clusters per group"
       x <- seq(2,25)
       N <- 2*input$n*x
       df <- 2*x-2
       ncp <- input$es*sqrt(input$n*x/2)*sqrt(1/(1+(input$n-1)*input$icc))
       if (input$cov == "Yes") {
         df <- df - 1
         ncp <- input$es*sqrt(input$n*x/2)*sqrt(1/(1+(input$n-1)*input$icc-(input$r2u+(input$n*input$r2c-input$r2u)*input$icc)))
       }
     }
     else if (input$type == "msrt") {
       whatisX <- "Total clusters"
       x <- seq(2,25)
       N <- 2*input$n*x
       df <- x-1
       ncp <- input$es*sqrt(input$n*x/2)*sqrt(1/(1+(input$upsilon*input$n/2-1)*input$icc))
       if (input$cov == "Yes") {
         df <- df - 1
         ncp <- input$es*sqrt(input$n*x/2)*sqrt(1/(1+(input$upsilon*input$n/2-1)*input$icc-(input$r2u+(input$upsilon*input$n/2*input$r2t-input$r2u)*input$icc)))
       }
     }
     power <- 1-(pt(qt(input$alpha/2,df, lower.tail = FALSE),df,ncp) - 
                   pt(-qt(input$alpha/2,df, lower.tail = FALSE),df,ncp))
     return(list(power = power, x = x, df = df, N = N, ncp = ncp, alpha = input$alpha, type = input$type, whatisX = whatisX))
   })
  
   computeMDES <- reactive({
     req(input$type)
     req(input$cov)
     req(input$es)
     req(input$power)
     if (input$type == "srs") {
       whatisX <- "Units per group"
       x <- seq(2,50)
       N <- 2*x
       df <- 2*x-2
       M <- qt(input$power,df)-qt(input$alpha/2,df)
       D <- 1
       sample <- x
       if (input$cov == "Yes") {
         df <- df - 1
         M <- qt(input$power,df)-qt(input$alpha/2,df)
         D <- sqrt(1-input$ryx^2)*sqrt((df+1)/(df))
       } 
     }
     else if (input$type == "crt") {
       whatisX <- "Clusters per group"
       x <- seq(2,25)
       N <- 2*input$n*x
       df <- 2*x-2
       M <- qt(input$power,df)-qt(input$alpha/2,df)
       sample <- x*input$n
       D <- sqrt(1+(input$n-1)*input$icc)

       if (input$cov == "Yes") {
         df <- df - 1
         M <- qt(input$power,df)-qt(input$alpha/2,df)
         D <- sqrt(1+(input$n-1)*input$icc-(input$r2u+(input$n*input$r2c-input$r2u)*input$icc))
       }
     }
     else if (input$type == "msrt") {
       whatisX <- "Total clusters"
       x <- seq(2,25)
       N <- 2*input$n*x
       df <- x-1
       M <- qt(input$power,df)-qt(input$alpha/2,df)
       sample <- x*input$n
       D <- sqrt(1+(input$upsilon*input$n/2-1)*input$icc)
       if (input$cov == "Yes") {
         df <- df - 1
         M <- qt(input$power,df)-qt(input$alpha/2,df)
         D <- sqrt(1+(input$upsilon*input$n/2-1)*input$icc-(input$r2u+(input$upsilon*input$n/2*input$r2t-input$r2u)*input$icc))
       }
     }
     mdes <- M*sqrt(2/sample)*D
     return(list(MDES = mdes, x = x, df = df, N = N, M = M, power = input$power, alpha = input$alpha, type = input$type, whatisX = whatisX))
   })
   
   plotMake <- reactive({
     req(input$out)
     req(input$type)
     if (input$out == "ep") {
       req(length(computeCurves()$t) == length(computeCurves()$cden))
       return (
         plot_ly(x = ~computeCurves()$t, 
                 y = ~computeCurves()$cden, 
                 type = 'scatter', 
                 mode = 'lines', 
                 fillcolor = 'white',
                 line = list(color='rgba(0, 0, 0, 1)'),
                 name = 'Central',
                 hoverinfo = 'text',
                 text = ~paste("alpha = ", as.character(input$alpha))) %>%
           add_trace(x = ~c(computeCurves()$alphaR.x, computeCurves()$alphaL.x), 
                     y = ~c(computeCurves()$alphaR.y, computeCurves()$alphaL.y), 
                     fill = 'tozeroy',
                     mode = 'none',
                     name = "alpha",
                     fillcolor = 'rgba(0, 0, 0, 1)',
                     line = list(color='rgba(0, 0, 0, 0)'),
                     hoverinfo = 'text',
                     text = ~paste("alpha = ",as.character(input$alpha))) %>%
           add_trace(x = ~computeCurves()$beta.x, 
                     y = ~computeCurves()$beta.y, 
                     fill = 'tozeroy',
                     name = 'beta',
                     mode = 'lines',
                     fillcolor = 'rgba(168, 216, 234, .5)',
                     line = list(color='rgba(168, 216, 234, 1)'),
                     hoverinfo = 'text',
                     text = ~paste("beta = ",as.character(round(computeCurves()$beta, digits = 3)))) %>%
           add_trace(x = ~computeCurves()$power.x, 
                     y = ~computeCurves()$power.y, 
                     fill = 'tozeroy',
                     name = "Power",
                     mode = 'lines',
                     fillcolor = 'rgba(100, 100, 234, .5)',
                     line = list(color='rgba(100, 100, 234, 1)'),
                     hoverinfo = 'text',
                     text = ~paste("power = ",as.character(round(computeCurves()$power, digits = 3)))) %>%
           add_trace(x = ~computeCurves()$ncp.x, 
                     y = ~computeCurves()$ncp.y, 
                     type = 'scatter',
                     name = "NCP",
                     hoverinfo = 'text',
                     marker = list(size = 12, color='rgba(80, 80, 80, .5)'),
                     line = list(color='rgba(0, 0, 0, 0)'),
                     text = ~paste("ncp = ", as.character(round(computeCurves()$ncp, digits = 3)))) %>%
           layout(
             xaxis = list(title = paste0("Critical = ",
                                         as.character(round(computeCurves()$ct, digits = 3)),
                                         ", NCP = ",
                                         as.character(round(computeCurves()$ncp, digits = 3)),
                                         ", beta = ",
                                         as.character(round(computeCurves()$beta, digits = 3)),
                                         ", power = ",
                                         as.character(round(computeCurves()$power, digits = 3))),
                          range = c(-5,5)),
             yaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE)
           )
       )
       
     }
    else if (input$out == "yPower" & length(computePower()$x) == length(computePower()$power)) {
      req(length(computeCurves()$x) == length(computeCurves()$power))
      return (
        plot_ly(
          x = ~computePower()$x, 
          y = ~computePower()$power,
          type = 'scatter', 
          mode = 'lines',
          hoverinfo = 'text',
          text = paste("Power=",as.character(round(computePower()$power, digits = 3)), ",",
                       computePower()$whatisX,"=", computePower()$x)) %>% 
          add_markers(marker = list(size = 10, color='rgba(80, 80, 80, .5)')) %>%
          layout (showlegend = FALSE,
            xaxis = list(title = computePower()$whatisX),
            yaxis = list(title = "Power", range = c(0,1))
          )
      )
    }
    else if (input$out == "yMDES" & length(computeMDES()$x) == length(computeMDES()$MDES)) {
      req(length(computeCurves()$x) == length(computeCurves()$power)) 
      return (
         plot_ly(
           x = ~computeMDES()$x, 
           y = ~computeMDES()$MDES,
           type = 'scatter', 
           mode = 'lines',
           hoverinfo = 'text',
           text = paste("MDES=",as.character(round(computeMDES()$MDES, digits = 3)), ",",
                        computeMDES()$whatisX,"=", computeMDES()$x)) %>% 
           add_markers(marker = list(size = 10, color='rgba(80, 80, 80, .5)')) %>%
           layout (showlegend = FALSE,
             xaxis = list(title = computeMDES()$whatisX),
             yaxis = list(title = "MDES", range = c(0,3))
           )
       )
     }
     else {
       return(
         plotly_empty(x = 0 , y = 0, text = "loading plot...") %>% add_text(textfont = list(size = 40), textposition = "center")
         )
     }
   })
   
   output$thePlot <- renderPlotly({
     plotMake()
   })
   
   output$showThePlot <- renderUI({
     plotlyOutput("thePlot", height = 400)
   })
   
   tableMake <- reactive({
     if (input$out == "ep") {
       leg <- rbind("es = effect size")
       d <- cbind(es = input$es, n = input$n)
       if (input$type == "srs") {
         leg <- rbind(leg, "n = units per group")
         d <- cbind(d, N = 2*input$n)
         leg <- rbind(leg, "N = total sample")
         if (is.null(input$ryx) == FALSE & input$cov == "Yes") {
           d <- cbind(d, ryx = input$ryx)
           leg <- rbind(leg, "ryx = corr(y,x)")
         }
         d <- cbind(d,
                    ncp = computeCurves()$ncp,
                    ct = computeCurves()$ct,
                    beta = computeCurves()$beta,
                    power = computeCurves()$power)
         leg <- rbind(leg, "ncp = non-centrality parameter")
         leg <- rbind(leg, "ct = critial t")
         leg <- rbind(leg, "beta = type 2 error")
         leg <- rbind(leg, "power = power")
       }
       else if (input$type == "crt") {
         leg <- rbind(leg, "n = units per cluster")
         d <- cbind(d, m = input$m, N = 2*input$n*input$m)
         leg <- rbind(leg, "m = clusters per group")
         leg <- rbind(leg, "N = total sample")
         if (is.null(input$r2u) == FALSE & input$cov == "Yes") {
           d <- cbind(d, r2u = input$r2u)
           leg <- rbind(leg, "r2u = R-square units")
         }
         if (is.null(input$r2c) == FALSE & input$cov == "Yes") {
           d <- cbind(d, r2c = input$r2c)
           leg <- rbind(leg, "r2c = R-square clusters")
         }
         d <- cbind(d,
                    ncp = computeCurves()$ncp,
                    ct = computeCurves()$ct,
                    beta = computeCurves()$beta,
                    power = computeCurves()$power)
         leg <- rbind(leg, "ncp = non-centrality parameter")
         leg <- rbind(leg, "ct = critial t")
         leg <- rbind(leg, "beta = type 2 error")
         leg <- rbind(leg, "power = power")
       }
       else if (input$type == "msrt") {
         leg <- rbind(leg, "n = units per group per cluster")
         d <- cbind(d, m = input$m, N = 2*input$n*input$m)
         leg <- rbind(leg, "m = total clusters")
         leg <- rbind(leg, "N = total sample")
         if (is.null(input$r2u) == FALSE & input$cov == "Yes") {
           d <- cbind(d, r2u = input$r2u)
           leg <- rbind(leg, "r2u = R-square units")
         }
         if (is.null(input$r2t) == FALSE & input$cov == "Yes") {
           d <- cbind(d, r2t = input$r2t)
           leg <- rbind(leg, "r2t = R-square treatment")
         }
         d <- cbind(d,
                    ncp = computeCurves()$ncp,
                    ct = computeCurves()$ct,
                    beta = computeCurves()$beta,
                    power = computeCurves()$power)
         leg <- rbind(leg, "ncp = non-centrality parameter")
         leg <- rbind(leg, "ct = critial t")
         leg <- rbind(leg, "beta = type 2 error")
         leg <- rbind(leg, "power = power")
       }
       return(list(d = d, leg = leg))
     }
     else if (input$out == "yPower") {
       return(computePower())
     }
     else if (input$out == "yMDES") {
       return(computeMDES())
     }
     else {
       return(c())
     }
   })
   
   output$theTable <- renderTable({
     if (input$out == "ep") {
       tableMake()$d
     }
     else {
       tableMake()
     }
     
   }, include.rownames=FALSE)
   
   output$theLegend <- renderTable({
     req(tableMake()$leg)
     tableMake()$leg
   }, include.colnames=FALSE)
   
   output$showTheTable <- renderUI({
     tableOutput("theTable")
   })
   
   output$showTheLegend <- renderUI({
     if (input$out == "ep") {
       tableOutput("theLegend")
     }
     
   })
   
   output$downloadTable <- downloadHandler(
     filename = function() {
       paste(input$out,input$type,Sys.Date(), ".csv", sep="")
     },
     content = function(file) {
       write.csv(tableMake()$d, file)
     }
   )
}

# Run the application 
shinyApp(ui = ui, server = server)

