library(shiny)
library(shinydashboard)
library(tidyverse)
library(rTensor)

source("AlphaNTF.R", local=T)
source("BetaNTF.R", local=T)
source("FastHALS.R", local=T)
source("BCD.R", local=T)

source("TensorGenerators.R", local=T)

set.seed(900707)

prob_type <- c(0.4,0.1,0.1,0.4)
prob_connect <- matrix(c(0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05,0.5,0.5,0.25,0.25,0.25,0.25,0.05,0.05),nrow=4)
xm_prob <- c(0.9,0.75,0.75,0.5) # Probability of Connection Remaining After Evaluation


ui <- dashboardPage(
  
  dashboardHeader(title = "Tensor Decomposition"),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Tensor Decomposition", tabName = "tensor", icon = icon("line-chart"))
    )
  ),
  
  dashboardBody(
    
    tabItems(
      
      tabItem(tabName = "tensor",
              
              box(title="Tensor Properties",
                  sliderInput("matchruns","# of Match Runs",min=8,max=30,value=10,step=1),
                  sliderInput("permatchrun","Nodes Per Match Run",min=1,max=10,value=4,step=1),
                  sliderInput("components","Number of Components",min=2,max=5,value=2,step=1)
              ),
              
              box(title="Parameters",
                  
                  selectInput("method","Method",c("CP","HALS","Alpha","Beta","BCD"), selected="CP"),
                  selectInput("tolerance","Tolerance",c(1e-4,1e-5,1e-6), selected=1e-4),
                  sliderInput("maxits","Max Iterations",min=1000,max=10000,value=1000,step=1000),
                  numericInput("tuning","Tuning Parameter",min=0,max=10,value=2,step=0.05)
              ),
              
              box(title="Components Plot",
                  plotOutput("componentplot", height = 400)
              ),
              
              box(title="Convergence Plot",
                  plotOutput("convergenceplot", height = 400)
              )
              
      )
      
    )
    
  )
)


server <- function(input, output){
  
  update.tensor <- reactive({
    
    per_mr <- as.numeric(input$permatchrun)
    num_mrs <- as.numeric(input$matchruns)
    
    return(generate.kpd(per_mr, num_mrs, prob_type, prob_connect, xm_prob))
    
  })
  
  update.initialization <- reactive({
    
    per_mr <- as.numeric(input$permatchrun)
    num_mrs <- as.numeric(input$matchruns)
    num_components <- as.numeric(input$components)
    
    init_matrix <- initialize.matrix(per_mr, num_mrs, num_components)
    
    return(init_matrix)
    
  })
  
  update.results <- reactive({
    
    kpd.tensor <- update.tensor()
    
    kpd <- kpd.tensor$kpd
    node_types <- kpd.tensor$node_types
    
    init_matrix <- update.initialization()
    
    maxits <- as.numeric(input$maxits)
    tol <- as.numeric(input$tolerance)
    tune <- input$tuning
    
    if (input$method=="HALS"){
      
      decomp <- fastHALS(kpd, init_matrix, max_k = maxits, tol = tol)
      
      B <- decomp$A
      norms <- decomp$objs
      value <- decomp$objs %>% last
      
    } else if (input$method=="Alpha"){
      
      decomp <- alphaNTF(kpd, init_matrix, alpha = tune, max_k = maxits, tol = tol)
      
      B <- decomp$A
      norms <- decomp$objs
      value <- decomp$objs %>% last
      
    } else if (input$method=="Beta"){
      
      decomp <- betaNTF(kpd, init_matrix, beta = tune, max_k = maxits, tol = tol)
      
      B <- decomp$A
      norms <- decomp$objs
      value <- decomp$objs %>% last
      
    } else if (input$method=="BCD"){
      
      decomp <- BCD(kpd, init_matrix, delta = tune, max_k = maxits, tol = tol)
      
      B <- decomp$A
      norms <- decomp$objs
      value <- decomp$objs %>% last
      
    } else {
      
      decomp <- cp(kpd, num_components=as.numeric(input$components), max_iter = maxits, tol = tol)
      
      B <- decomp$U
      norms <- decomp$all_resids
      value <- decomp$fnorm_resid
      
    }
    
    return(list(B=B,norms=norms,value=value))
    
  })
  
  output$componentplot <- renderPlot({
    
    kpd.tensor <- update.tensor()
    
    node_types <- kpd.tensor$node_types
    
    sim.results <- update.results()
    
    B <- sim.results$B
    norms <- sim.results$norms
    value <- sim.results$value
    
    node.frame <- rbind(data.frame(B[[1]],Node_Type=as.factor(node_types),u="Donor"),
                        data.frame(B[[2]],Node_Type=as.factor(node_types),u="Candidate"))
    
    
    p <- ggplot(data=node.frame, aes(x=X1,y=X2,color=Node_Type)) +
      facet_wrap(~u) +
      geom_point(size=4,alpha=0.5) +
      theme_bw() +
      xlab("Component 1") +
      ylab("Component 2") +
      ggtitle("Component Plot")
    
    return(p)
    
  })
  
  output$convergenceplot <- renderPlot({
    
    kpd.tensor <- update.tensor()
    
    node_types <- kpd.tensor$node_types
    
    sim.results <- update.results()
    
    norms <- sim.results$norms
    value <- sim.results$value
    
    norm.table <- data.frame(Iteration=1:length(norms),Objective=norms)
    
    p <- ggplot(data=norm.table,aes(x=Iteration,y=Objective)) +
      geom_path() +
      theme_bw() +
      ggtitle("Convergence Plot") +
      geom_hline(yintercept=value,color="red",linetype="dashed") +
      geom_label(aes(label=paste0(round(value,3)),x=10,y=value, color="red"),show.legend=FALSE)
    
    
    return(p)
    
  })
}

shinyApp(ui, server)