library(shiny)
library(shinythemes)
library(ggplot2)
library(scales)          
library(lubridate)
library(data.table)
library(shinyjs)
library(ncdf4)

Climate <- fread("Germany2.csv", sep = ',', stringsAsFactors = FALSE)
#options(shiny.port = ****)

ui <- fluidPage(
  tags$div(
    style = "
    position: absolute;
    top: 80px;       /* distance from top */
    bottom: auto;    /* set to 'auto' if not using */
    left: auto;      /* or a value like '30px' */
    right: 20px;     /* distance from right */
    z-index: 999;",
    img(src = "www/BNITM.png", height = "140px")
  ),
  theme = shinytheme("cerulean"),
  titlePanel("Welcome to the Zero-West-Nile-virus App for exploring different WNV control strategies in Germany."),
  tags$p("The app allows you to view simulated West Nile Virus infections in birds, humans and equids for a period of your choice. Select the period of your choice first and then click the Orange icon labeled GO!!! and the solver will dynamically solve the model based on your inputs."),
  tags$p("When using a mobile phone, it is best to rotate the screen for better visuals. Enjoy exploring and after exploring, please leave a review of the app based on how much you like it.!!!"),
  
  tags$p("Please enter your desired dates:"),
  sidebarLayout(
    sidebarPanel(width = 3,
                 dateInput("start_date", "Start Date", value = "2023-01-01", min = "2023-01-01", max = "2028-12-31"),
                 dateInput("end_date", "End Date", value = "2028-12-31", min = "2023-01-01", max = "2028-12-31"),
                 actionButton("apply", "GO!!!", class = "btn btn-primary btn-lg"),
                 sliderInput("u1", "Physical removal and destruction of potential breeding sites (BSR.)", min = 0, max = 0.9, value = 0, step = 0.25),
                 sliderInput("u2", "Larvicides (La.)", min = 0, max = 0.9, value = 0, step = 0.25),
                 sliderInput("u3", "Adulticides (Ad.)", min = 0, max = 0.9, value = 0, step = 0.25),
                 sliderInput("u4", "Individuals' sensitivity to WNV prevalence in birds", min = 0, max = 200, value = 0, step = 50),
                 sliderInput("u5", "Vaccination of equids (Vacc.)", min = 0, max = 0.9, value = 0, step = 0.25)),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Infectious birds", plotOutput("plot1"), tableOutput("error_output1")),
        tabPanel("Infected humans", plotOutput("plot2"), tableOutput("error_output2")),
        tabPanel("Infected equids", plotOutput("plot3"), tableOutput("error_output3")),
        tabPanel("Reviews",
                 column(10, offset = 1,
                        style = "margin-top: 50px;",
                        #h4("User Reviews"),
                        textInput("reviewer_name", "Enter a nickname (For privacy, please use a nickname or leave it blank to submit as 'Anonymous'):", placeholder = "Anonymous"),
                        sliderInput("rating", "Rate this app:", min = 1, max = 5, value = 3, step = 1),
                        textAreaInput("review_comment", "We would love to have your feedback!!!", 
                                      rows = 4, placeholder = "What do you think about this app?"),
                        actionButton("submit_rating", "Submit Review"),
                        br(), br(), br(),
                        tableOutput("ratings_table")))))),
  
  tags$hr(style = "margin-top: 10px; margin-bottom: 5px;"),
  
  tags$div(
    style = "text-align: center; font-size: 16px; color: black;",
    "This work was supported by the ",tags$a(href = "https://www.bmuv.de/en/", target = "_blank", "Federal Ministry for the Environment, Nature Conservation, Nuclear Safety and Consumer Protection"),
    " (Grant Number 3721484020), the", tags$a(href = "https://www.bmbf.de/EN/Home/home_node.html", target = "_blank", "Federal Ministry of Education and Research"),
    "(Grant Number 01Kl2022), and the", tags$a(href = "https://www.dfg.de/en", target = "_blank", "German Research Foundation"),
    "(Grant Number JO 1276/5-1)"
  ))

server <- function(input, output, session) {
  
  parameters <- reactiveVal()
  submitted <- reactiveVal(FALSE)
  
  observeEvent(input$apply, {
    parameters(list(
      input$start_date,
      input$end_date,
      input$u1,
      input$u2,
      input$u3,
      input$u4,
      input$u5
    ))
    
    german_states <- data.frame(
      State = c("Germany"),
      Latitude = c(48.79),
      Longitude = c(11.49))
    
    state_coords <- german_states
    
    # Define state coordinates
    state_lat <- state_coords$Latitude
    state_lon <- state_coords$Longitude
    
    # Open the NetCDF file
    file_path <- "Temperature_Data.nc"
    nc_data <- nc_open(file_path)
    lon <- ncvar_get(nc_data, "lon")
    lat <- ncvar_get(nc_data, "lat")
    time <- ncvar_get(nc_data, "time")
    t2m_daily <- ncvar_get(nc_data, "tas")
    time_units <- ncatt_get(nc_data, "time", "units")$value
    time_origin <- sub(".*since ", "", time_units)  
    dates <- as.Date(time, origin = time_origin)
    lat_index <- which.min(abs(lat - state_lat))
    lon_index <- which.min(abs(lon - state_lon))
    
    # Define time range
    start_date <- as.Date("2023-01-01")
    end_date <- as.Date("2028-12-31")
    
    # Filter indices for the chosen time period
    selected_indices <- which(dates >= start_date & dates <= end_date)
    
    # Extract daily temperature for the selected state and period
    T <- t2m_daily[lon_index, lat_index, selected_indices]
    
    selected_data <- data.frame(
      Date = dates[selected_indices],
      Temperature = T
    )
    
    Temperatur <- selected_data$Temperature
    
    nc_close(nc_data)
    
    head(selected_data)
    
    n <- length(T) - 1
    
    t<-seq(0,n)
    
    KM=3300000
    k=rep(0, n); k=0.344/(1+1.231*exp(-0.184*(T-20)))
    bL=rep(0, n); bL=2.325*k; bM=bL/10
    mL=rep(0, n); mL=0.0025*T^2-0.094*T+1.0257; mM=mL/10
    phi=45*pi/180
    day <- c(rep(seq(1,365),2),seq(1,366),rep(seq(1,365),3),seq(1,366),rep(seq(1,365)))
    epsilon=rep(0, n); epsilon=0.409*sin(2*pi*(day-80)/365)
    D=rep(0, n); D=7.639*asin(tan(epsilon)*tan(phi)+0.0146/(cos(epsilon)*cos(phi)))+12
    dM=rep(0, n); dM=1-1/(1+1775.7*exp(1.559*(D-18.177)))
    
    gammaM=rep(0, n+1)
    for (i in 1:n){if (T[i]>15) gammaM[i]=0.0093*T[i]-0.1352}
    
    # bird parameters
    KB=110000
    phiB=30
    bB=rep(0, n+1)
    bB=dgamma(day, scale=1.4, shape=86.4)
    mB=0.00034
    alphaB=0.4 
    gammaB=1.0 
    nuB=0.7
    
    # equid parameters
    
    phiE=300
    gammaE=0.05 
    alphaE=0.2 
    nuE=0.04
    bE=0.00016
    mE=0.00011
    omegaE=0.1; cE=0.95;
    thetaE1= 0.0; omegaE1=0.1; cE1=0.75;
    # human parameters
    
    phiH=0.03
    gammaH=0.25 
    alphaH=0.5
    nuH=0.004
    bH=0.000055
    mH=0.000034
    rH=bH-mH
    omegaH=0.1; cH=0.5;
    thetaH1= 0.0; omegaH1=0.1; cH1=0.5;
    
    # transmission parameters
    
    pM=1.0
    pB=0.125
    
    betaM=rep(0, n); betaM=k*pM
    betaB=rep(0, n); betaB=k*pB
    
    # initialisation
    
    LM=rep(0, n+1); SM=rep(500000, n+1); EM=rep(0, n+1); IM=rep(100, n+1); s=rep(0.01, n+1);
    NMmin=SM[1]
    NM=SM+EM+IM	
    
    SB=rep(110000, n+1); EB=rep(0, n+1); IB=rep(0, n+1); RB=rep(0, n+1); DB=rep(0, n+1)
    NB=SB+EB+IB+RB
    
    SE=rep(10700, n+1);VE=rep(1, n+1); EE=rep(0, n+1); IE=rep(0, n+1); RE=rep(0, n+1); DE=rep(0, n+1)
    NE=SE+VE+EE+IE+RE
    
    SH=rep(2700000, n+1); EH=rep(0, n+1); IH=rep(0, n+1); RH=rep(0, n+1); DH=rep(0, n+1)
    NH=SH+EH+IH+RH
    
    LM1=rep(0, n+1); SM1=rep(500000, n+1); EM1=rep(0, n+1); IM1=rep(100, n+1)
    NMmin1=SM1[1]
    NM1=SM1+EM1+IM1	
    
    SB1=rep(110000, n+1); EB1=rep(0, n+1); IB1=rep(0, n+1); RB1=rep(0, n+1); DB1=rep(0, n+1)
    NB1=SB1+EB1+IB1+RB1
    
    SE1=rep(10700, n+1);VE1=rep(1, n+1); EE1=rep(0, n+1); IE1=rep(0, n+1); RE1=rep(0, n+1); DE1=rep(0, n+1)
    NE1=SE1+VE1+EE1+IE1+RE1
    
    SH1=rep(2700000, n+1); EH1=rep(0, n+1); IH1=rep(0, n+1); RH1=rep(0, n+1); DH1=rep(0, n+1)
    NH1=SH1+EH1+IH1+RH1
    
    
    for(i in 1:n){
      
      lambdaBM1=dM[i]*betaB[i]*IB1[i]/KB
      lambdaMB1=dM[i]*betaM[i]*phiB*IM1[i]/KM
      lambdaME1=dM[i]*betaM[i]*phiE*IM1[i]/KM
      lambdaMH1=dM[i]*betaM[i]*phiH*IM1[i]/KM
      
      LM1[i+1]=LM1[i] + bL[i]*dM[i]*NM1[i]*(1-LM1[i]/KM)-bM[i]*LM1[i] - mL[i]*LM1[i]*(1-LM1[i]/KM)
      SM1[i+1]=SM1[i] - lambdaBM1*SM1[i]+bM[i]*LM1[i]-mM[i]*SM1[i]
      EM1[i+1]=EM1[i] + lambdaBM1*SM1[i]-gammaM[i]*EM1[i]-mM[i]*EM1[i]
      IM1[i+1]=IM1[i] + gammaM[i]*EM1[i]-mM[i]*IM1[i]
      NM1[i+1]=SM1[i+1]+EM1[i+1]+IM1[i+1]
      
      if (NM1[i+1]<NMmin1) {SM1[i+1]=SM1[i]; EM1[i+1]=EM1[i]; IM1[i+1]=IM1[i]}
      if (LM1[i+1]<0) {LM1[i+1]=0}
      
      SB1[i+1]=SB1[i] + (bB[i]-(bB[i]-mB)*NB1[i]/KB)*NB1[i]-lambdaMB1*SB1[i]-mB*SB1[i]
      EB1[i+1]=EB1[i] + lambdaMB1*SB1[i]-gammaB*EB1[i]-mB*EB1[i]
      IB1[i+1]=IB1[i] + gammaB*EB1[i]-alphaB*IB1[i]-mB*IB1[i]
      RB1[i+1]=RB1[i] + (1-nuB)*alphaB*IB1[i]-mB*RB1[i]
      DB1[i+1]=DB1[i] + nuB*alphaB*IB1[i]
      NB1[i+1]=SB1[i+1]+EB1[i+1]+IB1[i+1]+RB1[i+1]
      
      SE1[i+1]=SE1[i] + bE*NE1[i]-lambdaME1*SE1[i]-mE*SE1[i]-thetaE1*SE1[i]+omegaE1*VE1[i]
      VE1[i+1]=VE1[i] + thetaE1*SE1[i]-omegaE1*VE1[i]-mE*VE1[i] -cE1*lambdaME1*VE1[i]
      EE1[i+1]=EE1[i] + lambdaME1*SE1[i]+cE1*lambdaME1*VE1[i]-gammaE*EE1[i]-mE*EE1[i]
      IE1[i+1]=IE1[i] + gammaE*EE1[i]-alphaE*IE1[i]-mE*IE1[i]
      RE1[i+1]=RE1[i] + (1-nuE)*alphaE*IE1[i]-mE*RE1[i]
      DE1[i+1]=DE1[i] + nuE*alphaE*IE1[i]
      NE1[i+1]=SE1[i+1]+VE1[i+1]+EE1[i+1]+IE1[i+1]+RE1[i+1]
      
      SH1[i+1]=SH1[i] + bH*NH1[i]-lambdaMH1*SH1[i] -mH*SH[i]
      EH1[i+1]=EH1[i] + lambdaMH1*SH1[i]-gammaH*EH1[i]-mH*EH1[i]
      IH1[i+1]=IH1[i] + gammaH*EH1[i]-alphaH*IH1[i]-mH*IH1[i]
      RH1[i+1]=RH1[i] + (1-nuH)*alphaH*IH1[i]-mH*RH1[i]
      DH1[i+1]=DH1[i] + nuH*alphaH*IH1[i]
      NH1[i+1]=SH1[i+1]+EH1[i+1]+IH1[i+1]+RH1[i+1]
      
      lambdaBM=dM[i]*betaB[i]*IB[i]/KB
      lambdaMB=dM[i]*betaM[i]*phiB*IM[i]/KM
      lambdaME=dM[i]*betaM[i]*phiE*IM[i]/KM
      lambdaMH=dM[i]*betaM[i]*phiH*IM[i]/KM
      
      LM[i+1]=LM[i] + (1-input$u1/40)*bL[i]*dM[i]*NM[i]*(1-LM[i]/KM)-bM[i]*LM[i] - (mL[i]+input$u1/40+input$u2/10)*LM[i]*(1-LM[i]/KM)
      SM[i+1]=SM[i] - lambdaBM*SM[i]+bM[i]*LM[i]-(mM[i]+input$u3/400)*SM[i]
      EM[i+1]=EM[i] + lambdaBM*SM[i]-gammaM[i]*EM[i]-(mM[i]+input$u3/400)*EM[i]
      IM[i+1]=IM[i] + gammaM[i]*EM[i]-(mM[i]+input$u3/400)*IM[i]
      NM[i+1]=SM[i+1]+EM[i+1]+IM[i+1]
      
      if (NM[i+1]<NMmin) {SM[i+1]=SM[i]; EM[i+1]=EM[i]; IM[i+1]=IM[i]}
      if (LM[i+1]<0) {LM[i+1]=0}
      
      SB[i+1]=SB[i] + (bB[i]-(bB[i]-mB)*NB[i]/KB)*NB[i]-lambdaMB*SB[i]-mB*SB[i]
      EB[i+1]=EB[i] + lambdaMB*SB[i]-gammaB*EB[i]-mB*EB[i]
      IB[i+1]=IB[i] + gammaB*EB[i]-alphaB*IB[i]-mB*IB[i]
      RB[i+1]=RB[i] + (1-nuB)*alphaB*IB[i]-mB*RB[i]
      DB[i+1]=DB[i] + nuB*alphaB*IB[i]
      NB[i+1]=SB[i+1]+EB[i+1]+IB[i+1]+RB[i+1]
      
      SE[i+1]=SE[i] + bE*NE[i]-lambdaME*SE[i] - mE*SE[i] - SE[i]*input$u5/10 + omegaE*VE[i]
      VE[i+1]=VE[i] + SE[i]*input$u5/10-omegaE*VE[i]-mE*VE[i] -(1-cE)*lambdaME*VE[i]
      EE[i+1]=EE[i] + lambdaME*SE[i]+(1-cE)*lambdaME*VE[i]-gammaE*EE[i]-mE*EE[i]
      IE[i+1]=IE[i] + gammaE*EE[i]-alphaE*IE[i]-mE*IE[i]
      RE[i+1]=RE[i] + (1-nuE)*alphaE*IE[i]-mE*RE[i]
      DE[i+1]=DE[i] + nuE*alphaE*IE[i]
      NE[i+1]=SE[i+1]+VE[i+1]+EE[i+1]+IE[i+1]+RE[i+1]
      
      SH[i+1]=SH[i] + bH*NH[i]-(1-0.15*s[i])*lambdaMH*SH[i] -mH*SH[i]
      EH[i+1]=EH[i] + (1-0.15*s[i])*lambdaMH*SH[i]-gammaH*EH[i]-mH*EH[i]
      IH[i+1]=IH[i] + gammaH*EH[i]-alphaH*IH[i]-mH*IH[i]
      RH[i+1]=RH[i] + (1-nuH)*alphaH*IH[i]-mH*RH[i]
      DH[i+1]=DH[i] + nuH*alphaH*IH[i]
      s[i+1]=s[i] + 0.00003*s[i]*(1-s[i])*(-1+input$u4*IB1[i])
      NH[i+1]=SH[i+1]+EH[i+1]+IH[i+1]+RH[i+1]
    }
    
    start_date <- as.Date(input$start_date)
    end_date <- as.Date(input$end_date)
    
    output$plot1 <- renderPlot({
      y_breaks <- c(min(IB1), max(IB1))  
      y_labels <- c("Low", "High")  
      
      p <- ggplot(data = selected_data, aes(x = Date)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = y_labels) +
        xlab("") +
        ylab("Simulated population density") +
        scale_x_date(breaks = date_breaks("6 months"), 
                     labels = date_format("%b %y"), 
                     limits = c(start_date, end_date)) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1.5, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          legend.key.height = unit(1.5, 'cm'),
          axis.title.y = element_text(size = 20, face = 'bold',, color="black"),
          axis.text.y = element_text(size = 24, face = "bold", color="black"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 18, face = "bold")
        ) +
        labs(color = "Population Density") +  
        theme(legend.position = "top") +  
        scale_color_manual(
          breaks = c("No intervention", "10% Threshold", "With intervention"),
          values = c("No intervention" = "blue","10% Threshold"="green", "With intervention" = "red"))
      
      # Add lines based on user input
      if (input$u1 == 0 && input$u2 == 0 && input$u3 == 0&& input$u4 == 0 && input$u5 == 0) {
        p <- p + geom_line(data = selected_data, aes(x = Date, y = IB1, group = 1, color = "No intervention"), size = 1.8)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.1*max(IB1), group = 1, color = "10% Threshold"),linetype = "dashed", size = 1.8)
      } else {
        p <- p + geom_line(data = selected_data, aes(x = Date, y = IB, group = 1, color = "With intervention"), size = 1.2, linetype = "dashed") +
          geom_line(data = selected_data, aes(x = Date, y = IB1, group = 1, color = "No intervention"), size = 1.2)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.1*max(IB1), group = 1, color = "10% Threshold"),linetype = "dashed", size = 1.8)
      } 
      print(p)
    })
    
    output$error_output1 <- renderText({
      err1 <- max(abs(IB1 - IB)/(0.0000000001+IB1)*100) 
      formatted_err1 <- format(err1, scientific = FALSE)
      paste("Infections averted:", formatted_err1,"%")})
    
    
    ######## Second plot:###################################
    output$plot2 <- renderPlot({
      y_breaks <- c(min(IH1), max(IH1))  
      y_labels <- c("Low", "High")  
      
      p1 <- ggplot(data = selected_data, aes(x = Date)) +
        scale_y_continuous(breaks = y_breaks, # Set custom y-axis breaks
                           labels = y_labels) +
        xlab("") +
        ylab("Simulated population density") +
        scale_x_date(breaks = date_breaks("6 months"), 
                     labels = date_format("%b %y"), 
                     limits = c(start_date, end_date)) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1.5, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          legend.key.height = unit(1.5, 'cm'),
          axis.title.y = element_text(size = 20, face = 'bold',, color="black"),
          axis.text.y = element_text(size = 24, face = "bold", color="black"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 18, face = "bold")
        ) +
        labs(color = "Population Density") +  
        theme(legend.position = "top") +  
        scale_color_manual(
          breaks = c("No intervention", "10% Threshold", "With intervention"),
          values = c("No intervention" = "blue","10% Threshold"="green", "With intervention" = "red"))
      
      # Add lines based on user input
      if (input$u1 == 0 && input$u2 == 0 && input$u3 == 0&& input$u4 == 0 && input$u5 == 0) {
        p1 <- p1 + geom_line(data = selected_data, aes(x = Date, y = IH1, group = 1, color = "No intervention"), size = 1.8)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.1*max(IH1), group = 1, color = "10% Threshold"),linetype = "dashed", size = 1.8)
      } else {
        p1 <- p1 + geom_line(data = selected_data, aes(x = Date, y = IH, group = 1, color = "With intervention"), size = 1.2, linetype = "dashed") +
          geom_line(data = selected_data, aes(x = Date, y = IH1, group = 1, color = "No intervention"), size = 1.2)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.1*max(IH1), group = 1, color = "10% Threshold"),linetype = "dashed", size = 1.8)
      } 
      print(p1)
      
      
      output$error_output2 <- renderText({
        err2 <- max(abs(IH1 - IH) / (1e-18 + IH1)*100)
        rounded_err2 <- floor(err2)
        paste("Infections averted:", rounded_err2, "%")
      })
      
    })
    
    
    ## Third plot
    
    output$plot3 <- renderPlot({
      y_breaks <- c(min(IE1), max(IE1))  
      y_labels <- c("Low", "High")  
      
      p2 <- ggplot(data = selected_data, aes(x = Date)) +
        scale_y_continuous(breaks = y_breaks,
                           labels = y_labels) +
        xlab("") +
        ylab("Simulated population density") +
        scale_x_date(breaks = date_breaks("6 months"), 
                     labels = date_format("%b %y"), 
                     limits = c(start_date, end_date)) +
        theme_classic() +
        theme(
          plot.title = element_text(hjust = 0.5),
          legend.title = element_blank(),
          legend.text = element_text(size = 16),
          legend.key.size = unit(1.5, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          legend.key.height = unit(1.5, 'cm'),
          axis.title.y = element_text(size = 20, face = 'bold',, color="black"),
          axis.text.y = element_text(size = 24, face = "bold", color="black"),
          axis.text = element_text(size = 11, face = "bold"),
          axis.title = element_text(size = 18, face = "bold")
        ) +
        labs(color = "Population Density") +  
        theme(legend.position = "top") +  
        scale_color_manual(
          breaks = c("No intervention", "5% Threshold", "With intervention"),
          values = c("No intervention" = "blue","5% Threshold"="green", "With intervention" = "red"))
      
      # Add lines based on user input
      if (input$u1 == 0 && input$u2 == 0 && input$u3 == 0&& input$u4 == 0 && input$u5 == 0) {
        p2 <- p2 + geom_line(data = selected_data, aes(x = Date, y = IE1, group = 1, color = "No intervention"), size = 1.8)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.05*max(IE1), group = 1, color = "5% Threshold"),linetype = "dashed", size = 1.8)
      } else {
        p2 <- p2 + geom_line(data = selected_data, aes(x = Date, y = IE, group = 1, color = "With intervention"), size = 1.2, linetype = "dashed") +
          geom_line(data = selected_data, aes(x = Date, y = IE1, group = 1, color = "No intervention"), size = 1.2)+
          geom_hline(data = selected_data, aes(x = Date, yintercept = 0.05*max(IE1), group = 1, color = "5% Threshold"),linetype = "dashed", size = 1.8)
      } 
      print(p2)
      
      
      output$error_output3 <- renderText({
        err3 <- max(abs(IE1 - IE) / (1e-18 + IE1)*100)
        rounded_err3 <- floor(err3)
        paste("Infections averted:", rounded_err3, "%")
      })
    })
  })
  
  
  file_path <- "Reviews.csv"
  
  ratings <- reactiveValues(data = NULL)
  
  # Check if the CSV file exists and is not empty
  if (file.exists(file_path) && file.size(file_path) > 0) {
    ratings$data <- read.csv(file_path, stringsAsFactors = FALSE)
  } else {
    # Initialize an empty data frame with headers
    ratings$data <- data.frame(
      Reviewer = character(), 
      Rating = numeric(), 
      Date = character(), 
      Comments = character(),
      stringsAsFactors = FALSE
    )
    write.csv(ratings$data, file_path, row.names = FALSE)
  }
  
  # Observe event for when submitting a rating
  observeEvent(input$submit_rating, {
    reviewer <- ifelse(input$reviewer_name == "", "Anonymous", input$reviewer_name)
    review_date <- Sys.Date()
    review_comment <- input$review_comment
    
    new_rating <- data.frame(
      Reviewer = reviewer,
      Rating = input$rating,
      Date = as.character(review_date),
      Comments = review_comment,
      stringsAsFactors = FALSE
    )
    
    if (!any(apply(ratings$data, 1, function(row) all(row == unlist(new_rating))))) {
      ratings$data <- rbind(ratings$data, new_rating)
      write.csv(ratings$data, file_path, row.names = FALSE)
      submitted(TRUE)
      
      updateTextInput(session, "reviewer_name", value = "")
      updateSliderInput(session, "rating", value = 3)
      updateTextAreaInput(session, "review_comment", value = "")
    }
    
    
  })
  
  # Render the ratings table
  output$ratings_table <- renderTable({ratings$data})
}

shinyApp(ui = ui, server = server)
