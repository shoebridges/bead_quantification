#Install Packages-----------------
library(shiny)
library(matrixStats)
library(data.table)
library(tidyverse)
library(limma)
library(ggpubr)
library(limma)
library(rstatix)
library(nlme)
library(emmeans)
library(shinydashboard)
library(DT)
#custom functions
split_tibble <- function(tibble, col = 'col') tibble %>% split(., .[, col])

ui <- dashboardPage(
  dashboardHeader(title="Bead Quantification: Ver 1.1"),
  dashboardSidebar(width = 400, h2("Data Input"), 
                   p("Insert your csv files for analysis here"),
                     fileInput("files1", "Choose csv files from directory", multiple = T, accept = ".csv"), 
                   h2("Output"),
                   p("Merged Table"),
                     downloadButton("final_data", "Download merged table with exclusions", icon = icon("th-list", lib = "glyphicon")),
                   p("Statistics Table"), 
                      downloadButton("stat_test", "Download statistics", icon = icon(lib = "glyphicon", name = "search")),
                   p("Plots of Dataset"),
                     downloadButton("grouped_plot_raw", label = "Download grouped plot (raw)", icon = icon(lib = "glyphicon", name = "stats")),
                     downloadButton("grouped_plot_norm", label = "Download grouped plot (control-normalised)", icon = icon(lib = "glyphicon", name = "stats")),
                     downloadButton("histo_raw", label = "Download histogram of raw measurements", icon = icon(lib = "glyphicon", name = "stats")),
                     downloadButton("histo_corr", label = "Download histogram of normalised measurements", icon = icon(lib = "glyphicon", name = "stats")),
                   h2("Useful Documents"),
                   menuItem("Guide to Bead Quantification", href = "https://martens.apps.maxperutzlabs.ac.at/bead_quant_guide/A-guide-to-Bead-Quant.html", icon = icon("book", lib="glyphicon"), badgeLabel = "In development", badgeColor="red"),
                   tags$style(".skin-blue .sidebar .shiny-download-link { color: #444; }")),
  dashboardBody(
    fluidRow(box(title = "Introduction and method", height = 500,
                 h3("Introduction"),
                 p("The following programme analyses bead quantification data. Each csv file should have the following structure:"),
                 p(strong("User_Condition_BiologicalRepeat_Field_Channel_Date.csv")),
                 p("The fields channel and date are usually entered automatically. While the server can handle files that are not
                     .csv files, this will make the programme slower initially as it has to parse through the data. For this sake
                     please only submit .csv files to the work."),
                 h3("Methods"),
                 p("Significance is tested through a nested ANOVA. All statistics are performed on non-normalised/corrected data (raw values). Data that has been
                      corrected by a positive control is performed where the raw value is divided by the average value for the positive control with respect to
                      experiment. If the data is corrected with negative control the raw values are subtracted by the average negative control value for experiment.
                      If both controls are used then the data is subtracted by the negative control and then divided by the adjusted positive control. \n
                   To process this dataset enter the csv files into the panel on the right, select the settings you wish to use. Observe the warnings to see how many samples
                   have been exlcuded with the current settings. Observe the statistical test and plots produced. If you need to download the plots or statistics table,
                   the can be downloaded in the left panel."),
                 h3("Changelog"),
                 tagList("Statistical analysis can now handle expereiments with one biological replicate.", "Statistics have been improved to favour experiment replicate over bead replicate
                   within experiment.", "Some graph issues with one experimental replicate have been repaired."),),
             box(title="Settings", height = 500,
                 p("Please choose the exclusion criteria of your data (removing bad analysis points):"),
                 checkboxInput(inputId = "SD_Exclude", label = "SD is greater than 50% raw average", value = T),
                 checkboxInput(inputId = "Background_Exclude", label = "background is higher than raw average", value = T),
                 checkboxInput(inputId = "Neg_Exclude", label = "Raw Average is missing or below 0", value = T),
                 p("For relative (control adjusted) plots, please enter the label of negative control and choose adjust:"),
                 textInput("neg_control", label = "Insert the name of the negative control:"),
                 textInput("pos_control", label = "Insert the name of the positive control:"),
                actionButton("adjust", "Apply or adjust settings")),
    box(title = "Duplicates and Warnings", strong("SD is greater than 50% raw average"), textOutput("SD_Warn"),
        br(), strong("Background is higher than raw average"), textOutput("Neg_Warn"), 
        br(), strong("Duplicate entries"), tableOutput("duplicates"), br(), strong("Number of meauremnts excluded: "), textOutput("merge_final"), br(), strong("Number of Channels: "), textOutput("channels"), height= 400),#END OF BOX4
    box(title = "Files used", dataTableOutput("filenames", height=2), height = 400), #Files used table BOX2, One line
    box(title = "Preview", dataTableOutput("table1"), width = 12), #BOX3
    
    box(title = "Statistics (Nested ANOVA: Tukey Adjusted)", textOutput("test"),tableOutput("quant_pairs"),
            footer = "  Note that the Bonferroni adjustment is indeed an adjustment
             to the regular P values (p-value adjusted). But the Tukey-adjusted P values are an entirely separate calculation based on the t ratios
             but not on p-val-adj. So there is no default set of P values underlying the Tukey-adjusted P values."), #END OF BOX5
    box(title = "Graph", plotOutput("raw_plot"), plotOutput("cor_plot"), width = 12),
    box(title = "Histograms", plotOutput("raw_histo", height = "800px"), plotOutput("corr_histo", height = 500), width = 12)
        ) #Bracket for fluidRow
  ) #Bracket for dashboardBody
) # Bracket for Page
# Define server logic 
server <- function(input, output) {
  input_filelist <- reactive({
    filelist <- input$files1
    if (is.null(input$files1)) return()
    ext <- filelist$name
    stringr::str_subset(ext, ".csv")
  })
  

  output$filenames <- renderDataTable({as.data.frame(input_filelist())})
  
  preview_table <- reactive({
        filelist <- input$files1
        if (is.null(input$files1)) return()
        datanames <- filelist$datapath
        bead_data <- lapply(datanames, read_csv)
        names(bead_data) <- basename(filelist$name)
        merge_data <- rbindlist(bead_data, idcol = "ID",use.names = T)
        merge_data <- as.data.frame(merge_data)
        #Set columns to numeric
        ##Numeric
        colnames(merge_data)[colnames(merge_data)=="Profiles"] <- "Measurements"
        merge_data$Raw <- as.numeric(merge_data$Raw)
        merge_data$SD <- as.numeric(merge_data$SD)
        merge_data$Measurements <- as.numeric(merge_data$Measurements)
        merge_data$Corrected <- as.numeric(merge_data$Corrected)
        colnames(merge_data)[colnames(merge_data)=="#"] <- "Bead_Pos"
        
        #Calculate outliers and display warnings
        merge_data$SD_Exclude <-  FALSE
        merge_data$Neg_Exclude <- FALSE
        merge_data$Only_pos <- FALSE
        merge_data$SD_Exclude[merge_data$SD>merge_data$Raw*0.5] <-  TRUE
        merge_data$Neg_Exclude[merge_data$Corrected<0] <- TRUE
        merge_data$Only_pos[merge_data$Raw<0] <-  TRUE
        
        #Split label and make factors
        merge_data[c("User", "Condition","Experiment", "Image", "Channel", "Date")] <- str_split_fixed(str_remove(merge_data$ID, ".csv"), "_", 6)
        #Set Factors
        ##Factor
        merge_data$User <- factor(merge_data$User)
        merge_data$Bead_Pos <- factor(merge_data$Bead_Pos)
        merge_data$Condition <- factor(merge_data$Condition)
        merge_data$Image <- factor(merge_data$Image)
        merge_data$Experiment <- factor(merge_data$Experiment)
        merge_data$Channel <- factor(merge_data$Channel)
        merge_data
  })
  
  output$table1 <- renderDataTable({preview_table()})#Code line renders a datatable using DT package

  SD_Warn <- reactive({
    if (is.null(input$files1)) return()
    merge_data <- preview_table()
    paste("Percentage of measurements that have an SD that exceeds half of the mean corrected value", 
          signif(nrow(merge_data[merge_data$SD_Exclude==TRUE,])/nrow(merge_data[merge_data$SD_Exclude==FALSE,]), 3)*100, "%",
          "(", nrow(merge_data[merge_data$SD_Exclude==TRUE,]), "/", nrow(merge_data[merge_data$SD_Exclude==FALSE,]), ")")
  })
  Neg_Warn <- reactive({
    if (is.null(input$files1)) return()
    merge_data <- preview_table()
    paste("Percentage of measurements that have a negative mean corrected value", 
          signif(nrow(merge_data[merge_data$Neg_Exclude==TRUE,])/nrow(merge_data[merge_data$Neg_Exclude==FALSE,]), 3)*100, "%",
          "(", nrow(merge_data[merge_data$Neg_Exclude==TRUE,]), "/", nrow(merge_data[merge_data$Neg_Exclude==FALSE,]), ")"
    ) })
  output$SD_Warn <- renderText({SD_Warn()})
  output$Neg_Warn <- renderText({Neg_Warn()})
  output$duplicates <- renderTable({
    if (is.null(input$files1)) return("No files have been submitted")
    merge_data <- preview_table()
    table <- subset(merge_data[duplicated(paste(merge_data$ID, merge_data$Bead_Pos)),])
    if (nrow(table)==0) return()
    table
  })
  output$channels <- renderText({paste("There are", length(unique(final_data()$Channel)), "in this dataset.")})
  #Remove outliers to make final data
  final_data <- reactive({
    if (is.null(input$files1)) return()
    final_data <- preview_table()
    final_data <- if (input$SD_Exclude==T) subset(final_data, SD_Exclude==F) else return(final_data)
    final_data <- if (input$Background_Exclude==T) subset(final_data, Neg_Exclude==F) else return(final_data)
    final_data <- if (input$Neg_Exclude==T) subset(final_data, Only_pos==F) else return(final_data)
    return(final_data)
  })
  #How many measurements have been removed
  output$merge_final <- renderText({
    if (is.null(input$files1)) return()
    paste((nrow(preview_table())-nrow(final_data())), "measurements have been excluded from a total of", nrow(preview_table()))
  })
  
  #Statistical test
  quant_pairs <- reactive({
    if (is.null(input$files1)) return()
    if (length(unique(final_data()$Channel))==1) {
      quant_lme <- lme(Raw ~ Condition, random= ~1 |Experiment/ID, final_data())
      quant_emm <- emmeans(quant_lme, "Condition")
      quant_pairs <- data.frame(pairs(quant_emm, adjust="tukey"))
      quant_pairs <- quant_pairs[order(quant_pairs$p.value),]
      quant_pairs$channel <- paste(unique(final_data()$Channel))
      return(data.frame(quant_pairs)) }
    if (length(unique(final_data()$Channel))>1){
        split_data <- split(final_data(), f = final_data()$Channel)
        quant_pair_function <- function(x){
          quant_lme <- lme(Raw ~ Condition, random= ~1 |Experiment/ID, x)
          quant_emm <- emmeans(quant_lme, "Condition")
          quant_pairs <- data.frame(pairs(quant_emm, adjust="tukey"))
          quant_pairs <- quant_pairs[order(quant_pairs$p.value),]
          quant_pairs$Channel <- paste(unique(x$Channel))
          return(quant_pairs)
        }
        quant_pairs <- lapply(split_data, FUN = quant_pair_function)
        quant_pairs <- data.frame(rbindlist(quant_pairs))
        quant_pairs <- quant_pairs[order(quant_pairs$p.value),]
        return(as.data.frame(quant_pairs))
      }
  })
  output$quant_pairs <- renderTable({quant_pairs()})
  
  #Plot of Raw Data
  raw_plot <- reactive({
    if (is.null(input$files1)) return()
    final_data <- final_data()
    quant_pairs <- quant_pairs()
    quant_pairs[c("group1", "group2")] <- gsub(str_split_fixed(string = quant_pairs$contrast, " - ", 2), pattern = "[()]", replacement = "")
    max_fi <- aggregate(final_data$Raw, list(final_data$Condition), FUN=max)
    quant_pairs <- merge.data.frame(quant_pairs, max_fi, by.x = "group1", by.y = "Group.1");colnames(quant_pairs)[colnames(quant_pairs) == "x"] <- "group1_max"
    quant_pairs <- merge.data.frame(quant_pairs, max_fi, by.x = "group2", by.y = "Group.1");colnames(quant_pairs)[colnames(quant_pairs) == "x"] <- "group2_max"
    quant_pairs$max <- rowMaxs(as.matrix(quant_pairs[,c("group1_max", "group2_max")]))
    #Note that the Bonferroni adjustment is indeed an adjustment to the regular P values upv. 
    #But the Tukey-adjusted P values are an entirely separate calculation based on the t ratios but not on upv. 
    #So there is no "default" set of P values underlying the Tukey-adjusted P values.
    
    #Plot of raw values
    wsize <- max(as.numeric(final_data$Experiment))*0.1
    plot1 <- ggplot(final_data, aes(x=Condition, y=Raw))+
      geom_boxplot(width=wsize, size=2.5, color="black", outlier.alpha =0, alpha=0)+
      geom_point(position=position_dodge(width=wsize),shape=21, aes(fill=Experiment), size=2, alpha=0.5)+
      theme_bw()+
      facet_wrap(~Channel, scales = "free_x")+
      stat_pvalue_manual(data = subset(quant_pairs, p.value<0.05), label = "FDR = {signif(p.value, 3)} to 3.s.f.", y.position = subset(quant_pairs, p.value<0.05)$max*1.1, step.increase = 0.1)+
      labs(title = paste0("Quantification of Data Analysis (Grouped) (Raw Measurements): ", date()), y="Raw units of fluoresence", subtitle ="FDR: False discovery rate (p-adj, Tukey Adjusted)   |   3 s.f.: 3 significant figures (rounding)   |   Non-significant comparisons removed")
    plot1
  })
    output$raw_plot <- renderPlot({raw_plot()})
  
    #Plot a Histogram to shows sample spread
    raw_histo <- reactive({
      if (is.null(input$files1)) return()
      final_data <- final_data()
      plot1 <- ggplot(final_data, aes(x=Raw, fill=Experiment))+geom_histogram(alpha=0.5, position = "identity", bins = 20)+facet_grid(cols=vars(Channel), rows = vars(Condition), switch = "y")+
        theme(strip.text.y.left = element_text(angle = 0))+
        labs(title="Histogram of raw data sperated by experiment replicate", x="Raw Fluoresence Intensity", fill="Experiment/ Replicate No.")
      plot1
    }) 
        output$raw_histo <- renderPlot({raw_histo()})
        
    #Calculate corrected/Adjusted dataset
        corr_table <- eventReactive(eventExpr = input$adjust, {
          if (is.null(input$files1)) return()
          final_data <- final_data()
          #Corrected Data
          if(input$neg_control != "") {
            neg_control <- input$neg_control
            neg_data <- subset(final_data, Condition==neg_control)
            neg_average <- aggregate(neg_data$Raw, list(neg_data$Experiment), FUN=mean);rownames(neg_average) <- neg_average$Group.1
            final_data$contr_corrected <- final_data$Raw-neg_average[final_data$Experiment, "x"]
          }
          if(input$pos_control!="") { 
            if(input$neg_control=="") {
              pos_control <- input$pos_control
              pos_data <- subset(final_data, Condition==pos_control)
              pos_average <- aggregate(pos_data$Raw, list(pos_data$Experiment), FUN=mean);rownames(pos_average) <- pos_average$Group.1
              final_data$contr_corrected <-  final_data$Raw/pos_average[final_data$Experiment, "x"]
            } else {
              pos_control <- input$pos_control
              pos_data <- subset(final_data, Condition==pos_control)
              pos_average <- aggregate(pos_data$contr_corrected, list(pos_data$Experiment), FUN=mean);rownames(pos_average) <- pos_average$Group.1
              final_data$contr_corrected <-  final_data$neg_corrected/pos_average[final_data$Experiment, "x"]       
            }
          }
          final_data
        })   
  #Plot Corrected Graph
        cor_plot <- eventReactive(eventExpr = input$adjust, {
          if (is.null(input$files1)) return()
          if (input$pos_control=="" & input$neg_control=="") return()
          corr_data <- corr_table()
          quant_pairs <- quant_pairs()
          quant_pairs[c("group1", "group2")] <- gsub(str_split_fixed(string = quant_pairs$contrast, " - ", 2), pattern = "[()]", replacement = "")
          max_fi <- aggregate(corr_data$contr_corrected, list(corr_data$Condition), FUN=max)
          quant_pairs <- merge.data.frame(quant_pairs, max_fi, by.x = "group1", by.y = "Group.1");colnames(quant_pairs)[colnames(quant_pairs) == "x"] <- "group1_max"
          quant_pairs <- merge.data.frame(quant_pairs, max_fi, by.x = "group2", by.y = "Group.1");colnames(quant_pairs)[colnames(quant_pairs) == "x"] <- "group2_max"
          quant_pairs$max <- rowMaxs(as.matrix(quant_pairs[,c("group1_max", "group2_max")]))
          #Note that the Bonferroni adjustment is indeed an adjustment to the regular P values upv. 
          #But the Tukey-adjusted P values are an entirely separate calculation based on the t ratios but not on upv. 
          #So there is no "default" set of P values underlying the Tukey-adjusted P values.
          wsize <- max(as.numeric(corr_data$Experiment))*0.1
          plot1 <- ggplot(corr_data, aes(x=Condition, y=contr_corrected))+
            geom_boxplot(width=wsize, size=2.5, color="black", outlier.alpha =0, alpha=0)+
            geom_point(position=position_dodge(width=wsize),shape=21, aes(fill=Experiment), size=2, alpha=0.5)+
            theme_bw()+
            facet_wrap(~Channel)+
            stat_pvalue_manual(data = subset(quant_pairs, p.value<0.05), label = "FDR = {signif(p.value, 3)} to 3.s.f.", y.position = subset(quant_pairs, p.value<0.05)$max*1.1, step.increase = 0.1)+
            labs(title = paste0("Quantification of Data Analysis (Grouped): ", date()), y="Corrected units of fluoresence", subtitle ="FDR: False discovery rate (p-adj, Tukey Adjusted)   |   3 s.f.: 3 significant figures (rounding)   |   Non-significant comparisons removed")
          plot1
        })
        output$cor_plot <- renderPlot({cor_plot()})
#Plot Corrected Histogram
        corr_histo <- eventReactive(eventExpr = input$adjust, {
          if (is.null(input$files1)) return()
          if (input$pos_control=="" & input$neg_control =="") return()
          corr_table <- corr_table()
          plot1 <- ggplot(corr_table, aes(x=contr_corrected, fill=Experiment))+geom_histogram(alpha=0.5, position = "identity", bins = 20)+facet_grid(cols=vars(Channel), rows = vars(Condition), switch = "y")+
            theme(strip.text.y.left = element_text(angle = 0))+
            labs(title="Histogram of corrected data sperated by experiment replicate", x="Raw Fluoresence Intensity", fill="Experiment/ Replicate No.")
          plot1
        }) 
        output$corr_histo <- renderPlot({corr_histo()})
  #Save files
  output$grouped_plot_raw <- downloadHandler(filename = function() {paste("grouped_plot_raw_", final_data()[1,]$User, Sys.Date(), ".pdf", sep = "")}, content = function(file) {ggsave(filename = file, plot = raw_plot(), device = "pdf", units = "in", width = 8, height = 6)})
  output$grouped_plot_norm <- downloadHandler(filename = function() {paste("grouped_plot_cor_", final_data()[1,]$User, Sys.Date(), ".pdf", sep = "")}, content = function(file) {ggsave(filename = file, plot = cor_plot(), device = "pdf", units = "in", width = 8, height = 6)})
  output$histo_raw <- downloadHandler(filename = function() {paste("raw_histo_", final_data()[1,]$User, Sys.Date(), ".pdf", sep = "")}, content = function(file) {ggsave(filename = file, plot = raw_histo(), device = "pdf", units = "in", width = 8, height = 18)})
  output$histo_corr <- downloadHandler(filename = function() {paste("corr_histo_", final_data()[1,]$User, Sys.Date(), ".pdf", sep = "")}, content = function(file) {ggsave(filename = file, plot = corr_histo(), device = "pdf", units = "in", width = 8, height = 18)})
  output$final_data <- downloadHandler(filename = function() {paste("final_data_", final_data()[1,]$User, Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(final_data(), file)})
  output$stat_test <- downloadHandler(filename = function() {paste("statistics_", final_data()[1,]$User, Sys.Date(), ".csv", sep = "")}, content = function(file) {write.csv(quant_pairs(), file)})
  
}


# Run the application 
shinyApp(ui = ui, server = server)
