

library(shiny)
library(DT)
library(vroom)
library(caret)
library(ggplot2)
library(parallel)
library(doParallel)
library(dplyr)
library(plyr)
library(magrittr)
library(data.table)
library(MLeval)
library(limma)
library(purrr)
library(ggpubr)
library (gtools)
library(tidyverse) 
library(ggrepel)
library(tidyr)
library(e1071)
library(randomForest)
library(ggplot2)
library(tidyverse) 
library (readr)
library(PMCMRplus)
library (PMCMR)
library(ggpubr)
library (mosaic)
library (dplyr)
library (data.table)
library(reshape2)
library (gtools)
library(plyr)
library(limma)
library(ggrepel)
library(amap)
library(rstatix)
library(broom)
library(ggprism)
library(HDDesign)
library(caret)
library(rsample)
library(sandwich)
library(rpart)
library(rpart.plot)
library(randomForest)
library(RColorBrewer)
library(plotly)
library(purrr)
library(e1071)
library(ggraph)
library(igraph)
library(pscl)
library(parallel)
library(doParallel)
library(ROCR) 
library(qvalue)
library(corrr)
library(ggcorrplot)
library(ape)
library(forcats)
library(kernlab)
library(xgboost)
library(mlbench)
library(naivebayes)

#Helpful functions ------
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)
back_2_front <- function(df){
  df <- df[,c(ncol(df),1:(ncol(df)-1))]
  return(df)
}

# example file
met_file <- read.csv('20210507_AB_metabolites_baseline_matrix.csv')
comp_file <- read.csv('20210507_AB_composite_baseline_matrix.csv')

# deploy app to github using: rsconnect::deployApp('ML_taser/')
#options(repos = BiocManager::repositories()) 
thematic::thematic_shiny(font = "auto")

# Define UI -----
ui <- fluidPage( navbarPage(theme = bslib::bs_theme(bootswatch = "sandstone"),#sandstone, cerulean
                            
                            HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" class="active" href="#">Metabolic signature of MTX-responses in RA patients</a>'), id="nav",
                            windowTitle = "Metabolic signature of MTX-responses in RA patients"),
                 #h1("Metabolic signature of MTX-responses in RA patients"),
                 #tags$style(HTML("h1{color: navy;}")),
                 tabsetPanel(
                   tabPanel("Import File", 
                            sidebarPanel(
                              column(width=12),
                              # Importing file
                              tags$hr(),
                              h2('First Model Data'),
                              radioButtons('file_met', 'Select File', 'Baseline TaSER Metabolome File (Metabolites-Only)'),
                              tags$hr(),
                              h2('Second Model Data'),
                              radioButtons('file_comp', 'Select File', 'Baseline TaSER Metabolome File (Metabolites + Clinical Metadata'))),
                   
                   tabPanel("Exploratory Data Analysis",
                            h2("Histogram of DAS44 Outcomes*"),
                            h5('*vertical line denotes the EULAR-recommended cut-off for good/poor responses'),
                            tags$hr(),
                            h5("Patient responses based on change in DAS44 (ΔDAS44) and final DAS44 score after 3 months"),
                            h3('ΔDAS44'),
                            plotOutput("das44_hist"),
                            h3('DAS44 after 3 months'), 
                            plotOutput("das44_hist_end"),
                            tags$hr(),
                            h2("Class Balance: DAS44-Based Response Categories"),
                            plotOutput("class_balance"),
                            tags$hr(),
                            h2("PCA plot"),
                            plotOutput("PCA"),
                            tags$hr(),
                            h2("Volcano Plot: Differential Abundance of Metabolites Across DAS44-Based Responses"),
                            plotOutput("limma"),
                            tableOutput("limma_tab"),
                            
                   ),
                   tabPanel("Metabolites-Only Model Training",
                            
                            #h1("Training data"),
                            verticalLayout(
                              radioButtons('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", 10),
                              radioButtons('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", 50),
                              # numericInput('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", value= 10, min=1, step=1),
                              # numericInput('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", value= 100, min=1, step=1),
                              # numericInput('fs', 'Number of features used', value= 10, min=1, step=1),
                              h3('Please be patient, model training can take time...'),
                              tags$hr(),
                              h2("Feature Selection"),
                              img(src='20210516_met_FS.png'),
                              #plotOutput("gg_fs"),
                              tags$hr(),
                              h2("ROC Curves"),
                              img(src='20210516_met_mod_selection.png')
                              #plotOutput("model_selection")
                            ),
                            tags$hr(),
                            h2("Performance Metrics"),
                            img(src='20210516_met_table_perf.png')),
                            
                            #tableOutput("model_sel_tbl")),
                   
                   tabPanel("Metabolites-Only Model Testing",
                            sidebarPanel(
                              radioButtons('model', 'Select ML Algorithm', 'Naive Bayes')),
                              
                              #selectInput('model_from_train', 'Select ML algorithm', c('svmRadial', 'rf', 'xgbTree', 'knn', 'naive_bayes'))),
                            mainPanel(
                              h1("Evaluation of Final Metabolites-Only Model"),
                              tags$hr(),
                              #tableOutput("head_test"),
                              verticalLayout(
                                h2("ROC Curve"),
                                img(src='20210516_met_roc.png'),
                                tags$hr(),
                                #plotOutput("model_perf_roc"),
                                h2("Calibration Curve"),
                                img(src='20210516_met_cc.png'),
                                tags$hr(),
                                h2('Final model metrics'),
                                img(src='20210516_met_fin_perf.png')))),
                   
                   tabPanel("Composite Model Training",
                            
                            #h1("Training data"),
                            verticalLayout(
                              radioButtons('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", 10),
                              radioButtons('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", 50),
                              # numericInput('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", value= 10, min=1, step=1),
                              # numericInput('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", value= 100, min=1, step=1),
                              # numericInput('fs', 'Number of features used', value= 10, min=1, step=1),
                              h3('Please be patient, model training can take time...'),
                              tags$hr(),
                              h2("Feature Selection"),
                              img(src='20210516_comp_fs.png'),
                              tags$hr(),
                              #plotOutput("gg_fs"),
                              h2("ROC Curves"),
                              img(src='20210516_comp_mod_selection.png')
                              #plotOutput("model_selection")
                            ),
                            tags$hr(),
                            h2("Performance Metrics"),
                            img(src='20210516_comp_tab_perf.png')),
                   
                   tabPanel("Composite Model Testing",
                            sidebarPanel(
                              radioButtons('model', 'Select ML Algorithm', 'Naive Bayes')),
                            
                            #selectInput('model_from_train', 'Select ML algorithm', c('svmRadial', 'rf', 'xgbTree', 'knn', 'naive_bayes'))),
                            mainPanel(
                              h1("Evaluation of Final Composite Model"),
                              #tableOutput("head_test"),
                              verticalLayout(
                                tags$hr(),
                                h2("ROC Curve"),
                                img(src='20210516_comp_roc.png'),
                                tags$hr(),
                                #plotOutput("model_perf_roc"),
                                h2("Calibration Curve"),
                                img(src='20210516_comp_cc.png'),
                                tags$hr(),
                                h2('Final model metrics'),
                                img(src='20210516_comp_fin_perf.png')))),
                   
                   tabPanel("Comparison of Models",
                            
                            mainPanel(
                              h1("Comparing Predictive Performance of the Models"),
                              #tableOutput("head_test"),
                              verticalLayout(
                                tags$hr(),
                                h2("ROC Curves"),
                                img(src='20210516_comp_met_comparison.png'),
                                h2("Are the ROC curves statistically different?"),
                                img(src="20210516_comparison_delongs.png"),
                                tags$hr()))),
                   
                   tabPanel("Data",
                            h2("TaSER Metabolomic Data"),
                            h3('(Annotated metabolites included only)'),
                            mainPanel(
                              h2("Metabolites-Only"),
                              tableOutput("head_data"),
                              tags$hr(),
                              h2("Composite"),
                              tableOutput("head_data_comp")))))

# Define server ----- 
server <- function(input, output, session) ({
  thematic::thematic_shiny()
  load_file <- function(name, path) {
    ext <- tools::file_ext(name)
    switch(ext,
           csv = met_file,
           #tsv = vroom::vroom(path, delim = "\t"),
           validate("Invalid file; Please upload a .csv or .tsv file")
    )
  }
  
  
  data_df <- reactive({
    met_file
    })
  
  data_composite <- reactive({
    comp_file
  })
  
  output$head_data_comp <- renderTable({ # not included in the output right now
    head(data_composite(), 10)
  })
  
  data_train <- reactive({
    set.seed(42)
    index <- createDataPartition(data_df()$Response, p = 0.7, list = FALSE)
    train_data <- data_df()[index, ]
  })
  
  data_test <- reactive({
    set.seed(42)
    index <- createDataPartition(data_df()$Response, p = 0.7, list = FALSE)
    test_data <- data_df()[-index,]
  })
  
  output$das44_hist <- renderPlot({
    ggplot(data_df())+
      geom_histogram(aes(x=DAS44),
                     colour= 'black', fill='#56B4E9')+
      theme_minimal()+
      labs(x='ΔDAS44',
           y='Frequency',
           title='Distribution of DAS44 Changes Across 3 Months in TaSER Trial')+
      geom_vline(xintercept = -1.2, size=3, colour='red')
  })
  
  output$das44_hist_end <- renderPlot({
    ggplot(data_df())+
      geom_histogram(aes(x=End_DAS44),
                     colour= 'black', fill='#56B4E9')+
      theme_minimal()+
      labs(x='3-Months: DAS44',
           y='Frequency',
           title='Distribution of DAS44 Changes Across 3 Months in TaSER Trial')+
      geom_vline(xintercept = 2.4, size=3, colour='red')
  })
  
  
  
  PCA_data <- reactive({
    comb_ft_PCA <- data_df()[,-c(1,3:5)]
    comb_cut_t <- comb_ft_PCA[,-1]
    scaled_intensities <- scale((comb_cut_t))
    scaled_intensities[do.call(cbind, lapply(scaled_intensities, is.nan))] <- 0
    scaled_intensities<- as.data.frame(scaled_intensities)
    pca_data <- prcomp(scaled_intensities)
    pca_data
  })
  
  output$PCA <- renderPlot({
    pca_coord <- data.frame(PCA_data()$x)
    var_explained <- PCA_data()$sdev^2/sum(PCA_data()$sdev^2)
    var_explained[1:5]
    pca_coord$Response <-as.factor(data_df()$Response)
    ggplot(pca_coord) + 
      geom_point(size=2, alpha=0.7, 
                 aes(x=PC1,y=PC2, colour= Response, fill= Response))+
      labs(x=paste0("PC1: ",round(var_explained[1]*100,1),"%"),
           y=paste0("PC2: ",round(var_explained[2]*100,1),"%"))+
      geom_hline(yintercept = 0,
                 colour='navy',
                 linetype='dashed')+
      geom_vline(xintercept = 0,
                 colour='navy',
                 linetype='dashed')+
      theme_minimal()
  })
  
  output$head_data <- renderTable({ # not included in the output right now
    head(data_df(), 10)
  })
  
 
  
  
  output$head_test <- renderTable({ # not included in the output right now
    data_test()
  })
  
  output$class_balance <- renderPlot({
    ggplot(data_train())+
      geom_bar(aes(x=Response), fill= c('#E69F00', '#56B4E9'))+
      theme_minimal()
  })
  
  
  # Differential abundance
  volcano_df <- reactive({
    volc <- data_df()[,-c(1:5)]
    volc$Response <- as.factor(data_df()$Response)
    volc <- back_2_front(volc)
    rownames(volc) <- paste0(volc$Response, '_', rownames(volc))
    volc
  })
  
  limma_fun <- reactive({
    df1 <- volcano_df()
    df1 <- df1[,-1]
    df1 <- as.data.frame(t(df1))
    colnames(df1) <- volcano_df()$Response
    df1
    
    Group <- factor(colnames(df1), levels = c('Good', 'Poor'))
    design <- model.matrix (~Group)
    colnames(design) <- c('Good', 'GoodvsPoor')
    eset <- df1
    fit <- lmFit(eset, design)
    fit <- eBayes(fit)
    toptable <- topTable(fit, coef = 'GoodvsPoor', adjust = 'BH', number = 220)
    toptable <- as.data.frame(toptable)
    toptable$Putative_Metabolite <- rownames(toptable)
    toptable <- toptable[,c(ncol(toptable),1:(ncol(toptable)-1))]
    toptable$Sig <- 0
    toptable$Sig <- ifelse(toptable$adj.P.Val <0.05, '< 0.05', '> 0.05') 
    toptable$Sig_Names <-0
    toptable$Sig_Names <- ifelse(toptable$Sig =='< 0.05' ,toptable$Putative_Metabolite, '')
    toptable
    
  })
  
  output$limma_tab <- renderTable({
    head(limma_fun(), n=10)
  })
  
  output$limma <- renderPlot({
    ggplot(limma_fun(), aes(x=logFC, y=-log10(P.Value), 
                            colour=Sig, 
                            group=Sig)) +
      geom_point (alpha=0.7) +
      theme_minimal() +
      labs (x='LogFC',
            y='-Log p-value',
            colour='Adjusted \np-value')+
      geom_text_repel(aes(x = logFC, y = -log10(P.Value), label = Sig_Names),
                      box.padding =1,
                      size=2.5,
                      max.overlaps = Inf,
                      position = position_jitter(seed = 1),
                      arrow = arrow(length = unit(0.0015, "npc"))) +  
      theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
            axis.title = element_text(size = rel(1.25)))+
      scale_color_brewer(palette = "Set1",direction=-1)        })
  
  
  
  
  #Machine learning
  ft_sel <- reactive({

    options(warn=-1)
    lmProfile <- rfe(x=data_train()[,-c(1:5)], y=as.factor(data_train()$Response),
                     sizes = c(1:5, 10, 15, 20),
                     rfeControl = rfeControl(functions = rfFuncs,
                                             method = "repeatedcv",
                                             #repeats = 50,
                                             verbose = FALSE))
    var_imp_RFE <- lmProfile$variables
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$Overall > 0)
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'Galactonic/Gluconic acid')
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'LysoPC')
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'Slfasalazine')
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'Nicotine')
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'linolenic acid')
    var_imp_RFE <- subset(var_imp_RFE, var_imp_RFE$var != 'AMP-ubiquitin')
    
  })
  
  output$gg_fs <- renderPlot({
    ab <- head(ft_sel(), 10)
    ggplot(ab, aes(x=Overall, y=reorder(var, Overall)))+
      geom_col(aes(fill=Overall))+
      theme_minimal()+
      scale_fill_continuous(low='light green', high='navy')+
      theme(legend.position = 'none',
            axis.title.y = element_text(size = 12),
            axis.title.x = element_text(size = 12),
            axis.text.x = element_text(size = 10),
            axis.text.y = element_text(size = 10))+
      labs(x='Relative Importance', 
           y='Putative Metabolite')
  })
  
  model_df <- reactive({
    df <- data_train()[,-c(1:5)]
    df <- as.data.frame(t(df))
    df$Metabolite <- rownames(df)
    df <- back_2_front(df)
    #df$Metabolite <- as.factor(df$Metabolite)
    feat_sel <- as.data.frame(ft_sel())
    feat_sel <- head(feat_sel, 10)
    #feat_sel$var <- as.factor(feat_sel$var)
    df <- subset(df, df$Metabolite %in% feat_sel$var)
    rownames(df) <- paste0(df$Metabolite, '_', rownames(df))
    df <- df[,-1]
    df <- as.data.frame(t(df))
    colnames(df) <- gsub("(.*)_.*","\\1",colnames(df)) 
    df
  })
  

  resp_fs <- reactive({
    ab <- cbind.data.frame(data_train()$Response, model_df())
    names(ab)[1] <- 'Response'
    ab
  })
  
  output$resp_fs_tb <- renderTable({
    resp_fs()
  })
  
  model_sel <- reactive({

    control= trainControl(method="repeatedcv", 
                          number=10, 
                          summaryFunction = twoClassSummary,
                          savePredictions = TRUE, 
                          classProbs = TRUE, 
                          verboseIter = TRUE)
    
    # train the SVM model
    set.seed(42)
    SVM <- train(resp_fs()$Response~., data=resp_fs(), method="svmRadial", trControl=control)
    # train the RF model
    set.seed(42)
    RF <- train(resp_fs()$Response~., data=resp_fs(), method="rf", trControl=control)
    # train the XGBoost model
    set.seed(42)
    XGBoost <- train(resp_fs()$Response~., data=resp_fs(), method="xgbTree", trControl=control)
    # train the KNN model
    set.seed(42)
    KNN <- train(resp_fs()$Response~., data=resp_fs(), method="knn", trControl=control)
    # train the naive Bayes model
    set.seed(42)
    NB <- train(resp_fs()$Response~., data=resp_fs(), method="naive_bayes", trControl=control)
    
    comp_roc <- evalm(list(SVM, RF, XGBoost,  KNN, NB),
                      gnames=c('KNN', 'XGBoost', 'RF', 'SVM', 'NB'))
    
    
  }) 
  
  output$model_selection <- renderPlot({ # plot ROC for all models from training data
    model_sel()[[1]]
  })
  
  model_sel_tb <- reactive({
    ml_eval_output <- as.data.frame(model_sel()$stdres)
    ml_eval_output$Measure <- rownames(ml_eval_output)
    ml_eval_output <- back_2_front(ml_eval_output)
    ml_eval_output[c(1:3,8:13),]
  })
  
  ## Model generation
  model_gen <- reactive({

    set.seed(42)
    modell <- train(Response~., data=resp_fs(),
                    method=input$model_from_train, 
                    metric='Accuracy',
                    trControl= trainControl(method="repeatedcv", 
                                            number=10, 
                                            repeats=10,
                                            summaryFunction=twoClassSummary,
                                            savePredictions = TRUE, 
                                            classProbs = TRUE, 
                                            verboseIter = TRUE))
  })
  
  model_eval <- reactive({
    model_2_assess <- model_gen()
    performance <- evalm(model_2_assess, gnames='Model')
    performance
  })
  
  output$model_perf_roc <- renderPlot({ 
    model_eval()$roc
  })
  
  output$model_perf_cc <- renderPlot({ 
    model_eval()$cc
  })
  
  model_sel_tb <- reactive({
    ml_eval_output <- as.data.frame(model_sel()$stdres)
    ml_eval_output$Measure <- rownames(ml_eval_output)
    ml_eval_output <- back_2_front(ml_eval_output)
    ml_eval_output[c(1:3,8:13),]
  })
  
  output$model_sel_tbl <- renderTable({ 
    model_sel_tb()
  })
  
})

shinyApp(ui = ui, server = server)

