#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)
library(vroom)
library(caret)
library(ggplot2)
library(parallel)
library(doParallel)
library(dplyr)
library(magrittr)
library(data.table)
library(MLeval)
 
#Helpful functions ------
`%notin%` <- Negate(`%in%`)
`%notlike%` <- Negate(`%like%`)
back_2_front <- function(df){
    df <- df[,c(ncol(df),1:(ncol(df)-1))]
    return(df)
}

# Define UI -----
ui <- fluidPage( theme = bslib::bs_theme(bootswatch = "united"),
                 tabsetPanel(
    # Application title
    tabPanel("Import File", 
             sidebarPanel(
               column(width=12),
                 # Importing file
                 fileInput('data', 'Import file containing response-features matrix',accept = c(".csv", ".tsv")),
                 # Return table
                 numericInput('n', "Number of rows of file to show", value= 10, min=1, step=1),
                 numericInput('fs', 'Number of features used', value= 10, min=1, step=1)),
             mainPanel(
               tableOutput("head_data"))),
    
    tabPanel("Exploratory Data Analysis",
             h2("Histogram of DAS44 Outcomes"),
             h2("PCA plot"),
             h2("Volcano Plot"),
             h2("Baseline Correlations with DAS44 Outcomes")
             ),
    tabPanel("Model Training",
            
                 #h1("Training data"),
               verticalLayout(
                 numericInput('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", value= 10, min=1, step=1),
                 numericInput('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", value= 10, min=1, step=1),
                 
                 h2("Feature Selection"),
                 plotOutput("gg_fs"),
                 h2("ROC Curves"),
                 plotOutput("model_selection")
               ),
               h2("Performance Metrics"),
               tableOutput("model_sel_tbl")),
                
    tabPanel("Model Testing",
             sidebarPanel(
                 selectInput('model_from_train', 'Select ML algorithm', 
                             c('svmRadial', 'rf', 'xgbTree', 'knn', 'naive_bayes'))),
             mainPanel(
                 h1("Evaluation of Final Model"),
                 #tableOutput("head_test"),
                 verticalLayout(
                   h2("ROC Curve"),
                   plotOutput("model_perf_roc"),
                   h2("Calibration Curve"),
                   plotOutput("model_perf_cc"))))))

# Define server ----- 
server <- function(input, output, session) ({
  thematic::thematic_shiny()
    load_file <- function(name, path) {
        ext <- tools::file_ext(name)
        switch(ext,
               csv = vroom::vroom(path, delim = ","),
               tsv = vroom::vroom(path, delim = "\t"),
               validate("Invalid file; Please upload a .csv or .tsv file")
        )
    }
    
    data_df <- reactive({
      req(input$data)
      load_file(input$data$name, input$data$datapath)
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
    
    
    output$head_data <- renderTable({ # not included in the output right now
        head(data_df(), input$n)
    })
    
    output$head_test <- renderTable({ # not included in the output right now
      data_test()
    })
    
    output$gg <- renderPlot({
        ggplot(data_train())+
            geom_bar(aes(x=Response), fill= c('blue', 'red'))+
            theme_minimal()
    })
    
    ft_sel <- reactive({
        options(warn=-1)
        lmProfile <- rfe(x=data_train()[,-c(1:5)], y=as.factor(data_train()$Response),
                         sizes = c(1:5, 10, 15, 20),
                         rfeControl = rfeControl(functions = rfFuncs,
                                                 method = "repeatedcv",
                                                 repeats = 50,
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
    
    output$feats <- renderTable({
        head(ft_sel(), input$fs)
    })
    
    output$gg_fs <- renderPlot({
        ab <- head(ft_sel(), input$fs)
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
    
    output$DF_mod <- renderTable({
        head(model_df(),5)
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
    
    output$model_sel_tbl <- renderTable({ 
        model_sel_tb()
    })
    
    ## Model generation
    model_gen <- reactive({
      set.seed(42)
      registerDoParallel(4)
      getDoParWorkers()
      modell <- train(Response~., data=resp_fs(),
                      method=input$model_from_train, 
                      metric='Accuracy',
                      trControl= trainControl(method="repeatedcv", 
                                              number=input$repeatedcv_number, 
                                              repeats=input$mod_gen_repeats,
                                              summaryFunction=twoClassSummary,
                                              #summaryFunction = prSummary, # only needed with imbalanced classification
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
      model_eval()[1]
    })
    
    output$model_perf_cc <- renderPlot({ 
      model_eval()[4]
    })
    

})

shinyApp(ui = ui, server = server)
