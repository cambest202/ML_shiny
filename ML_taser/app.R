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
ui <- fluidPage(
    
    # Application title
    titlePanel("Predicting 3 Month DAS44 Outcomes Using Baseline Metabolome"),
    
    sidebarPanel(
        # Importing file
        fileInput('train', 'Import training file',accept = c(".csv", ".tsv")),
        fileInput('test', 'Import test file',accept = c(".csv", ".tsv")),
        
        # Return table
        numericInput('n', 'Rows', value= 10, min=1, step=1),
        
        numericInput('fs', 'Features', value= 10, min=1, step=1)
        
    ),
    
    mainPanel(
        h1("Training data"),
        #tableOutput("head_train"),
        #plotOutput("gg"),
        #tableOutput("feats"),
        #textOutput("resp_length"), 
        #textOutput("samp_count"),
        #tableOutput("DF_mod"),
        h2("Feature Selection: RFE"),
        plotOutput("gg_fs"),
        h2("Model Selection: ROC Curves"),
        plotOutput("model_selection"),
        h2("Model Selection: Performance Metrics"),
        tableOutput("model_sel_tbl"),
        #tableOutput("resp_fs_tb")
        #tableOutput("head_train_only"),
        #tableOutput("resp_only")
    )
)


# Define server ----- 
server <- function(input, output, session) ({
    load_file <- function(name, path) {
        ext <- tools::file_ext(name)
        switch(ext,
               csv = vroom::vroom(path, delim = ","),
               tsv = vroom::vroom(path, delim = "\t"),
               validate("Invalid file; Please upload a .csv or .tsv file")
        )
    }
    
    data_train <- reactive({
        req(input$train)
        load_file(input$train$name, input$train$datapath)
    })
    
    output$head_train <- renderTable({
        head(data_train(), input$n)
    })
    
    output$gg <- renderPlot({
        ggplot(data_train())+
            geom_bar(aes(x=Response), fill= c('blue', 'red'))+
            theme_minimal()
    })
    
    train_no_resp <- reactive({
        as.data.frame(data_train()[,-c(1:5)])
    })
    
    response <- reactive({
        data_train()$Response
    })
    
    
    output$head_train_only <- renderTable({
        train_no_resp()
    })
    
    output$resp_only <- renderTable({
        class(response())
    })
    
    resp_count <- reactive({
        length(response())
    })
    
    sample_count <- reactive({
        dim(train_no_resp())
    })
    
    output$resp_length <- renderText({
        resp_count()
    })
    
    output$samp_count <- renderText({
        sample_count()
    })
    
    
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
    
    output$feats <- renderTable({
        head(ft_sel(), input$fs)
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
        modelSvm <- train(resp_fs()$Response~., data=resp_fs(), method="svmRadial", trControl=control)
        # train the RF model
        set.seed(42)
        model_rf <- train(resp_fs()$Response~., data=resp_fs(), method="rf", trControl=control)
        # train the XGBoost model
        set.seed(42)
        model_xgb <- train(resp_fs()$Response~., data=resp_fs(), method="xgbTree", trControl=control)
        # train the KNN model
        set.seed(42)
        model_knn <- train(resp_fs()$Response~., data=resp_fs(), method="knn", trControl=control)
        # train the naive Bayes model
        set.seed(42)
        model_naivebayes <- train(resp_fs()$Response~., data=resp_fs(), method="naive_bayes", trControl=control)
        
        comp_roc <- evalm(list(modelSvm, model_rf, 
                               model_xgb,  model_knn, 
                               model_naivebayes),
                          gnames=c('SVM', 'RF', 'XGB', 'KNN', 'NB'))
        
    }) 
    
    output$model_selection <- renderPlot({ # plot ROC for all models from training data
        model_sel()
    })
    
    model_sel_tb <- reactive({
        ml_eval_output <- as.data.frame(model_sel()$stdres)
        ml_eval_output$Measure <- rownames(ml_eval_output)
        ml_eval_output <- back_2_front(ml_eval_output)
        ml_eval_output[c(1:3,8:13),]
    })
    
    output$model_sel_tbl <- renderTable({ # plot ROC for all models from training data
        model_sel_tb()
    })
    
    
})

# Run the application 
shinyApp(ui = ui, server = server)
