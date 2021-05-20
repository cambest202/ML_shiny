## Shiny app using R script from TaSER analysis
# refine to make dataset-agnostic
# imported file format: Column 1= Class/classes to predict; Columns 2:ncol= features (e.g. metabolites, genes, clinical data)
# additional columns including outcome-related data- e.g. DAS44 measures used to determine patient outcomes- should be included immediately after column 1 and reported in the appropriate input
library(shiny)
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
library (mosaic)
library(reshape2)
library(rsample)
library(sandwich)
library(rpart)
library(rpart.plot)
library(RColorBrewer)
library(ggraph)
library(igraph)
library(pscl)
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
# deploy app to github using: rsconnect::deployApp('ML_taser/')
# check repos: options(repos = BiocManager::repositories()) 

# match plots to shiny theme
thematic::thematic_shiny(font = "auto")

# Define UI -----
ui <- fluidPage( navbarPage(theme = bslib::bs_theme(bootswatch = "sandstone"),
                            
                            HTML('<a style="text-decoration:none;cursor:default;color:#FFFFFF;" 
                                 class="active" href="#">Developing Feature Profile for Sample Classes</a>'), 
                            id="nav",
                            windowTitle = "Developing Feature Profile for Sample Classes"),
                 tabsetPanel(
                     tabPanel("Import File", 
                              mainPanel(
                              column(width=12),
                              # Importing file
                              tags$hr(),
                              h2("Expected file layout*"),
                              h5('*additional measures can be included after Class column and before first Feature columns'),
                              h5('If additional measures included, make sure to note these below'),
                              tableOutput("demo_file"),
                              tags$hr(),
                              fileInput('data', 'Import file containing binary class + features',accept = c(".csv", ".tsv")),
                              tags$hr(),
                              numericInput('n', "Number of rows of file to show", value= 10, min=1, step=1),
                              tags$hr(),
                              selectInput('non_fts', 'Select if file includes additional non-feature, non-Class metrics e.g. Class-associated measures', 
                                          c('Yes', 'No')),
                              tags$hr(),
                              numericInput('col_ft', "Number of additional columns of non-feature data e.g. DAS44 measures", 
                                           value= 5, min=1, step=1),
                              tags$hr(),
                              textInput('class_one', 'Enter name of class one e.g. Good or Poor'),
                              tags$hr(),
                              textInput('class_two', 'Enter name of class two e.g. Good or Poor'),
                              # Horizontal line ----
                              tags$hr())),
                     
                     tabPanel("Data",
                              h2("Imported File Layout"),
                              h3('Check file matches expected layout:'),
                              h5("Column 1: Class, Columns 2:ncol- Features (additional measures can be included between these)"),
                              mainPanel(
                                  tableOutput("head_data"),
                                  tags$hr())),
                     
                     tabPanel("Exploratory Data Analysis",
                              tags$hr(),
                              h2("Class Balance"),
                              plotOutput("class_balance"),
                              # Horizontal line ----
                              tags$hr(),
                              h2("PCA plot"),
                              plotOutput("PCA"),
                              # Horizontal line ----
                              tags$hr(),
                              h2("Volcano Plot: Differential Levels of Features Across Classes*"),
                              h5("*note that features may include putatively identified molecules"),
                              plotOutput("limma"),
                              tags$hr(),
                              tableOutput("limma_tab")),
                     
                     tabPanel("Model Training",
                              
                              #h1("Training data"),
                              verticalLayout(
                                  numericInput('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", value= 10, min=1, step=1),
                                  radioButtons("repeats", "Number of folds in K-fold cross validation", 10),
                                  numericInput('fs', 'Number of features used', value= 10, min=1, step=1),
                                  h3('Please be patient, model training can take time'),
                                  # Horizontal line ----
                                  tags$hr(),
                                  h2("Feature Selection"),
                                  plotOutput("gg_fs"),
                                  # Horizontal line ----
                                  tags$hr(),
                                  h2("ROC Curves"),
                                  plotOutput("model_selection"),
                                  tags$hr()),
                              # Horizontal line ----
                              h2("Performance Metrics"),
                              tableOutput("model_sel_tbl")),
                     
                     
                     tabPanel("Model Testing",
                              sidebarPanel(
                                  selectInput('model_from_train', 'Select ML algorithm', c('Select Algorithm','svmRadial', 'rf', 'xgbTree', 'knn', 'naive_bayes')),
                                  numericInput('mod_gen_repeats', "Number of K-fold cross validations used in resampling: Enter number (e.g. 5-10)", value= 10, min=1, step=1),
                                  numericInput('repeatedcv_number', "Number of folds in K-fold cross-validation (e.g. 50-100)", value= 100, min=1, step=1)),
                              mainPanel(
                                  h1("Evaluation of Final Model"),
                                  #tableOutput("head_test"),
                                  verticalLayout(
                                      h2("Using the following features, the classes from the dataset can be predicted with the ROC-AUC metric providing a measure of the model performance:"),
                                      h2("ROC Curve"),
                                      plotOutput("model_perf_roc"),
                                      # Horizontal line ----
                                      tags$hr(),
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
    
    data_df_2 <- reactive({
        df <- data_df()
        names(df)[1] <- 'Class'
        df
    })
    
    demo_df <- reactive({
        Class <- c('Good', 'Poor')
        feats <- matrix(c(10,20,30,40,50,60,70,80,90,100), byrow = TRUE, nrow = 6, ncol = 10)
        feats <- as.data.frame(feats)                   
        names(feats) <- c(1:10)
        colnames(feats) <- paste0('Feature_', colnames(feats))
        df <- cbind.data.frame(Class, feats)
        df
    })
    
    output$demo_file <- renderTable({
        demo_df()
    })
    
    data_train <- reactive({
        set.seed(42)
        index <- createDataPartition(data_df_2()$Class, p = 0.7, list = FALSE)
        train_data <- data_df_2()[index, ]
    })
    
    data_test <- reactive({
        set.seed(42)
        index <- createDataPartition(data_df_2()$Class, p = 0.7, list = FALSE)
        test_data <- data_df_2()[-index,]
    })
    
    PCA_data <- reactive({
        comb_ft_PCA <- data_df_2()[,-c(2:input$col_ft)]
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
        pca_coord$Class <-as.factor(data_df_2()$Class)
        ggplot(pca_coord) + 
            geom_point(size=2, alpha=0.7, 
                       aes(x=PC1,y=PC2, colour= Class, fill= Class))+
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
    
    output$head_data <- renderTable({ 
        head(data_df_2(), input$n)
    })

    output$class_balance <- renderPlot({
        ggplot(data_train())+
            geom_bar(aes(x=Class), fill= c('#E69F00', '#56B4E9'))+
            theme_minimal()+
            theme(legend.position = 'none',
                  axis.title.y = element_text(size = 14),
                  axis.title.x = element_text(size = 14),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14))
    })
    
    
    # Differential abundance
    volcano_df <- reactive({
        volc <- data_df_2()[,-c(1:input$col_ft)]
        volc$Class <- as.factor(data_df_2()$Class)
        volc <- back_2_front(volc)
        rownames(volc) <- paste0(volc$Class, '_', rownames(volc))
        volc
    })
    
    output$Volc_table <- renderTable({
        head(volcano_df(),10)
    })
    
    limma_df <- reactive({
        df1 <- volcano_df()
        rownames(df1) <- paste(df1$Class,'_',rownames(df1))
        df1 <- df1[,-1]
        df1 <- as.data.frame(t(df1))
        colnames(df1) <- volcano_df()$Class
        df1
    })
    
    limma_fun <- reactive({
        df1 <- limma_df()
        colnames(df1) <- ifelse(colnames(df1) == input$class_one, 'Good', 'Poor')
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
            
            geom_hline(yintercept = -log(0.05), colour='navy', linetype='dashed')+
            theme(plot.title = element_text(size = rel(1.5), hjust = 0.5),
                  axis.title = element_text(size = rel(1.25)))+
            scale_color_brewer(palette = "Set1",direction=-1)        
    })
    
    #Machine learning
    ft_sel <- reactive({
        
        options(warn=-1)
        lmProfile <- rfe(x=data_train()[,-c(1:input$col_ft)], y=as.factor(data_train()$Class),
                         sizes = c(1:5, 10, 15, 20),
                         rfeControl = rfeControl(functions = rfFuncs,
                                                 method = "repeatedcv",
                                                 repeats = 10,
                                                 verbose = FALSE))
        var_imp_RFE <- lmProfile$variables
        var_imp_RFE <- head(var_imp_RFE, input$fs)
        var_imp_RFE
    })
    
    output$gg_fs <- renderPlot({
        ggplot(ft_sel(), aes(x=Overall, y=reorder(var, Overall)))+
            geom_col(aes(fill=Overall))+
            theme_minimal()+
            scale_fill_continuous(low='light green', high='navy')+
            theme(legend.position = 'none',
                  axis.title.y = element_text(size = 14),
                  axis.title.x = element_text(size = 14),
                  axis.text.x = element_text(size = 14),
                  axis.text.y = element_text(size = 14))+
            labs(x='Relative Importance', 
                 y='Putative Metabolite')
    })
    
    model_df <- reactive({
        df <- data_train()[,-c(1:input$col_ft)]
        df <- as.data.frame(t(df))
        df$Feature <- rownames(df)
        df <- back_2_front(df)
        feat_sel <- as.data.frame(ft_sel())
        feat_sel <- head(feat_sel, 10)
        #feat_sel$var <- as.factor(feat_sel$var)
        df <- subset(df, df$Feature %in% feat_sel$var)
        rownames(df) <- paste0(df$Feature, '_', rownames(df))
        df <- df[,-1]
        df <- as.data.frame(t(df))
        colnames(df) <- gsub("(.*)_.*","\\1",colnames(df)) 
        df
    })
    
    
    resp_fs <- reactive({
        ab <- cbind.data.frame(data_train()$Class, model_df())
        names(ab)[1] <- 'Class'
        ab
    })
    
    output$resp_fs_tb <- renderTable({
        resp_fs()
    })
    
    model_sel <- reactive({
        
        control= trainControl(method="repeatedcv", 
                              number=input$mod_gen_repeats, 
                              repeats=10,
                              summaryFunction = twoClassSummary,
                              savePredictions = TRUE, 
                              classProbs = TRUE, 
                              verboseIter = TRUE)
        # train the SVM model
        set.seed(42)
        SVM <- train(resp_fs()$Class~., data=resp_fs(), method="svmRadial", trControl=control)
        # train the RF model
        set.seed(42)
        RF <- train(resp_fs()$Class~., data=resp_fs(), method="rf", trControl=control)
        # train the XGBoost model
        set.seed(42)
        XGBoost <- train(resp_fs()$Class~., data=resp_fs(), method="xgbTree", trControl=control)
        # train the KNN model
        set.seed(42)
        KNN <- train(resp_fs()$Class~., data=resp_fs(), method="knn", trControl=control)
        # train the naive Bayes model
        set.seed(42)
        NB <- train(resp_fs()$Class~., data=resp_fs(), method="naive_bayes", trControl=control)
    
        evalm(list(SVM, RF, XGBoost,  KNN, NB),
              gnames=c('SVM', 'RF', 'XGBoost', 'KNN', 'NB'),optimise='INF')
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
        modell <- train(Class~., data=resp_fs(),
                        method=input$model_from_train, 
                        metric='ROC',
                        trControl= trainControl(method="repeatedcv", 
                                                number=input$mod_gen_repeats, 
                                                repeats=input$repeatedcv_number,
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

