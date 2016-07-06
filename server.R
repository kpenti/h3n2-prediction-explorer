library(shiny)
library(ggplot2)
library(caret)
library(xtable)
require(e1071);
library("party")
library(corrplot)

# compute Matthews correlation coefficient
mcc <- function(tab, pos = colnames(tab)[1])
{
  if(nrow(tab) != 2 | ncol(tab) != 2) stop("A 2x2 table is needed")
  neg <- colnames(tab)[colnames(tab) != pos]
  tp <- tab[pos, pos]
  tn <- tab[neg, neg]
  fp <- tab[pos,neg]
  fn <- tab[neg, pos]
  d1 <- as.double(tp + fp)
  d2 <- as.double(tp + fn)
  d3 <- as.double(tn + fp)
  d4 <- as.double(tn + fn)
  
  cat(tp, tn, fp, fn)
  
  if(d1 == 0 | d2 == 0 | d3 == 0 | d4 == 0) return(0)
  ((tp * tn) - (fp * fn))/sqrt(d1*d2*d3*d4)  
}

# compute stats
compute_stats <- function(predictions, results)
{
  
  predictions <- factor( ifelse(predictions > 0, 1, 0), levels=c(0,1) )
  print(33)
  #print(predictions)
  #print(results)
  #return (predictions)
  result = confusionMatrix(predictions, results, "1")
  cm <- as.matrix(result)
  print('cm-compute_stats')
  print(cm)
  mcc <- mcc(cm, '1')
  
  #print(results)
  #results <- factor( ifelse(results > 0, 1, 0), levels=c(0,1) )
  
  #print(results)
  #print(predictions)
  
  d <- data.frame(
    accuracy=c(result$overall[1]),
    sensitivity=c(sensitivity(predictions, results, "1")),
    specificity=c(sensitivity(predictions, results, "0")),
    mcc=c(mcc),
    predictions = length(predictions),
    results = length(results)
  )
  #print(d)
  return (d)
}

# compute stats
empty_stats <- function()
{
  
  d <- data.frame(
    accuracy=c(0),
    sensitivity=c(0),
    specificity=c(0),
    mcc=c(0),
    predictions = c(0),
    results = c(0),
    deviance = c(0),
    AIC = c(0),
    r2 = c(0)
  )
  #print(d)
  return (d)
}

show_table <- function(stats)
{
  #print(stats)
  
  d <- data.frame(
    accuracy=c(stats$accuracy),
    sensitivity=c(stats$sensitivity),
    specificity=c(stats$specificity),
    mcc=c(stats$mcc),
    predictions=c(stats$predictions),
    deviance=c(stats$deviance),
    AIC=c(stats$AIC),
    prediction_variance=c(stats$variance),
    r2=c(stats$r2)
  )
  
  return(d)
  
}

features <- c(
  #'cluster_0_mutations', "cluster_0_dist","cluster_1_mutations",'cluster_1_dist','external_mutations','ha1_reversions',
  #"cluster_2_mutations",'cluster_2_dist',"cluster_3_mutations",'cluster_3_dist',
  #"total_cluster_mutations", 'mutations_outside_main_cluster', 'nearest_cluster_to_rbs', 'mutations_outside_clusters',
  #"vol", "pI", "dist_to_rbs","rsa",'other_mutations', 'tail_mutations', 'buried_mutations',  'signal_peptide', 'fusion_peptide', 
  #'cluster_count', 'smith_important_residues', 'different_mutations', 'total_ha1_mutations', 'total_mutations', 
  #'patch_count', 'loops', 'strands', 'helixes' , 'patch_rank' ,
  #"bp_mutations", 'bp_different_mutations', 'bp_rsa', 'bp_loops', 'bp_strands', 'bp_helixes',
  'smith_important_residues', 'smith_145', 'smith_155', 'smith_156', 'smith_158', 'smith_159', 'smith_189', 'smith_193'
)
features <- c( 
  'constant','vol','dist_to_rbs','rsa','other_mutations','tail_mutations','buried_mutations','signal_peptide','fusion_peptide',
  'cluster_count','smith_important_residues','total_ha1_mutations','patch_count','loops','strands','helixes',
  'cluster_0_mutations','cluster_0_dist','cluster_1_mutations','cluster_1_dist','cluster_2_mutations','cluster_2_dist',
  'cluster_3_mutations','cluster_3_dist',
  'external_mutations','total_cluster_mutations','mutations_outside_main_cluster',
  #'nearest_cluster_to_rbs',
  'nearest_cluster_rbs_mutations',
  'nearest_cluster_rbs_distance',
  'nearest_cluster_core_mutations',
  'mutations_outside_clusters',
  'smith_145','smith_155','smith_156','smith_158','smith_159','smith_189','smith_193','pI','core_important_residues',
  'resi_1','resi_2','resi_3','resi_4','resi_5','resi_6','resi_7','resi_8','resi_9','resi_10','resi_11','resi_12','resi_13','resi_14','resi_15','resi_16','resi_17','resi_18','resi_19','resi_20','resi_21','resi_22','resi_23','resi_24','resi_25','resi_26','resi_27','resi_28','resi_29','resi_30','resi_31','resi_32','resi_33','resi_34','resi_35','resi_36','resi_37','resi_38','resi_39','resi_40','resi_41','resi_42','resi_43','resi_44','resi_45','resi_46','resi_47','resi_48','resi_49','resi_50','resi_51','resi_52','resi_53','resi_54','resi_55','resi_56','resi_57','resi_58','resi_59','resi_60','resi_61','resi_62','resi_63','resi_64','resi_65','resi_66','resi_67','resi_68','resi_69','resi_70','resi_71','resi_72','resi_73','resi_74','resi_75','resi_76','resi_77','resi_78','resi_79','resi_80','resi_81','resi_82','resi_83','resi_84','resi_85','resi_86','resi_87','resi_88','resi_89','resi_90','resi_91','resi_92','resi_93','resi_94','resi_95','resi_96','resi_97','resi_98','resi_99','resi_100','resi_101','resi_102','resi_103','resi_104','resi_105','resi_106','resi_107','resi_108','resi_109','resi_110','resi_111','resi_112','resi_113','resi_114','resi_115','resi_116','resi_117','resi_118','resi_119','resi_120','resi_121','resi_122','resi_123','resi_124','resi_125','resi_126','resi_127','resi_128','resi_129','resi_130','resi_131','resi_132','resi_133','resi_134','resi_135','resi_136','resi_137','resi_138','resi_139','resi_140','resi_141','resi_142','resi_143','resi_144','resi_145','resi_146','resi_147','resi_148','resi_149','resi_150','resi_151','resi_152','resi_153','resi_154','resi_155','resi_156','resi_157','resi_158','resi_159','resi_160','resi_161','resi_162','resi_163','resi_164','resi_165','resi_166','resi_167','resi_168','resi_169','resi_170','resi_171','resi_172','resi_173','resi_174','resi_175','resi_176','resi_177','resi_178','resi_179','resi_180','resi_181','resi_182','resi_183','resi_184','resi_185','resi_186','resi_187','resi_188','resi_189','resi_190','resi_191','resi_192','resi_193','resi_194','resi_195','resi_196','resi_197','resi_198','resi_199','resi_200','resi_201','resi_202','resi_203','resi_204','resi_205','resi_206','resi_207','resi_208','resi_209','resi_210','resi_211','resi_212','resi_213','resi_214','resi_215','resi_216','resi_217','resi_218','resi_219','resi_220','resi_221','resi_222','resi_223','resi_224','resi_225','resi_226','resi_227','resi_228','resi_229','resi_230','resi_231','resi_232','resi_233','resi_234','resi_235','resi_236','resi_237','resi_238','resi_239','resi_240','resi_241','resi_242','resi_243','resi_244','resi_245','resi_246','resi_247','resi_248','resi_249','resi_250','resi_251','resi_252','resi_253','resi_254','resi_255','resi_256','resi_257','resi_258','resi_259','resi_260','resi_261','resi_262','resi_263','resi_264','resi_265','resi_266','resi_267','resi_268','resi_269','resi_270','resi_271','resi_272','resi_273','resi_274','resi_275','resi_276','resi_277','resi_278','resi_279','resi_280','resi_281','resi_282','resi_283','resi_284','resi_285','resi_286','resi_287','resi_288','resi_289','resi_290','resi_291','resi_292','resi_293','resi_294','resi_295','resi_296','resi_297','resi_298','resi_299','resi_300','resi_301','resi_302','resi_303','resi_304','resi_305','resi_306','resi_307','resi_308','resi_309','resi_310','resi_311','resi_312','resi_313','resi_314','resi_315','resi_316','resi_317','resi_318','resi_319','resi_320','resi_321','resi_322','resi_323','resi_324','resi_325','resi_326','resi_327','resi_328'
)
training <- c()
testing <- c()
fit <- c()


current_dataset <<- ''
current_patch_size <<- ''

shinyServer(function(input, output, session) {
  
  values <- reactiveValues()
  
  runRegression <- reactive({
    
    print('selected_dataset')
    #loadDataset(input$selected_dataset, input$patch_size)
    my_dataset <<- read.csv( paste("datasets/", input$patch_size, "/", input$selected_dataset, ".csv", sep=""), header = TRUE, stringsAsFactors=FALSE)
    threshold <<- 2
    threshold <<- input$threshold
    print('threshold', threshold)
    print(threshold)
    print('threshold', input$threshold)
    print(input$threshold)
    my_dataset$escape <<- factor( ifelse(my_dataset$antigenic_distance >= threshold, 1, 0) )  
    my_dataset$constant <<- 1
    #print(my_dataset)
    #features <- c(    
    #'cluster_0_mutations', 'mutations_outside_main_cluster',  'tail_mutations'
    #)
    
    
    print(input$training_range[2])
    
    training <<- my_dataset[my_dataset$min_year>=as.integer(input$training_range[1]) & my_dataset$max_year<as.integer(input$training_range[2]),]
    testing <<- my_dataset[my_dataset$min_year>=as.integer(input$testing_range[1]) & my_dataset$max_year<=as.integer(input$testing_range[2]),]
    
    if(input$selected_dataset == 'liao_split') {
      
      training <<- my_dataset[my_dataset$split == 'T',]
      testing <<- my_dataset[my_dataset$split == 'V',]
      
    }
    
    print(training$max_year)
    print(testing$max_year)
    
    values$training = training
    values$testing = testing
    
    TrainData <<- training[,features]
    TestData <<- testing[,features]
    TrainData$escape <<- factor( ifelse(training$antigenic_distance >= threshold, 1, 0) )
    TestData$escape <<- factor( ifelse(testing$antigenic_distance >= threshold, 1, 0) )
    #print(TestResult)
    
    #features_str <- paste(input$selected_features, collapse="+")
    
    selected_features = c(input$s1, input$s2, input$s3, input$s4, input$s5, input$s6, input$s7)
    LL <- list()
    print("input$selected_features")
    print(input$selected_features)
    features_str <- paste(input$selected_features, collapse="+")
    for (resi in selected_features){
      resi_str <- paste("resi_", resi, sep = "")
      features_str <- paste(features_str, resi_str, sep="+")
    }
    #features_str <- substring(features_str, 2)
    if(features_str == "") {
      features_str = "constant"
    }
    
    values$feature_params <- unlist(strsplit(features_str, "\\+"))
    print(features_str)
    print(values$feature_params)
    features_str <- paste("antigenic_distance ~ ", features_str, sep = "")
    print(features_str)
    formula_str = as.formula(features_str)
    print(formula)
    fit$finalModel <<- lm( formula = formula_str, training)
    
    
    print(144)
    newpred <- predict(fit$finalModel, training, interval="confidence", level=input$level)
    print(newpred)
    #newpred <- predict(fit$finalModel, training)
    newpred <- as.data.frame(newpred)
    newpred$escape <- ifelse(newpred$lwr >= threshold, 1, 0)
    newpred$non_escape <- ifelse(newpred$upr <= threshold, 1, 0)
    newpred$remove <- ifelse(newpred$escape == newpred$non_escape, 1, 0)
    values$training_undecided = newpred$remove    
    TrainData$remove <- newpred$remove
    newpred$antigenic_distance <- training$antigenic_distance

    
    TrainData <- subset(TrainData, remove != 1)
    
    values$fitted_train_preds = newpred$fit
    values$fitted_train_observed = newpred$antigenic_distance
    
    newpred <- subset(newpred, remove != 1)
    new_pred_escape <- factor( ifelse(newpred$escape >= 1, 1, 0) )
    
    print(145)
    #print("new_pred_escape")
    #print(new_pred_escape)
    TrainData$holdoutpred <- as.numeric(as.character(new_pred_escape))
    #print(TrainData$holdoutpred)
    #TrainData$holdoutpred <- as.numeric(levels(new_pred_escape))[new_pred_escape]
    #TrainData$holdoutpred <- factor( ifelse(TrainData$holdoutpred >= threshold, 1, 0) )
    #TrainData$holdoutpred <- factor( ifelse(newpred >= threshold, 1, 0) )    
    #print("TrainData$holdoutpred")
    #print(TrainData$holdoutpred)

    #print(TrainData)
    print(234)
    print(newpred$fit)
    values$stats_training <- empty_stats()
    values$stats_testing <- empty_stats()
    values$testing_undecided <- 0
    stats <- empty_stats()    
    if(!identical(newpred$fit, numeric(0))) {
    #if(TRUE) {
      stats <- compute_stats(TrainData$holdoutpred,  TrainData$escape)
    }
    
    if(TRUE) {
      print(157)
      print(stats)
      stats$deviance = c(deviance(fit$finalModel))
      stats$AIC = c(AIC(fit$finalModel))
      print(158)
      print(newpred$fit == 0)
      print(158)
      
      print(newpred$antigenic_distance)
      print(newpred$fit)
      print(159)
      if(!identical(newpred$fit, numeric(0))) {
        variance_lm <- lm(newpred$antigenic_distance ~ newpred$fit)
        print(255)
        stats$variance = c(deviance(variance_lm))
        stats$r2 = summary(variance_lm)$r.squared
        print(stats)
      }
      values$stats_training <- stats
      
      print(168)
      newpred <- predict(fit$finalModel, testing, interval="confidence", level=input$level)
      newpred <- as.data.frame(newpred)
      newpred$escape <- ifelse(newpred$lwr >= threshold, 1, 0)
      newpred$non_escape <- ifelse(newpred$upr <= threshold, 1, 0)
      newpred$remove <- ifelse(newpred$escape == newpred$non_escape, 1, 0)
      TestData$remove <- newpred$remove
      newpred$antigenic_distance <- testing$antigenic_distance
      print(newpred)
      values$fitted_test_preds = newpred$fit
      values$fitted_test_observed = newpred$antigenic_distance 
      
      values$testing_undecided = newpred$remove
      newpred <- subset(newpred, remove != 1)
      new_pred_escape <- factor( ifelse(newpred$escape >= 1, 1, 0) )
      #print(newpred)
      TestData <- subset(TestData, remove != 1)
      
      
      #new_pred_escape <- factor( ifelse(newpred >= threshold, 1, 0) )
      TestData$holdoutpred <- as.numeric(as.character(new_pred_escape))
    }
    
    #TrainData$holdoutpred <- as.numeric(levels(new_pred_escape))[new_pred_escape]
    #TrainData$holdoutpred <- factor( ifelse(TrainData$holdoutpred >= threshold, 1, 0) )
    #TrainData$holdoutpred <- factor( ifelse(newpred >= threshold, 1, 0) )    
    if(!identical(newpred$fit, numeric(0))) {
      stats <- compute_stats(TestData$holdoutpred, TestData$escape)
      print(169)
      testLm <<- glm( formula = formula_str, testing, family=gaussian(link=identity))
      print(178)
      stats$deviance = c(deviance(testLm))
      stats$AIC = c(AIC(testLm))
      stats$r2 = summary(testLm)$r.squared
      
      
      #abline(lm(testing$antigenic_distance ~ TestData$holdoutpred))
      
      #print("abline")
      #print(testing$antigenic_distance)
      #print(newpred)
      
      variance_lm <- lm(newpred$antigenic_distance ~ newpred$fit)
      stats$variance = c(deviance(variance_lm))
      stats$r2 = summary(variance_lm)$r.squared                
      values$stats_testing <- stats
    }
    #par(mfrow = c(2,2))
    #dev.control(displaylist="enable") # enable display list 
    #plot(fit$finalModel)
    values$plotfinalModel <- fit$finalModel
    #values$plot <- recordPlot() # load displaylist into variable 
    values$summary = summary(fit$finalModel)
    #values$stats_training = stats_training
    #values$stats_testing = stats_collection
    
    #print("--- TESTING ---")
    #print(summary(testLm))
    #print(values$testing_results)
    
    print(formula)
    values$tree <- ctree(formula_str, data=training)
    trainpred <- predict(values$tree, newdata= training)
    trainpred <- factor( ifelse(trainpred >= threshold, 1, 0) )
    #print(trainpred)
    testpred <- predict(values$tree, newdata= testing)
    testpred <- factor( ifelse(testpred >= threshold, 1, 0) )
    #print(testpred)
    values$train_preds <- testpred
    values$test_preds <- testpred
    
    
    tree_stats <- data.frame(
      train_mcc = 0,
      test_mcc = 0
    )
    levels(trainpred) <- c(0,1)
    levels(testpred) <- c(0,1)
    result = confusionMatrix(trainpred, training$escape, "1")
    cm <- as.matrix(result)
    tree_stats$train_mcc <- mcc(cm, '1')
    
    result = confusionMatrix(testpred, testing$escape, "1")
    cm <- as.matrix(result)
    tree_stats$test_mcc <- mcc(cm, '1')
    #print( 'tree_stats' )
    #print( tree_stats )
    values$tree_stats <- tree_stats
    
    return (fit);
    #lm(as.formula(paste(input$dependent," ~ ",paste(input$independent,collapse="+"))),data=dat)
  })
  
  output$regTab <- renderTable({
    selected_features = c(input$selected_features, input$s1, input$s2, input$s3, input$s4, input$s5, input$s6, input$s7)
    if(!is.null(selected_features)){
      fit = runRegression();
      print(summary(fit$finalModel))
      #summary(fit$finalModel)$coefficients
      summary(fit$finalModel)
    } else {
      print(data.frame(Warning="Please select Model Parameters."))
    }
  })
  
  
  output$regCounts <- renderText({
    #loadDataset(input$selected_dataset, input$patch_size)
    print(values$training_undecided)
    paste("Split : ", nrow(values$training), "/", nrow(values$testing), "<br />training - antigenic: ", sum(values$training$escape == 1), "non-antigenic: ", sum(values$training$escape == 0), "<br />testing - antigenic: ", sum(values$testing$escape == 1), "non-antigenic: ", sum(values$testing$escape == 0), "<br />training undecided: ", sum(values$training_undecided == 1), "<br />testing undecided: ", sum(values$testing_undecided == 1))    
  })
  
  output$regPlot <- renderPlot({
    selected_features = c(input$selected_features, input$s1, input$s2, input$s3, input$s4, input$s5, input$s6, input$s7)
    if(!is.null(selected_features)){
      fit = runRegression();
      #print('regPlot')
      par(mfrow = c(2,2))
      return(plot(fit$finalModel))
    }
  })
  
  output$text <- renderTable({
    selected_features = c(input$selected_features, input$s1, input$s2, input$s3, input$s4, input$s5, input$s6, input$s7)
    #if(!is.null(selected_features)){
    if(TRUE){
      fit = runRegression();
      if(TRUE) {
        #print("predicting training runRegression...")
        #print(TrainData)
        probsTest <- predict(fit$finalModel, TrainData)
        #print(probsTest)
        TrainPred <- factor( ifelse(probsTest >= threshold, 1, 0), levels=c(0,1) )
        #print(TrainPred)
        #print(training$escape)
        result = confusionMatrix(TrainPred, training$escape, "1")
        
        cm <- as.matrix(result)
        training_mcc <- mcc(cm, '1')
        #print(result)
      }
      
      #print("predicting testing...")
      probsTest <- predict(fit$finalModel, TestData)
      TestPred <- factor( ifelse(probsTest >= threshold, 1, 0), levels=c(0,1) )
      TestResult <- factor( TestData$escape, levels=c(0,1) )
      #print(TestData)
      print("\n--- TESTING ---\n")
      print(length(TestPred))
      print(length(TestResult))
      print("\n--- TESTING ---\n")
      confusionTable <- table(TestPred, TestResult)
      #print(confusionTable)
      result = confusionMatrix(confusionTable, "1")
      
      cm <- as.matrix(result)
      testing_mcc <- mcc(cm, '1')
      #print( fit$finalModel[["Resid. Dev"]] )
      
      d <- data.frame(training_mcc=c(training_mcc),
                      testing_mcc=c(testing_mcc),
                      sensitivity=c(sensitivity(TestPred, TestResult, "1")),
                      specificity=c(sensitivity(TestPred, TestResult, "0")),
                      accuracy=c(result$overall[1]),
                      deviance=c(deviance(fit$finalModel)),
                      AIC=c(AIC(fit$finalModel))
      )
      
      return(d)
      
    }
    
    
  }, digits=2)
  output$plot <- renderPlot({
    par(mfrow = c(2,2))
    dev.control(displaylist="enable") # enable display list 
    #plot(values$plotfinalModel)
    #values$plotfinalModel <- fit$finalModel
    #values$plot <- recordPlot() # load displaylist into variable 
    return(plot(values$plotfinalModel))
  })
  output$plotTree <- renderPlot({
    #print(values$tree)
    return(plot(values$tree, 'simple'))
  })
  output$tableTest <- renderTable({
    return(show_table(values$stats_testing))
  })
  output$statsTree <- renderTable({
    return(values$tree_stats)
  })
  output$tableTrain <- renderTable({
    #print('values$stats_training')
    #print(values$stats_training)
    return(show_table(values$stats_training))
  })
  output$tableSummary <- renderTable({
    return(values$summary)
  })
  output$predictionsList <- renderText({
    return(as.numeric(as.character(values$testing_pred)))
  })
  output$actualList <- renderText({
    return(as.numeric(as.character(values$testing_results)))
  })
  
  output$plotFittedTraining <- renderPlot({
    #print("values$fitted_test")
    #print(values$fitted_train_preds)
    #print(values$fitted_train_observed)
    
    p <- plot(values$fitted_train_preds, values$fitted_train_observed, main="Training", 
              xlab="Predictions ", ylab="Observed")
    
    return(p)
  })
  
  output$plotFittedTesting <- renderPlot({
    #print("values$fitted_test")
    #print(values$fitted_test_preds)
    #print(values$fitted_test_observed)
    
    p <- plot(values$fitted_test_preds, values$fitted_test_observed, main="Testing", 
              xlab="Predictions ", ylab="Observed")
    
    return(p)
  })
  
  output$tableFittedTesting <- renderTable({
    return(values)
  })
  
  output$correlationMatrix  <- renderImage({
    
    print("correlationMatrix")    
    print(values$feature_params)
    print(values$training[,values$feature_params])
    
    outfile <- tempfile(fileext='.png')
    
    
    # Generate the PNG
    png(outfile, width=1200, height=800)
    coMatrix <- cor(values$training[,values$feature_params], use="pairwise.complete.obs")
    coMatrix[is.na(coMatrix)] <- 0
    #print(coMatrix)
    corrplot(coMatrix, method = "number")
    dev.off()
    
    # Return a list containing the filename
    list(src = outfile,
         contentType = 'image/png',
         width = 1200,
         height = 800,
         alt = "This is alternate text")
    }, deleteFile = TRUE)
  
})