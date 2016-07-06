library(shiny)
library(ggplot2)


#training <- read.csv("training.csv", header = TRUE, stringsAsFactors=FALSE)
#testing = read.csv("testing.csv", header = TRUE, stringsAsFactors=FALSE)

#threshold <- 2
#training$escape <- factor( ifelse(training$antigenic_distance >= threshold, 1, 0) )
#testing$escape <- factor( ifelse(testing$antigenic_distance >= threshold, 1, 0) )


features <- c(
  'cluster_0_mutations', 'mutations_outside_main_cluster',  'tail_mutations',  'ha1_reversions',  'smith_important_residues'
)

features <- c(
  'smith_important_residues', 'cluster_0_dist', 'cluster_1_dist',  'rsa',   'bp_mutations',  'bp_rsa',  'bp_strands'
)

features <- c(
  'smith_important_residues', 'constant'
)
feature_list <- c(
  'smith_important_residues', 
  'smith_145', 'smith_155', 'smith_156', 'smith_158', 'smith_159', 'smith_189', 'smith_193',
  'cluster_0_mutations', "cluster_0_dist","cluster_1_mutations",'cluster_1_dist','external_mutations','ha1_reversions',
  "cluster_2_mutations",'cluster_2_dist',"cluster_3_mutations",'cluster_3_dist',
  "total_cluster_mutations", 'mutations_outside_main_cluster', 'nearest_cluster_to_rbs', 'mutations_outside_clusters',
  "vol", "pI", "dist_to_rbs","rsa",'other_mutations', 'tail_mutations', 'buried_mutations',  'signal_peptide', 'fusion_peptide', 
  'cluster_count', 'different_mutations', 'total_ha1_mutations', 'total_mutations', 
  'patch_count', 'loops', 'strands', 'helixes', 'patch_rank' ,
  "bp_mutations", 'bp_different_mutations', 'bp_rsa', 'bp_loops', 'bp_strands', 'bp_helixes'
)
feature_list <- c(
  'smith_important_residues',
  'smith_145', 'smith_155', 'smith_156', 'smith_158', 'smith_159', 'smith_189', 'smith_193'
  #'cluster_0_mutations', "cluster_0_dist","cluster_1_mutations",'cluster_1_dist','external_mutations','ha1_reversions',
  #"cluster_2_mutations",'cluster_2_dist',"cluster_3_mutations",'cluster_3_dist',
  #"total_cluster_mutations", 'mutations_outside_main_cluster', 'nearest_cluster_to_rbs', 'mutations_outside_clusters',
  ##"vol", "pI", "dist_to_rbs","rsa",'other_mutations', 'tail_mutations', 'buried_mutations',  'signal_peptide', 'fusion_peptide', 
  #'cluster_count', 'different_mutations', 'total_ha1_mutations', 'total_mutations', 
  #'patch_count', 'loops', 'strands', 'helixes'
)
feature_list <- c( 
  'constant', 'smith_important_residues',
  #'smith_145', 'smith_155', 'smith_156', 'smith_158', 'smith_159', 'smith_189', 'smith_193',
  'vol','pI','dist_to_rbs','rsa','other_mutations','tail_mutations','buried_mutations','signal_peptide','fusion_peptide',
  'cluster_count','total_ha1_mutations','patch_count','loops','strands','helixes',
  'cluster_0_mutations','cluster_0_dist','cluster_1_mutations','cluster_1_dist','cluster_2_mutations','cluster_2_dist',
  'cluster_3_mutations','cluster_3_dist',
  'external_mutations','total_cluster_mutations','mutations_outside_main_cluster',
  #'nearest_cluster_to_rbs',
  'nearest_cluster_rbs_mutations',
  'nearest_cluster_rbs_distance',
  'nearest_cluster_core_mutations',
  'mutations_outside_clusters',
  'mutations_outside_clusters','core_important_residues'
)

selected_dataset = 'full'
dataset_list <- c(
  'liao',
  'bedford',
  'full'
)
dataset_list <- c(
  "Liao original split" = "liao_split",
  "Liao (1968-2004)" = "liao",
  "Bedford full (1968-2011)" = "bedford",
  "Bedford filtered (1968-2011)" = "bedford_min",
  "Everything  (1968-2014)" = "full"
)

selected_patch_size = 20

#TrainData <- training[,features]
#TestData <- testing[,features]

s1_selected <- c(
  #145,155,156,158,159,189,193
)

shinyUI(pageWithSidebar(
  
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "theme.css")
  ),
  
  sidebarPanel(
    selectInput('selected_dataset', 'Dataset', dataset_list, selected=selected_dataset),
    sliderInput('training_range','Training range', min=1968, max=2015, value=c(1968,2007),format="###0",sep = ""),
    sliderInput('testing_range','Testing range', min=1968, max=2015, value=c(2007,2015),format="###0",sep = ""),
    selectInput("patch_size", "Patch size", c(20,30), selected=selected_patch_size),
    numericInput("threshold", "Threshold:", 2, min = 1, max = 4, step=0.5),
    numericInput("level", "Level:", 0.95, min = 0.1, max = 1, step=0.01),
    #selectInput("method", "Method", c("glm","lm")),
    checkboxGroupInput("selected_features", "Features:", feature_list, selected=features),
    checkboxGroupInput("s1", "Smith residues:",
                       c(145,155,156,158,159,189,193), selected=s1_selected, inline = TRUE),
    
    checkboxGroupInput("s2", "RBS residues:",
                       c(98,134,135,136,137,138,153,155,183,190,194,195,224,225,226,227,228,229), inline = TRUE),
    
    
    checkboxGroupInput("s3", "Canonical site A:",
                       c(122, 124, 126, 130,131,132,133, 135, 137, 138, 140, 142,143,144,145,146, 150,152, 168), inline = TRUE),
    
    checkboxGroupInput("s4", "Canonical site B:",
                       c(128, 129, 155,156,157,158,159,160, 163,164,165, 186,187,188,189,190, 192,193,194, 196,197,198), inline = TRUE),
    
    checkboxGroupInput("s5", "Canonical site C:",
                       c(44, 45, 46, 47, 48, 50, 51, 53, 54, 273, 275, 276, 278, 279, 280, 294, 297, 299, 300, 304, 305, 307, 308, 309, 310, 311, 312), inline = TRUE),
    
    checkboxGroupInput("s6", "Canonical site D:",
                       c(96, 102, 103, 117, 121, 167,170, 171, 172, 173, 174, 175, 176, 177, 179, 182, 201, 203,207, 208, 209, 212, 213, 214, 215, 216, 217, 218, 219, 226, 227, 228, 229, 230,238, 240, 242, 244, 246, 247, 248), inline = TRUE),
    
    checkboxGroupInput("s7", "Canonical site E:",
                       c(57, 59, 62, 63, 67, 75, 78,80,81,82,83, 86,87,88, 91, 92, 94,109, 260,261,262, 265), inline = TRUE), 
    #br()
    submitButton("Submit")
  ),
  
  mainPanel(
    tableOutput("text"),
    hr(),
    #plotOutput('regPlot', width = "100%", height = "600px"),
    #br(),
    #hr(),
    plotOutput("plot", width = "100%", height = "600px"),
    br(),
    br(),
    "> training", tableOutput("tableTrain"),
    br(),
    "> testing", tableOutput("tableTest"),
    br(),
    "> summary", tableOutput("tableSummary"),
    hr(),
    htmlOutput("regCounts"),
    hr(),
    #"> predictions",
    #htmlOutput("predictionsList"),
    #"> actual",
    #htmlOutput("actualList"),
    "> fitted vs observed",
    plotOutput("plotFittedTraining", width = "100%", height = "600px"),
    plotOutput("plotFittedTesting", width = "100%", height = "600px"),
    imageOutput("correlationMatrix" )
    #tableOutput("tableFittedTesting")
    #hr(),
    #"> tree",
    #plotOutput("plotTree", width = "100%", height = "600px"),
    #br(),
    #"> stats", tableOutput("statsTree")
  )
))
