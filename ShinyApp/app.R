##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Shiny App
## Date: 25.07.2018
## Author: Mahmoud Hallal
##################################################

## Load libraries
#library(BiocManager)
#options(repos = BiocManager::repositories())
#source("http://bioconductor.org/biocLite.R")
#options(repos = BiocInstaller::biocinstallRepos())
#getOption("repos")
options(repos = BiocManager::repositories())
options("repos")
#biocLite("pathview")
#biocLite("GDCRNATools")

#library(GDCRNATools)
#library(XVector)
library(shiny)
library(pheatmap)
#install.packages("filesstrings", repos="http://cran.rstudio.com/")
#library(rstring)
library(filesstrings)
library(ggplot2)
library(pheatmap)
library(VennDiagram)
library(DT)
library(yaml)
library(plotly)
#library(bindrcpp)
library(reshape2)
library(dplyr)
#detach("package:dplyr", unload=TRUE)
library(pathview)
#data("gene.idtype.bods")
flog.threshold(ERROR)
library(org.Sc.sgd.db)

#detach("package::plotly",unload = TRUE)
#detach("package::dplyr",unload = TRUE)

## Create Shiny function
shinyPathview2 <- function () 
{  
  ## Define human pathways
  pathways_hsa <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                    'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                    'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
  
  ## Define yeast pathways
  pathways_sce <- c('Autophagy animal','MAPK Signaling Pathway','Hippo signaling pathway')
  
  ## Define user interface
  ui <- fluidPage(
    headerPanel("KAEA: Kinase Activity Enrichment Analysis"), 
    sidebarPanel(
      conditionalPanel(condition="input.tabselected==1",
                       h1("Data Upload"),
                       fileInput("file1", "Choose a Shiny output file",
                                 accept = c(".Rda")),
                       checkboxInput("checkbox1", label = "Human", value = TRUE),
                       checkboxInput("checkbox2", label = "Mouse", value = FALSE),
                       checkboxInput("checkbox3", label = "Yeast", value = FALSE),
                       

                       actionButton("go", "Go")
      ),
      conditionalPanel(#condition="input.checkbox2 == TRUE",
        condition="input.tabselected==2",
        uiOutput("dataset22"), 
        
        selectInput("xcol", "Pathway", pathways_hsa)),
      
      
      conditionalPanel(condition = "input.tabselected==3"),
      
      #topGO panel
      conditionalPanel(
        condition="input.tabselected==7",
        uiOutput("dataset23"), 
        selectInput("Direction", "Direction",c("less","greater")),
        selectInput("DB", "Database", c("MF","BP")))
    ),
    
    mainPanel(tabsetPanel(
      tabPanel("About", value = 1 ,helpText("Kinase Activity Enrichment Analysis (KAEA) is a tool to infer differential 
                                            kinase activity between tested conditions based on experimentally validated 
                                            databases and in-silico kinase-subsrate databse from NetworKIN")),
      tabPanel("Data Quality", value = 3, 
               tabsetPanel(
                 navbarMenu("Sites",
                            tabPanel("Histogram", plotOutput("hist_sites", width = "70%", height = "600px")),
                            tabPanel("Venn Diagram", plotOutput("venn_sites", width = "70%", height = "600px")),
                            tabPanel("Coefficient of Variation", plotOutput("CVs", width = "70%", height = "600px"))
                            ),
                 
                 navbarMenu("Peptides",
                            tabPanel("Histogram", plotOutput("hist_peps", width = "70%", height = "600px")),
                            tabPanel("Venn Diagram", plotOutput("venn_peps", width = "70%", height = "600px"))),
                 navbarMenu("Proteins",
                            tabPanel("Histogram", plotOutput("hist_prots", width = "70%", height = "600px")),
                            tabPanel("Venn Diagram", plotOutput("venn_prots", width = "70%", height = "600px")))
               )
               
      ),
      tabPanel("Phosphosites", value = 2, 
               tabPanel("Volcano plot", fluidRow(
                 column(6, plotlyOutput(outputId = "Volcano", height = "600px")),
                 column(6, plotlyOutput(outputId = "correlation", height = "600px"))
               )) ),
      tabPanel("Enrichment", value = 2,
               tabsetPanel(
                 tabPanel("Heatmap2",
                   fluidRow(
                    column(6, plotlyOutput("Heatmap", width = "100%", height = "600px")),
                    column(6, plotlyOutput(outputId = "correlation66", width = "100%", height = "600px"))
                          )),
                 
                 tabPanel("Barplot2", 
                          fluidRow(
                            column(6, plotlyOutput(outputId = "Barplot", width = "100%", height = "600px")),
                            column(6, plotlyOutput(outputId = "correlation3", width = "100%", height = "600px"))),
                          fluidRow(
                            column(10, offset=1, plotlyOutput(outputId = "correlation2", width = "100%", height = "500px"))
                          ),
                          
                          fluidRow(
                            column(10, offset=1, plotlyOutput(outputId = "correlation21", width = "100%", height = "500px"))
                            )
                          )
                 )
               ),
      
      tabPanel("KEGG pathways", value = 2,
               tabsetPanel(
                 tabPanel("Table",dataTableOutput("mytable1")),
                 tabPanel("KEGG pathway", plotOutput("plot1", width = "10%", height = "600px"))
               )),
      
      #tabPanel("GO terms", value = 7 ,dataTableOutput("mytable2")),
      id = "tabselected")
      )
    )
  
  ## Define server
  server <- function(input, output, session) {
    options(shiny.maxRequestSize=30*1024^2)
    
    ## Define the variables
    inputf <- reactive({input$file1$name})
    cell_line <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\2",inputf())})
    FDR_value <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\4",inputf())})
    P_value <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\3",inputf())})
    
    ## Define the directory
    directory <- eventReactive(input$go, {
      directory1 = reactive({paste0("../results/",cell_line(),'_',P_value(),'P_',FDR_value(),'FDR_imp',imp,"/pathview")})
      
      if (!dir.exists(directory1())) {
        dir.create(directory1())
        directory=directory1()
      }
      if (!endsWith(directory1(), "/")) {
        directory = paste(directory1(), "/", sep = "")
      }
      directory
    })
    
    ## Define the input file
    my_data <- eventReactive(input$go, {
      req(input$file1)
      inFile <- input$file1 
      if (is.null(inFile))
        return(NULL)
      data <- get(load(inFile$datapath,.GlobalEnv))
      data
    })
    
    ## Pathway names
    output$dataset22 <- renderUI({
      dd <- names(my_data()$table_rep)
      selectInput("DS", "Dataset", dd)
    })
    
    output$dataset23 <- renderUI({
      dd <- names(my_data()$table_rep)
      selectInput("DS", "Dataset", dd)
    })
    
    data_to_use <- reactive({req(input$DS)})
    
    ## Table
    output$mytable1 <- DT::renderDataTable({
      data_tb <- my_data()$table_rep
      data_tb2 <- as.data.frame(data_tb[[data_to_use()]])
      data_tb2
    }, options = list(searching=FALSE, paging = FALSE))
    
    ## KEGG pathways enrichment
    output$mytable2 <- renderDataTable({
      req(input$DB)
      
      enriched_pathways <- my_data()$pathways
      dbID <- strsplit(input$DB, "~", fixed = TRUE)[[1]][1]
      Dir <- strsplit(input$Direction, "~", fixed = TRUE)[[1]][1]
      
      enriched_pathways2 <- as.data.frame(enriched_pathways[[data_to_use()]][[Dir]])
      
      #choose the desired DB and filter out some columns
      enriched_pathways2 <- enriched_pathways2[enriched_pathways2$ontology == dbID,]
      
      enriched_pathways2
    }, options = list(searching=FALSE, paging = FALSE))
    
    #KEGG pathways
    #detach(package:plotly, unload=TRUE)
    #detach(package:dplyr, unload=TRUE)
    library(pathview)
    output$plot1 <- renderImage({
      
      material_for_waterfall <- my_data()$waterfall_rep
      datasets <- lapply(1:length(material_for_waterfall), function(x){
        df1 <- material_for_waterfall[[x]][!duplicated(material_for_waterfall[[x]]$Kinase),]
        
        #make dataframe
        df <- as.data.frame(df1$Enrichment.score)
        rownames(df) <- df1$Kinase
        colnames(df) <- "value"
        
        dd1 <- df
        dd2 <- as.data.frame(dd1)
        
        dd2
      })
      names(datasets) <- names(material_for_waterfall)
      
      ## KEGG Pathway IDs
      pathwayID <- strsplit(input$xcol, "~", fixed = TRUE)[[1]][1]
      if (input$checkbox2 == TRUE ){
        if (pathwayID == "Autophagy animal"){pathwayID <- "mmu04140"}
        if (pathwayID == "MAPK Signaling Pathway"){pathwayID <- "mmu04010"}
        if (pathwayID == "Hippo signaling pathway"){pathwayID <- "mmu04390"}
      } else if (input$checkbox3 == TRUE){
        if (pathwayID == "Autophagy animal"){pathwayID <- "sce04138"}
        if (pathwayID == "MAPK Signaling Pathway"){pathwayID <- "sce04011"}
        if (pathwayID == "Hippo signaling pathway"){pathwayID <- "sce04392"}
      } else if (input$checkbox1 == TRUE){
        if (pathwayID == "Autophagy animal"){pathwayID <- "hsa04140"}
        if (pathwayID == "mTOR Signaling Pathway"){pathwayID <- "hsa04150"}
        if (pathwayID == "Acute Myeloid Leukemia"){pathwayID <- "hsa05221"}
        if (pathwayID == "Pathways in Cancer"){pathwayID <- "hsa05200"}
        if (pathwayID == "RAS Signaling Pathway"){pathwayID <- "hsa04014"}
        if (pathwayID == "MAPK Signaling Pathway"){pathwayID <- "hsa04010"}
        if (pathwayID == "ERBB Signaling Pathway"){pathwayID <- "hsa04012"}
        if (pathwayID == "WNT Signaling Pathway"){pathwayID <- "hsa04310"}
        if (pathwayID == "TGF-BETA Signaling Pathway"){pathwayID <- "hsa04350"}
        if (pathwayID == "JAK-STAT Signaling Pathway"){pathwayID <- "hsa04630"}
        if (pathwayID == "PI3K-AKT Signaling Pathway"){pathwayID <- "hsa04151"}
        if (pathwayID == "Chronic Myeloid Leukemia"){pathwayID <- "hsa05220"}
      }
      
      ## Define WD
      setwd(paste0(getwd()))
      #print(sp_directory)
      
      ## Define output file
      outfile <- paste('./', pathwayID, ".pathview.png", sep = "")
      sp_directory <- c('./')
      #if (!file.exists(outfile)) {
      if (input$checkbox2 == TRUE){
        sp = "mmu"
        type = "SYMBOL"
      } else if (input$checkbox3 == TRUE){
        sp = "sce"
        type = "GENENAME"
      } else {
        sp = "hsa"
        type = "SYMBOL"
        }
      pathview(gene.data = datasets[[data_to_use()]], pathway.id = pathwayID,
               species = sp, gene.idtype = type, limit = list(gene = c(min(datasets[[data_to_use()]]), max(datasets[[data_to_use()]])), cpd = 1))

      list(src = outfile)
    }, deleteFile = FALSE)
    
    library(plotly)
    
    ## Heatmap
    output$Heatmap <- renderPlotly({
      
      data_hm <- as.data.frame(my_data()$Heatmap_rep[[data_to_use()]])
      
      data_hm1 <- reactive({as.data.frame(data_hm[,'Enrichment.score',drop=FALSE])})
      data_hm2 <- reactive({data_hm1()[data_hm1() != 0,,drop=FALSE]})
      #Define my_palette according to values given
      if (sum(data_hm2() >=0) == nrow(data_hm2())){
        my_palette <- my_palette <-colorRampPalette(c("white", "red"))(n = nrow(data_hm2()))
      } else if (sum(data_hm2()<=0) == nrow(data_hm2())) {
        my_palette <- colorRampPalette(c("blue", "white"))(n = nrow(data_hm2()))
      } else {
        my_palette <- colorRampPalette(c("blue", "white", "red"))(n = nrow(data_hm2()))
      }
      
      ## Plot heatmap
      data <- reactive({as.data.frame(data_hm2()[order(data_hm2()[,c('Enrichment.score')],decreasing =T),,drop=F])})
      
      rr <- ggplot(data = data(), mapping = aes(x=0,
                                                y = reorder(rownames(data()),Enrichment.score,order=TRUE),
                                                fill = Enrichment.score, 
                                                text=Enrichment.score)) +
        geom_tile(color = "black") +
        labs(fill='-log10 p-value')  +
        
        theme(axis.line=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks=element_blank(),
              axis.title.x=element_blank(),
              panel.background=element_blank(),
              panel.border=element_blank(),
              panel.grid.major=element_blank(),
              panel.grid.minor=element_blank(),
              plot.background=element_blank()) +

        scale_fill_gradient2(low = "royalblue4", mid = "white",
                             high = "red", midpoint = 0)+
      
        ylab("Kinases") 
      
      ggplotly(rr ,source="source66", tooltip=c("text")  )
    })
    
    ## Barplot
    data_bp <-  reactive({my_data()$waterfall_rep[[data_to_use()]]})
    output$Barplot <- renderPlotly({
      p <- plot_ly(data_bp(),source="source1", 
                   y = reorder(data_bp()$Kinase,data_bp()$Enrichment.score), 
                   x = data_bp()$Enrichment.score, 
                   type = 'bar', 
                   orientation = 'h', 
                   text = data_bp()$category, textposition = 'auto',
                   marker=list(color = data_bp()$Color)) %>%
        layout(showlegend = FALSE,
               xaxis=list(showline = FALSE,
                          showgrid = FALSE,
                          zeroline = FALSE
               ),
               yaxis=list(tickvals=data_bp()$Kinase, 
                          ticktext=data_bp()$Kinase,
                          automargin = TRUE
               ))
      p
    })
    
    ## Couple hover plot to Barplot 1
    output$correlation2 <- renderPlotly({
      all_dbs <- reactive({my_data()$DB})
      
      eventdata <- event_data("plotly_click", source = "source1")
      validate(need(!is.null(eventdata), "Press on a kinase bar to visualize its substrates"))
      
      #pointNumber describes the number of the point desired
      datapoint <- as.numeric(eventdata$pointNumber)[1]+1
      #curveNumber describes the direction of the bars here, 0 for negative and 1 for positive
      curveNum <- eventdata$curveNumber
      
      negative_length <- length(grep("-",data_bp()$Enrichment.score))
      positive_length <- length(grep("\\+",data_bp()$Enrichment.score))
      
      if (curveNum == 0){
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint])
      } else {
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint+negative_length])
      }
      exprs <- my_data()$ExpressSet
      kin_subs <- unique(all_dbs()[grep(kinase, all_dbs()$termName),]$geneID)
      hover_data <- as.data.frame(exprs)[rownames(exprs) %in% kin_subs,]
      hover_data <- as.data.frame(t(hover_data))
      hover_data$Replicate <- gsub("X","",rownames(hover_data))
      hover_data <- melt(hover_data)
      hover_data$Replicate <- factor(hover_data$Replicate, levels=unique(hover_data$Replicate))
      
      pl <- plot_ly(x = hover_data$Replicate, y=2^hover_data$value, name=hover_data$variable, type = "bar") %>%
        layout(title = paste0("Intensities of substrates of ",kinase,'_',FDR_value()), barmode = 'stack')
      pl
    })
    
    #############
    ## Couple hover plot to Barplot 2 (volcano plot)
    data_vlcs <-  reactive({my_data()$Volcano_special[[data_to_use()]]})
    
    output$correlation3 <- renderPlotly({
      all_dbs <- reactive({my_data()$DB})
      
      eventdata <- event_data("plotly_click", source = "source1")
      validate(need(!is.null(eventdata), "Press on a kinase bar to visualize its substrates"))
      
      #pointNumber describes the number of the point desired
      datapoint <- as.numeric(eventdata$pointNumber)[1]+1
      #curveNumber describes the direction of the bars here, 0 for negative and 1 for positive
      curveNum <- eventdata$curveNumber
      
      negative_length <- length(grep("-",data_bp()$Enrichment.score))
      positive_length <- length(grep("\\+",data_bp()$Enrichment.score))
      
      if (curveNum == 0){
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint])
      } else {
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint+negative_length])
      }
      exprs <- my_data()$ExpressSet
      kin_subs <- unique(all_dbs()[grep(paste0('^',kinase,'$'), all_dbs()$termName),]$geneID)
      hover_data <- as.data.frame(data_vlcs())[rownames(data_vlcs()) %in% kin_subs,]
      print(hover_data)
      print(data_vlcs())

      ## Plot volcano
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      pp <- plot_ly(type = 'scatter',source="source",
                    symbols = c('x','o' ,'circle'),
                    symbol = ~hover_data$share) %>%
        ##Overlay the symbols with legend
        add_trace(
          x = hover_data$topTable, 
          y = -log10(hover_data$padj),
          text = paste(hover_data$prots,' (',hover_data$gene,')'),
          color = I('grey'),
          mode = 'markers',
          
          hoverinfo = 'text',
          showlegend = T) %>%
        #overlay the colors with no legend
        add_trace(
          x = hover_data$topTable, 
          y = -log10(hover_data$padj),
          text = paste(hover_data$prots,' (',hover_data$gene,')'),
          color = hover_data$threshold,
          mode = 'markers',
          
          hoverinfo = 'text',
          #marker = list(color='green'),
          showlegend =F) %>%
        layout(xaxis = x, yaxis = y)
      pp
    })
    
    
    #############
    ## Couple hover plot to Barplot 2 (GSEA plot)
    data_gsea <-  reactive({my_data()$Volcano_special[[data_to_use()]]})
    
    output$correlation66 <- renderPlotly({
      all_dbs <- reactive({my_data()$DB})
      topTable <- reactive({my_data()$Volcano_special[[data_to_use()]]})
      
      ##define heatmap input
      data_hm <- as.data.frame(my_data()$Heatmap_rep[[data_to_use()]])
      
      data_hm1 <- reactive({as.data.frame(data_hm[,'Enrichment.score',drop=FALSE])})
      data_hm2 <- reactive({data_hm1()[data_hm1() != 0,,drop=FALSE]})
      
      data <- reactive({as.data.frame(data_hm2()[order(data_hm2()[,c('Enrichment.score')],decreasing =T),,drop=F])})
      
      #event data
      eventdata <- event_data("plotly_click", source = "source66")
      validate(need(!is.null(eventdata), "Press on a kinase bar to visualize its plot"))
      
      #pointNumber describes the number of the point desired
      datapoint <- as.numeric(eventdata$pointNumber[[1]][1])+1
      
      #curveNumber describes the direction of the bars here, 0 for negative and 1 for positive
      curveNum <- eventdata$curveNumber
      
      negative_length <- length(grep("-",data()[,c('Enrichment.score')]))
      positive_length <- length(data()[,c('Enrichment.score')]) - negative_length
      print(negative_length)
      print(positive_length)
      
      kinase <- rownames(data())[(nrow(data())+1)-datapoint]

      exprs <- my_data()$ExpressSet
      kin_subs <- unique(all_dbs()[grep(kinase, all_dbs()$termName),]$geneID)
      kin_subs_present <- rownames(exprs)[rownames(exprs) %in% kin_subs]
      toptable_1 <- topTable()$prots[order(topTable()$topTable, decreasing = T)]
      abline_indices <- which(toptable_1 %in% kin_subs_present)
      
      ## Define the axis of plots
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        rangemode = "tozero"
        
      )
      ax2 <- list(
        title = "",
        zeroline = TRUE,
        showline = FALSE,
        showticklabels = TRUE,
        showgrid = FALSE,
        rangemode = "tozero"
      )
      
      ## Get the topTable and order the values
      dd <- topTable()$topTable[order(topTable()$topTable, decreasing = T)]
      
      ## Create a color for every value
      pal <- colorRampPalette(c("red", "white", "blue"))(length(dd))
      rank.colors <- pal
       
      ## Create a df of color and value
      ll <- as.data.frame(cbind(as.numeric(dd),as.character(rank.colors)))
      
      ## Define x axis limit
      axx <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        range = c(0, nrow(ll))
      )
      
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      
      ## Plot1: red to blue heatmap
      s1 <- plot_ly(x=1:nrow(ll),y=rep(1,nrow(ll)), marker = list(color =  pal),type="bar") %>%
        layout(xaxis = axx,
               yaxis = ax,
               showlegend=F,
               title = kinase
        ) %>%
        #add phosphosites involved
        add_segments(x = ~abline_indices,
                     xend = ~abline_indices,
                     y = rep(0,length(abline_indices)),
                     yend=rep(1,length(abline_indices)),
                     color=I('black'),
                     mode="lines",
                     hoverinfo="text",
                     text=toptable_1[abline_indices])
      
      #Plot2: barplot of all phosphosites p-values
      s2 <- plot_ly(x=1:length(dd),y=dd,type="bar", color = I('grey')) %>%
        layout(xaxis = ax2, yaxis = ax2,showlegend=F)
      
      #Merge the 2 plots together
      subplot(
        s1, s2,
        nrows = 2,
        margin = 0,
        heights = c(0.5, 0.5)
      )
    })
    
    ## Volcano plots
    output$Volcano <- renderPlotly({
      inds <- data_to_use() ==  my_data()$table_rep
      
      data_vlc <-  my_data()$Volcano[[inds]]
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      pp <- plot_ly(type = 'scatter',source="source") %>%
        add_trace(
          x = data_vlc$topTable, 
          y = -log10(data_vlc$padj),
          text = paste(data_vlc$prots,' (',data_vlc$gene,')'),
          color = data_vlc$threshold,
          hoverinfo = 'text',
          #marker = list(color='green'),
          showlegend = F) %>%
        layout(xaxis = x, yaxis = y)
      pp
    })
    
    ## Couple hover plot
    output$correlation <- renderPlotly({
      
      ## Hover part
      exprs <- my_data()$ExpressSet
      
      eventdata <- event_data("plotly_hover", source = "source")
      validate(need(!is.null(eventdata), "Hover over the substrates points to get the quantifications"))
      
      datapoint <- as.numeric(eventdata$pointNumber)[1]+1
      
      site <- data_vlc()[datapoint,]$prots
      gene <- data_vlc()[datapoint,]$gene
      
      hover_data <- exprs[site,]
      #for the sake of silac data, we replace 0s by NAs before inverse log transformation
      hover_data[hover_data == 0] <- NA
      hover_data <- melt(hover_data)
      hover_data$Replicate <- rownames(hover_data)
      hover_data$value <- 2^hover_data$value
      hover_data$Replicate <- factor(hover_data$Replicate, levels=unique(hover_data$Replicate))
      
      # replace again NAs by 0
      hover_data[is.na(hover_data)] <- 0
      
      
      plot_ly(x = hover_data$Replicate, y=hover_data$value, type = "bar") %>%
        layout(title = paste0("Intensities of ",site,' (',gene,')')) #%>%
    })
    
    ## Quality report
    ## Q1: histograms
    ## Histogram sites
    output$hist_sites <- renderPlot({
      
      data1 <- my_data()$hist_sites
      
      ## plot
      ggplot(data1, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data1$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# sites") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    
    ## Histogram peps
    output$hist_peps <- renderPlot({
      data2 <- my_data()$hist_peps
      
      #plot
      ggplot(data2, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data2$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# peps") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    
    ## Histogram prots
    output$hist_prots <- renderPlot({
      data3 <- my_data()$hist_prots
      
      #plot
      ggplot(data3, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data3$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# proteins") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 16)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    
    ## Q2: Venn diagrams
    ## Venn sites
    output$venn_sites <- renderPlot({
      data_v1 <- my_data()$venn_sites
      
      #plot
      grid.draw(venn.diagram(data_v1,
                             fill = 2:(1+length(data_v1)),
                             filename = NULL, 
                             alpha = rep(0.5,length(data_v1)),
                             cat.fontface = 4,
                             lty =2, cex = 1.1, 
                             scaled = TRUE, euler.d = TRUE,
                             main.cex = 1.2))
    })
    
    ## Venn peps
    output$venn_peps <- renderPlot({
      data_v2 <- my_data()$venn_peps
      
      #plot
      grid.draw(venn.diagram(data_v2,
                             fill = 2:(1+length(data_v2)),
                             filename = NULL, 
                             alpha = rep(0.5,length(data_v2)),
                             cat.fontface = 4,
                             lty =2, cex = 1.1, 
                             scaled = TRUE, euler.d = TRUE,
                             main.cex = 1.2))
    })
    
    ## Venn proteins
    output$venn_prots <- renderPlot({
      data_v3 <- my_data()$venn_prots
      
      #plot
      grid.draw(venn.diagram(data_v3,
                             fill = 2:(1+length(data_v3)),
                             filename = NULL, 
                             alpha = rep(0.5,length(data_v3)),
                             cat.fontface = 4,
                             lty =2, cex = 1.1, 
                             scaled = TRUE, euler.d = TRUE,
                             main.cex = 1.2))
    })
    
    output$CVs <- renderPlot({
      CV_Data <- my_data()$CV
      g <- ggplot(CV_Data, aes(cell_line, CV, fill = cell_line))
      g <- g + geom_violin(alpha = 0.5,draw_quantiles = c(0.5)) +
        labs(x="Condition",
             y="CV(%)") +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(#axis.text=element_text(size=14),
          axis.title=element_text(size=24)) +
        guides(fill=FALSE)
      g
    })
    
    ## Additional GSEA with barplot, for check only
    ## Couple hover plot for Barplot 2 (GSEA plot)
    data_gsea <-  reactive({my_data()$Volcano_special[[data_to_use()]]})
    
    output$correlation21 <- renderPlotly({
      all_dbs <- reactive({my_data()$DB})
      topTable <- reactive({my_data()$Volcano_special[[data_to_use()]]})
      
      eventdata <- event_data("plotly_click", source = "source1")
      validate(need(!is.null(eventdata), "Press on a kinase bar to visualize its plot"))
      
      #pointNumber describes the number of the point desired
      datapoint <- as.numeric(eventdata$pointNumber)[1]+1
      #curveNumber describes the direction of the bars here, 0 for negative and 1 for positive
      curveNum <- eventdata$curveNumber
      
      negative_length <- length(grep("-",data_bp()$Enrichment.score))
      positive_length <- length(grep("\\+",data_bp()$Enrichment.score))
      
      if (curveNum == 0){
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint])
      } else {
        kinase <- gsub('(\\w+)([\\+-])','\\1',data_bp()$Kinase[datapoint+negative_length])
      }
      exprs <- my_data()$ExpressSet
      kin_subs <- unique(all_dbs()[grep(kinase, all_dbs()$termName),]$geneID)
      kin_subs_present <- rownames(exprs)[rownames(exprs) %in% kin_subs]
      toptable_1 <- topTable()$prots[order(topTable()$topTable, decreasing = T)]
      abline_indices <- which(toptable_1 %in% kin_subs_present)
      
      ## Define the axis of plots
      ax <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        rangemode = "tozero"
      )
      
      
      ax2 <- list(
        title = "",
        zeroline = TRUE,
        showline = FALSE,
        showticklabels = TRUE,
        showgrid = FALSE,
        rangemode = "tozero"
      )
      
      ## Get the topTable and order the values
      dd <- topTable()$topTable[order(topTable()$topTable, decreasing = T)]
      
      ## Create a color for every value
      pal <- colorRampPalette(c("red", "white", "blue"))(length(dd))
      rank.colors <- pal
      
      ## Create a df of color and value
      ll <- as.data.frame(cbind(as.numeric(dd),as.character(rank.colors)))
      
      ## Define x axis, we use this to extend the range of axis till the end
      axx <- list(
        title = "",
        zeroline = FALSE,
        showline = FALSE,
        showticklabels = FALSE,
        showgrid = FALSE,
        range = c(0, nrow(ll))
      )
      
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      
      ## Plot1: red to blue heatmap
      s1 <- plot_ly(x=0:nrow(ll),y=rep(1,(nrow(ll)+1)), marker = list(color =  pal),type="bar") %>%
        layout(xaxis = axx, 
               yaxis = ax,
               showlegend=F,
               title = kinase
               
        ) %>%
        #add phosphosites involved
        add_segments(x = ~abline_indices, 
                     xend = ~abline_indices, 
                     y = rep(0,length(abline_indices)), 
                     yend=rep(1,length(abline_indices)), 
                     color=I('black'), 
                     mode="lines",
                     hoverinfo="text", 
                     text=toptable_1[abline_indices])

      ## Plot2: barplot of all phosphosites p-values
      s2 <- plot_ly(x=1:length(dd), y=dd, type="bar", color = I('grey')) %>%
        layout(xaxis = ax2, yaxis = ax2, showlegend=F)
      
      #Merge the 2 plots together
      subplot(
        s1, s2, 
        nrows = 2, 
        margin = 0, 
        heights = c(0.5, 0.5)
      )
    })
    
    ## Volcano plots
    data_vlc <-  reactive({my_data()[['Volcano']][[data_to_use()]]})
    xx <- reactive({data_vlc()[['topTable']]})
    yy <-  reactive({data_vlc()[['p']]})
    output$Volcano <- renderPlotly({
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      pp <- plot_ly(type = 'scatter',source="source") %>%
        add_trace(
          x = xx(), 
          y = -log10(yy()),
          text = paste(my_data()[['Volcano']][[data_to_use()]][['prots']],' (',data_vlc()[['gene']],')'),
          color = data_vlc()$threshold,
          hoverinfo = 'text',
          showlegend = F
          ) %>%
        layout(xaxis = x, yaxis = y)
      pp
    })
    
    
    
  }
  shinyApp(ui, server)
}

## Run Shiny App
shinyPathview2()

