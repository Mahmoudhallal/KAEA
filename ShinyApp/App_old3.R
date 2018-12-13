##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Map enriched kinases to selected pathways
## Date: 25.07.2018
## Author: Mahmoud Hallal
##################################################

## Load libraries
#source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")
#biocLite("GDCRNATools")
#setwd("~/Desktop/PhD/Generalized_pipeline/snakemake")
#library(GDCRNATools)
library(shiny)
library(pheatmap)
#install.packages("filesstrings", repos="http://cran.rstudio.com/")
#library(rstring)
library(filesstrings)
#library(png)
library(ggplot2)
library(pheatmap)
library(VennDiagram)

library(yaml)
library(plotly)
#library(bindrcpp)
library(reshape2)
library(dplyr)
#detach("package:dplyr", unload=TRUE)
library(pathview)
data("gene.idtype.bods")
flog.threshold(ERROR)
#detach("package::plotly",unload = TRUE)
#detach("package::dplyr",unload = TRUE)

## Load parameters

## Create Shiny function
shinyPathview2 <- function () 
{  
  pathways <- c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia')
  
  ui <- fluidPage(
    headerPanel("KAEA: Kinase Activity Enrichment Analysis"), 
    sidebarPanel(
      conditionalPanel(condition="input.tabselected==1",
                       h1("Data Upload"),
                       fileInput("file1", "Choose a Shiny output file",
                                 accept = c(
                                   "text/csv",
                                   "text/comma-separated-values,text/plain",
                                   ".csv")),
                       checkboxInput("checkbox", label = "SILAC", value = TRUE),
                       checkboxInput("checkbox", label = "Impute missing values", value = TRUE),
                       
                       actionButton("go", "Go")
      ),
      conditionalPanel(condition="input.tabselected==2",
                       
                       uiOutput("cityControls"), 
                       
                       selectInput("xcol", "Pathway", pathways)),
      conditionalPanel(condition = "input.tabselected==3")
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
                            tabPanel("Coefficient of Variation", plotOutput("CVs", width = "70%", height = "600px"))),
                 
                 navbarMenu("Peptides",
                            tabPanel("Histogram", plotOutput("hist_peps", width = "70%", height = "600px")),
                            tabPanel("Venn Diagram", plotOutput("venn_peps", width = "70%", height = "600px"))),
                 navbarMenu("Proteins",
                            tabPanel("Histogram", plotOutput("hist_prots", width = "70%", height = "600px")),
                            tabPanel("Venn Diagram", plotOutput("venn_prots", width = "70%", height = "600px")))
               )
               
      ),
      
      tabPanel("Enrichment", value = 2,
               tabsetPanel(
                 tabPanel("Heatmap2", plotOutput("Heatmap", width = "50%", height = "600px")),
                 
                 tabPanel("Barplot2", fluidRow(
                   column(7, plotlyOutput(outputId = "Barplot", width = "120%", height = "700px")),
                   column(12, plotlyOutput(outputId = "correlation2", width = "100%", height = "500px")),id = "tabselected"
                 )),
                 
                 tabPanel("Volcano plot", fluidRow(
                   column(6, plotlyOutput(outputId = "Volcano", height = "600px")),
                   column(5, plotlyOutput(outputId = "correlation", height = "600px"))
                 )))
      ),
      
      tabPanel("KEGG pathways", value = 2,
               tabsetPanel(
                 tabPanel("Table", dataTableOutput("mytable1")),
                 tabPanel("KEGG pathway",  plotOutput("plot1", width = "10%", height = "600px"))
               )),id = "tabselected")))
  server <- function(input, output, session) {
    
    ## Define the variables
    inputf <- reactive({input$file1$name})
    cell_line <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\2",inputf())})
    FDR_value <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\4",inputf())})
    P_value <- reactive({gsub("(results_shiny_)(\\w+)_(0.0{1,2}\\d{1,2})P_(0.0{1,2}\\d{1,2})FDR(_\\w+)*(.Rda)","\\3",inputf())})
    
    ## Define the directory
    directory <- eventReactive(input$go, {
      directory1 = reactive({paste0("../results/",cell_line(),'_',P_value(),'P_',FDR_value(),'FDR/',"pathview")})
      
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
    ##
    output$cityControls <- renderUI({
      dd <- names(my_data()$table_rep)
      selectInput("DS", "Dataset", dd)
    })
    
    data_to_use <- reactive({data_to_use <- input$DS})
    
    ## Table
    output$mytable1 <- renderDataTable({
      data_tb2 <- as.data.frame(data_tb[[data_to_use()]])
      data_tb2
    }, options = list(searching=FALSE, paging = FALSE))
    
    #KEGG pathways
    #detach(package:plotly, unload=TRUE)
    #detach(package:dplyr, unload=TRUE)
    
    library(pathview)
    output$plot1 <- renderImage({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      material_for_waterfall <- my_data()$waterfall_rep
      datasets <- lapply(1:length(material_for_waterfall), function(x){
        df1 <- material_for_waterfall[[x]][!duplicated(material_for_waterfall[[x]]$Kinase),]
        #make dataframe
        df <- as.data.frame(df1$Enrichment.score)
        rownames(df) <- df1$Kinase
        colnames(df) <- "value"
        ##
        #dd <- df$value
        #names(dd) <- rownames(df)
        #names(dd) <- gsub('(\\w+)[+-]','\\1',rownames(df))
        #dd1 <- dd[!duplicated(names(dd))]
        dd1 <- df
        dd2 <- as.data.frame(dd1)
        #rownames(dd2) <- names(dd1)
        #colnames(dd2) <- "value"
        dd2
      })
      names(datasets) <- names(material_for_waterfall)
      ##
      
      
      pathwayID <- strsplit(input$xcol, "~", fixed = TRUE)[[1]][1]
      #sp_directory <- paste0(directory,'/',data_to_use(),'/')
      sp_directory <- paste0(directory(),'/',data_to_use(),'_',P_value(),'P_',FDR_value(),'FDR','/')
      print(sp_directory)
      outfile <- paste(sp_directory, pathwayID, ".pathview.png", sep = "")
      if (!file.exists(outfile)) {
        pathview(gene.data = datasets[[data_to_use()]], pathway.id = pathwayID,
                 species = "hsa", gene.idtype = "SYMBOL", limit = list(gene = max(abs(datasets[[data_to_use()]])), cpd = 1), 
                 kegg.dir = sp_directory, same.layer=T, kegg.native = T )
        file.move(paste(pathwayID, ".pathview.png", sep = ""), sp_directory, overwrite = TRUE)
      }
      list(src = outfile)
    }, deleteFile = FALSE)
    library(plotly)
    #heatmap
    #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    
    #print(colnames(data_hm2()))
    # rr <- reactive({sum(data_hm2()>=0)})
    # print(rr())
    # if (rr() >=0 ){
    #   print("OK")}
    # #
    
    output$Heatmap <- renderPlot({
      #reactive({get(load(input$datapath,.GlobalEnv))})
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data_hm <- as.data.frame(my_data()$Heatmap_rep)
      #palette
      
      # 
      data_hm2 <- reactive({as.data.frame(data_hm[,data_to_use(),drop=FALSE])})
      #Define my_palette according to values given
      if (sum(data_hm2() >=0) == nrow(data_hm2())){
        my_palette <- my_palette <-colorRampPalette(c("white", "red"))(n = 1000)
      } else if (sum(data_hm2()<=0) == nrow(data_hm2())) {
        my_palette <- colorRampPalette(c("blue", "white"))(n = 1000)
      } else {
        my_palette <- colorRampPalette(c("blue", "white", "red"))(n = 1000)
      }
      
      #plot heatmap
      rr <- pheatmap(data_hm2() ,color = my_palette, cluster_rows = T, cluster_cols = F, fontsize_row = 14,fontsize_col = 14)
      plot(rr$gtable)
    })
    
    ## Barplot
    #exprs <- shiny_results$ExpressSet
    
    #data <- my_data()$waterfall_rep
    #data_vlc <-  reactive({my_data()$Volcano[[data_to_use()]]})
    
    data_bp <-  reactive({my_data()$waterfall_rep[[data_to_use()]]})
    output$Barplot <- renderPlotly({
      #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      #reactive({get(load(input$datapath,.GlobalEnv))})
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      
      p <- plot_ly(data_bp(),source="source1", y = reorder(data_bp()$Kinase,data_bp()$Enrichment.score), x = data_bp()$Enrichment.score, type = 'bar', orientation = 'h', marker=list(color = data_bp()$Color)) %>%
        layout(showlegend = FALSE,xaxis=list(tickvals=data_bp()$Kinase, ticktext=data_bp()$Kinase))
      p
    })
    
    ## Couple hover plot for Barplot
    
    output$correlation2 <- renderPlotly({
      all_dbs <- read.csv(paste0("../results/",cell_line(),'_',P_value(),'P_',FDR_value(),'FDR/',"all_dbs.csv"))
      
      eventdata <- event_data("plotly_click", source = "source1")
      validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
      
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
      kin_subs <- unique(all_dbs[grep(kinase, all_dbs$termName),]$geneID)
      hover_data <- as.data.frame(exprs)[rownames(exprs) %in% kin_subs,]
      hover_data <- as.data.frame(t(hover_data))
      hover_data$Replicate <- gsub("X","",rownames(hover_data))
      hover_data <- melt(hover_data)
      hover_data$value <- 2^hover_data$value
      
      
      # pp <- ggplot(hover_data, aes(x = Replicate, y = value)) +
      #   geom_bar(stat="identity") +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   #ggtitle("Intensity levels of ABL1 (\"P00519\") phosphorylation sites in the PP") +
      #   theme(plot.title = element_text(size = 12))
      # pp <- ggplotly(pp)
      
      pl <- plot_ly(x = hover_data$Replicate, y=hover_data$value,name=hover_data$variable, type = "bar") %>%
        layout(title = paste0("Intensities of substrates of ",kinase,'_',FDR_value()), barmode = 'stack')
      pl
    })
    # output$Barplot <- renderPlot({
    #   ggplot(data=data_bp(), aes(x=reorder(Kinase,Enrichment.score), y=Enrichment.score, fill= Color)) +
    #     geom_bar(position="dodge",stat="identity", fill=data_bp()$Color) +
    #     coord_flip() +
    #     theme_bw() +
    #     theme(axis.text.x = element_text(angle = 90, hjust = 1, size=16)) +
    #     theme(axis.text.y = element_text(size=16)) +
    #     xlab("Kinase") +
    #     theme(#axis.text=element_text(size=14),
    #       axis.title=element_text(size=20))
    #   
    # })
    
    
    
    ## Volcano plots
    #reactive({get(load(input$datapath,.GlobalEnv))})
    data_vlc <-  reactive({my_data()$Volcano[[data_to_use()]]})
    
    output$Volcano <- renderPlotly({
      #data_v <- my_data()$Volcano
      
      x <- list(title = paste0("log2(FC)"))
      y <- list(title = "-log10(p-value)")
      pp <- plot_ly(type = 'scatter',source="source") %>%
        add_trace(
          x = data_vlc()$topTable, 
          y = -log10(data_vlc()$p),
          text = data_vlc()$prots,
          color = data_vlc()$threshold,
          hoverinfo = 'text',
          #marker = list(color='green'),
          showlegend = F) %>%
        layout(xaxis = x, yaxis = y)
      pp
      
      
      # p <- ggplot(data=data_vlc(),
      #        aes(x=topTable, y =-log10(p),colour=threshold, label=rownames(data_vlc()))) +
      #   geom_point(alpha=0.4, size=1.75) +
      #   xlab("log2 fold change") + ylab("-log10 p-value") +
      #   theme_bw() +
      #   theme(legend.position="none") +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
      #   theme(axis.text.y = element_text(size = 22)) +
      #   theme(axis.title=element_text(size=24))
      
      #ggplotly(p)
      
    })
    
    ## Coupler hover plot
    
    output$correlation <- renderPlotly({
      #get(load('../results/results_shiny_K562_0.05FDR_0.01P.Rda',.GlobalEnv))
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      exprs <- my_data()$ExpressSet
      
      eventdata <- event_data("plotly_hover", source = "source")
      validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
      
      datapoint <- as.numeric(eventdata$pointNumber)[1]+1
      
      site <- data_vlc()[datapoint,]$prots
      hover_data <- exprs[site,]
      hover_data <- melt(hover_data)
      hover_data$Replicate <- rownames(hover_data)
      hover_data$value <- 10^hover_data$value
      # pp <- ggplot(hover_data, aes(x = Replicate, y = value)) +
      #   geom_bar(stat="identity") +
      #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
      #   #ggtitle("Intensity levels of ABL1 (\"P00519\") phosphorylation sites in the PP") +
      #   theme(plot.title = element_text(size = 12))
      # pp <- ggplotly(pp)
      
      plot_ly(x = hover_data$Replicate, y=hover_data$value, type = "bar") %>%
        layout(title = paste0("Intensities of ",site))
      
    })
    
    ## Quality report
    ## Q1: histograms
    # Histogram sites
    #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    #data2 <- reactive ({data})
    output$hist_sites <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data1 <- my_data()$hist_sites
      
      ggplot(data1, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data1$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# sites") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    
    # Histogram peps
    #data2 <- reactive ({data})
    output$hist_peps <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data2 <- my_data()$hist_peps
      
      ggplot(data2, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data2$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# peps") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    
    # Histogram prots
    #data2 <- reactive ({data})
    output$hist_prots <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data3 <- my_data()$hist_prots
      
      ggplot(data3, aes(x=rep, y=prots, fill = Cell_line)) + 
        geom_bar(stat="identity") +
        geom_text(aes(label=data3$prots),vjust=-1, size=3) +
        xlab("Replicates") +
        ylab("# proteins") +
        theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
    })
    ## Q2: Venn diagrams
    # Venn sites
    #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    #data2 <- reactive ({data})
    output$venn_sites <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      #data_v1 <- shiny_results$venn_sites
      data_v1 <- my_data()$venn_sites
      
      grid.draw(venn.diagram(data_v1,
                             fill = 2:(1+length(data_v1)),
                             filename = NULL, 
                             alpha = rep(0.5,length(data_v1)),
                             cat.fontface = 4,
                             lty =2, cex = 1.1, 
                             scaled = TRUE, euler.d = TRUE,
                             main.cex = 1.2))
    })
    
    # Venn peps
    #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    #data2 <- reactive ({data})
    output$venn_peps <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data_v2 <- my_data()$venn_peps
      
      grid.draw(venn.diagram(data_v2,
                             fill = 2:(1+length(data_v2)),
                             filename = NULL, 
                             alpha = rep(0.5,length(data_v2)),
                             cat.fontface = 4,
                             lty =2, cex = 1.1, 
                             scaled = TRUE, euler.d = TRUE,
                             main.cex = 1.2))
    })
    
    # Venn proteins
    #load(paste0("./results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    #data2 <- reactive ({data})
    output$venn_prots <- renderPlot({
      #shiny_results <- get(load(input$file1$datapath,.GlobalEnv))
      
      data_v3 <- my_data()$venn_prots
      
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
  }
  shinyApp(ui, server)
}

#library(plotly)
# library(manhattanly)
# detach("package:manhattanly", unload=TRUE)
# # # 
# detach("package:plotly", unload=TRUE)
# library(dplyr)
# detach("package:dplyr", unload=TRUE)

shinyPathview2()

# library(rsconnect)
# rsconnect::setAccountInfo(name='mh-apps',
#                           token='FBA197B65BA6984F2E26F3FADE8F4956',
#                           secret='7AnJXNloBMKf7Zedt/a76q4ynK278AyFFWLq9vQl')
# deployApp()

