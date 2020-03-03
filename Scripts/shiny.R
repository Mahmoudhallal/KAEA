##################################################
## Project: Phosphoproteomics analysis of cell lines
## Script purpose: Map enriched kinases to selected pathways
## Date: 25.07.2018
## Author: Mahmoud Hallal
##################################################
## Load libraries
source("http://bioconductor.org/biocLite.R")
#biocLite("pathview")
#biocLite("GDCRNATools")

library(pathview)
#library(GDCRNATools)
library(shiny)
library(pheatmap)
#install.packages("filesstrings", repos="http://cran.rstudio.com/")
library(filesstrings)
#library(png)
library(ggplot2)
library(pheatmap)
data("gene.idtype.bods")
library(yaml)
## Load parameters
params <- read_yaml("./config.yaml")

## Load input file for SHINY
load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
material_for_waterfall <- shiny_results$waterfall_rep

to_sbmt <- lapply(1:length(material_for_waterfall), function(x){
  df1 <- material_for_waterfall[[x]][!duplicated(material_for_waterfall[[x]]$Kinase),]
  #make dataframe
  df <- as.data.frame(df1[,1])
  rownames(df) <- df1[,2]
  
  colnames(df) <- "value"
  ##
  dd <- df$value
  #names(dd) <- rownames(df)
  names(dd) <- gsub('(\\w+)[+-]','\\1',rownames(df))
  dd1 <- dd[!duplicated(names(dd))]
  dd2 <- as.data.frame(dd1)
  rownames(dd2) <- names(dd1)
  colnames(dd2) <- "value"
  dd2
})
names(to_sbmt) <- names(material_for_waterfall)


## Create Shiny function
shinyPathview2 <- function (datasets, pathways , directory) 
{
  if (!dir.exists(directory)) {
    dir.create(directory)
  }
  if (!endsWith(directory, "/")) {
    directory = paste(directory, "/", sep = "")
  }
  ui <- pageWithSidebar(headerPanel("KSEA"), 
                        sidebarPanel(
                          selectInput("DS", "Dataset", names(datasets)), 
                        
                          selectInput("xcol", "Pathway", pathways)), 

                        mainPanel(tabsetPanel(
                          tabPanel("Data Quality",
                                   tabsetPanel(
                                     navbarMenu("Sites",
                                                tabPanel("Histogram", plotOutput("hist_sites", width = "100%", height = "600px"),value=1),
                                                tabPanel("Venn Diagram", plotOutput("venn_sites", width = "100%", height = "600px"),value=1)),
                                     navbarMenu("Peptides",
                                                tabPanel("Histogram", plotOutput("hist_peps", width = "100%", height = "600px"),value=2),
                                                tabPanel("Venn Diagram", plotOutput("venn_peps", width = "100%", height = "600px"),value=2)),
                                     navbarMenu("Proteins",
                                                tabPanel("Histogram", plotOutput("hist_prots", width = "100%", height = "600px"),value=3),
                                                tabPanel("Venn Diagram", plotOutput("venn_prots", width = "100%", height = "600px"),value=3))
                                   )
                          
                          ),
                          
                          tabPanel("Enrichment",
                                   tabsetPanel(
                                     tabPanel("Heatmap2", plotOutput("Heatmap", width = "60%", height = "600px"),value=4),
                                     tabPanel("Barplot", plotOutput("Barplot", width = "80%", height = "600px"),value=5),
                                     tabPanel("Volcano plot", plotOutput("Volcano", width = "80%", height = "600px"),value=5)
                                     
                                   )
                                   
                          ),
                          
                          tabPanel("KEGG pathways",
                                   tabsetPanel(
                                      tabPanel("Table", dataTableOutput("mytable1"),value=6),
                                      tabPanel("KEGG pathway", plotOutput("plot1", width = "10%", height = "600px"),value=7)
                                   )
                                   
                          )

                          )))
  server <- function(input, output, session) {
    data_to_use <- reactive({data_to_use <- input$DS})
    
    ## Tbale
    output$mytable1 <- renderDataTable({
      load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      data_tb <- shiny_results$table_rep
      data_tb2 <- as.data.frame(data_tb[[data_to_use()]])
      data_tb2
    }, options = list(searching=FALSE, paging = FALSE))
    
    #KEGG pathways
    output$plot1 <- renderImage({
      pathwayID <- strsplit(input$xcol, "~", fixed = TRUE)[[1]][1]
      sp_directory <- paste0(directory,'/',data_to_use(),'/')
      outfile <- paste(sp_directory, pathwayID, ".pathview.png", sep = "")
      if (!file.exists(outfile)) {
        pathview(gene.data = datasets[[data_to_use()]], pathway.id = pathwayID,
                 species = "hsa", gene.idtype = "SYMBOL", limit = list(gene = max(abs(datasets[[data_to_use()]])), cpd = 1), 
                 kegg.dir = sp_directory, same.layer=T, kegg.native = T )
        file.move(paste(pathwayID, ".pathview.png", sep = ""), sp_directory, overwrite = TRUE)
        }
      list(src = outfile)
    }, deleteFile = FALSE)
    
    #heatmap
    load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    data_hm <- as.data.frame(shiny_results$Heatmap_rep)
    #palette

    # 
    data_hm2 <- reactive({as.data.frame(data_hm[,data_to_use(),drop=FALSE])})
    #print(colnames(data_hm2()))
    # rr <- reactive({sum(data_hm2()>=0)})
    # print(rr())
    # if (rr() >=0 ){
    #   print("OK")}
    # #

    output$Heatmap <- renderPlot({
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
    load(paste0(params$CWD,"/results/results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
    data <- shiny_results$waterfall_rep
    data_bp <- reactive ({data[[data_to_use()]]})
    output$Barplot <- renderPlot({
      ggplot(data=data_bp(), aes(x=reorder(Kinase,Enrichment.score), y=Enrichment.score, fill= Color)) +
        geom_bar(position="dodge",stat="identity", fill=data_bp()$Color) +
        coord_flip() +
        ggtitle("Kinase enrichments") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size=16)) +
        theme(axis.text.y = element_text(size=16)) +
        xlab("Kinase") +
        theme(#axis.text=element_text(size=14),
          axis.title=element_text(size=20))
      
    })
    
    ## Volcano plots
    data_v <- shiny_results$Volcano
    data_vlc <- reactive ({data_v[[data_to_use()]]})
   
    
    output$Volcano <- renderPlot({
      ggplot(data=data_vlc(),
             aes(x=topTable, y =-log10(p),colour=threshold, label=rownames(data_vlc()))) +
        geom_point(alpha=0.4, size=1.75) +
        xlab("log2 fold change") + ylab("-log10 p-value") +
        theme_bw() +
        theme(legend.position="none") +
        #ggtitle("Volcano plot",names(data)[[x]]) +
        theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
        theme(axis.text.y = element_text(size = 22)) +
        theme(axis.title=element_text(size=24))
      
      
      
    })
    
      ## Quality report
      ## Q1: histograms
      # Histogram sites
      load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      data1 <- shiny_results$hist_sites
      #data2 <- reactive ({data})
      output$hist_sites <- renderPlot({
        ggplot(data1, aes(x=rep, y=prots, fill = Cell_line)) + 
          geom_bar(stat="identity") +
          geom_text(aes(label=data1$prots),vjust=-1, size=3) +
          xlab("Replicates") +
          ylab("# sites") +
          labs(title="Phosphorylated sites") + 
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
          theme(axis.text.y = element_text(size = 22)) +
          theme(axis.title=element_text(size=24))
    })
      
      # Histogram peps
      data2 <- shiny_results$hist_peps
      #data2 <- reactive ({data})
      output$hist_peps <- renderPlot({
        ggplot(data2, aes(x=rep, y=prots, fill = Cell_line)) + 
          geom_bar(stat="identity") +
          geom_text(aes(label=data2$prots),vjust=-1, size=3) +
          xlab("Replicates") +
          ylab("# peps") +
          labs(title="Phosphorylated peps") + 
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
          theme(axis.text.y = element_text(size = 22)) +
          theme(axis.title=element_text(size=24))
      })
      
      # Histogram prots
      data3 <- shiny_results$hist_prots
      #data2 <- reactive ({data})
      output$hist_prots <- renderPlot({
        ggplot(data3, aes(x=rep, y=prots, fill = Cell_line)) + 
          geom_bar(stat="identity") +
          geom_text(aes(label=data3$prots),vjust=-1, size=3) +
          xlab("Replicates") +
          ylab("# proteins") +
          labs(title="Phosphorylated prots") + 
          theme_bw() + 
          theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 22)) +
          theme(axis.text.y = element_text(size = 22)) +
          theme(axis.title=element_text(size=24))
      })
      ## Q2: Venn diagrams
      # Venn sites
      load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      data_v1 <- shiny_results$venn_sites
      #data2 <- reactive ({data})
      output$venn_sites <- renderPlot({
        grid.draw(venn.diagram(data_v1,
                     fill = 2:1+length(data_v1),
                     filename = NULL, 
                     alpha = rep(0.5,length(data_v1)),
                     cat.fontface = 4,
                     lty =2, cex = 1.1, 
                     scaled = TRUE, euler.d = TRUE,
                     main = "Phosphorylated Peptides of K562 Ctrl and Drug", 
                     main.cex = 1.2))
      })
      
      # Venn peps
      load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      data_v2 <- shiny_results$venn_peps
      #data2 <- reactive ({data})
      output$venn_peps <- renderPlot({
        grid.draw(venn.diagram(data_v2,
                               fill = 2:1+length(data_v2),
                               filename = NULL, 
                               alpha = rep(0.5,length(data_v2)),
                               cat.fontface = 4,
                               lty =2, cex = 1.1, 
                               scaled = TRUE, euler.d = TRUE,
                               main = "Phosphorylated Peptides of K562 Ctrl and Drug", 
                               main.cex = 1.2))
      })
      
      # Venn proteins
      load(paste0(params$CWD,"/results/",params$cell_line,'_',params$pvalue_cutoff,'P_',params$fdr_cutoff,'FDR/',"results_shiny_",params$cell_line,"_",params$fdr_cutoff,"FDR_",params$pvalue_cutoff,"P.Rda"))
      data_v3 <- shiny_results$venn_prots
      #data2 <- reactive ({data})
      output$venn_prots <- renderPlot({
        grid.draw(venn.diagram(data_v3,
                               fill = 2:1+length(data_v3),
                               filename = NULL, 
                               alpha = rep(0.5,length(data_v3)),
                               cat.fontface = 4,
                               lty =2, cex = 1.1, 
                               scaled = TRUE, euler.d = TRUE,
                               main = "Phosphorylated Peptides of K562 Ctrl and Drug", 
                               main.cex = 1.2))
      })
  }
  shinyApp(ui, server)
}
#detach("package:dplyr", unload=TRUE)

shinyPathview2(to_sbmt, pathways
               = c('Autophagy animal','mTOR Signaling Pathway','Acute Myeloid Leukemia','Pathways in Cancer','RAS Signaling Pathway',
                   'MAPK Signaling Pathway','ERBB Signaling Pathway','WNT Signaling Pathway','TGF-BETA Signaling Pathway','JAK-STAT Signaling Pathway',
                   'PI3K-AKT Signaling Pathway','Chronic Myeloid Leukemia'), 
               directory = paste0(params$CWD,"/results/pathview/")
)

# library(rsconnect)
# rsconnect::setAccountInfo(name='mh-apps',
#                           token='FBA197B65BA6984F2E26F3FADE8F4956',
#                           secret='7AnJXNloBMKf7Zedt/a76q4ynK278AyFFWLq9vQl')
#deployApp()
