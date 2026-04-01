library("shiny")
library("bslib")
library("bsicons")
library("devtools")
library("ComplexHeatmap")
library("dplyr")
library("tidyverse")
library("RColorBrewer")
library("ggrepel")
library("edgeR")
library("limma")
library("factoextra")
library("FactoMineR")


ui <- page_sidebar(
  title = "My Differential Gene Expression App",
  sidebar = sidebar(helpText("Only .txt files accepted"),
                    fileInput("gexpr", 
                              label = "Upload the gene expression data below:",
                              accept = c(".txt")
                              ),
                    fileInput("long",
                              label = "Upload the reference gene length data below:",
                              accept = c(".txt")
                              ),
                    
                    actionButton("go", "Upload data"),
                    
                    selectInput("num", 
                                label = h3("How many genes would you like to visualize?"), 
                                choices = list(10, 50, 100, 250, 500, 1000),
                                selected = 100),
                    
                    width = 400
  ),
  
  mainPanel(
    tabsetPanel(
      tabPanel(
        "Summary",
        tableOutput("table")
      ),
      
      
      tabPanel(
        "Heatmap", 
        plotOutput("heatmap", height = "900px")
      ),
      
      tabPanel(
        "DGE Table",
        "Differential Gene Expression Table",
        tableOutput("dge_table")
      ),
      
      tabPanel(
        "Volcano Plot",
        plotOutput("volcano", height = "900px")
      ),
      
      tabPanel(
        "MA Plot",
        plotOutput("ma", height = "900px")
      ),
      
      tabPanel(
        "Long Gene", 
        plotOutput("longplot", height = "900px")
      )
    )
    
  )
)

server <- function(input, output, session) {
  
  expr_data <- eventReactive(input$go, {
    req(input$gexpr)
    req(input$long)
    
    withProgress(message = "Generating a table...",{
      for (i in 1:5){
        incProgress(1/5)
      }
    
    read.table(input$gexpr$datapath, header = TRUE)
    
    })
  })
  
  ## Gene Expression Data Table
  output$table <- renderTable({
    expr_data()
  })
  
  ## Rank and Normalize Most Variable Genes
  data_pipeline <- reactive({
    req(input$gexpr)
    req(input$num)
    numgenes <- as.numeric(input$num)
    
    df <- expr_data()
    
    # Identify numeric columns
    df1 <- sapply(df, is.numeric)
    
    # Apply Log-Transformation
    df[, df1] <- lapply(df[, df1], log2)
    
    # Find the variation of genes
    df$Variation <- NA
    df2 <- numeric(ncol(df)-2)
    for (r in 1:nrow(df)){
      for (c in 2:(ncol(df)-1)){
        if(is.infinite(df[r,c]) || is.na(df[r,c])){
          df2[c-1] <- 0
          df[r,c] <- 0 
        }
        else{
          df2[c-1] <- df[r,c]
        }
      }
      df$Variation[r] <- mad(df2)
    }
    
    # Rank the top most variable genes
    mostvar <- df %>% slice_max(order_by = Variation, n = numgenes)
    rownames(mostvar) <- mostvar[,1]
    mostvar <- subset(mostvar, select = -Gene)
    
    # Normalize, Z-scale function per gene across samples 
    mat <- as.matrix(subset(mostvar, select = -Variation))
    scalevar <- t(scale(t(mat)))
    
  })
  
  ## Generate a heatmap
  output$heatmap <- renderPlot({
    print(ComplexHeatmap::pheatmap(data_pipeline(),
                             fontsize = 10, 
                             heatmap_legend_param = list(
                               title = "Gene Expression", 
                               position = "center"
                             )
                            )
    )
  })
  
  ## Calculate Differential Gene Expression
  dge_data <- reactive({
    req(expr_data())
    
    # Import data and rename
    gexpr <- expr_data()
    rownames(gexpr) <- gexpr$Gene
    gexpr$Gene <- NULL
    
    # Filter data on lower count rate, FPKM < 1
    cutoff <- 1
    keep <- apply(gexpr, 1, max) >= cutoff
    gexpr_filtered <- gexpr[keep,]
    dim(gexpr_filtered) # number of genes left
    
    # Log-transform
    log2_gexpr <- log2(gexpr + 1)
    
    # Set sample names 
    snames <- colnames(gexpr)
    s1 <- substr(snames, 4, 7)
    s2 <- substr(snames, nchar(snames)-2, nchar(snames))
    group <- ifelse(substr(s1, 1, 1) == "S", "USP7wt", "USP7down")
    group <- factor(group, levels = c("USP7wt", "USP7down"))
    labels <- paste(s1, s2, sep = "_")
    
    ##--Data Normalization and Making Models--##
    mm <- model.matrix(~ 0 + group)
    voom.y.expr <- voom(gexpr_filtered, mm, plot = F)
    
    # Fit data into lm model
    fit <- lmFit(voom.y.expr, mm)
    coef.fit <- fit$coefficients
    
    ##--Establish sample group for DEGs analysis--##
    contr <- makeContrasts(groupUSP7wt - groupUSP7down, levels = colnames(coef(fit)))
    
    # Extract a table of the top-ranked genes from a linear model fit and run Empirical Bayes Statistics
    tmp <- contrasts.fit(fit, contr)
    tmp <- eBayes(tmp)
    top.table <- topTable(tmp, sort.by = "P", n = Inf)
    
    # Compare to coefficients
    coef_gexpr <- coef.fit[rownames(coef.fit) %in% rownames(top.table)[1:5],]
    
    # Filter based on p-value < 0.05
    gexpr_f <- top.table %>% arrange(logFC) %>% filter(adj.P.Val < 0.05)
    
    # Export gene list to txt
    top.table$Gene <- rownames(top.table)
    top.table <- top.table[,c("Gene", names(top.table)[1:6])]
    
    top.table
  })
  
  ## Differential Gene Expression Table
  output$dge_table <- renderTable({
    dge_data()
  })
  
  ## Prepare for Volcano Plot
  volcano_data <- reactive({
    req(dge_data())
    
    df <- dge_data()
    
    # Set threshold, t
    t = 0
    df$diffexpressed <- "NO"
    df$diffexpressed[df$logFC > t & df$adj.P.Val < 0.05] <- "UP"
    df$diffexpressed[df$logFC < -t & df$adj.P.Val < 0.05] <- "DOWN"
    
    # Highlight top 10 up/downregulated genes
    
    top <- head(df[order(df$adj.P.Val), "Gene"], 30)
    df$delabel <- ifelse(df$Gene %in% top, df$Gene, NA)
    
    return(df)
  })
  
  ## Generate a Volcano Plot
  output$volcano <- renderPlot({
    t <- 0
    
    ggplot(data = volcano_data(), aes(x = logFC, y = -log10(adj.P.Val), col = diffexpressed, label = delabel)) +
      geom_vline(xintercept = c(-t, t), col = 'gray', linetype = 'dashed') +
      geom_hline(yintercept = c(1.3), col = 'gray', linetype = 'dashed') +
      geom_point() +
      geom_point(data = volcano_data()[volcano_data()$Gene == "USP7", ],
                 color = "darkgoldenrod1", 
                 size = 4) + 
      scale_color_manual(values = c("turquoise3", "grey", "red2"), 
                         labels = c("Downregulated", "Not significant", "Upregulated")) +
      coord_cartesian(ylim = c(0, 10), xlim = c(-4, 4)) +
      scale_x_continuous(breaks = seq(-4, 4, 2)) +
      labs(color = "",
           x = expression("log"[2]*"FC"),
           y = expression("-log"[10]*"p-value")) +
      ggtitle("USP7 Wild-Type vs. Knockdown Jurkat Cells") +
      geom_text_repel(max.overlaps = Inf) +
      geom_label_repel(data = volcano_data()[volcano_data()$Gene == "USP7", ],
                       aes(label = Gene),
                       fill = "darkgoldenrod1",
                       color = "black")
  })
  
  ## Prepare MA plot data
  ma_data <- reactive({
    req(dge_data())
    
    df <- dge_data()
    
    # Create a data frame
    df$MA_DE <- abs(df$logFC) > 0
    
    # Highlight most differentially expressed genes
    genes_to_label <- df %>% 
      dplyr::filter(MA_DE == TRUE) %>%
      dplyr::arrange(desc(abs(df$logFC))) %>%
      dplyr::slice(1:30)
    
    # Add significance categories
    df$significance <- "NO"
    df$significance[df$adj.P.Val < 0.05] <- "YES"
    
    list(
      df = df, 
      genes_to_label = genes_to_label
    )
    
  })
  
  ## Generate an MA Plot
  output$ma <- renderPlot({
    df <- ma_data()$df
    genes_to_label <- ma_data()$genes_to_label
    
    ggplot(data = df, aes(x = df$AveExpr, y = df$logFC)) +
      geom_point()
    
    ggplot(data = df, aes(x = AveExpr, y = logFC, color = significance)) +
      geom_point() +
      geom_point(data = df[df$Gene == "USP7", ],
                 color = "dodgerblue", 
                 size = 2) + 
      scale_color_manual(values = c("NO" = "grey", "YES" = "skyblue")) + 
      geom_hline(yintercept = 0, color = "midnightblue", linetype = "dashed", linewidth = 1) +
      labs(title = "MA Plot: USP7 Wild-Type vs. Knockdown",
           x = "A (Average Expression)",
           y = expression(M~"(log"[2]*"FC)")) +
      ggrepel::geom_text_repel(data = genes_to_label, 
                               aes(label = Gene),
                               color = "skyblue",
                               max.overlaps = Inf, 
                               box.padding = 0.4, 
                               point.padding = 0.3,
                               segment.color = "grey") +
      ggrepel::geom_label_repel(data = df[df$Gene == "USP7", ],
                                aes(label = Gene),
                                fill = "dodgerblue",
                                color = "black", 
                                point.padding = 0.3) + 
      theme(legend.position = "top")
    
  })

  ## Prepare data for Long Gene Analysis
  long_data <- reactive({
    req(input$long)
    req(dge_data())

    len <- read.table(input$long$datapath, header = TRUE)
    df <- dge_data()

     # Find common genes and make a data frame
    df1 <- df %>%
      inner_join(len, by = c("Gene" = "GeneID"))

    df1$LogAvgLength <- log(df1$AvgLength)

    # Highlight most differentially expressed genes
    genes_to_label <- df1 %>%
      dplyr::arrange(desc(abs(logFC))) %>%
      dplyr::slice(1:30)

    # Add significance categories
    df1$significance <- "Not significant at alpha = 0.05"
    df1$significance[df1$adj.P.Val < 0.05] <- "Significant at alpha = 0.05"

    list(
      df1 = df1,
      genes_to_label = genes_to_label
    )
  })

  ## Plot Gene Length Plot with Log x-Axis Scale
  output$longplot <- renderPlot({
    req(long_data())

    df <- long_data()$df1
    genes_to_label <- long_data()$genes_to_label

    ggplot(data = df, aes(x = AvgLength, y = logFC, color = significance)) +
      geom_point(size = 1.2) +
      scale_x_log10() +
      geom_point(data = df[df$Gene == "USP7", , drop = FALSE],
                 color = "dodgerblue",
                 size = 2) +
      scale_color_manual(values = c("Not significant at alpha = 0.05" = "grey", "Significant at alpha = 0.05" = "skyblue")) +
      geom_hline(yintercept = 0, color = "midnightblue", linetype = "dashed", linewidth = 1) +
      labs(title = "Differential Expression and Gene Length: USP7 Wild-Type vs. Knockdown",
           x = "Gene Length",
           y = expression(M~"(log"[2]*"FC)")) +
      ggrepel::geom_text_repel(data = genes_to_label,
                               aes(label = Gene),
                               color = "skyblue",
                               max.overlaps = Inf,
                               box.padding = 0.4,
                               point.padding = 0.3,
                               segment.color = "grey") +
      ggrepel::geom_label_repel(data = df[df$Gene == "USP7", , drop = FALSE],
                                aes(label = "USP7"),
                                fill = "dodgerblue",
                                color = "black",
                                point.padding = 0.3) +
      theme(legend.position = "top")
  })
}

shinyApp(ui = ui, server = server)

