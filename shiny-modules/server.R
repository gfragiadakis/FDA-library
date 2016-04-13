library(shiny)
library(ggplot2)

compiled_df <- read.csv("~/Documents/FDAlibrary/saved_data_structures/human_signaling_fold_asinh_0.2.csv", row.names = 1)
df <- compiled_df
df <- dplyr::select(df, Donor, Condition_Feature, value)
df <- reshape2::dcast(df, Donor ~ Condition_Feature)
df <- data.frame(df[,-1], row.names = df$Donor)

correlations <- cor(df)

# function to cluster using correlation between variables as distance (add to FDA library)
reorder_cormat <- function(cormat){
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

cormat_reordered <- reorder_cormat(correlations)
adj <- make_adjacency(cormat_reordered, cor_threshold = 0.5)
adj_melted <- reshape2::melt(adj)
p <- ggplot(data = adj_melted, aes(Var1, Var2, fill = value)) + geom_tile() +
  scale_fill_gradient2(low = "red", high = "blue", mid = "white", midpoint = 0, limits=c(-1, 1)) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + coord_fixed(ratio = 1)

shinyServer(function(input, output) {

  df <- reactive({
    brushedPoints(adj_melted, input$plot_brush, xvar = "Var1", yvar = "Var2")
  })

  vals <- reactiveValues(pdata = p)

  output$info <- renderPrint({
    df()
  })
  observeEvent(input$save, {
      vals$pdata <- vals$pdata + annotate("rect", xmin = input$plot_brush$xmin, xmax = input$plot_brush$xmax,
                        ymin = input$plot_brush$ymin, ymax = input$plot_brush$ymax,
                        fill = "green", color = "black", alpha = 0.5 )
      print(paste("Module", input$save, "saved", sep = " "))
      # save module file
      write.csv(df(), file = paste("output/module_", input$save, ".csv", sep = ""))
      # print(input$plot_brush$xmin)
      # print(input$plot_brush$xmax)
      # print(input$plot_brush$ymin)
      # print(input$plot_brush$ymax)

    })
  observeEvent(vals$pdata, {
    output$plot <- renderPlot({
      isolate(vals$pdata)
    }, height = 1200, width = 1200)
  })
  })



