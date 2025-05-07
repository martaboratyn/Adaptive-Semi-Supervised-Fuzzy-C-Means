symulacja <- function(n,t_k,scenariusz=c(1,2)){
  symu <- data.frame(x1 = numeric(n), x2 = numeric(n),y=numeric(n))
  t<-1
  while (t <= n) {
    rand <- runif(1)
    if (t <= t_k) {
      if (rand < 1/2) {
        symu$x1[t] <- rnorm(1, mean = 0, sd = 0.5)
        symu$x2[t] <- rnorm(1, mean = 1, sd = 0.5)
        symu$y[t] <- '1'
      } else {
        symu$x1[t] <- rnorm(1, mean = 1, sd = 0.5)
        symu$x2[t] <- rnorm(1, mean = 2, sd = 0.5)
        symu$y[t] <- '2' 
      }
    } else {
      if (scenariusz==1){
        if (rand < 1/3) {
          symu$x1[t] <- rnorm(1, mean = 0, sd = 0.5)
          symu$x2[t] <- rnorm(1, mean = 1, sd = 0.5)
          symu$y[t] <- '1'
        } else if (rand >= 1/3 & rand < 2/3) {
          symu$x1[t] <- rnorm(1, mean = 1, sd = 0.5)
          symu$x2[t] <- rnorm(1, mean = 2, sd = 0.5)
          symu$y[t] <- '2' 
        } else {
          symu$x1[t] <- rnorm(1, mean = 1.5, sd = 0.5)
          symu$x2[t] <- rnorm(1, mean = 2.5, sd = 0.5)
          symu$y[t] <- '2' 
        }
      }
      if (scenariusz==2){
        if (rand < 1/3) {
          symu$x1[t] <- rnorm(1, mean = 0, sd = 1)
          symu$x2[t] <- rnorm(1, mean = 1, sd = 1)
          symu$y[t] <- '1'
        } else if (rand >= 1/3 & rand < 2/3) {
          symu$x1[t] <- rnorm(1, mean = 1, sd = 1)
          symu$x2[t] <- rnorm(1, mean = 2, sd = 1)
          symu$y[t] <- '2' 
        } else {
          symu$x1[t] <- rnorm(1, mean = 0.5, sd = 1)
          symu$x2[t] <- rnorm(1, mean = 1.5, sd = 1)
          symu$y[t] <- '2' 
        }
      }
    }
    t <- t + 1
  }
  return(symu)
}

chunki <- function(data,n,p,label_column){
  chunks_list <- list()  
  obserwacje_w_czesci <- ceiling(nrow(data) / n)
  indeksy <- rep(1:n, each = obserwacje_w_czesci)
  chunk <- split(data, indeksy)
  for ( i in 1:n){
    chunk1 <- chunk[[i]]
    chunk1$cgic_1 <- chunk1[,label_column]
    etyk <- chunk1[chunk1[,label_column]!='lack',]
    do_zmiany <- sample(nrow(etyk), nrow(etyk)*p)
    etyk$cgic_1 <- as.character(etyk$cgic_1)
    chunk1$cgic_1 <- as.character(chunk1$cgic_1)
    etyk$cgic_1[do_zmiany] <- 'lack'
    chunk1[chunk1[, label_column] != 'lack', 'cgic_1'] <- etyk$cgic_1
    chunks_list[[i]] <- chunk1
  }
  
  return(chunks_list)
}



szereg <- function(data, kolumna, y1) {
  unikalne_wartosci <- unique(data[,y1])
  kolory <- c("red","gray", "blue", "green", "orange", "purple")
  
  plot <- ggplot(data, aes(x = czas, y = !!rlang::sym(kolumna), color = data[,y1])) +
    geom_point(size = 2, shape = 19) +
    geom_line( color = "black") +
    
    theme_minimal() +
    labs(title = "Szereg czasowy", x = "Czas", y = "Zmienna", color = "Condition")
  
  return(plot)
}





model_kmeans <- function(chunk, k, selected_columns, x_limit = c(-1, 1), y_limit = c(-1, 1), chunk_number, label,new_label) {
  X <- chunk[, selected_columns]
  model_klasyczny <- kmeans(X, centers = k)
  colors <- as.factor(model_klasyczny$cluster)
  
  if (label == 'euthymia') {
    colory <- ifelse(colors == '1', 'disease', 'euthymia') 
    palette <- c("disease" = "blue", "euthymia" = "green")
    shapes <- c("disease" = 16, "euthymia" = 17, "lack"=3)
  } else {
    present_colors <- unique(colors)
    colory <- colors
    all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
    palette <- all_colors[present_colors]
    
    shapes <- setNames(c(16, 17, 18, 19, 20), c("1", "2", "3", "4", "5"))
    shapes["lack"] <- 3 
    shapes <- shapes[unique(chunk[[new_label]])]
  }
  
  data_plot <- as.data.frame(model_klasyczny$centers)
  colnames(data_plot) <- c("V1", "V2")
  
  plot <- ggplot(chunk, aes(x = chunk[, selected_columns[1]], y = chunk[, selected_columns[2]], 
                            shape = as.factor(chunk[[new_label]]), color = colory)) +
    geom_point(size = 3, stroke = 0.5) +
    labs(title = paste(chunk_number), 
         x = selected_columns[1], 
         y = selected_columns[2],
         color = "Predicted", 
         shape = "Observed") +
    geom_point(data = data_plot, aes(x = V1, y = V2), size = 6, 
               color = "red", shape = 18, alpha = 0.8) +
    xlim(x_limit) +  
    ylim(y_limit) +  
    scale_color_manual(values = palette) +
    scale_shape_manual(values = shapes) + 
    theme_minimal() +
    theme(
      legend.text = element_text(size = 22),        
      legend.title = element_text(size = 24),       
      plot.title = element_text(size = 26, face = "bold"),  
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 22),       
      axis.title.y = element_text(size = 22)           
    )  +
    guides(color = guide_legend(override.aes = list(shape = 15)))
  
  model_klasyczny$plot <- plot
  model_klasyczny$colors <- colors
  model_klasyczny$V <- model_klasyczny$centers
  return(model_klasyczny)
}


model_dissfcm <- function(chunks, C, F_s, selected_columns, label_column, new_label, x_limit = c(-1, 1), y_limit = c(-1, 1)) {
  model <- DISSFCM(chunks = chunks, C = C, F_s = F_s, selected_columns = selected_columns)
  plots <- list()
  colors <- list()
  reconstruction_error <- model$reconstruction_error
  
  for (i in 1:length(chunks)) {
    chunk <- chunks[[i]]
    color <- as.factor(model$U[[i]])
    str(chunk)
    
    if (label_column == 'euthymia') {
      color <- ifelse(color == '1', 'disease', 'euthymia')
      palette <- c("disease" = "blue", "euthymia" = "green")
      shapes <- c("disease" = 16, "euthymia" = 17, "lack"=3)
    } else {
      present_colors <- unique(color)
      present_colors <- as.character(present_colors)
      color <- color
      all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
      palette <- all_colors[present_colors]
      
      shapes <- setNames(c(16, 17, 18, 19, 20, 3), c("1", "2", "3", "4", "5", "lack"))
      shapes <- shapes[unique(chunk[[new_label]])]
      chunk[[new_label]] <- as.factor(chunk[[new_label]])
    }
    
    data_plot <- as.data.frame(model$V[[i]])
    colnames(data_plot) <- c("V1", "V2")
    length(chunk[,selected_columns[1]])
    length(chunk[,selected_columns[2]])
    plot <- ggplot(chunk, aes(x = chunk[,selected_columns[1]], y = chunk[,selected_columns[2]], 
                              shape = chunk[,new_label], color = color)) +
      geom_point(size = 3, stroke = 0.5) +
      labs(title = paste("Chunk", i), 
           x = selected_columns[1], 
           y = selected_columns[2],
           color = "Predicted", 
           shape = "Observed") +
      geom_point(data = data_plot, aes(x = V1, y = V2), size = 6, color='red',shape = 18, alpha = 0.8) +
      xlim(x_limit) +  
      ylim(y_limit) +  
      scale_color_manual(values = palette) +
      scale_shape_manual(values = shapes) + 
      theme_minimal() +
      theme(
        legend.text = element_text(size = 22),        
        legend.title = element_text(size = 24),       
        plot.title = element_text(size = 26, face = "bold"),  
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),       
        axis.title.y = element_text(size = 22)           
      )         +
      guides(color = guide_legend(override.aes = list(shape = 15)))
    
    plot_gtable <- ggplot_gtable(ggplot_build(plot))
    plots[[i]] <- plot_gtable
    colors[[i]] <- color
  }
  
  model$plots <- plots
  model$colors <- colors
  model$reconstruction_error <- reconstruction_error
  model$F_s <- F_s
  model$V <- model$V
  return(model)
}




model_ssfcm <- function(chunk, C, F_, selected_columns, label_column, new_label, x_limit=c(-1,1), y_limit=c(-1,1), chunk_number) {
  X <- chunk[, selected_columns]
  model <- ssfclust::SSFCM(X = X, C = C, alpha = 1, F_ = F_)
  colors <- as.factor(max.col(model$U))
  if (label_column == 'euthymia') {
    colory <- ifelse(colors == '1', 'disease', 'euthymia') 
    palette <- c("disease" = "blue", "euthymia" = "green")
    shapes <- c("disease" = 16, "euthymia" = 17, "lack"=3)  
  } else {
    present_colors <- unique(colors)
    colory <- colors
    all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
    palette <- all_colors[present_colors]
    
    shapes <- setNames(c(16, 17, 18, 19, 20), c("1", "2", "3", "4", "5"))
    shapes["lack"] <- 3 
    shapes <- shapes[unique(chunk[[new_label]])]
  }
  
  chunk[[label_column]] <- as.factor(chunk[[label_column]])
  
  data_plot <- as.data.frame(model$V)
  colnames(data_plot) <- c("V1", "V2")
  
  plot <- ggplot(chunk, aes(x = chunk[, selected_columns[1]], y = chunk[, selected_columns[2]], 
                            shape = chunk[[new_label]], color = colory)) +
    geom_point(size = 3, stroke = 0.5) +
    labs(title = paste( chunk_number), 
         x = selected_columns[1], 
         y = selected_columns[2],
         color = "Predicted", 
         shape = "Observed") +
    geom_point(data = data_plot, aes(x = V1, y = V2), size = 6, 
               color = "red", shape = 18, alpha = 0.8) +
    xlim(x_limit) +  
    ylim(y_limit) +  
    scale_color_manual(values = palette) +
    scale_shape_manual(values = shapes) + 
    theme_minimal() +
    theme(
      legend.text = element_text(size = 22),        
      legend.title = element_text(size = 24),       
      plot.title = element_text(size = 26, face = "bold"),  
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 22),       
      axis.title.y = element_text(size = 22)           
    )  +
    guides(color = guide_legend(override.aes = list(shape = 15)))
  
  model$plot <- plot
  model$colors <- colors
  model$V <- model$V
  
  return(model)
}

model_assfcm<- function(chunks, C, F_s, selected_columns, label_column, new_label, threshold, x_limit = c(-1, 1), y_limit = c(-1, 1)) {
  model <- assfcm(chunks = chunks, C = C, F_s = F_s, selected_columns = selected_columns, threshold = threshold)
  plots <- list()
  colors <- list()
  reconstruction_error <- model$reconstruction_error
  
  for (i in 1:length(chunks)) {
    chunk <- chunks[[i]]
    color <- as.factor(model$U[[i]])
    
    
    if (label_column == 'euthymia') {
      color <- ifelse(color == '1', 'disease', 'euthymia')
      palette <- c("disease" = "blue", "euthymia" = "green")
      shapes <- c("disease" = 16, "euthymia" = 17, "lack"=3)
    } else {
      present_colors <- unique(color)
      present_colors <- as.character(present_colors)
      color <- color
      all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
      palette <- all_colors[present_colors]
      
      shapes <- setNames(c(16, 17, 18, 19, 20, 3), c("1", "2", "3", "4", "5", "lack"))
      shapes <- shapes[unique(chunk[[new_label]])]
      chunk[[label_column]] <- as.factor(chunk[[label_column]])
    }
    chunk[[label_column]] <- as.factor(chunk[[label_column]])
    data_plot <- as.data.frame(model$V[[i]])
    colnames(data_plot) <- c("V1", "V2")
    
    plot <- ggplot(chunk, aes(x = chunk[, selected_columns[1]], y = chunk[, selected_columns[2]], 
                              shape = chunk[[new_label]], color = color)) +
      geom_point(size = 3, stroke = 0.5) +
      labs(title = paste("Chunk", i), 
           x = selected_columns[1], 
           y = selected_columns[2],
           color = "Predicted", 
           shape = "Observed") +
      geom_point(data = data_plot, aes(x = V1, y = V2), size = 6, shape = 18, color='red',alpha = 0.8) +
      xlim(x_limit) +  
      ylim(y_limit) +  
      scale_color_manual(values = palette) +
      scale_shape_manual(values = shapes) + 
      theme_minimal() +
      theme(
        legend.text = element_text(size = 22),        
        legend.title = element_text(size = 24),       
        plot.title = element_text(size = 26, face = "bold"),  
        axis.text.x = element_text(size = 20),
        axis.text.y = element_text(size = 20),
        axis.title.x = element_text(size = 22),       
        axis.title.y = element_text(size = 22)           
      ) +
      guides(color = guide_legend(override.aes = list(shape = 15)))
    
    plot_gtable <- ggplot_gtable(ggplot_build(plot))
    plots[[i]] <- plot_gtable
    colors[[i]] <- color
  }
  model$reconstruction_error <- reconstruction_error
  model$plots <- plots
  model$colors <- colors
  model$V <- model$V
  model$F_s <- F_s
  return(model)
}





acf_fun <- function(series, ...) {
  lapply(series, function(x) { as.numeric(acf(x, lag.max = 50L,plot=FALSE)$acf) })
}



model_szeregi <- function(chunk, k, columns, col) {
  X <- chunk
  model_szereg <- tsclust(X[, columns], type = "fuzzy", k = k,
                          preproc = acf_fun, distance = "L2",
                          constrained = list(labels = X$eu))
  colors <- model_szereg@cluster
  colors <- ifelse(colors=='1', 'disease','euthymia') 
  plot1<- ggplot(data = X, aes(x = czas, y = .data[[col]])) +
    geom_line() +
    geom_point(data = subset(X, colors == 'disease'), aes(color = "Cluster 1"), shape = 19) +
    geom_point(data = subset(X, colors == 'euthymia'), aes(color = "Cluster 2"), shape = 19) +
    scale_color_manual(values = c("red", "blue")) +
    labs(color = "Cluster") +
    theme_minimal()
  
  return(list(model_szereg = model_szereg, plot = plot1, colors = colors))
  
}


model_dissfcm_modyfikacja <- function(chunks, C, F_s, selected_columns, label_column, threshold, x_limit = c(-1, 1), y_limit = c(-1, 1)) {
  model <- DISSFCM_modyfikacja(chunks = chunks, C = C, F_s = F_s, selected_columns = selected_columns, threshold = threshold)
  plots <- list()
  colors <- list()
  max_distance <- model$max_distance
  
  for (i in 1:length(chunks)) {
    chunk <- chunks[[i]]
    color <- as.factor(model$U[[i]])
    chunk[[label_column]] <- as.factor(chunk[[label_column]])
    
    if (label_column == 'euthymia') {
      color <- ifelse(color == '1', 'disease', 'euthymia')
      palette <- c("disease" = "blue", "euthymia" = "green")
      shapes <- c("disease" = 16, "euthymia" = 17, "bral"=3)
    } else {
      present_colors <- unique(color)
      present_colors <- as.character(present_colors)
      all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
      palette <- all_colors[present_colors]
      
      shapes <- setNames(c(16, 17, 18, 19, 20), c("1", "2", "3", "4", "5"))
      shapes["lack"] <- 3  
      shapes <- shapes[unique(chunk[[label_column]])]
    }
    
    data_plot <- as.data.frame(model$V[[i]])
    colnames(data_plot) <- c("V1", "V2")
    
    plot <- ggplot(chunk, aes(x = chunk[, selected_columns[1]], y = chunk[, selected_columns[2]], 
                              shape = chunk[[new_label]], color = color)) +
      geom_point(size = 3, stroke = 0.5) +
      labs(title = paste("Chunk", i), 
           x = selected_columns[1], 
           y = selected_columns[2],
           color = "Predicted", 
           shape = "Observed") +
      geom_point(data = data_plot, aes(x = V1, y = V2), size = 6, 
                 color = "red", shape = 18, alpha = 0.8) +
      xlim(x_limit) +  
      ylim(y_limit) +  
      scale_color_manual(values = palette) +
      scale_shape_manual(values = shapes) + 
      theme_minimal() +
      theme(
        legend.text = element_text(size = 14),        
        legend.title = element_text(size = 16),       
        plot.title = element_text(size = 18, face = "bold"),  
        axis.title.x = element_text(size = 16),       
        axis.title.y = element_text(size = 16)        
      ) +
      guides(color = guide_legend(override.aes = list(shape = 15)))
    
    plot_gtable <- ggplot_gtable(ggplot_build(plot))
    plots[[i]] <- plot_gtable
    colors[[i]] <- color
  }
  
  model$plots <- plots
  model$colors <- colors
  model$V <- model$V
  model$max_distance <- max_distance
  model$F_s <- F_s
  return(model)
}



macierz_pomylek<-function(w,u){
  m <- matrix(0,nrow=2,ncol=2)
  m[1,1] <- length(which((w=='1' | w=='disease')  & (u =='1' | u =='disease')))
  m[2,1] <- length(which((w=='1'| w=='disease')  & (u =='2' | u =='euthymia'))) 
  m[1,2] <- length(which((w=='2' | w == 'euthymia' ) & (u =='1' | u == 'disease')))
  m[2,2] <- length(which((w=='2' | w=='euthymia') & (u =='2' | u == 'euthymia' ))) 
  return(m)
}

miary <- function(confusion_matrix, miara = c('accuracy', 'recall', 'precision', 'F1')) {
  TP <- confusion_matrix[1, 1]
  FP <- confusion_matrix[1, 2]
  TN <- confusion_matrix[2, 2]
  FN <- confusion_matrix[2, 1]
  total <- sum(confusion_matrix)
  
  if (total == 0) {
    warning("Pusta macierz konfuzji — brak danych")
    return(NA)
  }
  
  if (miara == 'accuracy') {
    return((TP + TN) / total)
  }
  
  if (miara == 'recall') {
    if ((TP + FN) == 0) {
      warning("Brak rzeczywistych pozytywnych przypadków — recall nieokreślony")
      return(NA)
    }
    return(TP / (TP + FN))
  }
  
  if (miara == 'precision') {
    if ((TP + FP) == 0) {
      warning("Brak przewidzianych pozytywnych przypadków — precision nieokreślona")
      return(NA)
    }
    return(TP / (TP + FP))
  }
  
  if (miara == 'F1') {
    if ((TP + FP) == 0 || (TP + FN) == 0) {
      warning("Nie można obliczyć F1 — precision lub recall nieokreślone")
      return(NA)
    }
    precision <- TP / (TP + FP)
    recall <- TP / (TP + FN)
    return(2 * (precision * recall) / (precision + recall))
  }
  
  warning("Nieznana miara: ", miara)
  return(NA)
}


models <- function(data, n, p, C, method, selected_columns, label_column, new_label, plot = FALSE, x_limit, y_limit, iteration) {
  results <- character()
  chunk_accuracies <- list()
  chunk_precision <- list()
  chunk_recall <- list()
  chunk_f1 <- list()
  chunk_list <- chunki(data, n, p, label_column)
  F_s <- vector("list", length(chunk_list))
  chunks <- list()
  all_unique_labels <- character()
  v_matrices <- list()  # Lista do przechowywania macierzy 'v'
  
  for (chunk_number in 1:n) {
    chunk <- chunk_list[[chunk_number]]
    chunks[[chunk_number]] <- chunk
    etykieta <- which(chunk[,new_label] != chunk[[label_column]])
    true_labels <- chunk[[label_column]]
    chunk_accuracy <- c()
    chunk_unique_labels <- sort(unique(chunk[[label_column]][chunk[[label_column]] != "lack"]))
    all_unique_labels <- unique(c(all_unique_labels, chunk_unique_labels))
    new_C <- length(all_unique_labels)
    
    if (method == "ssfcm") {
      if (label_column == 'euthymia') {
        if (any(chunk[[new_label]] %in% c('disease', 'euthymia'))) {
          F_s[[chunk_number]] <- matrix(0, nrow = nrow(chunk), ncol = C)
          F_s[[chunk_number]][which(chunk[,new_label] == 'disease'), 1] <- 1
          F_s[[chunk_number]][which(chunk[,new_label] == 'euthymia'), 2] <- 1
        }
        else if (label_column == 'cgic') {
            #F <- matrix(0, nrow = nrow(chunks[[i]]), ncol = length(all_unique_labels))
          F_s[[chunk_number]] <- matrix(0, nrow = nrow(chunk), ncol = new_C)
          for (l in 1:new_C) {
            label <- all_unique_labels[l]
            F_s[[chunk_number]][chunk[[new_label]] == label, l] <- 1
          }
        }
      } else if (label_column == 'y') {
        F_s[[chunk_number]] <- matrix(0, nrow = nrow(chunk), ncol = new_C)
        for (l in 1:new_C) {
          label <- all_unique_labels[l]
          F_s[[chunk_number]][chunk[[new_label]] == label, l] <- 1
        }
      } 
    }
    F_s[[chunk_number]]
    for (i in 1:(length(selected_columns) - 1)) {
      for (j in (i + 1):length(selected_columns)) {
        col1 <- selected_columns[i]
        col2 <- selected_columns[j]
        column_pair <- paste("columns:", col1, "and", col2)
        chunk_number1 <- paste("Chunk", chunk_number)
        if (method == "ssfcm") {
          modelik <- model_ssfcm(chunk, C = ncol(F_s[[chunk_number]]), 
                                 F_ = F_s[[chunk_number]], 
                                 selected_columns = c(col1, col2),
                                 label_column=label_column, new_label=new_label,
                                 chunk_number = chunk_number1, x_limit,y_limit)
          if (plot) {
            # Zapis wykresu
            file_name <- paste0("ssfcm_chunk_", chunk_number, "_columns_", col1, "_", col2, new_label, iteration, ".png")
            ggsave(file_name, plot = modelik$plot, bg = "white")
          }
          v_matrices[[paste0("chunk_", chunk_number, "_columns_", col1, "_", col2)]] <- modelik$V
          predicted_labels <- modelik$colors
        } 
        else if (method == "kmeans") {
          model <- model_kmeans(chunk = chunk, k = new_C, selected_columns = c(col1, col2), chunk_number = chunk_number1, label = label_column,new_label=new_label,x_limit,y_limit)
          if (plot && !is.null(model$plot)) {
            # Zapis wykresu
            file_name <- paste0("plot_kmeans_chunk_", chunk_number, "_columns_", col1, "_", col2, new_label, iteration,".png")
            ggsave(file_name, plot = model$plot, bg = "white")
          }
          v_matrices[[paste0("chunk_", chunk_number, "_columns_", col1, "_", col2)]] <- model$V
          predicted_labels <- model$colors
        } else {
          return(chunk_list)
        }
        if (label_column=='euthymia'){
          confusion_matrix <- macierz_pomylek(true_labels[etykieta], predicted_labels[etykieta])
          accuracy <- miary(confusion_matrix, miara = 'accuracy')
          precision <- miary(confusion_matrix, miara = 'precision')
          recall <- miary(confusion_matrix, miara = 'recall')
          f1_score <- miary(confusion_matrix, miara = 'F1')
        }
        else if (label_column=='y' | label_column=='cgic'){
          confusion_matrix <- table(Predicted = as.character(predicted_labels[etykieta]), Actual = as.character(true_labels[etykieta]))
          confusion_matrix[is.na(confusion_matrix)] <- 0
          if (nrow(as.matrix(confusion_matrix)) < 2 | ncol(as.matrix(confusion_matrix)) < 2){
            accuracy <- NA
            precision <- NA
            recall <- NA
            f1_score <- NA
          } else {
            accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
            TP <- diag(confusion_matrix)
            FP <- colSums(confusion_matrix) - TP
            FN <- rowSums(confusion_matrix) - TP  
            TP_total <- sum(TP)
            FP_total <- sum(FP)
            FN_total <- sum(FN)
            precision <- TP_total / (TP_total + FP_total)
            recall <- TP_total / (TP_total + FN_total)
            f1_score <- ifelse((precision + recall) == 0, 0, 
                               2 * (precision * recall) / (precision + recall))
          }
        }
        chunk_accuracy <- c(chunk_accuracy, accuracy)
        chunk_precision <- c(chunk_precision, precision)  # Średnia dla wieloklasowego
        chunk_recall <- c(chunk_recall, recall)
        chunk_f1 <- c(chunk_f1, f1_score)
        result <- paste("chunk:", chunk_number, ",", column_pair, ", accuracy:", accuracy)
        results <- c(results, result)
        #print(paste("Result added:", result)) 
      }
    }
    chunk_accuracies[[chunk_number]] <- chunk_accuracy
    chunk_precision[[chunk_number]] <- chunk_precision
    chunk_recall[[chunk_number]] <- chunk_recall
    chunk_f1[[chunk_number]] <- chunk_f1
  }
  accuracy_matrix <- do.call(rbind, chunk_accuracies)
  
  
  all_accuracies <- unlist(chunk_accuracies)
  all_precision <- unlist(chunk_precision)
  all_recall <- unlist(chunk_recall)
  all_f1 <- unlist(chunk_f1)
  mean_accuracy <- mean(all_accuracies, na.rm = TRUE)
  std_accuracy <- sd(all_accuracies, na.rm = TRUE)
  min_accuracy <- min(all_accuracies, na.rm = TRUE)
  max_accuracy <- max(all_accuracies, na.rm = TRUE)
  mean_precision <- mean(all_precision, na.rm = TRUE)
  mean_recall <- mean(all_recall, na.rm = TRUE)
  mean_f1 <- mean(all_f1, na.rm = TRUE)
  models <- list()
  models$mean_accuracy <- mean_accuracy
  models$std_accuracy <- std_accuracy 
  models$min_accuracy <- min_accuracy
  models$max_accuracy <- max_accuracy
  models$all_accuracies <- all_accuracies
  models$Fs <- F_s
  models$chunki <- chunks
  models$v_matrices <- v_matrices  # Dodanie macierzy 'v'
  models$accuracy <- c(results, mean_accuracy)
  models$precision <- c(results, mean_precision)
  models$recall <- c(results, mean_recall)
  models$f1 <- c(results, mean_f1)
  return(models)
}


models_dynamic <- function(data, n, p, C, method, selected_columns, label_column, new_label, threshold, plot = FALSE, x_limit, y_limit,iteration) {
  results <- character()
  chunks <- list()
  etykiets <- list()
  true_labels <- list()
  F_s <- list()
  v_matrices <- list()  # Lista do przechowywania macierzy 'v'
  all_plots <- list()  # Lista do przechowywania wykresów
  reconstruction_error_list <- list()
  chunks <- chunki(data, n, p, label_column)
  num_columns <- length(selected_columns)
  accuracy_s <- lapply(1:length(chunks), function(x) matrix(0, num_columns - 1, num_columns))
  precision_s <- lapply(1:length(chunks), function(x) matrix(0, num_columns - 1, num_columns))
  recall_s <- lapply(1:length(chunks), function(x) matrix(0, num_columns - 1, num_columns))
  f1_s <- lapply(1:length(chunks), function(x) matrix(0, num_columns - 1, num_columns))
  all_unique_labels <- character()
  accuracy <- matrix(0, nrow = length(selected_columns), ncol = length(selected_columns))
  
  for (i in 1:length(chunks)) {
    etykiets[[i]] <- which(chunks[[i]][, new_label] != chunks[[i]][[label_column]])
    true_labels[[i]] <- chunks[[i]][[label_column]]
    chunk_unique_labels <- as.character(sort(unique(chunks[[i]][[label_column]][chunks[[i]][[label_column]] != "lack"])))
    all_unique_labels <- unique(c(all_unique_labels, chunk_unique_labels))
    
    if (label_column == 'euthymia') {
      if (any(chunks[[i]][[new_label]] %in% c('disease', 'euthymia'))) {
        F <- matrix(0, nrow = nrow(chunks[[i]]), ncol = length(all_unique_labels))
        F[which(chunks[[i]][, new_label] == 'disease'), 1] <- 1
        F[which(chunks[[i]][, new_label] == 'euthymia'), 2] <- 1
        true_labels[[i]] <- ifelse(true_labels[[i]] == 'disease', 1, 2)
      }
    else if (label_column == 'cgic') {
      F <- matrix(0, nrow = nrow(chunks[[i]]), ncol = length(all_unique_labels))
      for (l in 1:length(all_unique_labels)) {
        label <- all_unique_labels[l]
        F[chunks[[i]][[new_label]] == label, l] <- 1
        }  
      }
    } else if (label_column == 'y') {
      F <- matrix(0, nrow = nrow(chunks[[i]]), ncol = length(all_unique_labels))
      for (l in 1:length(all_unique_labels)) {
        label <- all_unique_labels[l]
        F[chunks[[i]][[new_label]] == label, l] <- 1
      }
    }
    F_s[[i]] <- F
  }
  
  for (i in 1:(length(selected_columns) - 1)) {
    for (j in (i + 1):length(selected_columns)) {
      col1 <- selected_columns[i]
      col2 <- selected_columns[j]
      column_pair <- paste("columns:", col1, "and", col2)
      if (method == "dissfcm") {
        modelik <- model_dissfcm(chunks, C = C, F_s = F_s, selected_columns = c(col1, col2), label_column = label_column, new_label=new_label,x_limit, y_limit)
        predicted_labels <- modelik$colors
        v_matrices[[paste0("cols_", col1, "_", col2)]] <- modelik$V
        reconstruction_error_list[[paste0("cols_", col1, "_", col2)]] <- modelik$reconstruction_error
        all_plots[[paste0("cols_", col1, "_", col2)]] <- modelik$plots  # Dodanie listy wykresów do wyników
      } else if (method == 'assfcm') {
        modelik <- model_assfcm(chunks, C = C, F_s = F_s, selected_columns = c(col1, col2), label_column = label_column, new_label=new_label, threshold, x_limit, y_limit)
        predicted_labels <- modelik$colors
        v_matrices[[paste0("cols_", col1, "_", col2)]] <- modelik$V
        reconstruction_error_list[[paste0("cols_", col1, "_", col2)]] <- modelik$reconstruction_error
        all_plots[[paste0("cols_", col1, "_", col2)]] <- modelik$plots  # Dodanie listy wykresów do wyników
      }
      
      if (plot) {
        # Zapisywanie wszystkich wykresów w modelik$plots
        for (l in seq_along(all_plots[[paste0("cols_", col1, "_", col2)]])) {
          plot_name <- paste0(deparse(substitute(data)),'_', method, "_chunk_", l, "_cols_", col1, "_", col2, new_label,"_", iteration,".png")
          ggsave(plot_name, plot = all_plots[[paste0("cols_", col1, "_", col2)]][[l]], bg = "white")
        }
      }
      
      for (l in 1:length(chunks)) {
        if (label_column == 'euthymia') {
          confusion_matrix <- macierz_pomylek(true_labels[[l]][etykiets[[l]]], predicted_labels[[l]][etykiets[[l]]])
          accuracy_s[[l]][i, j] <- miary(confusion_matrix, miara = 'accuracy')
          precision_s[[l]][i, j] <- miary(confusion_matrix, miara = 'precision')
          recall_s[[l]][i, j] <- miary(confusion_matrix, miara = 'recall')
          f1_s[[l]][i, j] <- miary(confusion_matrix, miara = 'F1')
        } else if (label_column == 'y' | label_column == 'cgic') {
          confusion_matrix <- table(Predicted = predicted_labels[[l]][etykiets[[l]]], Actual = true_labels[[l]][etykiets[[l]]])
          confusion_matrix <- confusion_matrix[rowSums(confusion_matrix) > 0, colSums(confusion_matrix) > 0]
          #print(confusion_matrix)
          #print(str(confusion_matrix))
          if (nrow(as.matrix(confusion_matrix)) < 2 | ncol(as.matrix(confusion_matrix)) < 2) {
            accuracy_s[[l]][i, j] <- sum(as.character(true_labels[[l]][etykiets[[l]]]) == as.character(predicted_labels[[l]][etykiets[[l]]])) / length(as.character(true_labels[[l]][etykiets[[l]]]))
            precision_s[[l]][i, j] <- NA
            recall_s[[l]][i, j] <- NA
            f1_s[[l]][i, j] <- NA
            } else {
            accuracy_s[[l]][i, j] <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
            TP <- diag(confusion_matrix)
            FP <- colSums(confusion_matrix) - TP
            FN <- rowSums(confusion_matrix) - TP
            precision_s[[l]][i, j] <- mean(TP / (TP + FP))
            recall_s[[l]][i, j] <- mean(TP / (TP + FN))
            f1_s[[l]][i, j] <- 2 * (precision_s[[l]][i, j] * recall_s[[l]][i, j]) / (precision_s[[l]][i, j] + recall_s[[l]][i, j])
          }
        }
      }
    }
  }
  
  all_accuracies <- unlist(lapply(accuracy_s, function(x) na.omit(as.vector(x))))
  all_accuracies <- all_accuracies[all_accuracies != 0]
  mean_accuracy_all <- mean(all_accuracies)
  std_accuracy_all <- sd(all_accuracies)
  min_accuracy_all <- min(all_accuracies)
  max_accuracy_all <- max(all_accuracies)
  
  all_precisions <- unlist(lapply(precision_s, function(x) na.omit(as.vector(x))))
  all_recalls <- unlist(lapply(recall_s, function(x) na.omit(as.vector(x))))
  all_f1s <- unlist(lapply(f1_s, function(x) na.omit(as.vector(x))))
  
  mean_precision_all <- mean(all_precisions, na.rm = TRUE)
  mean_recall_all <- mean(all_recalls, na.rm = TRUE)
  mean_f1_all <- mean(all_f1s, na.rm = TRUE)
  
  models <- list()
  #models$F <- F_s
  #models$chunki <- chunks
  models$accuracy <- accuracy_s
  models$precision <- precision_s
  models$recall <- recall_s
  models$f1 <- f1_s
  models$mean_accuracy <- mean_accuracy_all
  models$std_accuracy <- std_accuracy_all
  models$min_accuracy <- min_accuracy_all
  models$max_accuracy <- max_accuracy_all
  models$all_accuracies <- all_accuracies
  models$v_matrices <- v_matrices 
  models$reconstruction_error <- reconstruction_error_list
  models$all_plots <- all_plots 
  return(models)
}






mean_without_zeros <- function(x) {
  values <- unlist(x)
  values <- values[values != 0 & !is.na(values)]
  if (length(values) > 0) {
    return(mean(values))
  } else {
    return(NA)
  }
}


generate_plot <- function(data, label_column, number) {
  
  # Definicja palety kolorów
  my_palette <- c("blue", "green", "red", "purple", "grey")
  
  # Ustalenie kolorów i kształtów w zależności od wartości label_column
  if (label_column == 'euthymia') {
    data$color <- ifelse(data[[label_column]] == '1', 'disease', 'euthymia')
    palette <- c("disease" = "blue", "euthymia" = "green")
    shapes <- c("disease" = 16, "euthymia" = 17, "lack" = 3)
  } else {
    present_colors <- unique(as.character(data[[label_column]]))
    all_colors <- c("1" = "blue", "2" = "green", "3" = "brown", "4" = "purple", "5" = "grey")
    palette <- all_colors[present_colors]
    shapes <- setNames(c(16, 17, 18, 19, 20, 3), c("1", "2", "3", "4", "5", "lack"))
    shapes <- shapes[unique(data[[label_column]])]
    data[[label_column]][is.na(data[[label_column]])] <- 'lack'
    data[[label_column]] <- as.factor(data[[label_column]])
  }
  
  # Tworzenie mapowania kolorów
  color_mapping <- scale_color_manual(values = palette, name = "Class Y")
  shape_mapping <- scale_shape_manual(values = shapes, guide = "none")
  
  # Generowanie wykresu
  ggplot(data, aes(x = x1, y = x2, color = factor(!!sym(label_column)), shape = factor(!!sym(label_column)))) +
    geom_point(size = 3) +
    labs(title = paste("Chunk", number)) +
    color_mapping +
    shape_mapping +
    xlim(c(-1, 4)) +  
    ylim(c(-0.5, 4.5)) +
    theme_minimal() +
    theme(
      legend.text = element_text(size = 22),        
      legend.title = element_text(size = 24),       
      plot.title = element_text(size = 26, face = "bold"),  
      axis.text.x = element_text(size = 20),
      axis.text.y = element_text(size = 20),
      axis.title.x = element_text(size = 22),       
      axis.title.y = element_text(size = 22)           
    ) +
    guides(color = guide_legend(override.aes = list(size = 5)))
}


