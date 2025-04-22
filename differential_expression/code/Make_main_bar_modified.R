Make_main_bar_modified <- function (Main_bar_data, Q, show_num, ratios, customQ, number_angles, 
          ebar, ylabel, ymax, scale_intersections, text_scale, attribute_plots) 
{
  library(ggbreak)
  bottom_margin <- (-1) * 0.65
  if (is.null(attribute_plots) == FALSE) {
    bottom_margin <- (-1) * 0.45
  }
  if (length(text_scale) > 1 && length(text_scale) <= 6) {
    y_axis_title_scale <- text_scale[1]
    y_axis_tick_label_scale <- text_scale[2]
    intersection_size_number_scale <- text_scale[6]
  }
  else {
    y_axis_title_scale <- text_scale
    y_axis_tick_label_scale <- text_scale
    intersection_size_number_scale <- text_scale
  }
  if (is.null(Q) == F) {
    inter_data <- Q
    if (nrow(inter_data) != 0) {
      inter_data <- inter_data[order(inter_data$x), ]
    }
    else {
      inter_data <- NULL
    }
  }
  else {
    inter_data <- NULL
  }
  if (is.null(ebar) == F) {
    elem_data <- ebar
    if (nrow(elem_data) != 0) {
      elem_data <- elem_data[order(elem_data$x), ]
    }
    else {
      elem_data <- NULL
    }
  }
  else {
    elem_data <- NULL
  }
  if (is.null(ymax) == T) {
    ten_perc <- ((max(Main_bar_data$freq)) * 0.1)
    ymax <- max(Main_bar_data$freq) + ten_perc
  }
  if (ylabel == "Intersection Size" && scale_intersections != 
      "identity") {
    ylabel <- paste("Intersection Size", paste0("( ", scale_intersections, 
                                                " )"))
  }
  if (scale_intersections == "log2") {
    Main_bar_data$freq <- round(log2(Main_bar_data$freq), 
                                2)
    ymax <- log2(ymax)
  }
  if (scale_intersections == "log10") {
    Main_bar_data$freq <- round(log10(Main_bar_data$freq), 
                                2)
    ymax <- log10(ymax)
  }
  break_1 <- 800
  break_2 <- 10000
  my_squish <- trans_new(
    name      = "my_squish_1k_10k",
    transform = function(x) {
      ifelse(
        x <= break_1, 
        x,
        ifelse(
          x <= break_2,
          break_1 + 0.01 * (x - break_1),
          break_1 + 0.01 * (break_2 - break_1) + (x - break_2)
        )
      )
    },
    inverse = function(y) {
      cut_point <- break_1 + 0.01 * (break_2 - break_1)
      ifelse(
        y <= break_1,
        y,
        ifelse(
          y <= cut_point,
          break_1 + (y - break_1) / 0.01,
          break_2 + (y - cut_point)
        )
      )
    },
    domain = c(0, Inf)
  )
  
  Main_bar_plot <- (ggplot(data = Main_bar_data, aes_string(x = "x", 
                                                            y = "freq")) +
                      scale_y_continuous(
                        trans      = my_squish,
                        limits     = c(0, 10200),
                        breaks     = c(0, 100, 200, 300, 400, 500, 600, 700, 800, 10000, 10100, 10200),                   # optional: nicer labels
                        minor_breaks = NULL                            # or seq(0, ymax, by = some_small)
                      ) 
                    
                    +
                      geom_bar(stat = "identity", width = 0.6,
                                               fill = Main_bar_data$color) + scale_x_continuous(limits = c(0, 
                                                                                                           (nrow(Main_bar_data) + 1)), expand = c(0, 0), breaks = NULL) + 
                      xlab(NULL) + ylab(ylabel) + labs(title = NULL) + theme(panel.background = element_rect(fill = "white"), 
                                                                             plot.margin = unit(c(0.5, 0.5, bottom_margin, 0.5), 
                                                                                                "lines"), panel.border = element_blank(), axis.title.y = element_text(vjust = -0.8, 
                                                                                                                                                                      size = 8.3 * y_axis_title_scale), axis.text.y = element_text(vjust = 0.3, 
                                                                                                                                                                                                                                   size = 7 * y_axis_tick_label_scale)))
  if ((show_num == "yes") || (show_num == "Yes")) {
    Main_bar_plot <- (Main_bar_plot + geom_text(aes_string(label = "freq"), 
                                                size = 2.2 * intersection_size_number_scale, vjust = -1, 
                                                angle = number_angles, colour = Main_bar_data$color))
  }
  bInterDat <- NULL
  pInterDat <- NULL
  bCustomDat <- NULL
  pCustomDat <- NULL
  bElemDat <- NULL
  pElemDat <- NULL
  if (is.null(elem_data) == F) {
    bElemDat <- elem_data[which(elem_data$act == T), ]
    bElemDat <- bElemDat[order(bElemDat$x), ]
    pElemDat <- elem_data[which(elem_data$act == F), ]
  }
  if (is.null(inter_data) == F) {
    bInterDat <- inter_data[which(inter_data$act == T), 
    ]
    bInterDat <- bInterDat[order(bInterDat$x), ]
    pInterDat <- inter_data[which(inter_data$act == F), 
    ]
  }
  if (length(customQ) != 0) {
    pCustomDat <- customQ[which(customQ$act == F), ]
    bCustomDat <- customQ[which(customQ$act == T), ]
    bCustomDat <- bCustomDat[order(bCustomDat$x), ]
  }
  if (length(bInterDat) != 0) {
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bInterDat, 
                                              aes_string(x = "x", y = "freq"), fill = bInterDat$color, 
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if (length(bElemDat) != 0) {
    Main_bar_plot <- Main_bar_plot + geom_bar(data = bElemDat, 
                                              aes_string(x = "x", y = "freq"), fill = bElemDat$color, 
                                              stat = "identity", position = "identity", width = 0.6)
  }
  if (length(bCustomDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_bar(data = bCustomDat, 
                                               aes_string(x = "x", y = "freq2"), fill = bCustomDat$color2, 
                                               stat = "identity", position = "identity", width = 0.6))
  }
  if (length(pCustomDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pCustomDat, 
                                                 aes_string(x = "x", y = "freq2"), colour = pCustomDat$color2, 
                                                 size = 2, shape = 17, position = position_jitter(width = 0.2, 
                                                                                                  height = 0.2)))
  }
  if (length(pInterDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pInterDat, 
                                                 aes_string(x = "x", y = "freq"), position = position_jitter(width = 0.2, 
                                                                                                             height = 0.2), colour = pInterDat$color, size = 2, 
                                                 shape = 17))
  }
  if (length(pElemDat) != 0) {
    Main_bar_plot <- (Main_bar_plot + geom_point(data = pElemDat, 
                                                 aes_string(x = "x", y = "freq"), position = position_jitter(width = 0.2, 
                                                                                                             height = 0.2), colour = pElemDat$color, size = 2, 
                                                 shape = 17))
  }
  Main_bar_plot <- (Main_bar_plot + geom_vline(xintercept = 0, 
                                               color = "gray0") + geom_hline(yintercept = 0, color = "gray0"))
  Main_bar_plot <- ggplotGrob(Main_bar_plot)
  return(Main_bar_plot)
}
