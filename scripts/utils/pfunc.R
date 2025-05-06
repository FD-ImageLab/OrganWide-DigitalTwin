source("scripts/utils/Func.R")
library("mediation")
library(ggplot2)
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)
library(RColorBrewer)
library(data.table)
# Function to generate phenotype data based on x and y, including residuals and distances
gene_new_pheno <- function(x, y, x_name, y_name, col_seq = "_") {
    # Find common subjects between x and y
    enroll_id <- intersect(rownames(na.omit(x)), rownames(na.omit(y)))
    x <- x[enroll_id, ]
    y <- y[enroll_id, ]

    x <- data.frame(scale(x), check.names = FALSE)
    y <- data.frame(scale(y), check.names = FALSE)

    health_x <- x[intersect(enroll_id, health_id), ]
    health_y <- y[intersect(enroll_id, health_id), ]

    results <- PMA::CCA(health_x, health_y, typex = "standard", typez = "standard")
    u <- results$u
    v <- results$v
    x_prime <- as.matrix(x) %*% u
    y_prime <- as.matrix(y) %*% v

    x_prime <- scale(x_prime)
    y_prime <- scale(y_prime)

    # Print correlation and p-value
    print(paste0(
        x_name, " and ", y_name, " correlation is ", results$cors,
        ", p-value is ", cor.test(x_prime, y_prime)$p.value
    ))

    data <- data.frame(x_prime = x_prime, y_prime = y_prime)

    line <- lm(y_prime ~ x_prime, data = data)
    y_res <- line$residuals %>% as.data.frame()
    x_res <- -y_res / line$coefficients[2]
    dis <- abs(x_res * y_res / sqrt(x_res^2 + y_res^2))
    colnames(y_res) <- c(paste0("y_res_", x_name, col_seq, y_name))
    colnames(x_res) <- c(paste0("x_res_", x_name, col_seq, y_name))
    colnames(dis) <- c(paste0("dis_", x_name, col_seq, y_name))
    new_pheno <- cbind(y_res, x_res, dis)
    rownames(new_pheno) <- rownames(y_res)
    new_pheno <- as.data.frame(new_pheno)
    return(y_res)
}

# Function to compute CCA-transformed x and y variables (without further analysis)
gene_cca_xy_prime <- function(x, y, x_name, y_name) {
    enroll_id <- intersect(rownames(na.omit(x)), rownames(na.omit(y)))
    x <- x[enroll_id, ]
    y <- y[enroll_id, ]

    x <- data.frame(scale(x), check.names = FALSE)
    y <- data.frame(scale(y), check.names = FALSE)

    health_x <- x[intersect(enroll_id, health_id), ]
    health_y <- y[intersect(enroll_id, health_id), ]

    results <- PMA::CCA(health_x, health_y, typex = "standard", typez = "standard")
    u <- results$u
    v <- results$v

    x_prime <- as.matrix(x) %*% u
    y_prime <- as.matrix(y) %*% v

    x_prime <- scale(x_prime)
    y_prime <- scale(y_prime)
    return(list(x_prime = x_prime, y_prime = y_prime))
}


read_ukb <- function(raw_file_path, ref_file_path, select_column, select_name) {
    ref_df <- read.csv(ref_file_path)
    UDI <- ref_df$UDI[ref_df[[select_column]] == select_name]
    df <- fread(raw_file_path, select = UDI, sep = ",")
}











# Function to categorize phenotype data into quantiles
divide_pheno_by_quantile <- function(data = data,
                                     probs = c(0, 0.2, 0.8, 1),
                                     labels = c("lower", "middle", "upper")) {
    # 计算分位数
    quantiles <- quantile(data, probs = probs, na.rm = TRUE)

    # 使用cut函数分组
    groups <- cut(data, breaks = quantiles, labels = labels, include.lowest = TRUE)
    if (length(data) != length(groups)) {
        error("输入和输出的array长度不一致")
    }
    return(groups)
}




####################### Without rescale #######################
###### use the 3th col as the edge ##########
draw_organ_complexnet <- function(data, output_name = "plot/complexnet/ukb_organ_axis_cca_complexnet.png",
                                  limits = c(0, 5), node_size = c(0.2, 0.6),
                                  col = "#1f78b4", show_title = FALSE) {
    smallest_size <- 0.0000000001
    colnames(data)[3] <- "r"
    data$organ1 <- rename_based_on_df(data$organ1, color_map, "organ", "short_organ")
    data$organ2 <- rename_based_on_df(data$organ2, color_map, "organ", "short_organ")
    # build network
    exports <- data %>%
        distinct(organ1) %>%
        rename(label = organ1)
    imports <- data %>%
        distinct(organ2) %>%
        rename(label = organ2)
    nAttr <- nrow(full_join(exports, imports, by = "label"))

    nodes <- full_join(exports, imports, by = "label")
    nodes$label <- nodes$label[order(match(nodes$label, color_map$organ))]
    nodes <- nodes %>%
        mutate(id = 1:nAttr) %>%
        dplyr::select(id, everything())


    # change node color
    nodes$label <- factor(nodes$label, levels = nodes$label)
    node_actual_name <- as.character(nodes$label)
    node_color <- rename_based_on_df(node_actual_name, color_map,
        from = "short_organ", to = "color"
    )


    edges <- data %>%
        left_join(nodes, by = c("organ1" = "label")) %>%
        rename(from = id)

    edges <- edges %>%
        left_join(nodes, by = c("organ2" = "label")) %>%
        rename(to = id)

    edges <- dplyr::select(edges, from, to, r)

    network <- tbl_graph(
        nodes = nodes, edges = edges, directed = TRUE
    )

    x0 <- seq(0, nAttr - 1)
    coords1 <- data.frame(
        x = 2 * sin(x0 * 360 / nAttr * pi / 180),
        y = 2 * cos(x0 * 360 / nAttr * pi / 180)
    )

    r_alpha <- c(edges$r)
    point_size <- c()
    for (i in 1:nAttr) {
        point_size[i] <- (c(mean(edges$r[which(edges$from == i | edges$to == i)], na.rm = TRUE)) - node_size[1]) * node_size[2]
    }
    point_size[is.nan(point_size)] <- smallest_size
    print(point_size)
    g <- ggraph(network, layout = coords1) +
        geom_edge_fan(aes(edge_width = (r)),
            color = col, show.legend = F
        ) +
        geom_node_circle(aes(color = label, fill = label, r = point_size),
            show.legend = F
        ) +
        coord_fixed() +
        scale_edge_width(range = c(1, 8), limits = limits) +
        scale_size(range = c(0, 2), limits = limits) +

        scale_edge_alpha(range = c(0.1, 1)) +
        geom_node_text(aes(label = label), fontface = "bold", repel = F, size = 15, family = FONTFAMILY) +
        # scale_color_identity(str_replace_all(, "_", "\n"),aesthetics ="label")+
        scale_color_manual(values = node_color) +
        scale_fill_manual(values = node_color) +
        theme(plot.margin = margin(t = 0, r = 1, b = 0, l = 1, unit = "cm")) +
        theme_graph()
    if (show_title) {
        g <- g +
            ggtitle(basename(output_name))
    }


    print(g)
    check_path(paste0(output_name))
    ggsave(plot = g, filename = paste0(output_name), width = 7, height = 7, units = "in", dpi = 300)  # Save the plot
}



# 表型分组函数，输入为array
divide_pheno_by_quantile <- function(data = data,
                                     probs = c(0, 0.2, 0.8, 1),
                                     labels = c("lower", "middle", "upper")) {
    # 计算分位数
    quantiles <- quantile(data, probs = probs, na.rm = TRUE)

    # 使用cut函数分组
    groups <- cut(data, breaks = quantiles, labels = labels, include.lowest = TRUE)
    if (length(data) != length(groups)) {
        error("输入和输出的array长度不一致")
    }
    return(groups)
}



# Function for performing mediation analysis
# Input: df - dataframe containing three columns: X, Y, and M (no missing values)
# Output: Mediation results as a dataframe with XYM, mediation proportion, p-values, etc.
my_mediation <- function(df = df, boot = TRUE, sims = 500) {
    # mediator ~ x
    formula1 <- my_formula(
        y_variable = colnames(df)[3],
        x_variables = colnames(df)[1]
    )
    fit.mediator <- lm(formula1, data = df)
    fit.mediator$call$formula <- formula1
    a <- fit.mediator$coefficients[2]
    # y ~ x + mediator
    formula2 <- my_formula(
        y_variable = colnames(df)[2],
        x_variables = c(colnames(df)[1], colnames(df)[3])
    )
    fit.dv <- lm(formula2, data = df)
    fit.dv$call$formula <- formula2
    b <- fit.dv$coefficients[3]
    # mediate
    if (boot) {
        results <- mediation::mediate(fit.mediator, fit.dv,
            treat = colnames(df)[1],
            mediator = colnames(df)[3], boot = TRUE, sims = sims
        )
    } else {
        results <- mediation::mediate(fit.mediator, fit.dv,
            treat = colnames(df)[1],
            mediator = colnames(df)[3], boot = FALSE
        )
    }

    results <- summary(results)
    acme <- results$d0
    acme.ci.l <- results$d0.ci[1]
    acme.ci.u <- results$d0.ci[2]
    acme.p <- results$d0.p
    ade <- results$z0
    ade.ci.l <- results$z0.ci[1]
    ade.ci.u <- results$z0.ci[2]
    ade.p <- results$z0.p
    total.effect <- results$tau.coef
    total.effect.ci.l <- results$tau.ci[1]
    total.effect.ci.u <- results$tau.ci[2]
    total.effect.p <- results$tau.p
    prop.mediated <- results$n0
    prop.mediated.ci.l <- results$n0.ci[1]
    prop.mediated.ci.u <- results$n0.ci[2]
    prop.mediated.p <- results$n0.p
    temp_csv <- data.frame(
        x = colnames(df)[1], y = colnames(df)[2],
        mediator = colnames(df)[3],
        acme = acme, acme.p = acme.p, acme.ci.l = acme.ci.l, acme.ci.u = acme.ci.u,
        ade = ade, ade.p = ade.p, ade.ci.l = ade.ci.l, ade.ci.u = ade.ci.u,
        total.effect = total.effect, total.effect.p = total.effect.p,
        total.effect.ci.l = total.effect.ci.l, total.effect.ci.u = total.effect.ci.u,
        prop.mediated = prop.mediated, prop.mediated.p = prop.mediated.p,
        prop.mediated.ci.l = prop.mediated.ci.l, prop.mediated.ci.u = prop.mediated.ci.u,
        a = a, b = b, N = nrow(df), row.names = NULL, check.names = FALSE
    )
    return(temp_csv)
}

# Simplify joining dataframes using rownames as ID columns
# Input: df1, df2 - dataframes to join
# Output: Joined dataframe with rownames as ID columns
my_join <- function(df1, df2, .f = "full_join") {
    if (.f == "full_join") {
        df1 <- df1 %>% rownames_to_column("ID")
        df2 <- df2 %>% rownames_to_column("ID")
        df <- full_join(x = df1, y = df2, by = "ID") %>%
            column_to_rownames("ID")
        return(df)
    } else if (.f == "inner_join") {
        df1 <- df1 %>% rownames_to_column("ID")
        df2 <- df2 %>% rownames_to_column("ID")
        df <- inner_join(x = df1, y = df2, by = "ID") %>%
            column_to_rownames("ID")
        return(df)
    } else if (.f == "left_join") {
        df1 <- df1 %>% rownames_to_column("ID")
        df2 <- df2 %>% rownames_to_column("ID")
        df <- left_join(x = df1, y = df2, by = "ID") %>%
            column_to_rownames("ID")
        return(df)
    } else {
        stop("Input '.f' must be one of {inner_join/full_join/left_join}")  # Error for invalid join type
    }
}





generate_uniform_points_in_circle <- function(n, center = c(1, 1), radius = 1) {
    # Generate n uniformly distributed points within a circle.
    #
    # :param n: Number of points to generate
    # :param center: Numeric vector representing the center of the circle (x, y)
    # :param radius: Radius of the circle
    # :return: Matrix of points (x, y)

    angles <- seq(-1 / 2 * pi, 3 / 2 * pi, length.out = n + 1)[1:n]
    x <- center[1] + radius * cos(angles)
    y <- center[2] + radius * sin(angles)
    return(cbind(x, y))
}


# Function to map categorical labels to coordinates in a circle with varying radii
# Input: 
# - category: A vector of category labels (e.g., node labels).
# - r_n2: The radius for points within each category (controls spread of points inside each node group).
# - coef: A scaling factor for overall positioning.
# Output: A matrix of coordinates for each point, corresponding to their category positions in the circle.

get_position <- function(category = node$node_organ, r_n2 = 0.2, coef = 1) {
    # Order the categories to maintain consistent positioning
    order_category <- order(order(category))
    
    # Get the number of unique categories (n1_n)
    n1_n <- length(unique(category))
    
    # Generate initial positions for the unique categories (n1_n) in a circle
    n1_pos <- generate_uniform_points_in_circle(n = n1_n, center = c(coef, coef), radius = coef)
    
    # Count the number of occurrences of each category (n2_n)
    n2_n <- table(category)
    
    # Initialize an empty matrix for storing positions of points within each category
    n2_pos <- NULL
    
    # Loop through each unique category and generate positions for points inside each category
    for (i in c(1:n1_n)) {
        n1_ipos <- n1_pos[i, ]  # Get the center position for the current category
        n2_in <- n2_n[i]  # Get the number of points in the current category
        # Generate positions for the points within the current category (using smaller radius)
        n2_ipos <- generate_uniform_points_in_circle(n = n2_in, center = n1_ipos, radius = r_n2 * coef)
        n2_pos <- rbind(n2_pos, n2_ipos)  # Append the generated positions to the overall list
    }
    
    # Reorder the points to match the original category order
    n2_pos <- n2_pos[order_category, ]
    
    # Return the matrix of positions for all points
    return(n2_pos)
}




devtools::install_github("fbreitwieser/sankeyD3")
library(sankeyD3)
draw_sankey <- function(df = df, source = "source", target = "target", value = "value",
                        node_color_map = NULL, from = NULL, to = NULL,
                        output_path = "") {
    df <- as.data.frame(df)
    node <- c(unique(df[[source]]), unique(df[[target]]))
    node <- data.frame(node = node, ID = c(0:(length(node) - 1)))

    node_color <- NULL
    if (!is.null(node_color_map)) {
        node$color <- rename_based_on_df(node$node,
            nmapdf = node_color_map,
            from = from, to = to
        )
    }
    df[[source]] <- rename_based_on_df(df[[source]], node, from = "node", to = "ID")
    df[[target]] <- rename_based_on_df(df[[target]], node, from = "node", to = "ID")
    df[[source]] <- as.double(df[[source]])
    df[[target]] <- as.double(df[[target]])



    sankey <- sankeyD3::sankeyNetwork(df, node,
        Source = source, NodeColor = "color", showNodeValues = FALSE,
        Target = target, Value = value, NodeID = "node", NodeGroup = "node", orderByPath = TRUE,
        units = "TWh", nodeWidth = 80, fontSize = 50, align = "none"
    )

    check_path(output_path)
    saveWidget(sankey, output_path)
}



draw_pie <- function(df) {
    col <- c("#EE924F", "#BDBBD7", "#7DBFA6", "#F2B670", "#BDDD78", "#EB8677", "#BDB0A5")
    ggplot(df, aes(x = 2.2, y = value, fill = category)) + 
        geom_bar(stat = "identity", width = 1, color = "white") +
        coord_polar(theta = "y") +
        xlim(c(1, 3)) +
        theme_void() + 
        scale_fill_manual(values = col) +
        scale_color_manual(values = col) + 
        labs(fill = "Category")
    # geom_text(aes(label = label), position = position_stack(vjust = 0.5),
    #           size=15) # 添加标签
}






draw_sankey_3layer <- function(df = df, source = "source", target = "target", value = "value",
                               node_color_map = NULL, from = NULL, to = NULL,
                               output_path = "", sort_by_node_value = TRUE) {
    browser()
    df <- as.data.frame(df)

    x <- setdiff(df[[source]], df[[target]])
    y <- intersect(df[[source]], df[[target]])
    z <- setdiff(df[[target]], df[[source]])
    node <- c(x, y, z)
    node <- data.frame(node = node, ID = c(0:(length(node) - 1)))

    if (sort_by_node_value) {
        x_sort <- df[df[[source]] %in% x, ] %>%
            group_by(source) %>%
            dplyr::summarise(node_value = sum(mean))
        new_x <- x_sort$source[rev(order(x_sort$node_value))]
        y_sort <- df[df[[source]] %in% y, ] %>%
            group_by(source) %>%
            dplyr::summarise(node_value = sum(mean))
        new_y <- y_sort$source[rev(order(y_sort$node_value))]
        z_sort <- df[df[[target]] %in% z, ] %>%
            group_by(target) %>%
            dplyr::summarise(node_value = sum(mean))
        new_z <- z_sort[[target]][rev(order(z_sort$node_value))]
        node$node <- c(new_x, new_y, new_z)
        compareFunction <- function(a, b) {
            return(a.ID - b.ID)
        }
    }


    node_color <- NULL
    if (!is.null(node_color_map)) {
        node$color <- rename_based_on_df(node$node,
            nmapdf = node_color_map,
            from = from, to = to
        )
    }


    df[[source]] <- rename_based_on_df(df[[source]], node, from = "node", to = "ID")
    df[[target]] <- rename_based_on_df(df[[target]], node, from = "node", to = "ID")
    df[[source]] <- as.double(df[[source]])
    df[[target]] <- as.double(df[[target]])
    df <- df %>% arrange(source, target)


    sankey <- sankeyD3::sankeyNetwork(df, node,
        Source = source, NodeColor = "color", showNodeValues = FALSE,
        Target = target, Value = value, NodeID = "node", NodeGroup = "node", orderByPath = TRUE,
        units = "TWh", nodeWidth = 80, fontSize = 50, align = "none", dragY = TRUE,
    )

    check_path(output_path)
    saveWidget(sankey, output_path)
}
