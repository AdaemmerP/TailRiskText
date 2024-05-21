library(tidyverse)
library(conflicted) 
conflict_prefer("filter", "dplyr")

# Set path to folder 'Results'
# setwd("")

# Color palette
cbbPalette <- c("#999999", "#3a5795", "#000000", "#56B4E9", "#009E73", "#D55E00", "#CC79A7")

sheets <- readxl::excel_sheets("Results_linear_all.xlsx")
df_all <- purrr::map_df(
  sheets,
  ~ mutate(readxl::read_excel("Results_linear_all.xlsx", sheet = .x),
           sheetname = .x
  )
)

# Create tibble
df_plots_var <- df_all |>
  select(Row, QLmean, sheetname, QL_DM) |>
  mutate(model = (str_extract(Row, "\\w+_|^\\w+\\s\\w+") |>
                    str_remove("fred") |>
                    str_remove("_"))) |>
  mutate(fred = parse_number(str_extract(Row, "fred=\\d"))) |>
  mutate(text = parse_number(str_extract(Row, "text=\\d+"))) |>
  mutate(
    Predictors = str_c(fred, text),
    Predictors = case_when(
      Predictors == "16" ~ "K = 80",
      Predictors == "17" ~ "K = 68",
      Predictors == "18" ~ "K = 100",
      Predictors == "19" ~ "K = 80, full sample",
      Predictors == "111" ~ "K = 80, last week",
      Predictors == "112" ~ "K = 80, 15k words",
      Predictors == "113" ~ "K = 80, STM model",
    )
  ) |>
  mutate(Predictors = fct_relevel(Predictors, c("K = 68", "K = 80", "K = 100", "K = 80, 15k words", "K = 80, last week"))) |>
  mutate(
    q = parse_number(str_extract(sheetname, "q=\\d+")),
    q = fct_relevel(as.factor(q), c("5", "10"))
  ) |>
  mutate(h = parse_number(str_extract(sheetname, "h=\\d+"))) |>
  mutate(ylag = parse_number(str_extract(Row, "ylag=\\d+"))) |>
  mutate(
    response = str_extract(sheetname, "CPI|INDPRO|UMCSENT|CE16"),
    response = case_when(
      response == "CE16" ~ "Employment",
      response == "CPI" ~ "Inflation",
      response == "INDPRO" ~ "Production",
      response == "UMCSENT" ~ "Sentiment"
    )
  ) |>
  mutate(model = str_replace(model, "Random Forest", "QR Forest")) |>
  drop_na()

# -----------------------------------------------------------------------------#
# ---                               Upper Panel                        --------#
# -----------------------------------------------------------------------------#
df_plots_var |>
  filter(h == 1 & ylag == 12) |>
  filter(model %in% c("Horseshoe", "Lasso", "Ridge")) |>
  mutate(model = fct_relevel(model, "Horseshoe", "Lasso", "Ridge")) |>
  mutate(shape_point = if_else(QL_DM < 0.1, "p < 0.1", "p >= 0.1")) |>
  # make ggplot
  ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = q, y = QLmean, color = Predictors, group = Predictors)) +
  geom_point(aes(x = q, y = QLmean, color = Predictors, shape = shape_point)) +
  facet_grid(model ~ response, scales = "free_y") +
  labs(
    y = "Relative Quantile Loss",
    x = "Quantile"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -10),
    strip.background = element_rect(fill = alpha("grey", 0.1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 11),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  ) +
  scale_color_manual(values = cbbPalette) +
  scale_shape_manual(values = c(19, 1)) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

# -----------------------------------------------------------------------------#
# ---                               Lower Panel                        --------#
# -----------------------------------------------------------------------------#
df_plots_var |>
  filter(h == 2 & ylag == 12) |>
  filter(model %in% c("Horseshoe", "Lasso", "Ridge")) |>
  mutate(model = fct_relevel(model, "Horseshoe", "Lasso", "Ridge")) |>
  mutate(shape_point = if_else(QL_DM < 0.1, "p < 0.1", "p >= 0.1")) |>
  # make ggplot
  ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = q, y = QLmean, color = Predictors, group = Predictors)) +
  geom_point(aes(x = q, y = QLmean, color = Predictors, shape = shape_point)) +
  facet_grid(model ~ response, scales = "free_y") +
  labs(
    y = "Relative Quantile Loss",
    x = "Quantile"
  ) +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.margin = margin(t = -10),
    strip.background = element_rect(fill = alpha("grey", 0.1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 12),
    axis.title = element_text(size = 9),
    axis.text = element_text(size = 11),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12)
  ) +
  scale_color_manual(values = cbbPalette) +
  scale_shape_manual(values = c(19, 1)) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

