library(tidyverse)

# Set path to folder 'Results'
# setwd("")

# Colors
r_colors <- c(
  "#000000", "#E69F00", "#3a5795", "#56B4E9", "#009E73", "#0072B2",
  "#D55E00", "#CC79A7"
)

# -----------------------------------------------------------------------------#
# --------- Figure 5, upper panel (data in folder 'Results') ------------------#
# -----------------------------------------------------------------------------#
sheets <- readxl::excel_sheets("Results_nonlinear_all.xlsx")
df_all <- purrr::map_df(
  sheets,
  ~ mutate(readxl::read_excel("Results_nonlinear_all.xlsx", sheet = .x),
           sheetname = .x
  )
)

# Modify Tibble to extract responses
df_plots_var <- df_all |>
  select(Row, QLmean, sheetname, QL_DM) |>
  mutate(model = (str_extract(Row, "\\w+_|^\\w+\\s\\w+") |>
                    str_remove("fred") |>
                    str_remove("_"))) |>
  mutate(fred = parse_number(str_extract(Row, "fred=\\d"))) |>
  mutate(text = parse_number(str_extract(Row, "text=0y|text=1y|text=2y|text=6y"))) |>
  mutate(
    Predictors = str_c(fred, text),
    Predictors = case_when(
      Predictors == "10" ~ "FRED only",
      Predictors == "11" ~ "FRED & Topics",
      Predictors == "16" ~ "FRED & Topics & SentTopics"
    )
  ) |>
  mutate(Predictors = fct_relevel(Predictors, c("FRED only", "FRED & Topics"))) |> # sort baseline
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

# Make plot
df_plots_var |>
  filter(h == 7 & ylag == 12) |>
  filter(model %in% c("Gaussian Processes", "QR Forest")) |>
  mutate(shape_point = if_else(QL_DM < 0.1, T, F)) |>
  mutate(shape_point = if_else(QL_DM < 0.1, "p < 0.1", "p >= 0.1")) |>
  # make plot
  ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = q, y = QLmean, color = Predictors, group = Predictors)) +
  geom_point(aes(x = q, y = QLmean, color = Predictors, shape = shape_point)) +
  facet_grid(model ~ response, scales = "free_y") +
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
  scale_color_manual(values = c(r_colors[2], r_colors[5], r_colors[3])) +
  scale_shape_manual(values = c(19, 1)) +
  labs(
    y = "Relative Quantile Score",
    x = "Quantile"
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

# -----------------------------------------------------------------------------#
# ---------           Figure 5, lower panel                  ------------------#
# -----------------------------------------------------------------------------#

# Make plot
df_plots_var |>
  filter(h == 13 & ylag == 12) |>
  filter(model %in% c("Gaussian Processes", "QR Forest")) |>
  mutate(shape_point = if_else(QL_DM < 0.1, T, F)) |>
  mutate(shape_point = if_else(QL_DM < 0.1, "p < 0.1", "p >= 0.1")) |>
  # make plot
  ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = q, y = QLmean, color = Predictors, group = Predictors)) +
  geom_point(aes(x = q, y = QLmean, color = Predictors, shape = shape_point)) +
  facet_grid(model ~ response, scales = "free_y") +
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
  scale_color_manual(values = c(r_colors[2], r_colors[5], r_colors[3])) +
  scale_shape_manual(values = c(19, 1)) +
  labs(
    y = "Relative Quantile Score",
    x = "Quantile"
  ) +
  guides(
    colour = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )
