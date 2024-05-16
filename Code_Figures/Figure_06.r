library(tidyverse)

# Set path to folder 'Results'
# setwd("")

# Colors
r_colors <- c(
  "#000000", "#E69F00", "#3a5795", "#56B4E9", "#009E73", "#0072B2",
  "#D55E00", "#CC79A7"
)

# Data in 'Results'
sheets <- readxl::excel_sheets("beta_results_80_10.xlsx")
df_betas <- purrr::map_df(sheets, ~ dplyr::mutate(
  readxl::read_excel("beta_results_80_10.xlsx",
    sheet = .x
  ),
  sheetname = .x
)) |>
  relocate(Row, sheetname, CPIAUCSL, INDPRO, UMCSENTx, CE16OV)

# Create beta tibble
df_betas <- df_betas |>
  mutate(
    q = parse_number(str_extract(sheetname, "q=\\d+")),
    h = parse_number(str_extract(sheetname, "hfore\\d+")),
    response = str_extract(sheetname, "INDPRO|CPI|UMCSENT|CE16"),
    response = case_when(
      response == "CE16" ~ "Employment",
      response == "CPI" ~ "Inflation",
      response == "INDPRO" ~ "Production",
      response == "UMCSENT" ~ "Sentiment"
    )
  ) |>
  select(-sheetname) |>
  mutate(
    Model = str_extract(Row, "Model=\\w+_"),
    Model = str_extract(Model, "=\\w+"),
    Model = str_remove_all(Model, "=|_|Text|Fred|Both"),
    Model = case_when(
      Model == "HS" ~ "Horseshoe",
      Model == "LASSO" ~ "Lasso",
      Model == "RIDGE" ~ "Ridge",
      Model == "BayesianKernel" ~ "Gaussian Processes",
      Model == "Tree" ~ "QR Forest"
    )
  ) |>
  relocate(Model) |>
  select(-Row) |>
  pivot_longer(!c(Model, q, response, h),
    names_to  = "var",
    values_to = "value"
  ) |>
  mutate(value = abs(value)) |>
  group_by(Model, q, response) |> # pred_set,
  arrange(desc(value), .by_group = T) |>
  ungroup() |>
  mutate(Model = fct_relevel(Model, "Horseshoe", "Lasso", "Ridge")) |>
  mutate(var = str_replace(var, "-", "")) |>
  mutate(var = if_else(str_detect(var, "ylag"), str_c("Lag", str_extract(var, "\\d+")), var))

# Make ggplot
df_betas |>
  mutate(Model == str_replace(Model, "Random Forest", "QR Forest")) |>
  filter(q %in% 10) |>
  mutate(var = str_squish(var)) |>
  mutate(var_type = case_when(
    str_detect(var, "^V\\d+") ~ "Topic",
    str_detect(var, "sent_") ~ "SentTopic",
    TRUE ~ "FRED"
  )) |>
  mutate(var_type = fct_relevel(var_type, c("FRED", "Topic", "SentTopic"))) |>
  mutate(value = if_else(value != 0, 1, 0)) |>
  filter(value == 1 & !str_detect(var, "Lag")) |>
  group_by(Model, response, var_type, q, h) |>
  count() |>
  ungroup() |>
  mutate(q = as.factor(q)) |>
  mutate(h = if_else(h == 1, "h=0", "h=1")) |>
  ggplot() +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = alpha("gray", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1, "lines"),
    text = element_text(size = 18),
    axis.text = element_text(size = 16),
    strip.text.x = element_text(size = 18),
    strip.text.y = element_text(size = 18),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    legend.margin = margin(t = -10),
    legend.key.size = unit(2, "lines")
  ) +
  geom_col(aes(h, n, fill = var_type), alpha = 1, width = 0.5) +
  ylim(c(0, 250)) +
  facet_grid(Model ~ response, scales = "free") +
  labs(
    y = "Number of nonzero coefficients",
    x = ""
  ) +
  scale_fill_manual(values = c(r_colors[2], r_colors[5], r_colors[3]))

