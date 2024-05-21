library(tidyverse)
library(conflicted) 
conflict_prefer("filter", "dplyr")

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

# Make tibble for ggplot
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
  group_by(Model, q, response) |>
  arrange(desc(value), .by_group = T) |>
  ungroup() |>
  mutate(Model = fct_relevel(Model, "Horseshoe", "Lasso", "Ridge")) |>
  mutate(var = str_replace(var, "-", "")) |>
  mutate(var = if_else(str_detect(var, "ylag"), str_c("Lag", str_extract(var, "\\d+")), var))

# Set order for factors
plot_b_df <- df_betas |>
  filter(
    q == 90,
    h == 1
  ) |>
  group_by(Model, response) |>
  drop_na() |>
  mutate(value = value / max(value, na.rm = T)) |>
  slice_max(n = 5, order_by = value, with_ties = T, na_rm = T) |>
  mutate(value = 1) |>
  mutate(
    fred_news = case_when(
      str_detect(var, "Lag") ~ "Lags of y",
      str_detect(var, "^V\\d+$") ~ "Topic",
      str_detect(var, "sent_") ~ "SentTopic",
      !str_detect(var, "Lag|Topic") ~ "FRED"
    ),
    fred_news = as.factor(fred_news)
  ) |>
  mutate(fred_news = fct_relevel(as_factor(fred_news), "FRED", "Topic", "SentTopic")) |>
  mutate(
    var = str_replace(var, "V", "Topic "),
    var = str_replace(var, "sent_Topic ", "SentTopic ")
  ) |>
  mutate(var = as_factor(var)) |>
  ungroup()

# relevel such that fred first, then topics and then senttopics
all_levels <- plot_b_df$var |>
  levels() |>
  sort()
level_no_topics <- all_levels[!str_starts(all_levels, "Topic|SentTopic")] |> str_sort(numeric = T)
level_topics <- all_levels[str_starts(all_levels, "Topic")] |> str_sort(numeric = T)
level_sent_topics <- all_levels[str_starts(all_levels, "SentTopic")] |> str_sort(numeric = T)
plot_b_df$var <- fct_relevel(plot_b_df$var, c(level_no_topics, level_topics, level_sent_topics))

# Make ggplot
plot_b_df |>
  ggplot() +
  theme_bw() +
  geom_col(aes(value, var, fill = fred_news), alpha = 1, width = 0.5) +
  facet_grid(Model ~ response, scales = "free_y") +
  theme(
    strip.background = element_rect(fill = alpha("gray", 0.1)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.spacing.x = unit(0.5, "lines"),
    panel.spacing.y = unit(1.5, "lines"),
    text = element_text(size = 12),
    axis.text = element_text(size = 8),
    strip.text.x = element_text(size = 12),
    strip.text.y = element_text(size = 12),
    legend.position = "bottom",
    legend.title = element_blank(),
    legend.text = element_text(size = 12),
    legend.margin = margin(t = -10),
    legend.key.size = unit(2, "lines")
  ) +
  labs(
    x = "",
    y = ""
  ) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = rev) +
  scale_fill_manual(values = c(r_colors[2], r_colors[5], r_colors[3], alpha("gray", 1)))

