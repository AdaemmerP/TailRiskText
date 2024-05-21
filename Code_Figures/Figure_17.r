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

# Recession tibble
recessions_df <- tribble(
  ~Peak, ~Trough,
  "2007-12-01", "2009-06-01",
  "2020-02-01", "2020-04-01"
) |> mutate(across(everything(), as.Date))


sheets <- readxl::excel_sheets("Results_linear_Overtime.xlsx")
df_all <- purrr::map_df(sheets, ~ dplyr::mutate(readxl::read_excel("Results_linear_Overtime.xlsx", sheet = .x),
                                                sheetname = .x
)) |> filter(str_detect(Row, "Ridge|Lasso|Horseshoe|Gaussian Processes|Random Forest"))

# Modify tibble
df_plots_time <- df_all |>
  pivot_longer(cols = starts_with("ql_time_")) |>
  rename(QLmean = value) |>
  select(-name) |>
  mutate(model = str_extract(Row, "\\w+_|^\\w+\\s\\w+") |>
           str_remove("fred") |>
           str_remove("_")) |>
  filter(!str_detect(model, "AR1")) |>
  mutate(fred = parse_number(str_extract(Row, "fred=\\d"))) |>
  mutate(text = parse_number(str_extract(Row, "text=\\d"))) |>
  mutate(Predictors = str_c(fred, text)) |>
  mutate(Predictors = case_when(
    Predictors == "10" ~ "FRED only",
    Predictors == "11" ~ "FRED & Topics",
    Predictors == str_c("16") ~ "FRED & Topics & SentTopics"
  )) |>
  mutate(Predictors = fct_relevel(Predictors, "FRED only", "FRED & Topics", "FRED & Topics & SentTopics")) |>
  mutate(q = parse_number(str_extract(sheetname, "q=\\d+"))) |>
  mutate(q = fct_relevel(as.factor(q), c("5", "10"))) |>
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
  mutate(
    model = case_when(model == "Random Forest" ~ "QR Forest",
                      .default = model
    ),
    model = fct_relevel(model, c("Horseshoe", "Lasso", "Ridge"))
  )

# Date sequence for forecasts
date_seq <- seq(as.Date("1999-09-01"), as.Date("2021-12-01"), by = "month")

# Show ggplot
df_plots_time |>
  filter(h == 2 & q == 10 & ylag == 12) |>
  filter(Predictors != "FRED only") |>
  select(-Row, -sheetname, -text, -fred) |>
  drop_na() |>
  group_by(model, Predictors, q, response, h) |>
  mutate(date = rep_len(date_seq, length.out = n())) |>
  filter(date >= "2005-01-01") |>
  ungroup() |>
  # Make plot
  ggplot() +
  theme_bw() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_line(aes(x = date, y = QLmean, color = Predictors)) +
  geom_rect(data = recessions_df, aes(xmin = Peak, xmax = Trough, ymin = -Inf, ymax = +Inf), fill = "gray", alpha = 0.5) +
  facet_grid(model ~ response, scales = "free_y") +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = alpha("gray", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  scale_color_manual(values = c(r_colors[5], r_colors[3])) + # r_colors[2],
  labs(
    x = "Year",
    y = "Relative Cumulative Quantile Score",
    color = NULL
  )
