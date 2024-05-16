library(tidyverse)

# Set path to folder 'Data/SentTopics'
# setwd("")

# Load data (in folder 'Data/SentTopics')
sent_daily <- readRDS("sentiment_daily.RData")

# NBER recessions
usrec <- tribble(
  ~peak, ~trough,
  "1980-01-01", "1980-07-01",
  "1981-07-01", "1982-11-01",
  "1990-07-01", "1991-03-01",
  "2001-04-01", "2001-11-01",
  "2007-12-12", "2009-06-01",
  "2020-02-01", "2020-04-01"
) |>
  mutate(across(everything(), as.Date))

# Make plot
sent_daily |>
  group_by(year_month) |>
  summarise(sent_bignomics = mean(sent_bignomics, na.rm = T)) |>
  ggplot() +
  geom_line(aes(x = year_month, y = sent_bignomics)) +
  theme_bw() +
  geom_hline(aes(yintercept = 0), linetype = "dashed") +
  theme_classic() +
  geom_rect(data = usrec, aes(xmin = peak, xmax = trough, ymin = -Inf, ymax = Inf), fill = "grey", alpha = 0.5) +
  labs(
    x = "Year",
    y = "Sentiment value"
  )

