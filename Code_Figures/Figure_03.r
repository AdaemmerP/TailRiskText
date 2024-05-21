library(tidyverse)
library(conflicted) 
conflict_prefer("filter", "dplyr")

# Set path to folder 'Data/SentTopics/K_80'
# setwd("")

# Load topics
topics <- read_delim("ctm_80_10_topicsonly_monthly.csv",
  delim = ";",
  escape_double = FALSE,
  locale = locale(decimal_mark = ","),
  trim_ws = TRUE
) |>
  mutate(year_month = as.Date(year_month, "%d.%m.%Y"))

# Load SentTopics
sent_topics <- read_delim("ctm_80_10_senttopics_bignomics_monthly.csv",
  delim = ";",
  escape_double = FALSE,
  locale = locale(decimal_mark = ","),
  trim_ws = TRUE
) |>
  mutate(year_month = as.Date(year_month, "%d.%m.%Y"))

topics_df <- topics |>
  left_join(sent_topics, join_by(year_month)) |>
  pivot_longer(!year_month, names_to = "topic", values_to = "value") |>
  filter(str_detect(topic, "V9|V28|V46")) |>
  mutate(topic = case_when(
    topic == "V9" ~ "Topic 9: Debt Crisis",
    topic == "V28" ~ "Topic 28: Oil and War",
    topic == "V46" ~ "Topic 46: Housing",
    topic == "sent_V9" ~ "SentTopic 9: Debt Crisis",
    topic == "sent_V28" ~ "SentTopic 28: Oil and War",
    topic == "sent_V46" ~ "SentTopic 46: Housing",
    TRUE ~ topic
  )) |>
  mutate(topic = fct_relevel(topic, c(
    "Topic 9: Debt Crisis",
    "Topic 28: Oil and War",
    "Topic 46: Housing",
    "SentTopic 9: Debt Crisis",
    "SentTopic 28: Oil and War",
    "SentTopic 46: Housing"
  )))

topics_df |>
  ggplot() +
  geom_line(aes(x = year_month, y = value)) +
  geom_hline(
    data = filter(topics_df, str_detect(topic, "SentTopic")), aes(yintercept = 0), linetype = "dashed"
  ) +
  facet_wrap(~topic, scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = alpha("gray", 0.2)),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    text = element_text(size = 18),
    strip.text.x = element_text(size = 16),
    strip.text.y = element_text(size = 16)
  ) +
  labs(x = "Year", y = "(Sent)Topic Proportion")
