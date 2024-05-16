# Load packages
library(tidyverse)
library(quanteda)
library(stm)
library(zoo)

#-------------------------------------------------------------------------------#
#-----------------------       Training Data      ------------------------------#
#-------------------------------------------------------------------------------#

# Load prepared dfm
dfm_all <- readRDS("")

# Start and end dates (docvars of dfm contain a column with dates for each document)
start_train <- 1
end_train   <- tail(which(docvars(dfm_all)$date == "1999-08-31"), 1)

# Prepare training dfm
dfm_train <- dfm_all[start_train:end_train, ]

# Get feature names to keep via tf-idf
feat_names <- featnames(dfm_train)[order(colSums(dfm_tfidf(dfm_train)), decreasing = T)]
feat_keep <- feat_names[1:10000] # The main results are based on the top 10k words

# Keep only top features
dfm_train <- dfm_train |>
  dfm_select(feat_keep,
    selection = "keep",
    verbose = T
  )

# Remove empty documents
dfm_train <- dfm_train[rowSums(dfm_train) > 0, ]

# Convert to stm object
dfm_train_ctm <- dfm_train |>
  quanteda::convert(
    "stm",
    docvars = docvars(dfm_train)
  )

# Estimate ctm with training data
ctm_train <- stm(
  documents = dfm_train_ctm$documents,
  vocab = dfm_train_ctm$vocab,
   #data      = dfm_train_ctm$meta,  # For STM
   #prevalence = ~type,              # For STM
  init.type = "Spectral",
  K = 80,
  verbose = T
)

#-------------------------------------------------------------------------------#
#-----------------------       Testing Data       ------------------------------#
#-------------------------------------------------------------------------------#

# Test documents
dfm_test <- dfm_all[docvars(dfm_all)$date >= "1999-09-01", ]
dfm_test <- dfm_test |>
  dfm_select(
    ctm_train$vocab,
    selection = "keep",
    valuetype = "fixed",
    verbose   = T
  )

# Remove empty documents 
dfm_test <- dfm_test[rowSums(dfm_test) > 0, ]

# Convert and align testing data
dfm_test_ctm <- convert(
  dfm_test,
  to = "stm",
  docvars = docvars(dfm_test)
)
align_corpus <- alignCorpus(new = dfm_test_ctm, old.vocab = ctm_train$vocab)

# Fit new documents
fit_test <- fitNewDocuments(
  model = ctm_train,
  documents = align_corpus$documents,
   #newData          = dfm_test_ctm$meta, # For STM
   #prevalence       =~type,              # For STM
  prevalencePrior = "None"
)

#-------------------------------------------------------------------------------#
#-----------------------   Compute topics                -----------------------#
#-------------------------------------------------------------------------------#

# Merge meta data
meta_data <- dfm_train_ctm$meta |>
  as_tibble() |>
  bind_rows(as_tibble(dfm_test_ctm$meta)) |>
  select(id, date, year_month)

# Get training and test topic proportions and join with metadata
topics_all <- ctm_train$theta |>
  as_tibble() |>
  bind_rows(as_tibble(fit_test$theta)) |>
  bind_cols(meta_data) |>
  relocate(id, date, year_month) |>
  mutate(year_month = as.Date(as.yearmon(year_month)))

# Compute monthly topic proportions
topics_monthly <- topics_all |>
  select(-id, -date) |>
  group_by(year_month) |>
  summarize_all(mean)

#-------------------------------------------------------------------------------#
#-----------------------      Compute SentTopics         -----------------------#
#-------------------------------------------------------------------------------#
# Load daily sentiment tibble
sent_daily <- readRDS("sentiment_daily.RData")

# Multiply each topic with sentiment value
topics_sent_monthly <- topics_all |>
  left_join(sent_daily, join_by(id, date, year_month)) |>
  mutate(across(starts_with("V"), ~ . * sent_bignomics)) |>
  rename_with(~ str_c("sent_", .), starts_with("V")) |>
  group_by(year_month) |>
  summarise_at(vars(starts_with("sent_")), mean, na.rm = T)


