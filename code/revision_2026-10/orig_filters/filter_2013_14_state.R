# Filter function in the state that produced the published 2013-2014 SRI values (file z_stuff/old_filtering_func.R of the earlier study's project: no 1440-minute rule, partial days of placement and removal kept, participant-days with 120 min or more of non-wear dropped).
# Copied verbatim for the reproduction of the published values (repro_published_sri.R); do not edit.
filter_accelo_data_NHANES_11_to_14 <- function(df, n_sub) {
  
  # Step 1: Filter the data based on n_sub
  if (n_sub == "all") {
    df_filtered <- df
  } else {
    df_filtered <- df %>%
      filter(SEQN %in% unique(df$SEQN)[1:n_sub])
  }
  
  # Step 2: Apply the remaining filters
  df_filtered <- df_filtered %>%
    filter(PAXTSM == 60) %>% # only 60s observations for the current minute
    group_by(SEQN, PAXDAYM) %>%
    #_old_filter(sum(PAXPREDM == "non-wear") < 480) %>% # at least 16 hours of wearing; < 60*8=480 of non-wear time
    filter(sum(PAXPREDM == "non-wear") < 2*60) %>% # at least 22 hours of wearing; See Validation Paper, Lunsford-Avery
    ungroup() %>%
    group_by(SEQN) %>%
    #mutate(number_of_valid_days = n_distinct(PAXDAYM)) %>%
    #ungroup() %>%
    #filter(number_of_valid_days >= 5) %>% # See Validation Paper, Lunsford-Avery
    filter(PAXQFM == 0) %>% # data quality flag count
    mutate(binary_sleep = case_when(
      PAXPREDM == "sleep" ~ 1,
      PAXPREDM == "wake" ~ 0,
      TRUE ~ NA_real_
    ))
  
  # New Step: Filter for at least 5 consecutive days including a Saturday or Sunday
  #df_filtered <- df_filtered %>%
  #  group_by(SEQN) %>%
  #  filter({
  #    day_seq <- sort(as.numeric(unique(PAXDAYWM))) # Convert PAXDAYWM to numeric
  #    consecutive <- max(rle(diff(day_seq))$lengths) >= 4 # at least 4 consecutive differences of 1 (5 consecutive days)
  #    weekend_day <- any(day_seq %in% c(1, 7)) # check if there is a Saturday (7) or Sunday (1)
  #    consecutive && weekend_day
  #  }) %>%
  #  ungroup()
  
  df_filtered <- df_filtered %>%
    group_by(SEQN) %>%
    filter({
      day_seq <- sort(as.numeric(unique(PAXDAYWM))) # Convert PAXDAYWM to numeric
      if (length(day_seq) >= 5 && !all(is.na(day_seq))) {
        consecutive <- max(rle(diff(day_seq))$lengths) >= 4 # at least 4 consecutive differences of 1 (5 consecutive days)
        weekend_day <- any(day_seq %in% c(1, 7)) # check if there is a Saturday (7) or Sunday (1)
        consecutive && weekend_day
      } else {
        FALSE # Skip this group if day_seq is empty or invalid
      }
    }) %>%
    ungroup()
  
  # Step 2: Calculate missingness for each SEQN
  missingness <- df_filtered %>%
    group_by(SEQN) %>%
    dplyr::summarize(
      #SEQN = unique(SEQN),
      total_observations = n(),
      missing_count = sum(is.na(binary_sleep)),
      missingness_rate = missing_count / total_observations
    ) %>%
    ungroup() %>%
    as.data.table()
  
  # Step 3: Filter out SEQNs with missingness rate >= 40%
  df_filtered <- df_filtered %>%
    filter(SEQN %in% missingness[missingness_rate < 0.4,]$SEQN)
  
  # Step 4: Calculate the percentage of rows lost to filtering
  if (n_sub == "all") {
    rows_lost_percentage <- 100 - (dim(df_filtered)[1] / dim(df)[1] * 100)
  } else {
    rows_lost_percentage <- 100 - (dim(df_filtered)[1] / dim(df %>%
                                                               filter(SEQN %in% unique(df$SEQN)[1:n_sub]))[1] * 100) 
  }
  
  # Step 5: Check missingness in binary_sleep
  missingness_binary_sleep <- sum(is.na(df_filtered$binary_sleep)) / length(df_filtered$binary_sleep)
  
  # Step 6: Generate plots
  #barplot(table(df_filtered$binary_sleep), main = "Binary Sleep Distribution")
  #barplot(table(df_filtered$PAXPREDM), main = "PAXPREDM Distribution")
  
  # Step 7: Return the filtered data and summary statistics
  return(list(
    filtered_data = df_filtered,
    rows_lost_percentage = rows_lost_percentage,
    missingness_binary_sleep = missingness_binary_sleep
  ))
}
