# Filter function in the state that produced the published 2011-2012 SRI values (file z_stuff/2.1_Calculate_SRI_manually.R of the earlier study's project; the state of the filter when the 2011-2012 values were computed and saved).
# Copied verbatim for the reproduction of the published values (repro_published_sri.R); do not edit.
filter_accelo_data_NHANES_11_to_14 <- function(df, n_sub) {
  
  # Filter the data based on n_sub
  if (n_sub == "all") {
    df_filtered <- df
  } else {
    df_filtered <- df %>%
      filter(SEQN %in% unique(df$SEQN)[1:n_sub])
  }
  
  # Initialize lists to store the number of participants and excluded SEQNs after each step
  participants_tracker <- list()
  excluded_seqn_tracker <- list()
  participants_tracker$initial <- length(unique(df$SEQN))
  
  # Step 1: Filter for only 60s observations for the current minute (necessary if one wants to use a time stamp)
  initial_seqn <- unique(df_filtered$SEQN)
  df_filtered <- df_filtered %>%
    filter(PAXTSM == 60)
  participants_tracker$after_60s_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_60s_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  # Step 3: Filter for data quality flag count
  initial_seqn <- unique(df_filtered$SEQN)
  df_filtered <- df_filtered %>%
    filter(PAXQFM == 0)
  participants_tracker$after_quality_flag_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_quality_flag_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  # Define binary_sleep variable
  df_filtered <- df_filtered %>%
    mutate(binary_sleep = case_when(
      PAXPREDM == "sleep" ~ 1,
      PAXPREDM == "wake" ~ 0,
      TRUE ~ NA_real_
    )) 
  
  # Step 2: Filter for at least 22 hours of wearing (non-wear filter) (see Validation paper)
  # filter out participants with less than 5 days of <120min non-wear time
  initial_seqn <- unique(df_filtered$SEQN)
  valid_non_wear_time_days <- df_filtered %>%
    group_by(SEQN, PAXDAYM) %>% 
    summarise(non_wear_time = sum(PAXPREDM == "non-wear")) %>% 
    ungroup() %>%
    filter(non_wear_time < 2*60) %>%  # Keeping days with less than 2 hours (120 minutes) of non-wear time
    group_by(SEQN) %>% 
    summarise(valid_non_wear_days = n_distinct(PAXDAYM)) %>% 
    ungroup() %>%
    filter(valid_non_wear_days >= 5)
  df_filtered <- df_filtered %>%
    filter(SEQN %in% valid_non_wear_time_days$SEQN)
  participants_tracker$after_non_wear_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_non_wear_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  # Step 4: filter for 1440 time steps per day
  # filter out participants with less than 5 days of 1440 time steps days.
  initial_seqn <- unique(df_filtered$SEQN)
  
  df_filtered <- df_filtered %>%
    dplyr::group_by(SEQN, PAXDAYM) %>%
    dplyr::mutate(time_steps_for_PAXDAYM = n()) %>%
    dplyr::filter(time_steps_for_PAXDAYM == 1440) %>%
    dplyr::ungroup()
  
  # ... this leads to dimension problems later...
  #valid_days_count <- df_filtered %>%
  #     group_by(SEQN, PAXDAYM) %>% 
  #     summarise(time_steps_for_PAXDAYM = n()) %>% 
  #     ungroup() %>%
  #     filter(time_steps_for_PAXDAYM == 1440) %>% 
  #     group_by(SEQN) %>%
  #     summarise(valid_days = n_distinct(PAXDAYM)) %>% 
  #     ungroup() %>% 
  #     filter(valid_days >= 5)
  #df_filtered <- df_filtered %>% filter(SEQN %in% valid_days_count$SEQN)
  
  participants_tracker$after_time_steps_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_time_steps_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  # Step 5: Filter for at least 5 consecutive days including a Saturday or Sunday
  initial_seqn <- unique(df_filtered$SEQN)
  df_filtered <- df_filtered %>%
    group_by(SEQN) %>%
    filter({
      day_seq <- sort(as.numeric(unique(PAXDAYWM))) # Convert day of week PAXDAYWM to numeric
      if (sum(!is.na(day_seq)) < 5) {
        consecutive <- FALSE
      } else {
        #consecutive <- max(rle(diff(day_seq))$lengths, na.rm = TRUE) >= 4 # At least 5 consecutive days
        consecutive <- max(rle(diff(day_seq))$lengths, na.rm = TRUE) >= 2 # At least 3 consecutive days
      }
      weekend_day <- any(day_seq %in% c(1, 7), na.rm = TRUE) # Includes Saturday or Sunday
      consecutive && weekend_day
    }) %>%
    ungroup()
  participants_tracker$after_consecutive_days_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_consecutive_days_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  # Step 6: Calculate missingness for each SEQN
  initial_seqn <- unique(df_filtered$SEQN)
  
  missingness <- df_filtered %>%
    group_by(SEQN) %>%
    dplyr::summarize(
      total_observations = n(),
      missing_count = sum(is.na(binary_sleep)),
      missingness_rate = missing_count / total_observations
    ) %>%
    ungroup() %>%
    as.data.table()
  
  df_filtered <- df_filtered %>% 
    filter(SEQN %in% missingness$SEQN[missingness$missingness_rate < 0.3])
  
  participants_tracker$after_missingness_filter <- length(unique(df_filtered$SEQN))
  excluded_seqn_tracker$after_missingness_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
  
  
  # Step 7: Calculate the percentage of rows and participants lost to filtering
  if (n_sub == "all") {
    rows_lost_percentage <- 100 - (dim(df_filtered)[1] / dim(df)[1] * 100)
  } else {
    rows_lost_percentage <- 100 - (dim(df_filtered)[1] / dim(df %>%
                                                               filter(SEQN %in% unique(df$SEQN)[1:n_sub]))[1] * 100) 
  }
  percent_participants_lost <- 100 - length(unique(df_filtered$SEQN)) / participants_tracker$initial * 100
  
  # Step 8: Check missingness in binary_sleep
  missingness_binary_sleep <- sum(is.na(df_filtered$binary_sleep)) / length(df_filtered$binary_sleep)
  
  # Step 9: Return the filtered data, summary statistics, participants tracker, and excluded SEQNs
  return(list(
    filtered_data = df_filtered,
    rows_lost_percentage = rows_lost_percentage,
    missingness_binary_sleep = missingness_binary_sleep,
    participants_tracker = participants_tracker,
    percent_participants_lost = percent_participants_lost,
    excluded_seqn_tracker = excluded_seqn_tracker
  ))
}
