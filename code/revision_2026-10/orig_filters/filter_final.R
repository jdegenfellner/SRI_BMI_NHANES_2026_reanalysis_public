# Final archived filter function of the earlier study (R_Codes/2_Read_Clean_Prepare_for_Analysis.R), the filter described in its Methods; with the final pairing rule it gives the rebuilt series of this re-analysis.
# Copied verbatim for the reproduction of the published values (repro_published_sri.R); do not edit.
filter_accelo_data_NHANES_11_to_14 <- function(df, n_sub) {
    
    # Filter the data based on n_sub
    if (n_sub == "all") {
      df_filtered <- df
    } else {
      df_filtered <- df %>%
        dplyr::filter(SEQN %in% unique(df$SEQN)[1:n_sub])
    }
    
    # Initialize lists to store the number of participants and excluded SEQNs after each step
    participants_tracker <- list()
    excluded_seqn_tracker <- list()
    participants_tracker$initial <- length(unique(df$SEQN))
    
    # Step 0: Define binary_sleep variable
    df_filtered <- df_filtered %>%
      dplyr::mutate(binary_sleep = case_when(
        PAXPREDM == "sleep" ~ 1,
        PAXPREDM == "wake" ~ 0,
        TRUE ~ NA_real_
      )) 
    
    # Step 1: Filter for only 60s observations for the current minute
    initial_seqn <- unique(df_filtered$SEQN)
    df_filtered <- df_filtered %>%
      dplyr::filter(PAXTSM == 60)
    participants_tracker$after_60s_filter <- length(unique(df_filtered$SEQN))
    excluded_seqn_tracker$after_60s_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
    
    # Step 2: Filter for data quality flag count
    initial_seqn <- unique(df_filtered$SEQN)
    df_filtered <- df_filtered %>%
      dplyr::filter(PAXQFM == 0)
    participants_tracker$after_quality_flag_filter <- length(unique(df_filtered$SEQN))
    excluded_seqn_tracker$after_quality_flag_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
    
    # Step 3: Filter for at least 22 hours of wearing (non-wear filter) (see Validation paper)
    # filter out participants with less than 3 days of <120min non-wear time
    initial_seqn <- unique(df_filtered$SEQN)
    valid_non_wear_time_days <- df_filtered %>% # defined as: days with less than 2 hours of non-wear time
      dplyr::group_by(SEQN, PAXDAYM) %>%  # unique(df_PAX1112$PAXDAYM) : "1" "2" "3" "4" "5" "6" "7" "8" "9"
      dplyr::summarise(non_wear_time = sum(PAXPREDM == "non-wear"), .groups = 'drop') %>% 
      dplyr::filter(non_wear_time < 2*60) %>%  # Keeping days with less than 2 hours (120 minutes) of non-wear time
      dplyr::ungroup() %>%
      dplyr::group_by(SEQN) %>% 
      dplyr::summarise(valid_non_wear_days = n_distinct(PAXDAYM), .groups = 'drop') %>% 
      dplyr::filter(valid_non_wear_days >= 3) %>% # valid non wear day is a days with less than 2 hours of non-wear time
      dplyr::ungroup()
    df_filtered <- df_filtered %>%
      dplyr::filter(SEQN %in% valid_non_wear_time_days$SEQN)
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
    
    participants_tracker$after_time_steps_filter <- length(unique(df_filtered$SEQN))
    excluded_seqn_tracker$after_time_steps_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
    
    # Step 5: Filter for at least X consecutive days including a Saturday or Sunday
    initial_seqn <- unique(df_filtered$SEQN)
    df_filtered <- df_filtered %>%
      dplyr::group_by(SEQN) %>%
      dplyr::filter({
        day_seq <- sort(as.numeric(unique(PAXDAYWM))) # Convert day of week PAXDAYWM to numeric
        if (sum(!is.na(day_seq)) < 3) { # At least 3 days of data per participant
          consecutive <- FALSE # formally not entirely correct: there are just not 3 days, but the 2 can be consecutive nevertheless, but should work.
        } else {
          #consecutive <- max(rle(diff(day_seq))$lengths, na.rm = TRUE) >= 4 # At least 5 consecutive days
          consecutive <- max(rle(diff(day_seq))$lengths, na.rm = TRUE) >= 2 # At least 3 consecutive days
        }
        weekend_day <- any(day_seq %in% c(1, 7), na.rm = TRUE) # Includes Saturday or Sunday; see for example: https://wwwn.cdc.gov/Nchs/Nhanes/2013-2014/PAXDAY_H.htm#PAXDAYWD
        consecutive && weekend_day
      }) %>%
      dplyr::ungroup()
    participants_tracker$after_consecutive_days_filter <- length(unique(df_filtered$SEQN))
    excluded_seqn_tracker$after_consecutive_days_filter <- setdiff(initial_seqn, unique(df_filtered$SEQN))
    
    # Step 6: Calculate missingness for each SEQN and exclude participants
    # with missingness rate > 30%
    initial_seqn <- unique(df_filtered$SEQN)
    
    missingness <- df_filtered %>%
      dplyr::group_by(SEQN) %>%
      dplyr::summarize(
        total_observations = n(),
        missing_count = sum(is.na(binary_sleep)),
        missingness_rate = missing_count / total_observations
      ) %>%
      dplyr::ungroup() %>%
      as.data.table()
    
    df_filtered <- df_filtered %>% 
      dplyr::filter(SEQN %in% missingness$SEQN[missingness$missingness_rate < 0.3]) # missingness per SEQN < 30%
    
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
