library(dplyr)
library(tidyr)
df1 <- data.frame(date = seq(as.Date("2021-01-01"), as.Date("2023-01-25"), by = "1 month"),
                  station = sample(100:200, 25),
                  value = c(value = c(sample(1:10, 10), rep(Inf, 5), sample(1:10, 10))),
                  valueTwo = c(value = c(sample(50:75, 5), rep(Inf, 5), sample(50:75, 15)))) %>% 
  mutate(month = format(date, format = "%b"))

# Want max per month, replace Inf with that
dfMax <- df1 %>% 
  group_by(month) %>% 
  summarise(across(c(value, valueTwo), ~max(.x[which(is.finite(.x))], na.rm = T)),
            .groups = "drop") %>%
  # summarise(across(contains("value"), ~max(.x[which(is.finite(.x))], na.rm = T)), 
  #           .groups = "drop") %>% 
  # rename_with(.cols = -month, ~paste0(.x, "_replace"))
  pivot_longer(-month, names_to = "variable", values_to = "replaceValue")

# Now, join and replace
what <- df1 %>% 
  pivot_longer(contains("value"), names_to = "variable", values_to = "values") %>% 
  left_join(dfMax, 
            by = c("month", "variable")) %>% 
  mutate(valuesTest = ifelse(is.infinite(values), replaceValue, values)) %>% 
  select(-c(replaceValue, values)) %>%
  pivot_wider(names_from = variable, values_from = valuesTest) 
  
  
tryThis <- c("value", "valueTwo")
# or
# tryThis <- names(df1)[which(grepl("value", names(df1)))]

what2 <- lapply(tryThis, function(x) {
  x <- sym(x)
  
  df <- df1 %>% 
    select(date, station, month, !!x)
  
  dfMax <- df %>% group_by(month) %>% 
    summarise(replace = max({{x}}[which(is.finite({{x}}))], na.rm = T))

  dfFin <- df %>% 
    left_join(dfMax, by = "month") %>% 
    mutate(!!x := ifelse(is.infinite(!!x), replace, !!x)) %>% 
    select(-replace)
}) %>% 
  purrr::reduce(full_join, by = c("date", "station", "month"))

all.equal(data.frame(what), what2)


