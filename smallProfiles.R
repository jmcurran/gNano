gnano.df = readData()

library(tidyverse)
smallPeaks = gnano.df %>%
  arrange(prof,loc) %>%
  filter(obs < 500) %>%
  select(prof) %>%
  pull() %>%
  as.numeric() %>%
  unique() %>%
  sort()

smallProf = gnano.df %>%
  filter(prof %in% smallPeaks) %>%
  group_by(prof) %>%
  summarise(min = min(obs), q05 = quantile(obs, 0.05), med = median(obs), mean = mean(obs), q95 = quantile(obs, 0.95), max = max(obs))

smallProf

smallProf = gnano.df %>%
  filter(prof %in% smallPeaks) %>%
  group_by(prof) %>%
  filter(obs < 500)

smallProf
