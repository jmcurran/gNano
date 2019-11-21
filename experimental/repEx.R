library(dplyr)
library(magrittr)

example.df = data.frame(
  person = rep(1:10, rep(8, 10)),
  locus = rep(rep(1:4, rep(2, 4)), 10),
  dye = rep(rep(1:2, 4), 10),
  obs = round(runif(80, 100, 500))
)

example.df = example.df %>% 
  group_by(person, locus) %>% 
  summarise(d = diff(obs))
            
            