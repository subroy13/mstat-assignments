setwd("D:/Assignments/Quantitative_Finance_DM")

library(tidyverse)

# SESAME SEEDS
dat <- read_csv('./FAOSTAT_data_10-15-2020.csv')

dat %>% 
    ggplot(aes(x = Year, y = Value, fill = Area)) +
    geom_col(stat='identity', position='dodge') +
    theme_bw() + 
    ylab("Total Production of Sesame Seeds") +
    ggsave("sesame-total.pdf")


# SUGAR CANE
dat <- read_csv('./FAOSTAT_data_sugar.csv')
dat %>% 
    dplyr::filter(Year > 2000, Item == "Sugar cane") %>%
    ggplot(aes(x = Year, y = Value, fill = Area)) +
    geom_col(stat='identity', position='dodge') +
    theme_bw() + 
    ylab("Total Production of Sugar Crops") +
    ggsave("sugar-total.pdf")

    
















