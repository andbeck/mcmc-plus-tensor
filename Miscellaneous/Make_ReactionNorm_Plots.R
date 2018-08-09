library(tidyverse)

load('EcolLettSupp.Rdata')

# choose data
df <- data1[[1]]
glimpse(df)

# stack correctly
df_e1 <- df %>% select(sire, trait1, trait2, trait3, trait4, trait5) %>% 
  mutate(environment = rep("E1", dim(df)[1])) %>% 
  mutate(population = rep("P1", dim(df)[1]))

df_e2 <- df %>% select(sire, trait6, trait7, trait8, trait9, trait10) %>% 
  mutate(environment = rep("E2", dim(df)[1])) %>% 
  mutate(population = rep("P2", dim(df)[1])) %>% 
  rename(trait1 = trait6, trait2 = trait7, trait3 = trait8, trait4 = trait9, trait5 = trait10)

df_use <- bind_rows(df_e1, df_e2) %>% 
  group_by(sire, population) %>% 
  summarise(trait1_mean = mean(trait1),
            trait2_mean = mean(trait2),
            trait3_mean = mean(trait3))

# plot and annotate
ggplot(df_use, aes(x = population, y = trait1_mean, group = sire))+
  geom_line()+
  ylab("Trait Value")+
  annotate('text', x = 2.2, y = 0, label = "VARIANCE", angle = 90) +
  theme_bw(base_size = 15)
