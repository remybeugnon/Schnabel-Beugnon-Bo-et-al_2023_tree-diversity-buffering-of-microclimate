#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# R version 4.1.3 (2022-03-10)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.6
#
# Packages: 
# - tidyverse_1.3.1
# - nlme_3.1-157
# - emmeans_1.7.3 
# - mgcv_1.8-39
# - piecewiseSEM_2.1.2
# - ggplot2_3.3.5 
# - ggeffects_1.1.1 
# - ggpubr_0.4.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Cleaning workspace
rm(list = ls())
set.seed(1603)

#### > 1. Packages ####
libs <- c(
  'tidyverse',
  'nlme', "mgcv", 'piecewiseSEM',
  'ggplot2', 'ggpubr', 'ggeffects'
)

invisible(lapply(libs, library, character.only = T))

#### > 2. Data ####
load("workspace_data.RData")

df.phe = read.csv("leaf-phe.csv")

df.biomass <- 
  readxl::read_xlsx("Stand-volume_2015-2020.xlsx") %>% 
  select(site = SITE, plot = PLOT, 
         `2020` = VOL.a.2020, `2019` = VOL.a.2019, `2018` = VOL.a.2018,
         `2017` = VOL.a.2017, `2016` = VOL.a.2016, `2015` = VOL.a.2015) %>%
  pivot_longer(cols = 3:8, 
               names_to = 'year', 
               values_to = 'volume', 
               names_transform = as.numeric, 
               values_transform = as.numeric)

df.macro.T = read.csv("macro_1901_2021.csv") %>%
  filter(year %in% 2015:2020) %>% 
  select(year, month, 
         avg.t = tmp_celsius, 
         min.t = tmn_celsius, 
         max.t = tmx_celsius)

df.macro.spei = read.csv("spei_1901_2021.csv") %>%
  filter(year %in% 2015:2020) %>% 
  select(date, year,spei = SPEI12) %>%
  filter(str_detect(date,'-12-')) %>% 
  select(-date)

d.1 = df.clim %>% 
  filter(TreeDiv > 0 & year != 2014) %>%
  group_by(site, plot, year, month) %>%
  summarise(Buff = mean(T, na.rm = T)/ sd(T, na.rm = T) , 
            TreeDiv = mean(TreeDiv)) %>%
  left_join(., 
            df.phe, 
            by = c('site', 'plot', 'month' = 'time')) %>%
  left_join(., 
            df.biomass |> 
              select(site, plot, year, volume), 
            by = c('site', 'plot', 'year')) %>%
  left_join(., 
                df.design, 
                by = c('site', 'plot')) %>%
  left_join(df.macro.T, 
            by = c('year', 'month')) %>% 
  filter(!is.na(phe.asyn)) %>%
  mutate(year.month = year + (month-1)/12)

d.1$phe.asyn = scale(d.1$phe.asyn)
d.1$Buff = scale(d.1$Buff)
d.1$phe.cwm = scale(d.1$phe.cwm)
d.1$phe.asyn.y = scale(d.1$phe.asyn.y)
d.1$volume = scale(as.numeric(d.1$volume))
d.1$avg.t = scale(d.1$avg.t)
d.1$log.TreeDiv = log(d.1$TreeDiv, base = 2)

#### > 3. Statistical analyses ####
#### >> 3.1. Monthly buffering drivers ####
month.mod = 
  psem(
    lme(phe.asyn ~ log.TreeDiv + avg.t,
        data = d.1, 
        random = ~ 1|site/plot,
        correlation=corCAR1(form = ~year.month)),
    lme(phe.cwm ~ log.TreeDiv + avg.t,
        data = d.1, 
        random = ~ 1|site/plot,
        correlation=corCAR1(form = ~year.month)),
    lme(volume ~ log.TreeDiv + avg.t,
        data = d.1, 
        random = ~ 1|site/plot, 
        correlation=corCAR1(form = ~year.month)),
    lme(Buff ~ avg.t + (phe.cwm + phe.asyn + volume + log.TreeDiv),
        data = d.1, 
        random = ~ 1|site/plot,
        correlation=corCAR1(form = ~year.month)),
    phe.cwm %~~% phe.asyn,
    volume %~~% phe.asyn,
    data = d.1 %>% data.frame()
  )

summary(month.mod)
plot(month.mod)

#### >> 3.2. Monthly fits ####
mod.sem = function(m, d.1){
  dd = d.1 |> 
    filter(month == m)
  
  month.mod = psem(
    lme(phe.asyn ~ log.TreeDiv + avg.t,
        data = dd, 
        random = ~ 1|site/plot),
    lme(phe.cwm ~ log.TreeDiv + avg.t,
        random = ~ 1|site/plot,
        data = dd),
    lme(volume ~ log.TreeDiv + avg.t,
        data = dd, 
        random = ~ 1|site/plot),
    lme(Buff ~ avg.t + (phe.cwm + phe.asyn + volume + log.TreeDiv),
        data = dd, 
        random = ~ 1|site/plot),
    phe.cwm %~~% phe.asyn,
    volume %~~% phe.asyn,
    data = dd %>% data.frame())
  a = summary(month.mod)$coefficients |> 
    select(Response, Predictor, Std.Estimate, Std.Error, P.Value) |> 
    mutate(month = m) %>%
    mutate(R2 = summary(month.mod)$R2[4,5])
  rm(dd,month.mod)
  return(a)
}

model.SEM = 
  1:12 |> 
  map_df(~mod.sem(.x, d.1 = d.1))

model.SEM = model.SEM |> 
  mutate(signif = if_else(P.Value<0.1,1,0))

model.SEM$Predictor =
  model.SEM$Predictor %>%
  str_replace_all('log.TreeDiv', 'Tree species richness') %>%
  str_replace_all('avg.t', 'Monthly\nmacro-temperature') %>%
  str_replace_all('phe.cwm', 'Leaf phenological\npotential') %>%
  str_replace_all('phe.asyn', 'Leaf phenological\nasynchrony') %>%
  str_replace_all('volume', 'Wood volume') %>%
  factor(levels = c('Monthly\nmacro-temperature',
                    'Tree species richness',
                    'Leaf phenological\npotential',
                    'Leaf phenological\nasynchrony',
                    'Wood volume'))

#### > 4. Plot ####
p = ggplot(data = model.SEM %>% filter(Response == 'Buff'), 
       aes(x = month, y = Std.Estimate)) + 
  geom_point(data = model.SEM %>% filter(Response == 'Buff'), 
             aes(x = month, y = Std.Estimate)) + 
  geom_smooth(method = 'loess',
              data = model.SEM %>% filter(Response == 'Buff'), 
              aes(x = month, y = Std.Estimate), color = 'black', 
              se = T) +
  facet_grid(rows = vars(Predictor), scales = 'free') +
  scale_x_continuous(breaks = 1:12, minor_breaks = T) +
  scale_y_continuous( minor_breaks = T) +
  labs(x = 'Month', y = 'Estimate') + 
  theme_bw() + 
  theme()

p

# All plots were made manually using Inkscape