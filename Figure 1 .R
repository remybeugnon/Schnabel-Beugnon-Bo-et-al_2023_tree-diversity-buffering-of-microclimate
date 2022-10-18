#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# R version 4.1.3 (2022-03-10)
# Platform: aarch64-apple-darwin20 (64-bit)
# Running under: macOS Monterey 12.6
#
# Packages: 
# - dplyr_1.0.9 
# - stringr_1.4.0 
# - nlme_3.1-157
# - emmeans_1.7.3 
# - mgcv_1.8-39 
# - car_3.0-12 
# - ggplot2_3.3.5 
# - ggeffects_1.1.1 
# - ggpubr_0.4.0
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~# 

# Cleaning workspace
rm(list = ls())
set.seed(1603)

#### > 1. Packages ####
libs <- c(
  'dplyr', 'stringr',
  'nlme','emmeans',"mgcv",'car',
  'ggplot2', 'ggpubr', 'ggeffects'
)
invisible(lapply(libs, library, character.only = T))

#### > 2. Data ####
load("workspace_data.RData")

df.design = df.clim.month.sum %>%
  filter(TreeDiv > 0) %>% 
  mutate(year.month = year + (month -1 )/12) %>%
  mutate(month.f = factor(month)) %>%
  select(site, plot, 
         year, month, month.f, year.month, 
         TreeDiv)

d.1 = df.clim %>% 
  filter(!is.na(T) & 
           TreeDiv > 0 & 
           year == 2017) %>% # focusing on 2017 which represent an average year 
  mutate(date = paste0(year, month, day)) %>%
  mutate(hour.f = factor(hour))

d.2 = df.clim %>% 
  filter(TreeDiv > 0) %>%
  group_by(site, plot, year, month) %>%
  summarise(T.max = median(T[T>quantile(T,probs=c(.95),na.rm = T)], 
                           na.rm = T),
            T.med = median(T, na.rm = T),
            T.min = median(T[T<quantile(T,probs=c(.05),na.rm = T)], 
                           na.rm = T)) %>%
  left_join(.,
            df.design,
            by = c('site','plot','year','month'))

#### > 3. Statistical analyses ####
#### >> 3.1 Daily buffering ####
mod.1 = lme(T ~ log(TreeDiv) * hour.f, 
            random = ~ 1|site/plot/date,
            data = d.1,
            correlation=corCAR1())
summary(mod.1)

mod.1.res = lme(T ~ 1,
                random = ~ 1|site/plot/date,
                data = d.1,
                correlation=corCAR1())

d.1$res = residuals(mod.1.res) + mod.1.res$coefficients$fixed
d.1$group = d.1$hour.f

pred.1 = ggpredict(model = mod.2,
                   terms = c("TreeDiv", 'hour.f'))

t.lab = c('0-1','1-2','2-3','3-4','4-5','5-6',
          '6-7','7-8','8-9','9-10','10-11','11-12',
          '12-13','13-14','14-15','15-16','16-17',
          '17-18','18-19','19-20','20-21','21-22',
          '22-23','23-0')
names(t.lab) = 0:23

#### >> 3.2 Monthly extrema ####
#### >> 3.2.1 temperature maximum ####
mod.2.tmax = lme(T.max ~ log(TreeDiv, base = 2) * month.f,
                 random = ~ 1|site/plot/year,
                 data = d.2,
                 correlation=corCAR1())
summary(mod.2.tmax)

pred.2.max   = ggpredict(model = mod.2.tmax,   
                         terms = c("TreeDiv", 'month.f'))

mod.2.tmax.res = lme(T.max ~ 1,
                     random = ~ 1|site/plot/year,
                     data = d.2)

summary(mod.2.tmax)
Anova(mod.2.tmax, type = "II")

d.2$res.max = residuals(mod.2.tmax.res) + 
  mod.2.tmax.res$coefficients$fixed

df.out.max <- data.frame(matrix(nrow = 1, ncol = 6))
colnames(df.out.max) <- c("Value", "Std.Error", "DF", 
                          "t-value", "p-value", "month")

for (m in 1:12) {
  mod.month <- lme(T.max ~ log(TreeDiv, base = 2),
                   random = ~ 1|site/plot/year,
                   data = d.2 %>% filter(month == m),
                   correlation=corCAR1())
  
  df <- summary(mod.month)$tTable[2,] %>% 
    matrix(nrow=1) %>% 
    as.data.frame()
  
  colnames(df) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  df$month <- m
  df.out.max <- rbind(df.out.max,df)
}
df.out.max <- df.out.max[-1,] # remove intercept

#### >> 3.2.2 temperature median ####
mod.2.tmed = lme(T.med ~ log(TreeDiv, base = 2) * month.f,
                 random = ~ 1|site/plot/year,
                 data = d.2,
                 correlation=corCAR1())
summary(mod.2.tmed)
Anova(mod.2.tmed, type = "II")

pred.2.med   = ggpredict(model = mod.2.tmed,   
                         terms = c("TreeDiv", 'month.f'))

mod.2.tmed.res = lme(T.med ~ 1,
                     random = ~ 1|site/plot/year,
                     data = d.2)

d.2$res.med = residuals(mod.2.tmed.res) + 
  mod.2.tmed.res$coefficients$fixed

df.out.med <- data.frame(matrix(nrow = 1, ncol = 6)) 
colnames(df.out.med) <- c("Value", "Std.Error", "DF", 
                          "t-value", "p-value", "month")

for (m in 1:12) {
  mod.month <- lme(T.med ~ log(TreeDiv, base = 2),
                   random = ~ 1|site/plot/year,
                   data = d.2 %>% filter(month == m),
                   correlation=corCAR1())
  
  df <- summary(mod.month)$tTable[2,] %>% 
    matrix(nrow=1) %>% 
    as.data.frame()
  
  colnames(df) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  df$month <- m
  df.out.med <- rbind(df.out.med,df)
}
df.out.med <- df.out.med[-1,] # remove intercept

#### >> 3.2.3 temperature minimum ####
mod.2.tmin = lme(T.min ~ log(TreeDiv, base = 2) * month.f,
                 random = ~ 1|site/plot/year,
                 data = d.2,
                 correlation=corCAR1())

summary(mod.2.tmin)
Anova(mod.2.tmin, type = "II")

pred.2.min = ggpredict(model = mod.2.tmin,
                       terms = c("TreeDiv", 'month.f'))

mod.2.tmin.res = lme(T.min ~ 1,
                     random = ~ 1|site/plot/year,
                     data = d.2)

d.2$res.min = residuals(mod.2.tmin.res) +
  mod.2.tmin.res$coefficients$fixed

df.out.min <- data.frame(matrix(nrow = 1, ncol = 6)) 
colnames(df.out.min) <- c("Value", "Std.Error", "DF", 
                          "t-value", "p-value", "month")

for (m in 1:12) {
  mod.month <- lme(T.min ~ log(TreeDiv, base = 2),
                   random = ~ 1|site/plot/year,
                   data = d.2 %>% filter(month == m),
                   correlation=corCAR1())
  
  df <- summary(mod.month)$tTable[2,] %>% 
    matrix(nrow=1) %>% 
    as.data.frame()
  
  colnames(df) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  df$month <- m
  df.out.min <- rbind(df.out.min,df)
}
df.out.min <- df.out.min[-1,] # remove intercept

#### > 4. Figure 1 ####
#### >> 4.1 Figure 1.A. ####
p.day =
  ggplot(data = pred.1,
         aes(x, predicted)) +
  geom_line(data = pred.1,
            aes(x, predicted),
            col = 'black') +
  geom_ribbon(data = pred.1, 
              aes(ymin = conf.low, ymax = conf.high), 
              alpha = .1) +
  scale_x_continuous(trans = 'log2', breaks = c(1,4,24)) +
  labs(title = "Daily temperature modulation", 
       subtitle = "(Sp. Rich.: p = 0.14, 
       Hour: p < 0.001, 
       Sp. Rich. x Hour: p < 0.001)",
       x = 'Tree species richness x hour', 
       y = expression(paste('Hourly temperature [',
                            ~degree,'C]',sep=''))) +
  facet_grid(cols = vars(group),
             labeller = labeller(group = t.lab)) +
  lims(y = c(min(pred.2$conf.low, pred.2$conf.high), 
             max(pred.2$conf.low, pred.2$conf.high))) + 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        panel.spacing.x = unit(0.6, "lines"), 
        strip.text.x = element_text(size = 8.5)) 
p.day

#### >> 4.2 Figure 1.B. ####
d.2$group = d.2$month.f

text_max <- 
  data.frame(
    lab = c(paste0("p ", 
                   ifelse(round(df.out.max$`p-value` <.001), 
                          "< 0.001", 
                          paste0("= ", round(df.out.max$`p-value`, 3))))
    ),
    group = factor(1:12, levels=c(1:12))
  )

text_med <- 
  data.frame(
    lab = c(paste0("p ", 
                   ifelse(round(df.out.med$`p-value` <.001), 
                          "< 0.001", 
                          paste0("= ", round(df.out.med$`p-value`, 3))))
    ),
    group = factor(1:12, levels=c(1:12))
  )

text_min <- 
  data.frame(
    lab = c(paste0("p ", 
                   ifelse(round(df.out.min$`p-value` <.001), 
                          "< 0.001", 
                          paste0("= ", round(df.out.min$`p-value`, 3))))
    ),
    group = factor(1:12, levels=c(1:12))
  )

p.month =
  ggplot(data = pred.2.max, 
         aes(x, predicted)) + 
  geom_jitter(data = d.2, 
              aes(x = TreeDiv, y = res.max), 
              alpha = .02, color = 'red',
              size = .5) + 
  geom_line(data = pred.2.max, 
            aes(x, predicted), 
            col = 'red') +
  geom_ribbon(data = pred.2.max, 
              aes(ymin = conf.low, ymax = conf.high), 
              fill = 'red', alpha = .2) +
  geom_jitter(data = d.2, 
              aes(x = TreeDiv, y = res.med), 
              alpha = .02, color = 'black',
              size = .5) + 
  geom_line(data = pred.2.med, 
            aes(x, predicted), 
            col = 'black') +
  geom_ribbon(data = pred.2.med, 
              aes(ymin = conf.low, ymax = conf.high), 
              fill = 'black', alpha = .2) +
  geom_jitter(data = d.2, 
              aes(x = TreeDiv, y = res.min), 
              alpha = .02, color = 'blue',
              size = .5) + 
  geom_line(data = pred.2.min, 
            aes(x, predicted), 
            col = 'blue') +
  geom_ribbon(data = pred.2.min, 
              aes(ymin = conf.low, ymax = conf.high), 
              fill = 'blue', alpha = .2) +
  scale_x_continuous(trans = 'log2', breaks = c(1,2,4,8,24)) +
  labs(x = 'Tree species richness x month', 
       y = expression(paste('Monthly temperature [',
                            ~degree,'C]',sep='')), 
       title = "Monthly temperature patterns") +
  facet_grid(cols = vars(group)) +
  theme_bw() + 
  theme(axis.ticks.y.right = element_line(color = 'blue'),
        axis.text.y.right  = element_text(color = 'blue'),
        axis.title.y.right =  element_text(color = 'blue'),
        panel.grid = element_blank()) + 
  # addition p-values
  geom_text(aes(x=4,y=55, label = lab),
            data = text_max, 
            size = 3, 
            color = "red") +
  geom_text(aes(x=4,y=50, label = lab),
            data = text_med, 
            size = 3, 
            color = "black") +
  geom_text(aes(x=4,y=45, label = lab),
            data = text_min, 
            size = 3, 
            color = "blue")
p.month

#### >> 4.3 Figure ####
p = ggarrange(p.day,
              p.month, 
              nrow = 2, 
              labels = paste0(LETTERS[1:2],'.')
              ) 
p

ggsave(filename = 'Figure1.png', 
       height = 17, width = 28,
       units = 'cm')