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
  'nlme','emmeans',"mgcv",
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
  filter(TreeDiv > 0) %>%
  group_by(site, plot, year, month) %>%
  summarise(T.buff = mean(T, na.rm = T) / sd(T, na.rm = T)) %>%
  left_join(.,
            df.design ,
            by = c('site','plot','year','month'))

d.2 = df.clim %>% 
  filter(TreeDiv > 0) %>%
  group_by(site, plot, year) %>%
  summarise(T.buff = mean(T, na.rm = T) / sd(T, na.rm = T) ,
            TreeDiv = mean(TreeDiv)) %>% 
  mutate(year = factor(year))

df.macro.T = 
  read.csv("macro_1901_2021.csv") %>%
  filter(year %in% 2015:2020) %>% 
  select(year, month, 
         avg.t = tmp_celsius, 
         min.t = tmn_celsius, 
         max.t = tmx_celsius)

df.macro.T.year = df.macro.T %>%
  group_by(year) %>%
  summarise(mean.avg.t = mean(avg.t),
            sd.avg.t = sd(avg.t), 
            min.t = mean(min.t),
            max.t = mean(max.t)
            )

df.macro.plot =
  df.macro.T %>%
  filter(year %in% 2015:2020) %>%
  group_by(month) %>%
  summarise(
    avg = mean(avg.t),
    int.pos = mean(avg.t) + 1.96 * sd(avg.t),
    int.neg = mean(avg.t) - 1.96 * sd(avg.t)
  ) %>%
  mutate(group = as.factor(month))

df.macro.spei = 
  read.csv("spei_1901_2021.csv") %>%
  filter(year %in% 2015:2020) %>% 
  select(date, year, spei = SPEI12) %>%
  filter(str_detect(date,'-12-')) %>% 
  select(-date)

#### > 3. Statistical analyses ####
#### >> 3.1 Monthly buffering ####
mod.1 = lme(T.buff ~ log(TreeDiv, base = 2) * month.f,
              random = ~ 1|site/plot/year,
              data = d.1,
              correlation=corCAR1(form = ~year.month))
summary(mod.1)
Anova(mod.1)

df.out <- data.frame(matrix(nrow = 1, ncol = 6)) 
colnames(df.out) <- c("Value", "Std.Error", "DF", 
                      "t-value", "p-value", "month")

for (m in 1:12) {
  mod.month <- lme(T.buff ~ log(TreeDiv, base = 2),
                   random = ~ 1|site/plot/year,
                   data = d.1 %>% filter(month == m),
                   correlation=corCAR1(form = ~year))
  
  df <- summary(mod.month)$tTable[2,] %>% 
    matrix(nrow=1) %>% 
    as.data.frame()
  
  colnames(df) <- c("Value", "Std.Error", "DF", "t-value", "p-value")
  df$month <- m
  df.out <- rbind(df.out,df)
}
df.out$month <- as.factor(df.out$month)
df.out <- df.out[-1,]

mod.1.2.res = lme(T.buff ~ 1,
                  random = ~ 1|site/plot/year,
                  data = d.1,
                  correlation=corCAR1(form = ~year.month))

d.1$res = residuals(mod.1.res) + mod.1.res$coefficients$fixed
d.1$group = d.1$month.f
pred.1 = ggpredict(model = mod.1,
                   terms = c("TreeDiv", 'month.f'))

ann_text <- data.frame(
  lab = c(paste0("p ", 
                 ifelse(round(df.out$`p-value` <.001), 
                        "< 0.001", paste0("= ", 
                                          round(df.out$`p-value`, 3))
                 ))
  ),
  group = factor(1:12, levels=c(1:12)))

#### >> 3.2 Yearly buffering ####
mod.2 = lme(T.buff ~ log(TreeDiv, base = 2) * year, 
              random = ~ 1|site/plot,
              data = d.2,
              correlation=corCAR1(form = ~year))
summary(mod.2)
Anova(mod.2)

mod.2.2 = lme(T.buff ~ log(TreeDiv) , 
                random = ~ 1|site/plot/year,
                data = d.2,
                correlation=corCAR1())

summary(mod.2.2)

mod.2.2.res = lme(T.buff ~ 1, 
                    random = ~ 1|site/plot,
                    data = d.2,
                    correlation=corCAR1())
mod.2.2.res
d.2$res = residuals(mod.2.2.res) + mod.2.2.res$coefficients$fixed 
d.2$group = d.2$year

pred.2   = ggpredict(model = mod.2,   terms = c("TreeDiv", 'year'))
pred.2.2 = ggpredict(model = mod.2.2, terms = c("TreeDiv", 'year'))

#### > 4. Figure 2 ####
#### >> 4.1 Figure 2.A. ####
# Second y-axis coefficients
coef = (9/30)
int = 0

p.month =
  ggplot(data = NULL) + 
  geom_point(data = df.macro.plot, 
             aes(x = 4, y = (avg-int) * coef),
             fill = 'blue', alpha = 1) + 
  geom_linerange(data = df.macro.plot, 
                 aes(x = 4, ymin = (int.neg-int)*coef, 
                     ymax = (int.pos-int)*coef),
                 color = 'blue', alpha = 1) + 
  geom_jitter(data = d.1.2, 
              aes(x = TreeDiv, y = res), 
              alpha = .05) + 
  geom_line(data = pred.1, 
            aes(x, predicted),
            col = 'black') +
  geom_ribbon(data = pred.1, 
              aes(x = x,
                  ymin = conf.low, 
                  ymax = conf.high), 
              alpha = .2) +
  scale_x_continuous(trans = 'log2', 
                     breaks = c(1,2,4,8,24)) + 
  scale_y_continuous(
    sec.axis = sec_axis(~(./ coef) + int ,
                        name = expression(
                          paste('Mean monthly macroclimatic temperature [',
                                ~degree,'C]',
                                sep='')))) +
  labs(x = 'Tree species richness x month', 
       y = "Monthly temperature\nstability (1/CV)", 
       title = "Monthly temperature buffering", 
       subtitle = "n = 4'728") +
  facet_grid(cols = vars(group)) +
  theme_bw() + 
  theme(axis.ticks.y.right = element_line(color = 'blue'),
        axis.text.y.right  = element_text(color = 'blue'),
        axis.title.y.right =  element_text(color = 'blue'),
        panel.grid = element_blank()) + 
  geom_text(data = ann_text,
            aes(x = 1.5, y = 11.2, 
                label = lab, angle = 90),
            size = 3)

p.month

#### >> 4.2 Figure 2.B. ####
df.macro.T.year$group = df.macro.T.year$year %>% factor
df.macro.spei$group = df.macro.spei$year %>% factor
pred.2.p = left_join(pred.2  %>% data.frame(), 
                     df.macro.spei %>% select(group, spei),
                     by = c('group'))
df.macro.spei$year = df.macro.spei$year %>% factor
d.2 = left_join(d.2  %>% data.frame(), 
                       df.macro.spei %>% select(year, spei),
                       by = c('year'))
coef.1 = 6
int.1 = 16.5

p.year =
  ggplot(data = pred.2.2, 
         aes(x, predicted)) +
  geom_jitter(data = d.2, aes(x = TreeDiv, y = res, 
                                color = spei), 
              alpha = .2, width = 0.05, size = .5) +
  geom_line(data = pred.2.2,
            aes(x, predicted, group = group),
            color = 'black') +
  geom_line(data = pred.2.p,
            aes(x, predicted, group = group, color = spei),
            lty = 2) +
  geom_label(data = data.frame(
                      x = 27.5,
                      y = c(2.22, 2.03, 2.095, 1.97, 2.155, 2.28),
                      label = seq(2015,2020)
                      ),
             aes(x = x, y = y, label = label),
             size = 3,
             color = 'gray50') +
  annotate(geom = 'text', x = 20 , y = 1.75, label = 'p < 0.001') + 
  scale_x_continuous(trans = 'log2', breaks = c(1,2,4,8,16, 24)) + 
  scale_color_gradient(low = 'red', high = 'blue') + 
  labs(x = 'Tree species richness', y = "Annual temperature\nstability (1/CV)", 
       title = "Yearly temperature buffering", 
       subtitle = "n = 375",
       color = "SPEI 12") +
  lims(y = c(1.7,2.3)) + 
  theme_bw() + 
  theme(axis.ticks.y.right = element_line(color = 'blue'),
        axis.text.y.right  = element_text(color = 'blue'),
        axis.title.y.right =  element_text(color = 'blue'),
        panel.grid = element_blank())
p.year

#### >> 4.3 Figure ####
p = ggarrange(p.month, p.year, 
              nrow = 2, 
              heights = c(.6,.4),
              labels = paste0(LETTERS[1:2],'.')) 
p

ggsave(filename = 'Figure2.png', height = 17, width = 20, units = 'cm')

