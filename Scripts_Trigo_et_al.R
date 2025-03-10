# Scripts for data analysis and plots of the manuscript by Trigo et al.

# Date: August, 2024

# Journal: Applied Vegetation Science

# Libraries ----

library(tidyverse)
library(readxl)
library(patchwork)
library(ggrepel)
library(iNEXT)
library(patchwork)
library(lme4)
library(nlme)

# Data loading and manipulation ----

data <- read_excel("data_Trigo_et_al.xlsx", sheet = "renovales15")

# Convert columns to factors

data <- data %>%
  mutate(across(c(Bloque, Clausura, Parcela, Especie), as.factor))

# Filter saplings and remove shrub species (Molle and Abreboca)
# which do not have a tree form in the study area.

saplings <- data %>% 
  select(-c("Orden", "Fecha")) %>% 
  filter(Renoval == "Si",
         Especie != "Abreboca",
         Especie != "Molle") %>% 
  select(-"Renoval")

# Sapling rinchness and diversity ----

## Fig. 2, Panel a) Rank-abundance curves ----

RA <- read_excel("data_Trigo_et_al.xlsx", sheet = "rank_ab")

RA <- RA %>% 
  mutate(
    sp = case_when(
      Especies == "Quebracho_colorado" ~ "Sch.lor",
      Especies == "Quebracho_blanco" ~ "Asp.que",
      Especies == "Mistol" ~ "Sar.mis",
      Especies == "Sombra_de_toro" ~ "Jod.rho",
      Especies == "Pata" ~ "Xim.ame",
      Especies == "Guayacan" ~ "Lib.par",
      TRUE ~ Especies
    )
  )

(a <- ggplot(RA, aes(x = Orden, y = pi)) +
    geom_point(aes(color = Clausura, shape = Clausura), size = 3) +
    geom_line(aes(color = Clausura), linewidth = 0.7) +
    geom_text_repel(
      aes(label = sp),
      nudge_x = 0.7, # Horizontal adjustment to move the text to the right
      size = 3,
      box.padding = 0.15,
      point.padding = 0.3,
      segment.color = "grey50",
      fontface = "italic"
    ) +
    labs(
      y = "Proportional abundance",
      x = "Species rank"
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.75, 0.9),
      legend.direction = "horizontal", # coordinates between 0 and 1
      legend.title = element_blank(), # remove legend title
      axis.text.x = element_blank(), # remove x-axis numbers
      legend.background = element_rect(fill = "transparent")) # transparent legend background
)

## Fig. 2, Panel b) Rarefaction curves ----

raref <- read_excel("data_Trigo_et_al.xlsx", sheet = "rarefaction") # each group in separate columns
raref <- raref[,-1] # remove the column with species names
raref <- data.frame(raref) # convert to data.frame
out <- iNEXT(raref, q = c(0, 1, 2), datatype = "abundance", endpoint = 300)

# Sample-size-based R/E curves, separating plots by "site"

# As it is difficult to edit the size of points and lines with the ggiNEXT() function,
# we extract the information from the 'out' object containing the rarefaction analysis results
# and use this dataframe to create the plot with ggplot:

df <- fortify(out)

df <- df %>% 
  mutate(
    Order.q = case_when(
      Order.q == 0 ~ "q = 0: Species richness",
      Order.q == 1 ~ "q = 1: Shannon",
      TRUE ~ "q = 2: Simpson"),
  )

df.point <- df[which(df$Method == "Observed"), ]
df.line <- df[which(df$Method != "Observed"), ]
df.line$Method <- factor(
  df.line$Method,
  c("Rarefaction", "Extrapolation"),
  c("Rarefaction", "Extrapolation")
)

(b <- ggplot(df, aes(x = x, y = y, colour = Assemblage)) +
    geom_point(aes(shape = Assemblage), size = 3, data = df.point) +
    geom_line(aes(linetype = Method), linewidth = 1, data = df.line) +
    geom_ribbon(aes(ymin = y.lwr, ymax = y.upr, fill = Assemblage, colour = NULL), alpha = 0.2) +
    scale_y_continuous(limits = c(0, 7), breaks = 0:7) +
    labs(x = "Number of individuals", y = "Species diversity") +
    facet_wrap(~Order.q) +
    theme_bw() +
    theme(
      legend.position = "none",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
)

## Fig. 2, with two panels ----

(a / b) + plot_annotation(tag_levels = 'a', tag_suffix = ")") & 
  theme(plot.tag = element_text(size = 10))

# ggsave("Fig. 2.png",
#        width = 5, height = 6, dpi = 600)


# Sapling density ----

## Total density ----

total <- saplings %>%
  group_by(Bloque, Clausura, Parcela) %>%
  dplyr::summarise(Abundancia = sum(Abundancia))

total$Clausura <- factor(total$Clausura, levels = c("Af", "Ad"))

total %>% 
  group_by(Clausura) %>% 
  dplyr::summarise(
    Abu_media = mean(Abundancia),
    Abu_sd = sd(Abundancia),
    n = n()
  )

m.abund <- lme(Abundancia ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = total)
anova(m.abund)

# Replace NAs with zero

m <- saplings %>% 
  pivot_wider(names_from = Especie, values_from = Abundancia, values_fill = 0)

## Schinopsis lorentzii density ----

m.sch <- lme(Q_colorado ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.sch)
summary(m.sch)

## Aspidosperma quebracho-blanco density ----

m.asp <- lme(Q_blanco ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.asp)
summary(m.asp)

## Sarcoramphus mistol density ----

m.mistol <- lme(Mistol ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.mistol)
summary(m.mistol)

## Libidibia paraguariensis density ----

m.guaya <- lme(Guayacan ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.guaya)
summary(m.guaya)

## Xymenia americana density ----

m.pata <- lme(Pata ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.pata)
summary(m.pata)

## Jodinia rhombifolia density ----

m.sombra <- lme(Sombra_de_toro ~ Clausura, random = ~ 1|Bloque/Clausura/Parcela, data = m)
anova(m.sombra)
summary(m.sombra)

## Fig. 3 ----

# Plotting effects with confidence intervals

t <- qt(0.025, 4, lower.tail = FALSE)

sps <- factor(c("Total", "Asp.que", "Sch.lor", "Sar.mis", "Lib.par", "Xim.ame", "Jod.rho"))
coef <- c(2.9, -4.6, 1.5, 0.3, -0.7, 0.3, 0.3) # changed signs to reflect the exclusion effect
se <- c(7.306104, 1.639688, 5.993316, 0.5700876, 0.8062258, 0.3000004, 0.5042679)
ee <- se * t

effects <- data.frame(sps, coef, ee)

effects$sps <- factor(effects$sps,
                      levels = c("Total", "Asp.que", "Sch.lor", "Sar.mis", "Lib.par", "Xim.ame", "Jod.rho"))

ggplot(efectos, aes(x = sps, y = coef)) +
  geom_errorbar(aes(ymin = coef - ee, ymax = coef + ee), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_vline(xintercept = 1.5, color = "grey") +
  geom_point(size = 5, show.legend = F) +
  coord_flip() +
  labs(
    x = "",
    y = "Exclusion effect size"
  ) +
  ylim(-25, 25) +
  # scale_y_continuous(
  #   # limits = c(0, 10),
  #   # expand = expansion(mult = c(0.1, 0.5))  # Ajustar la expansión del eje y
  # ) +
  annotate("text",
           label = "Lower densities\ninside the exclosure",
           x = 7.3, y = -14, size = 3, colour = "black"
  ) +
  annotate("text",
           label = "Higher densities\noutside the exclosure",
           x = 7.3, y = 14, size = 3, colour = "black"
  ) +
  theme_bw() +
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.y = element_text(size = 14, face = "italic"),
    axis.text.x = element_text(size = 14)
  )

# ggsave("Fig. Efectos.png",
#        width = 5,
#        height = 5,
#        dpi = 600)
  

# Height and diameter of "quebrachos" saplings ----

# Load the data
queb <- read_excel("data_Trigo_et_al.xlsx", sheet = "Quebrachos")

# Convert columns to factors
queb <- queb %>%
  mutate(across(c(Bloque, Clausura, Individuo, Especie), as.factor))

## Models for Aspidosperma quebracho-blanco ----

asp <- queb %>% filter(Especie == "Asp.que")

### Height model for Aspidosperma quebracho-blanco ----
m.altura.asp <- lme(Altura ~ Clausura,
                    random = ~ 1 | Bloque / Clausura / Individuo,
                    data = asp)
anova(m.altura.asp)
summary(m.altura.asp)

### Diameter model for Aspidosperma quebracho-blanco ----
m.diam.asp <- lme(Diámetro ~ Clausura, random = ~ 1 | Bloque / Clausura / Individuo, data = asp)
anova(m.diam.asp)
summary(m.diam.asp)

## Models for Schinopsis lorentzii ----

sch <- queb %>% filter(Especie == "Sch.lor")

### Height model for Schinopsis lorentzii ----
m.altura.sch <- lme(Altura ~ Clausura,
                    random = ~ 1 | Bloque / Clausura / Individuo,
                    data = sch)
anova(m.altura.sch)
summary(m.altura.sch)

### Diameter model for Schinopsis lorentzii ----
m.diam.sch <- lme(Diámetro ~ Clausura, random = ~ 1 | Bloque / Clausura / Individuo, data = sch)
anova(m.diam.sch)
summary(m.diam.sch)

## Effect plots for both quebrachos ----

# Prepare data for plots
sp <- factor(c("Asp.que", "Sch.lor"))
coef.alt <- c(13.4, 49.28)
coef.diam <- c(0.1, 0.472) # Changed signs to reflect the exclusion effect
# Diameters are converted to cm to keep units consistent
se.alt <- c(16.79941, 7.603008)
se.diam <- c(0.1773122, 0.5875981) 
ee.alt <- se.alt * t
ee.diam <- se.diam * t

df.alt <- data.frame(sp, coef.alt, ee.alt)
df.diam <- data.frame(sp, coef.diam, ee.diam)

# Plot for sapling height
ef.altura <- ggplot(df.alt, aes(x = sp, y = coef.alt)) +
  geom_errorbar(aes(ymin = coef.alt - ee.alt, ymax = coef.alt + ee.alt), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5, show.legend = FALSE) +
  coord_flip() +
  labs(title = "Sapling height (cm)",
       x = "",
       y = "Exclusion effect size") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.y = element_text(size = 12, face = "italic"),
        axis.text.x = element_text(size = 12)) +
  ylim(-75, 75)

# Plot for sapling diameter
ef.diam <- ggplot(df.diam, aes(x = sp, y = coef.diam)) +
  geom_errorbar(aes(ymin = coef.diam - ee.diam, ymax = coef.diam + ee.diam), width = 0.1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_point(size = 5, show.legend = FALSE) +
  coord_flip() +
  labs(title = "Sapling diameter (cm)",
       x = "",
       y = "Exclusion effect size") +
  theme_bw() +
  theme(axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_blank()) +
  ylim(-3, 3)

# Combine the plots
ef.altura + ef.diam +
  plot_layout(axis_titles = "collect")

# Save the figure
# ggsave("Fig. Effects_height_and_diameter.png",
#        width = 6, height = 4, dpi = 600)



