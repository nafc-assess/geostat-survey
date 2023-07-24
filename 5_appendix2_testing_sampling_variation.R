#Appendix 2 --- Testing sampling variation

library(here)
load(here("data2", "index_cod_all_scenarios_10x100.Rdata"))
load(here("data2", "index_cod_all_scenarios_100x10.Rdata"))

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggpubr)

### Calculating statistics for 100 populations x 10 surveys [Process error]

head(index_cod_all_scenarios_100x10)

index_cod_base_acc_samp_myear_100x10 <- # mean accuracy calculations for all years
  index_cod_all_scenarios_100x10 |>
  ungroup() |>
  group_by(pop, sim, type, scenario, species) |>
  summarise(MRE = mean((N - true) / true),
            RSME =  sqrt(mean((log(N) - log(true))^2)))

### Calculating statistics for 10 populations x 100 surveys [Sampling error]

head(index_cod_all_scenarios_10x100)

index_cod_base_acc_samp_myear_10x100 <- # mean accuracy calculations for all years
  index_cod_all_scenarios_10x100 |>
  ungroup() |>
  group_by(pop, sim, type, scenario, species) |>
  summarise(MRE = mean((N - true) / true),
            RSME =  sqrt(mean((log(N) - log(true))^2)))


### combining IQR plots

IQR_table <- rbind(
  index_cod_base_acc_samp_myear_100x10 |>
    group_by(sim, type) |>
    summarise(IQR_MRE = IQR(MRE),
              IQR_RSME = IQR(RSME)) |>
    group_by(type) |>
    summarise(mean_IQR_MRE = mean(IQR_MRE),
              mean_IQR_RSME = mean(IQR_RSME)) |>
    mutate(run = "Process + Observarion error"), #100 populations x 10 surveys
  index_cod_base_acc_samp_myear_10x100 |>
    group_by(pop, type) |>
    summarise(IQR_MRE = IQR(MRE),
              IQR_RSME = IQR(RSME))|>
    group_by(type) |>
    summarise(mean_IQR_MRE = mean(IQR_MRE),
              mean_IQR_RSME = mean(IQR_RSME))|>
    mutate(run = "Observation error")) #10 populations x 100 surveys


a <- ggplot(IQR_table) +
  geom_point(aes(mean_IQR_RSME, type, colour = run)) +
  theme_bw()+
  labs(x = "Mean IQR for RMSE", y = "", colour = "Error type", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs("")+
  xlim(0,0.1) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

b <- ggplot(IQR_table) +
  geom_point(aes(mean_IQR_MRE, type, colour = run)) +
  theme_bw()+
  labs(x = "Mean IQR for MRE", y = "", colour = "Error type", fill = "Estimator") +
  theme(text = element_text(size = 14))+
  theme(strip.background = element_rect(fill = "grey97")) +
  labs("")+
  xlim(0,NA)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

figure_apx_iqr <- ggarrange(a, b, ncol = 2, nrow = 1, legend = "right", common.legend = TRUE)
figure_apx_iqr

ggsave("data2/figure_apx_iqr.pdf", plot = figure_apx_iqr, width = 10, height = 6, units = "in", dpi = 500, bg = "white")


#### numeric


IQR_table_cbind <- merge(
  index_cod_base_acc_samp_myear_100x10 |>
    group_by(sim, type) |>
    summarise(IQR_MRE = IQR(MRE),
              IQR_RSME = IQR(RSME)) |>
    group_by(type) |>
    summarise(mean_IQR_MRE_PO = mean(IQR_MRE),
              mean_IQR_RSME_PO = mean(IQR_RSME)),
  index_cod_base_acc_samp_myear_10x100 |>
    group_by(pop, type) |>
    summarise(IQR_MRE = IQR(MRE),
              IQR_RSME = IQR(RSME))|>
    group_by(type) |>
    summarise(mean_IQR_MRE_O = mean(IQR_MRE),
              mean_IQR_RSME_O = mean(IQR_RSME)), by="type")


IQR_table_cbind |>
  group_by(type) |>
  mutate(MRE_PO_O = mean_IQR_MRE_PO - mean_IQR_MRE_O,
         RMSE_PO_O = mean_IQR_RSME_PO  - mean_IQR_RSME_O) |>
  ungroup()|>
  summarise(mean(MRE_PO_O), mean(RMSE_PO_O))

