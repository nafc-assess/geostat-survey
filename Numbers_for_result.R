library(tidyr)
library(dplyr)
library(here)

load(here("Data", "index_all_scenarios.Rdata"))

### Calculating MRE, RMSLE, and pulling AIC
index_all_scenarios_accuracy <-
  index_all_scenarios |>
  filter(year > 10) |>
  ungroup() |>
  filter(type != "Bootstrapped") |>
  group_by(pop, type, scenario, species) |>
  summarise(MRE = mean((log(N) - log(true)) / log(true)),
            RMSE= sqrt(mean((log(N) - log(true))^2)),
            AIC = mean(AIC))

#### Base scenario

# maximum MRE of the model-based estimators for cod-like
index_all_scenarios_accuracy |>
  filter(species == "Cod-like") |>
  filter(type != "Design-based") |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  group_by(type) |>
  summarise(mMRE=mean(MRE)) |>
  slice_max(abs(mMRE))

# For the yellowtail-like species, compared to the design-based approach,
#models without a depth covariate had slightly worse accuracy (mean XX\% higher RMSE)
#and were slightly biased (mean MRE = YY), while models with a depth covariate had
#improved accuracy and were effectively unbiased

m_mRMSE <- index_all_scenarios_accuracy |> # mean for models without depth
  filter(species == "Yellowtail-like") |>
  filter(type != "Design-based") |>
  filter(type != "Bootstrapped") |>
  filter(!grepl("Depth", type)) |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  #group_by(species, type, scenario) |>
  summarise(mRMSE=mean(RMSE))
m_mRMSE

d_mRMSE <-  index_all_scenarios_accuracy |> # mean for design
  filter(species == "Yellowtail-like") |>
  filter(type == "Design-based") |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  #group_by(species, type, scenario) |>
  summarise(mRMSE=mean(RMSE))
d_mRMSE

#result_table |> filter(scenario == "Base") |> filter(type == "Design-based")

m_mRMSE/d_mRMSE * 100 - 100 ## XX
m_mRMSE * 100 /d_mRMSE  - 100


index_all_scenarios_accuracy |>
  filter(species == "Yellowtail-like") |>
  filter(type == "Design-based") |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  summarise(mean(abs(MRE)))

index_all_scenarios_accuracy |>
  filter(species == "Yellowtail-like") |>
  filter(type != "Design-based") |>
  filter(type != "Bootstrapped") |>
  filter(!grepl("Depth", type)) |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  summarise(mean(MRE)) ## YY


#### Set reduction

#The 30\% set reduction scenario showed a similar pattern to the base scenario for both species,
#although with slightly reduced accuracy (mean 23\% and 18\% higher RMSE for cod-like and yellowtail-like, respectively

r30_numbers <- merge(index_all_scenarios_accuracy |>
  filter(type != "Bootstrapped") |>
  filter(scenario ==  "30% set reduction") |>
  ungroup() |>
  group_by(species, type) |>
  summarise(mRMSE_r30 = mean(RMSE)),

index_all_scenarios_accuracy |>
  filter(type != "Bootstrapped") |>
  filter(scenario ==  "Base") |>
  ungroup() |>
  group_by(species, type) |>
  summarise(mRMSE_base = mean(RMSE)))

r30_numbers <- r30_numbers |>
  mutate(diff = mRMSE_r30 - mRMSE_base,
         perc_change = diff * 100 / mRMSE_base)

r30_numbers |>
  #group_by(species) |>
  summarise(mean(perc_change))

#Confidence interval coverage was slightly reduced (mean 4.5\% reduction)
#for the design-based bootstrap confidence intervals compared to the base case

CI_coverage <- index_all_scenarios |>
  filter(type !=  "Design-based") |>
  filter(type !=  "Bootstrapped" | scenario !=  "Strata removal") |>
  mutate(covered = lwr < true & upr > true) |>
  group_by(type, scenario, species) |>
  summarise(mc = mean(covered)) |>
  group_by(scenario) |>
  arrange(mc)

CI_coverage |>
  filter(scenario == "Base") |>
  filter(species == "Cod-like")

CI_coverage |>
  filter(scenario == "Base") |>
  filter(species == "Yellowtail-like")

CI_coverage |>
  filter(scenario == "30% set reduction") |>
  filter(species == "Cod-like")

CI_coverage |>
  filter(scenario == "30% set reduction") |>
  filter(species == "Yellowtail-like")

CI_coverage |>
  filter(type == "Bootstrapped") |>
  filter(scenario == "Base" | scenario == "30% set reduction")

r30_CIs <- merge(CI_coverage |>
                       filter(type == "Bootstrapped") |>
                       filter(scenario ==  "30% set reduction") |>
                       ungroup() |>
                       group_by(species, type) |>
                       summarise(mc_r30 = mean(mc)),

                 CI_coverage |>
                   filter(type == "Bootstrapped") |>
                   filter(scenario ==  "Base") |>
                   ungroup() |>
                   group_by(species, type) |>
                   summarise(mc_base = mean(mc)))

r30_CIs <- r30_CIs |>
  mutate(diff = mc_r30  - mc_base,
         perc_change = diff * 100 / mc_base)

mean(r30_CIs$perc_change)

#### Area blocked reduction scenario

#With a blocked area removed from the survey domain,
#the design-based estimator remained unbiased although was slightly less accurate than
#the best geostatistical models (0.4 vs. 0.3 for cod-like and 0.3 vs. 0.17 for yellowtail-like)
#in this scenario for both species

## best model

b30_RMSE <- rbind(index_all_scenarios_accuracy |>
                  group_by(type, species, scenario) |>
                  summarise(mMRE = mean(RMSE)) |>
                  filter(scenario ==  "30% area blocked") |>
                   group_by(species, type) |>
                   summarise(mRMSE_b30 = mean(mMRE)) |>
                  filter(type != "Design-based")|>
                  group_by(species) |>
                  slice_min(abs(mRMSE_b30)),

                index_all_scenarios_accuracy |>
                  group_by(type, species, scenario) |>
                  summarise(mMRE = mean(RMSE)) |>
                  filter(scenario ==  "30% area blocked") |>
                  group_by(species, type) |>
                  summarise(mRMSE_b30 = mean(mMRE)) |>
                  filter(type == "Design-based"))
b30_RMSE


