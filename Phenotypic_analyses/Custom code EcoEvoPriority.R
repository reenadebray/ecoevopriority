## Code for calculating strength of priority effects
# Reena Debray
# 2022

# Pantoea dispersa
## invasion success - Pantoea - ancestral
Di_ji<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_focalP0D3","focal_CFU"]+1)
Di_0i<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP0D0_MgD0","focal_CFU"]+1)
Di_ij<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP0D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP0D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_Pan_anc<-d$Var1-d$Var2

## invasion success - Pantoea - evolved
Di_ji<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_focalP6D3","focal_CFU"]+1)
Di_0i<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP6D0_MgD0","focal_CFU"]+1)
Di_ij<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP6D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP6D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_Pan_evo<-d$Var1-d$Var2

## resistance - Pantoea - ancestral
Di_ji<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP0D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_focalP0D3","competitor_CFU"]+1)
Di_i0<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_Pan_anc<-d$Var1-d$Var2

## resistance - Pantoea - evolved
Di_ji<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="focalP6D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_focalP6D3","competitor_CFU"]+1)
Di_i0<-log(Pantoea_dispersa[Pantoea_dispersa$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_Pan_evo<-d$Var1-d$Var2

# Pseudomonas protegens
## invasion success - protegens - ancestral
Di_ji<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_focalP0D3","focal_CFU"]+1)
Di_0i<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP0D0_MgD0","focal_CFU"]+1)
Di_ij<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP0D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP0D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_prot_anc<-d$Var1-d$Var2

## invasion success - protegens - evolved
Di_ji<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_focalP6D3","focal_CFU"]+1)
Di_0i<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP6D0_MgD0","focal_CFU"]+1)
Di_ij<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP6D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP6D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_prot_evo<-d$Var1-d$Var2

## resistance - protegens - ancestral
Di_ji<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP0D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_focalP0D3","competitor_CFU"]+1)
Di_i0<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_prot_anc<-d$Var1-d$Var2

## resistance - protegens - evolved
Di_ji<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="focalP6D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_focalP6D3","competitor_CFU"]+1)
Di_i0<-log(Pseudomonas_protegens[Pseudomonas_protegens$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_prot_evo<-d$Var1-d$Var2

# Pseudomonas syringae

## invasion success - syringae - ancestral
Di_ji<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_focalP0D3","focal_CFU"]+1)
Di_0i<-log(EEPE_syr[EEPE_syr$treatment=="focalP0D0_MgD0","focal_CFU"]+1)
Di_ij<-log(EEPE_syr[EEPE_syr$treatment=="focalP0D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(EEPE_syr[EEPE_syr$treatment=="focalP0D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_syr_anc<-d$Var1-d$Var2

## invasion success - syringae - evolved
Di_ji<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_focalP6D3","focal_CFU"]+1)
Di_0i<-log(EEPE_syr[EEPE_syr$treatment=="focalP6D0_MgD0","focal_CFU"]+1)
Di_ij<-log(EEPE_syr[EEPE_syr$treatment=="focalP6D0_competitorD3","focal_CFU"]+1)
Di_i0<-log(EEPE_syr[EEPE_syr$treatment=="focalP6D0_MgD3","focal_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
inv_syr_evo<-d$Var1-d$Var2

## resistance - syringae - ancestral
Di_ji<-log(EEPE_syr[EEPE_syr$treatment=="focalP0D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_focalP0D3","competitor_CFU"]+1)
Di_i0<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_syr_anc<-d$Var1-d$Var2

## resistance - syringae - evolved
Di_ji<-log(EEPE_syr[EEPE_syr$treatment=="focalP6D0_competitorD3","competitor_CFU"]+1)
Di_0i<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_MgD0","competitor_CFU"]+1)
Di_ij<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_focalP6D3","competitor_CFU"]+1)
Di_i0<-log(EEPE_syr[EEPE_syr$treatment=="competitorD0_MgD3","competitor_CFU"]+1)
d<-expand.grid(Di_ji,Di_ij)
res_syr_evo<-d$Var1-d$Var2
