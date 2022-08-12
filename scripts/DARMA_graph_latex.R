library(ggplot2)
library(dplyr)
library(latex2exp)

# Carregando simulações
miceadds::load.Rdata("./simulacoes/simun100CG.RData", 'simu_100')
miceadds::load.Rdata("./simulacoes/simun250CG.RData", 'simu_250')
miceadds::load.Rdata("./simulacoes/simun500CG.RData", 'simu_500')

miceadds::load.Rdata("./simulacoes/simu_ARMA(2,1)100CGp2.RData", 'simu_100_ar2')
miceadds::load.Rdata("./simulacoes/simu_ARMA(2,1)250CGp2.RData", 'simu_250_ar2')
miceadds::load.Rdata("./simulacoes/simu_ARMA(2,1)500CGp2.RData", 'simu_500_ar2')

simu_darma = vector(mode='list', 6)

#Unindo as simulações
simu_darma[[1]] = simu_100 %>% mutate(param1 = 0)
simu_darma[[2]] = simu_250 %>% mutate(param1 = 0)
simu_darma[[3]] = simu_500 %>% mutate(param1 = 0)
simu_darma[[4]] = simu_100_ar2
simu_darma[[5]] = simu_250_ar2
simu_darma[[6]] = simu_500_ar2


#Valores dos parâmetros
theta0 <- c(0.5, 0.8, 3.5, 11.2, 3, 5, 3, 0)
theta_ar2 <- c(0.5,-0.4, 0.8, 3.5, 11.2, 3, 5, 3)

n_ = c(100,250,500,100,250,500)

# Métricas para avaliar as simulações
resultados_darma = data.frame(matrix(ncol = 8, nrow  = 0))
colnames(resultados_darma) = c('mean', 'sd' ,'parametros', 'Vies', 'EQM',
                               'Assimetria', 'Curtose', 'N')

j = 1
for(df in simu_darma){
  if(j<=3){
    theta = theta0
  }else{
    theta = theta_ar2
  }
  df_ = df |> 
    colMeans() %>%
    cbind(apply(df, 2, sd)) %>%
    as.data.frame() |> 
    rename(mean = `.`, sd = V2) |>
    mutate(parametros = theta) |>
    mutate(Vies = mean - parametros) |>
    mutate(EQM = Vies^2 + sd^2,
           Assimetria = moments::skewness(df),
           Curtose = moments::kurtosis(df),
           N = rep(n_[j], length(theta)))
  resultados_darma = resultados_darma %>%
    add_row(df_)
  j = j +1
}

resultados_darma = resultados_darma %>% 
  mutate(estrutura = rep(c("ARMA(1,1)", "ARMA(2,1)"), c(24,24))) %>%
  rename(Média = mean, Viés = Vies, DP = sd, EQM = EQM, CA = Assimetria, 
         CC = Curtose, N = N)

#Construindo latex
resultados_darma %>% filter(estrutura == "ARMA(2,1)", N == 100) %>%
  select(-parametros, -N, - estrutura) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
resultados_darma %>% filter(estrutura == "ARMA(2,1)", N == 250) %>% 
  select(-parametros, -N, - estrutura) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
resultados_darma %>% filter(estrutura == "ARMA(2,1)", N == 500) %>% 
  select(-parametros, -N, - estrutura) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
#Construindo gráfico

# resultados_darma %>%
#   # group_by(parametros, N) %>% # Also group_by(across(-comment)) would work with the example
#   # slice_tail() %>%
#   # ungroup() %>%
# filter(parametros == 5) %>%
resultados_darma %>% 
  filter(parametros == 0.5 | parametros == 0.8 | parametros == -0.4) %>%
  mutate(estrutura_parametros = paste0(estrutura, parametros, sep = ' ')) %>%
  ggplot() +
  geom_line(aes(x=N, y=CC, color = estrutura_parametros), 
            size = 2) +
  geom_point(aes(x=N, y=CC, color = estrutura_parametros), size = 3) +
  scale_color_brewer(labels = unname(TeX(c("$\\widehat{\\phi}$ - C1", 
                                             "$\\widehat{\\theta}$  - C1",
                                             "$\\widehat{\\phi_2}$ - C2", 
                                             "$\\widehat{\\phi_1}$ - C2",
                                             "$\\widehat{\\theta}$ - C2"))),
                     palette = 'Paired') +
  theme_minimal() +
  labs(x = 'Tamanho amostral', y= "Coeficiente de Curtose", 
       color = "Estimadores")

# resultados_darma %>%
#   group_by(parametros, N) %>% # Also group_by(across(-comment)) would work with the example
#   slice_head() %>%
#   ungroup()


