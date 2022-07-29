library(ggplot2)
library(dplyr)

# Carregando simulações
miceadds::load.Rdata("./simulacoes/simun100CG.RData", 'simu_100')
miceadds::load.Rdata("./simulacoes/simun250CG.RData", 'simu_250')
miceadds::load.Rdata("./simulacoes/simun500CG.RData", 'simu_500')

simu_darma = vector(mode='list', 3)

#Unindo as simulações
simu_darma[[1]] = simu_100
simu_darma[[2]] = simu_250
simu_darma[[3]] = simu_500

#Valores dos parâmetros
theta0 <- c(0.5, 0.8, 3.5, 11.2, 3, 5, 3)
n_ = c(100,250,500)

# Métricas para avaliar as simulações
resultados_darma = data.frame(matrix(ncol = 8, nrow  = 0))
colnames(resultados_darma) = c('mean', 'sd' ,'parametros', 'Vies', 'EQM',
                               'Assimetria', 'Curtose', 'N')

j = 1
for(df in simu_darma){
  df_ = df |> 
    colMeans() %>%
    cbind(apply(df, 2, sd)) %>%
    as.data.frame() |> 
    rename(mean = `.`, sd = V2) |>
    mutate(parametros = theta0) |>
    mutate(Vies = mean - parametros) |>
    mutate(EQM = Vies^2 + sd^2,
           Assimetria = moments::skewness(df),
           Curtose = moments::kurtosis(df),
           N = rep(n_[j], length(theta0))) 
  resultados_darma = resultados_darma %>%
    add_row(df_)
  j = j +1
}

resultados_darma = resultados_darma %>% rename(Média = mean, Viés = Vies, 
                                               DP = sd, EQM = EQM, 
                                               CA = Assimetria, CC = Curtose, 
                                               N = N)

#Construindo latex
resultados_darma %>% distinct(parametros,N) #%>% filter(N == 100) 
  select(-parametros, -N) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
resultados_darma %>% filter(N == 250) %>% 
  select(-parametros, -N) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
resultados_darma %>% filter(N == 500) %>% 
  select(-parametros, -N) %>% round(3) %>%
  mypdf1::pdf1_tbl(format='latex')
#Construindo gráfico

resultados_darma %>%
  # group_by(parametros, N) %>% # Also group_by(across(-comment)) would work with the example
  # slice_tail() %>%
  # ungroup() %>%
filter(parametros == 5) %>%
  ggplot() +
  geom_line(aes(x=N, y=Viés, color = "Viés"), 
            size = 2) +
  geom_line(aes(x=N, y=CA, color = "Assimetria"), 
            size = 2) +
  geom_line(aes(x=N, y=CC, color = "Curtose"), 
            size = 2) +
  geom_point(aes(x=N, y=Viés, color = 'Viés'), size = 3) + 
  geom_point(aes(x=N, y=CA, color = 'Assimetria'), size = 3) + 
  geom_point(aes(x=N, y=CC, color = 'Curtose'), size = 3) + 
  scale_colour_manual("Medidas", 
                      values = c("Viés" = "#003f5c",
                                 "Assimetria" = "#933494", 
                                 "Curtose" = "#ff0d0d")) + 
  theme_minimal() +
  labs(x = 'Tamanho amostral', y= "Valores das Medidas", colour = "Quantis")

# resultados_darma %>%
#   group_by(parametros, N) %>% # Also group_by(across(-comment)) would work with the example
#   slice_head() %>%
#   ungroup()


