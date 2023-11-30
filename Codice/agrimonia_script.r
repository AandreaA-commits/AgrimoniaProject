install.packages('tseries')
install.packages('imputeTS')
install.packages('GGPlot')
library(tseries)
library(imputeTS)
library(ggplot2)


load("Agrimonia_Dataset_v_3_0_0.Rdata")
subset <- AgrImOnIA_Dataset_v_3_0_0[AgrImOnIA_Dataset_v_3_0_0$IDStations=="677",]

#AQ_nh3 variabile risposta, si può usare l'altra come regressore
AQ_nh3 <- subset$AQ_nh3

#numero di valori nan al''interno della response variable AQ_NH3
apply ( X = is.na(subset), MARGIN = 2, FUN = sum)

#kalman smoother per missing values
Y <- na_kalman(AQ_nh3, model = "StructTS", smooth = TRUE)

#plotting della serie storica
plot(AQ_nh3)
plot(Y)

#istogramma di AQ_nh3
hist(Y)

#istogramma del logaritmo di AQ_nh3
hist(log(Y))

#test normalità della Y
jarque.bera.test(Y)

#test normalità della Y
jarque.bera.test(log(Y))

#autocorrelation function -> per ordine dell'MA
acf(Y)

#partial autocorrelation function -> perordine dell'AR
pacf(Y)

#test di stazionarietà della response variable
adf.test(Y) #stazionario

#calcolo della varianza e media di AQ_nh3
var_AQ_nh3 <- var(Y)
mean_AQ_nh3 <- mean(Y)


#plot
# Creazione del grafico ggplot2
grafico <- ggplot(subset, aes(x = subset$Time, y = Y)) +
  geom_line() +
  geom_hline(yintercept = mean_AQ_nh3, linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept = mean_AQ_nh3 + 2*sqrt(var_AQ_nh3), linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = mean_AQ_nh3 - 2*sqrt(var_AQ_nh3), linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = mean_AQ_nh3 + 3*sqrt(var_AQ_nh3), linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept = mean_AQ_nh3 - 3*sqrt(var_AQ_nh3), linetype = "dashed", color = "blue", size = 1) +
  labs(title = "Serie Storica con Media e Varianza",
       x = "Data",
       y = "Valore") +
  theme_minimal()

# Visualizzazione del grafico
print(grafico)

#standardizzazione e logaritmo dei dati
Y <- log(Y)
Y_std <- (Y-mean_AQ_nh3)/sqrt(var_AQ_nh3)

#matrici dei regressori
X <- subset[, c(13:28)]
names <- colnames(X)
X <- model.matrix (reformulate(names) , subset)


#stima del modello con regressori del tempo

mdl <- list()

for (i in 1:4){
  for (j in 1:4){
    mdl <- arima(Y, order = c(i, 0, j))
  }
}



