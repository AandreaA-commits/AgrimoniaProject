#install.packages('tseries')
#install.packages('imputeTS')
#install.packages('GGPlot')
#install.packages('Matrix')
install.packages('clusterSim')
install.packages("lmtest")
install.packages('forecast')
library(tseries)
library(imputeTS)
library(ggplot2)
library(Matrix)
library(clusterSim)
library(lmtest)
library(forecast)

load("Agrimonia_Dataset_v_3_0_0.Rdata")
subset <- AgrImOnIA_Dataset_v_3_0_0[AgrImOnIA_Dataset_v_3_0_0$IDStations=="677",]

#AQ_nh3 variabile risposta, si può usare l'altra come regressore
AQ_nh3 <- subset$AQ_nh3

#logaritmo di AQ_nh3
AQ_nh3_log <- log(AQ_nh3)

#numero di valori nan al''interno della response variable AQ_NH3
apply ( X = is.na(subset), MARGIN = 2, FUN = sum)

#hist di log aq_nh3
hist(AQ_nh3_log)

#Eliminazione degli outliers (fuori da IC al 95%)
q_995 <- quantile(AQ_nh3_log, probs = 0.995, na.rm = TRUE)
q_005 <- quantile(AQ_nh3_log, probs = 0.005, na.rm = TRUE)

# metto Na i valori che non sono all'interno del 99o percentile
sum(is.na(AQ_nh3_log))
AQ_nh3_log[AQ_nh3_log > q_995 | AQ_nh3_log < q_005] <- NA
sum(is.na(AQ_nh3_log))

#kalman smoother per missing values del log
Y <- na_kalman(AQ_nh3_log, model = "StructTS", smooth = TRUE)

grafico <- ggplot(subset, aes(x = subset$Time, y = Y)) +
  geom_line() +
  geom_hline(yintercept = mean(Y), linetype = "dashed", color = "red", size = 1) +
  geom_hline(yintercept =  + q_995, linetype = "dashed", color = "blue", size = 1) +
  geom_hline(yintercept =  + q_005, linetype = "dashed", color = "blue", size = 1) +
  labs(title = "Serie Storica con Media e Varianza",
       x = "Data",
       y = "Valore") +
  theme_minimal()

# Visualizzazione del grafico
print(grafico)

Y_norm <- (Y-mean(Y))/sd(Y)

#plotting della serie storica
plot(Y_norm, type = 'l')

#istogramma di AQ_nh3
hist(Y_norm)

#test normalità della Y
jarque.bera.test(Y_norm)


#autocorrelation function -> per ordine dell'MA
acf(Y_norm)
acf(Y_norm, na.action = na.pass, lag.max = length(Y))

#partial autocorrelation function -> perordine dell'AR
pacf(Y_norm)
pacf(Y_norm, na.action = na.pass, lag.max = length(Y))
#test di stazionarietà della response variable
adf.test(Y_norm) #rifiuta ipotesi con p-value = 0 => stazionario

#matrici dei regressori
X <- subset[, c(13:28)];
names <- colnames(X);
X <- model.matrix (reformulate(names) , subset);

X_norm <- data.Normalization(X, type = "n1");
X_norm <- X_norm[, -1];


#stima del modello con regressori del tempo

p <- 2;
q <- 2;
P <- 1;
Q <- 1;

mdl <- arima(Y_norm, order = c(p, 0, q), xreg = X_norm, seasonal = list(order = c(P, 0, Q), period = 7))
Box.test(mdl$residuals, lag = 20, type = "Ljung-Box")
#verifica stazionarietà dei residui
res <- mdl$residuals

pacf(mdl$residuals, lag.max = length(Y_norm))
acf(mdl$residuals, lag.max = length(Y_norm))

R2 <- 1 - var(res)/var(Y_norm)

plot(res)
media <- mean(res, na.rm = TRUE)
adf.test(res)
hist(res)

#verifica normalità dei residui
jarque.bera.test(res) #non è normale

# calcolo dei p value
coeftest(mdl)

idx_da_usare <- 1:31;
idx_usati <- list();

models_sw <- list();
aic_sw <- list();


for (i in 1:length(idx_da_usare)){
  aic_temp <- list();
  for (j in 1:length(idx_da_usare)){
    idx <- c(idx_usati, idx_da_usare[j])
    print(as.numeric(idx))
    mdl <- arima(Y_std, order = c(1, 0, 1), xreg = X_norm[, as.numeric(idx)], )
    aic_temp <- c(aic_temp, mdl$aic)
  }
  
  indice <- which.min(aic_temp)
  print("ciao")
  print(indice)
  idx_usati <<- c(idx_usati, idx_da_usare[indice])
  idx_da_usare <<- idx_da_usare[-indice]
  aic_sw <<- c(aic_sw, aic_temp[indice])
}

best <- which.min(aic_sw)

#modello migliore dopo sw
best_mdl_sw <- arima(Y_std, order = c(1,0,1), xreg = X_norm[, as.numeric(idx[1:best])])

res_sw <- best_mdl_sw$residuals

res_sw_prova <- na_kalman(res_sw, model = "StructTS", smooth = TRUE)

adf.test(res_sw_prova)
jarque.bera.test(res_sw_prova)
plot(res_sw_prova)
mean(res_sw_prova)
acf(res_sw_prova)
pacf(res_sw_prova)

#staz residui iniziale
#analizzare modello migliore sw
#bootstrap epr validazione (vedere normalità dei residui)

mdl <- list()

max_p <- 4
max_q <- 4
max_P <- 2  # Massimo ordine stagionale per P
max_Q <- 2  # Massimo ordine stagionale per Q
m <- 12     # Frequenza stagionale

# Loop per esplorare le combinazioni di ordini AR, MA e stagionali
mdl <- list()
for (i in 1:max_p) {
  for (j in 1:max_q) {
    for (k in 0:max_P) {
      for (l in 0:max_Q) {
        stimato <- arima(Y_std, order = c(i, 0, j), seasonal = list(order = c(k, 0, l), period = 12), xreg = X_norm)
        mdl <- c(mdl, stimato$aic);
      }
    }
  }
}

#utilizzo auto.arimo
auto <- auto.arima(Y_std, max.p=5, max.q = 5, max.P = 2, max.Q = 2, start.P = 0, start.Q = 0, start.p = 1, start.q = 1, 
           stationary = TRUE, seasonal = TRUE, ic=c("aicc", "aic", "bic"), stepwise=TRUE, trace=FALSE,
           xreg = X_norm, seasonal.test="ch", test=c("kpss","adf","pp"),
           allowdrift=TRUE, allowmean=TRUE, lambda=NULL, biasadj=FALSE)


arima(Y_std, )

