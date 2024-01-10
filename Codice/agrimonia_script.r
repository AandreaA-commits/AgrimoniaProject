#install.packages('tseries')
#install.packages('imputeTS')
#install.packages('GGPlot')
#install.packages('Matrix')
#install.packages('clusterSim')
#install.packages("lmtest")
#install.packages('forecast')
#install.packages("https://cran.r-project.org/src/contrib/Archive/rlang/rlang_1.1.1.tar.gz", repos = NULL, type="source")
#install.packages("devtools")
#devtools::install_github("lamferzon/block-bootstrap-for-R")
installed.packages("forecast")
install.packages("RSNNS")
install.packages("e1071")
library(e1071)
library(RSNNS)
library(forecast)
library(bboot)
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
P <- 0;
Q <- 0;

mdl <- arima(Y_norm, order = c(p, 0, q), xreg = X_norm, seasonal = list(order = c(P, 0, Q), period = 7), include.mean = FALSE)
Box.test(mdl$residuals, lag = 20, type = "Ljung-Box") #non correlati
#verifica stazionarietà dei residui
res <- mdl$residuals

pacf(mdl$residuals, lag.max = length(Y_norm))
acf(mdl$residuals, lag.max = length(Y_norm))

R2 <- 1 - var(res)/var(Y_norm)

plot(res)
media <- mean(res, na.rm = TRUE)
adf.test(res) #stazionario
hist(res)

#verifica normalità dei residui
jarque.bera.test(res) #non è normale

# calcolo dei p value
coeftest(mdl)

idx_da_usare <- 1:31;
idx_usati <- list();

models_sw <- list();
R2_sw <- list();


for (i in 1:length(idx_da_usare)){
  R2_temp <- list();
  for (j in 1:length(idx_da_usare)){
    idx <- c(idx_usati, idx_da_usare[j])
    #print(as.numeric(idx))
    mdl <- arima(Y_norm, order = c(p, 0, q), xreg = X_norm[, as.numeric(idx)], )
    R2 <- 1- var(mdl$residuals)
    R2_temp <- c(R2_temp, R2)
  }
  
  indice <- which.max(R2_temp)
  idx_usati <<- c(idx_usati, idx_da_usare[indice])
  idx_da_usare <<- idx_da_usare[-indice]
  R2_sw <<- c(R2_sw, R2_temp[indice])
}

aic_best_sw <- list()

for (i in 1:length(idx_usati)){
  mdl <- arima(Y_norm, order = c(p, 0, q), xreg = X_norm[, as.numeric(idx_usati[1:i])], include.mean = FALSE)
  aic_best_sw <<- c(aic_best_sw, mdl$aic)
}

idx_aic <- which.min(aic_best_sw)

best_mdl_sw <- arima(Y_norm, order = c(p, 0, q), xreg = X_norm[, as.numeric(idx_usati[1:idx_aic])], include.mean = FALSE)
R2 <- 1- var(best_mdl_sw$residuals)
coeftest(best_mdl_sw)


#modello migliore dopo sw
res_sw <- best_mdl_sw$residuals

adf.test(res_sw)
jarque.bera.test(res_sw)
plot(res_sw)
mean(res_sw)
acf(res_sw)
pacf(res_sw)

Box.test(res_sw, lag = 20, type = "Ljung-Box") #non correlati

#staz residui iniziale
#analizzare modello migliore sw
#bootstrap per validazione (vedere normalità dei residui)

#SELEZIONE DEI COEFFICIENTI SIGNIFICATIVI IN TRAIN
N <- 1000 #N è la lunghezza di ciascuna simulazione
K <- 100 #numero di ripetizioni
L <- 100 #lunghezza media dei blocchi

idx_boot <- blockboot(Y_norm, N, K, L, l_gen = "poisson")

col_name_reg <- colnames(X_norm[, as.numeric(idx_usati[1:idx_aic])])
col_name <- names(best_mdl_sw$coef)

matrix_beta <- matrix(ncol = length(col_name), nrow = K)

for(i in 1:K){
  mdl <- arima(Y_norm[idx_boot[,i]], order = c(p, 0, q), xreg = X_norm[idx_boot[,i], as.numeric(idx_usati[1:idx_aic])], include.mean = FALSE)
  matrix_beta[i, 1:length(col_name)] = as.vector(mdl$coef)
  print(i)
}

#distribuzioni e analisi significatività
#calcolo degli intervalli di confidenza
IC_up <- numeric(length(col_name))
IC_down <- numeric(length(col_name))
medie_coeff <- setNames(numeric(length(col_name)), col_name)

isSignificativo<- setNames(numeric(length(col_name)), col_name)

for (k in 1:length(col_name)){
  IC_up[k] <- quantile(matrix_beta[,k], 0.975)
  IC_down[k] <- quantile(matrix_beta[,k], 0.025)
  medie_coeff[k] <- mean(matrix_beta[,k])
  isSignificativo[k] <- ifelse(IC_down[k]*IC_up[k] >=0, 1, 0)
}

tabella_pred <- data.frame(medie_coeff,  IC_down, IC_up, isSignificativo)
sum(isSignificativo)
vettore_sel <- as.numeric(isSignificativo[5:length(isSignificativo)])
indici_uno <- which(vettore_sel == 1)
X_norm_significativi <- X_norm[, as.numeric(idx_usati[1:idx_aic])]
X_norm_significativi <- X_norm_significativi[,indici_uno]

final_mdl <- arima(Y_norm, xreg = X_norm_significativi, order = c(p, 0, q), include.mean = FALSE)
final_mdl$coef
coeftest(final_mdl$coef)
mse_final <- var(final_mdl$residuals)
R2_final <- 1- mse/var(Y_norm)

#in isSignificativo o tutti e soli i coefficienti singificativi

mse_list <- list()
R2_list <- list()

for (i in 1:K){
  idx_train <- idx_boot[1:(N*0.75),i];
  idx_test <- idx_boot[(N*0.75+1):N,i];
  
  mdl_train <- Arima(Y_norm[idx_train], order = c(p, 0, q), xreg = X_norm_significativi[idx_train,], include.mean = FALSE)
  mdl_test <- Arima(Y_norm[idx_test], model = mdl_train, xreg = X_norm_significativi[idx_test,], include.mean = FALSE )
  
  res <- mdl_test$residuals
  #temp <- forecast(mdl_train, h = 250, newxreg = X_norm_significativi[idx_train,])
  #temp <- as.vector(forecast(mdl_train, h = 250, xreg = X_norm_significativi[idx_test,])$mean)
  #y_hat_test <- as.vector(temp$pred)
  #y_hat_test <- as.vector(temp$mean)
  
  
  #mse <- mean((Y_norm[idx_test] - y_hat_test)^2)
  mse <- var(res)
  R2 <- 1- mse/var(Y_norm[idx_test])
  R2_list <- c(R2_list, R2)
  mse_list <- c(mse_list, mse)
  
}

mse_medio <- mean(as.numeric(mse_list))
R2_medio <- mean(as.numeric(R2_list))

#we have a model

#Support Vector Machine per time forecasting

#aggiunta
ar <- 3
Y_svm <- Y_norm[(ar+1):length(Y_norm)]
X_svm <- X_norm[(ar+1):length(Y_norm),]

for (i in 1:ar){
  X_svm <- cbind(Y_norm[i:(length(Y_norm)-ar+i-1)], X_svm)
  n = -ar+i-1
  colnames(X_svm)[1] <- paste("Y(t", n, ")", sep = "")
}

train_control = trainControl(method = "cv", number = 5)
train(x = X_svm, y = Y_svm, method = "svmRadialSigma", trControl = train_control)

#altro package
mse_list_svm <- list()
R2_list_svm <- list()

for (i in 1:K){
  idx_train <- idx_boot[1:(N*0.75),i];
  idx_test <- idx_boot[(N*0.75+1):N,i];
  
  mdl_train_svm <- svm(x = X_svm[idx_train,], y = Y_svm[idx_train], kernel = 'linear', scale = FALSE, shrinkage = FALSE, epsilon =0.5, gamma = 100)
  y_hat_svm <- predict(mdl_train_svm, X_svm[idx_test,])
  res <- Y_svm[idx_test] - y_hat_svm
 
  mse <- var(res)
  R2 <- 1- mse/var(Y_svm[idx_test])
  R2_list_svm <- c(R2_list_svm, R2)
  mse_list_svm <- c(mse_list_svm, mse)
}

mse_medio_svm <- mean(as.numeric(mse_list_svm))
R2_medio_svm <- mean(as.numeric(R2_list))
