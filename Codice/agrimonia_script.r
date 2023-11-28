install.packages('tseries')
library(tseries)


load("Agrimonia_Dataset_v_3_0_0.Rdata")
subset <- AgrImOnIA_Dataset_v_3_0_0[AgrImOnIA_Dataset_v_3_0_0$IDStations=="627",]

#AQ_nh3 variabile risposta, si può usare l'altra come regressore
AQ_nh3 <- subset$AQ_nh3

#numero di valori nan al''interno della response variable AQ_NH3
apply ( X = is.na(subset), MARGIN = 2, FUN = sum)

#kalman smoother per missing values

#autocorrelation function -> per ordine dell'MA
acf(AQ_nh3)

#partial autocorrelation function -> perordine dell'AR
pacf(AQ_nh3)

#test di stazionarietà della response variable
adf.test(AQ_nh3)

#