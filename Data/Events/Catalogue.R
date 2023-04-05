cat <-  read.csv(file= "C:/Users/murphyc4/OneDrive - Lancaster University/STOR-i/PhD/Projects/Spatio-temporal/Data/Events/2022-04-12_15-09-25_cat.csv")
cat$date <- ymd(cat$date)
plot(cat$mag)
h<-c(rep(1.05,750), rep(1.65, 750))
lines(h,col="red", lwd=2)

Cat.df <- data.frame(Date=cat$date, M=cat$mag, X=cat$RD_X, Y=cat$RD_Y, Year=cat$Year)
Cat.df <- Cat.df[Cat.df$Date >= as.Date("1995-01-01"),]
save(Cat.df, file="STOR-i/PhD/Projects/Spatio-temporal/Data/Catalogue.Rda")
rm(Cat.df)

