
KNMI_list <- read.csv("Data/Geophones/Geophone data/KNMI Station list.csv", header=F) #KNMI list sent by Steve
iris_all <- read.csv("Data/Geophones/Geophone data/query.csv") #Full list from iris database
iris_all$End_date <- stringr::str_sub(iris_all$EndTime,1,10) #Converting endtime to date 
iris_all_current <- iris_all[iris_all$End_date > "2022-12-13",] #Current geophones

ind <- KNMI_list[,2] %in% unique(iris_all_current$Station) #Which KNMI station lie within iris dataset 
(KNMI_missed <- KNMI_list[!ind,]) #Finding missing stations
KNMI_missed[,2] #Missing station IDs

write.csv(iris_all_current, "Data/Geophones/Geophone data/iris_current.csv")
iris_all_current$End_date

ind_2 <- KNMI_missed[,4] %in% round(iris_all_current$Latitude,4)
KNMI_missed[!ind_2,]

KNMI_list[,4]
