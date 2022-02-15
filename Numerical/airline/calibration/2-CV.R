library(data.table)
library(plyr)
library(maxLik)
library(openxlsx)
library(wsrf)

rm(list = ls())
gc()

skip_header = 0

Sys.setlocale("LC_TIME", "C")
setwd("C:/Code/MUProject/Estimation")

# market
orig_selected = "ORG"
dest_selected = "DST"
year_selected = 2019
data_selected = "MON"

for(month_iter in seq(1,12)){ 
  print(month_iter)

if(!skip_header){
  month_selected = month_iter
  if(month_selected==12){
    dept_date_start = as.Date(paste0(year_selected,"-",month_selected,"-1")) 
    dept_date_end = as.Date(paste0(year_selected,"-",month_selected,"-31"))
  }else{
    dept_date_start = as.Date(paste0(year_selected,"-",month_selected,"-1")) 
    dept_date_end = as.Date(paste0(year_selected,"-",month_selected+1,"-1"))
  }
  dept_date_set = seq.Date(from = dept_date_start, 
                           to = dept_date_end,
                           by = "day")
  # should exclude holidays
  holiday_set = seq.Date(from = as.Date(paste0(year_selected,"-10-1")), 
                         to = as.Date(paste0(year_selected,"-10-7")), 
                         by = "day")
  # specify the holiday as announced by the government
  if(year_selected %in% c(2017)){
    holiday_set = c(holiday_set,as.Date(paste0(year_selected,"-10-08")))
  }
  # take set difference
  dept_date_set = as.Date(setdiff(as.character(dept_date_set),
                                  as.character(holiday_set)))
  dept_date_set = dept_date_set[weekdays(dept_date_set)=="Monday"]
  
  # OA
  OA_set = c("OA1","OA2")
  OA1 = "OA1"
  OA1_set = c("OA1")
  OA2 = "OA2" 
  OA2_set = c("OA2")
  
  # model specification
  enable_Y_cabin_only = 1
  enable_OA = 0
  enable_price_combo = 1 # e.g., DtD x channel 
  add_DtD_to_price_combo = 1
  add_bookhour_to_price_combo = 0
  add_bookday_to_price_combo = 0
  add_channel_to_price_combo = 1
  enable_refundable_combo = 0 # e.g., carrier x channel
  add_carrier_to_refundable_combo = 0
  add_channel_to_refundable_combo = 0
  enable_spike_combo = 0 # carrier x subclass x channel
  add_carrier_to_spike_combo = 0
  add_subclass_to_spike_combo = 0
  add_channel_to_spike_combo = 0
  
  filename = paste0("data/",orig_selected,"-",dest_selected,
                    "_Choice_Set_",year_selected,"_",
                    data_selected,"_",month_selected,".csv")
  choice_set = read.csv(filename, stringsAsFactors = FALSE)
  source("load_subclass_mapping.R")
  source("define_price_col.R")
}

# function specification
enable_my_implementation = 1
enable_result_output = 1
enable_xlsx_output = 1

postfix = "test"

#############################################################################

result.file.name = paste0("output/",orig_selected,"-",dest_selected,
                          "_Result_",year_selected,"_",
                          data_selected,"_",month_selected,".txt")
result.xlsx.name = paste0("output/",orig_selected,"-",dest_selected,
                          "_Result_",year_selected,"_",
                          data_selected,"_",month_selected,".xlsx")

setDT(choice_set)
choice_set = choice_set[weekdays(as.Date(Dept_Dt))=="Monday",]

if(enable_Y_cabin_only){
  choice_set = choice_set[which(Cabin=="Y"),]
}

if(!enable_OA){
  choice_set = choice_set[which(Carrier=="Host"),]
}


set.seed(19260817)
# split choice set to get training and testing data
pnr_list = unique(choice_set$PNR)
temp = sample(pnr_list, round(length(pnr_list)/2))
data.training = choice_set[ -which(PNR %in% temp),]
data.testing = choice_set[ which(PNR %in% temp),]

# clear output file and write data info
if(enable_result_output){
  write(paste0("****************** ",as.character(Sys.time()),
               " ******************"),file = result.file.name)
  write(paste(c("Market:",orig_selected,"-",dest_selected,
                year_selected,data_selected),collapse = " "),
        file = result.file.name, append = TRUE)
  write(paste0(c("Duration:", as.character(dept_date_start),
                 "-",as.character(dept_date_end)),collapse=" "),
        file = result.file.name, append = TRUE)
  if(enable_OA){
    write("With OA",file = result.file.name, append = TRUE)
  }else{
    write("Without OA",file = result.file.name, append = TRUE)
  }
  write(paste0("Training data: ", length(pnr_list)-length(temp), " PNRs, ",
               dim(data.training)[1]," rows"), 
        file = result.file.name, append = TRUE)
  write(paste0("Testing data: ", length(temp), " PNRs, ",
               dim(data.testing)[1]," rows"), 
        file = result.file.name, append = TRUE)
}

# create workbook
if(enable_xlsx_output){
  workbook <- createWorkbook()
}

###########################################################################

show("Calibrating choice model...")

# use complete choice set
# data.training = choice_set

if(month_selected==12){
  filename = paste0("data/",orig_selected,"-",dest_selected,
                    "_Choice_Set_",year_selected+1,"_",
                    data_selected,"_",1,".csv")
}else{
  filename = paste0("data/",orig_selected,"-",dest_selected,
                    "_Choice_Set_",year_selected,"_",
                    data_selected,"_",month_selected+1,".csv")
}

data.testing = read.csv(filename, stringsAsFactors = FALSE)
setDT(data.testing)[,Flt_No := as.character(Flt_No)]

# MNL model
mnl.model.name = "MNL"
if(enable_price_combo){
  mnl.x.names = c(col_price_combo,"Mileage_Gain")
}else{
  mnl.x.names = c("Tkt_Price","Mileage_Gain")
}
if(enable_OA){
  mnl.x.names = c(mnl.x.names,OA_set)
  if(enable_refundable_combo){
    mnl.x.names = c(mnl.x.names,col_refundable_combo)
  }else{
    mnl.x.names = c(mnl.x.names,"is_refundable")
  }
}else{
  mnl.x.names = c(mnl.x.names,"is_refundable")
}
# mnl.x.names = c("Tkt_Price")
source("calibrate_my_MNL_model.R")
if(!is.null(mnl.model)){ source("calculate_MNL_accuracy.R") }
# SMNL model
mnl.model.name = "SMNL"
if(enable_price_combo){
  mnl.x.names = c(col_price_combo,"Mileage_Gain")
}else{
  mnl.x.names = c("Tkt_Price","Mileage_Gain")
}
if(enable_OA){
  mnl.x.names = c(mnl.x.names,OA_set)
  if(enable_refundable_combo){
    mnl.x.names = c(mnl.x.names,col_refundable_combo)
  }else{
    mnl.x.names = c(mnl.x.names,"is_refundable")
  }
  if(enable_spike_combo){
    mnl.x.names = c(mnl.x.names,col_spike_combo)
  }else{
    mnl.x.names = c(mnl.x.names,"is_cheapest")
  }
}else{
  mnl.x.names = c(mnl.x.names,"is_refundable","is_cheapest")
}
source("calibrate_my_MNL_model.R")
if(!is.null(mnl.model)){ source("calculate_MNL_accuracy.R") }

###########################################################################################
data.training[,Carrier := as.factor(Carrier)]
data.testing[,Carrier := as.factor(Carrier)]
data.training[,num_item := .N, by=c("PNR")]
data.testing[,num_item := .N, by=c("PNR")]
data.training[,avg_price := mean(Tkt_Price), by=c("PNR")]
data.testing[,avg_price := mean(Tkt_Price), by=c("PNR")]

rf.model.name = "RF"
rf.x.names = c("Tkt_Price","Change_Fee","Mileage_Gain",
               "book_DtD","book_hour","book_weekday_idx","book_channel_idx",
               "is_refundable")
rf.formula = as.formula(paste0("Booked_Itin~",paste(rf.x.names,collapse = "+")))
rf.x.names = c("Booked_Itin",rf.x.names)
ptm <- proc.time()
rf.model <- wsrf(rf.formula,data = data.training[,..rf.x.names],ntree=200)
elapsed_time = proc.time() - ptm
show(elapsed_time[3])
rf.model.predict <- predict(rf.model,data.testing[,..rf.x.names],type="waprob");
rf.model.prob <- rf.model.predict$waprob[,2];
data.testing[,RF := rf.model.prob]
data.testing[,RF_sum := sum(RF), by=c("PNR")]
data.testing[,RF := RF/RF_sum]
rf.model.testing.loglik = sum(log(data.testing$RF)*data.testing$Booked_Itin*data.testing$pnr_weight)
show(rf.model.testing.loglik)

if(enable_result_output){
  write(paste0("RF-loglik: ",rf.model.testing.loglik),
        file = result.file.name, append = TRUE)
  write(paste0("Time: ",elapsed_time[3]),
        file = result.file.name, append = TRUE)
}

########################################################################################

data.testing[, flight_is_booked := sum(Booked_Itin), by=c("PNR","Carrier","Flt_No")]
temp = data.testing[Carrier == "Host" & flight_is_booked,]
temp[, cheapest_M := sum((Subclass == "M")*is_cheapest), by=c("PNR")]
temp[, cheapest_N := sum((Subclass == "N")*is_cheapest), by=c("PNR")]
temp[, cheapest_R := sum((Subclass == "R")*is_cheapest), by=c("PNR")]
#temp[, cheapest_H := sum((Subclass == "H")*is_cheapest), by=c("PNR")]
temp[, cheapest_K := sum((Subclass == "K")*is_cheapest), by=c("PNR")]
temp[, cheapest_L := sum((Subclass == "L")*is_cheapest), by=c("PNR")]

enable_by_date = 0
if(enable_by_date){
# open 8 classes
data2 = temp[cheapest_M==1,]
M_base = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "Booked_Itin")
setDT(M_base)[,base := freq/sum(freq),by=c("Dept_Dt")]
setorder(M_base,Dept_Dt,-Tkt_Price)
M_mnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "MNL")
setDT(M_mnl)[,mnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(M_mnl,Dept_Dt,-Tkt_Price)
M_smnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "SMNL")
setDT(M_smnl)[,smnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(M_smnl,Dept_Dt,-Tkt_Price)

M_rf = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "RF")
setDT(M_rf)[,rf := freq/sum(freq),by=c("Dept_Dt")]
setorder(M_rf,Dept_Dt,-Tkt_Price)

data2 = temp[cheapest_N==1,]
N_base = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "Booked_Itin")
setDT(N_base)[,base := freq/sum(freq),by=c("Dept_Dt")]
setorder(N_base,Dept_Dt,-Tkt_Price)
N_mnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "MNL")
setDT(N_mnl)[,mnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(N_mnl,Dept_Dt,-Tkt_Price)
N_smnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "SMNL")
setDT(N_smnl)[,smnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(N_smnl,Dept_Dt,-Tkt_Price)
N_rf = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "RF")
setDT(N_rf)[,rf := freq/sum(freq),by=c("Dept_Dt")]
setorder(N_rf,Dept_Dt,-Tkt_Price)

data2 = temp[cheapest_R==1,]
R_base = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "Booked_Itin")
setDT(R_base)[,base := freq/sum(freq),by=c("Dept_Dt")]
setorder(R_base,Dept_Dt,-Tkt_Price)
R_mnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "MNL")
setDT(R_mnl)[,mnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(R_mnl,Dept_Dt,-Tkt_Price)
R_smnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "SMNL")
setDT(R_smnl)[,smnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(R_smnl,Dept_Dt,-Tkt_Price)
R_rf = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "RF")
setDT(R_rf)[,rf := freq/sum(freq),by=c("Dept_Dt")]
setorder(R_rf,Dept_Dt,-Tkt_Price)

data2 = temp[cheapest_K==1,]
K_base = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "Booked_Itin")
setDT(K_base)[,base := freq/sum(freq),by=c("Dept_Dt")]
setorder(K_base,Dept_Dt,-Tkt_Price)
K_mnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "MNL")
setDT(K_mnl)[,mnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(K_mnl,Dept_Dt,-Tkt_Price)
K_smnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "SMNL")
setDT(K_smnl)[,smnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(K_smnl,Dept_Dt,-Tkt_Price)
K_rf = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "RF")
setDT(K_rf)[,rf := freq/sum(freq),by=c("Dept_Dt")]
setorder(K_rf,Dept_Dt,-Tkt_Price)

data2 = temp[cheapest_L==1,]
L_base = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "Booked_Itin")
setDT(L_base)[,base := freq/sum(freq),by=c("Dept_Dt")]
setorder(L_base,Dept_Dt,-Tkt_Price)
L_mnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "MNL")
setDT(L_mnl)[,mnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(L_mnl,Dept_Dt,-Tkt_Price)
L_smnl = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "SMNL")
setDT(L_smnl)[,smnl := freq/sum(freq),by=c("Dept_Dt")]
setorder(L_smnl,Dept_Dt,-Tkt_Price)
L_rf = count(data2,vars= c("Tkt_Price","Dept_Dt"),wt_var = "RF")
setDT(L_rf)[,rf := freq/sum(freq),by=c("Dept_Dt")]
setorder(L_rf,Dept_Dt,-Tkt_Price)

M_stat = merge(M_base,M_mnl,by=c("Tkt_Price","Dept_Dt"))
M_stat = merge(M_stat,M_smnl,by=c("Tkt_Price","Dept_Dt"))
M_stat = merge(M_stat,M_rf,by=c("Tkt_Price","Dept_Dt"))
M_stat = M_stat[,c("Tkt_Price","Dept_Dt","base","mnl","smnl","rf")]
setorder(M_stat,Dept_Dt,-Tkt_Price)

N_stat = merge(N_base,N_mnl,by=c("Tkt_Price","Dept_Dt"))
N_stat = merge(N_stat,N_smnl,by=c("Tkt_Price","Dept_Dt"))
N_stat = merge(N_stat,N_rf,by=c("Tkt_Price","Dept_Dt"))
N_stat = N_stat[,c("Tkt_Price","Dept_Dt","base","mnl","smnl","rf")]
setorder(N_stat,Dept_Dt,-Tkt_Price)

R_stat = merge(R_base,R_mnl,by=c("Tkt_Price","Dept_Dt"))
R_stat = merge(R_stat,R_smnl,by=c("Tkt_Price","Dept_Dt"))
R_stat = merge(R_stat,R_rf,by=c("Tkt_Price","Dept_Dt"))
R_stat = R_stat[,c("Tkt_Price","Dept_Dt","base","mnl","smnl","rf")]
setorder(R_stat,Dept_Dt,-Tkt_Price)

K_stat = merge(K_base,K_mnl,by=c("Tkt_Price","Dept_Dt"))
K_stat = merge(K_stat,K_smnl,by=c("Tkt_Price","Dept_Dt"))
K_stat = merge(K_stat,K_rf,by=c("Tkt_Price","Dept_Dt"))
K_stat = K_stat[,c("Tkt_Price","Dept_Dt","base","mnl","smnl","rf")]
setorder(K_stat,Dept_Dt,-Tkt_Price)

L_stat = merge(L_base,L_mnl,by=c("Tkt_Price","Dept_Dt"))
L_stat = merge(L_stat,L_smnl,by=c("Tkt_Price","Dept_Dt"))
L_stat = merge(L_stat,L_rf,by=c("Tkt_Price","Dept_Dt"))
L_stat = L_stat[,c("Tkt_Price","Dept_Dt","base","mnl","smnl")]
setorder(L_stat,Dept_Dt,-Tkt_Price)

}else{
  # open 8 classes
  data2 = temp[cheapest_M==1,]
  M_base = count(data2,vars= c("Tkt_Price"),wt_var = "Booked_Itin")
  setDT(M_base)[,base := freq/sum(freq)]
  setorder(M_base,-Tkt_Price)
  M_mnl = count(data2,vars= c("Tkt_Price"),wt_var = "MNL")
  setDT(M_mnl)[,mnl := freq/sum(freq)]
  setorder(M_mnl,-Tkt_Price)
  M_smnl = count(data2,vars= c("Tkt_Price"),wt_var = "SMNL")
  setDT(M_smnl)[,smnl := freq/sum(freq)]
  setorder(M_smnl,-Tkt_Price)
  M_rf = count(data2,vars= c("Tkt_Price"),wt_var = "RF")
  setDT(M_rf)[,rf := freq/sum(freq)]
  setorder(M_rf,-Tkt_Price)

  data2 = temp[cheapest_N==1,]
  N_base = count(data2,vars= c("Tkt_Price"),wt_var = "Booked_Itin")
  setDT(N_base)[,base := freq/sum(freq)]
  setorder(N_base,-Tkt_Price)
  N_mnl = count(data2,vars= c("Tkt_Price"),wt_var = "MNL")
  setDT(N_mnl)[,mnl := freq/sum(freq)]
  setorder(N_mnl,-Tkt_Price)
  N_smnl = count(data2,vars= c("Tkt_Price"),wt_var = "SMNL")
  setDT(N_smnl)[,smnl := freq/sum(freq)]
  setorder(N_smnl,-Tkt_Price)
  N_rf = count(data2,vars= c("Tkt_Price"),wt_var = "RF")
  setDT(N_rf)[,rf := freq/sum(freq)]
  setorder(N_rf,-Tkt_Price)

  data2 = temp[cheapest_R==1,]
  R_base = count(data2,vars= c("Tkt_Price"),wt_var = "Booked_Itin")
  setDT(R_base)[,base := freq/sum(freq)]
  setorder(R_base,-Tkt_Price)
  R_mnl = count(data2,vars= c("Tkt_Price"),wt_var = "MNL")
  setDT(R_mnl)[,mnl := freq/sum(freq)]
  setorder(R_mnl,-Tkt_Price)
  R_smnl = count(data2,vars= c("Tkt_Price"),wt_var = "SMNL")
  setDT(R_smnl)[,smnl := freq/sum(freq)]
  setorder(R_smnl,-Tkt_Price)
  R_rf = count(data2,vars= c("Tkt_Price"),wt_var = "RF")
  setDT(R_rf)[,rf := freq/sum(freq)]
  setorder(R_rf,-Tkt_Price)

  data2 = temp[cheapest_K==1,]
  K_base = count(data2,vars= c("Tkt_Price"),wt_var = "Booked_Itin")
  setDT(K_base)[,base := freq/sum(freq)]
  setorder(K_base,-Tkt_Price)
  K_mnl = count(data2,vars= c("Tkt_Price"),wt_var = "MNL")
  setDT(K_mnl)[,mnl := freq/sum(freq)]
  setorder(K_mnl,-Tkt_Price)
  K_smnl = count(data2,vars= c("Tkt_Price"),wt_var = "SMNL")
  setDT(K_smnl)[,smnl := freq/sum(freq)]
  setorder(K_smnl,-Tkt_Price)
  K_rf = count(data2,vars= c("Tkt_Price"),wt_var = "RF")
  setDT(K_rf)[,rf := freq/sum(freq)]
  setorder(K_rf,-Tkt_Price)

  data2 = temp[cheapest_L==1,]
  L_base = count(data2,vars= c("Tkt_Price"),wt_var = "Booked_Itin")
  setDT(L_base)[,base := freq/sum(freq)]
  setorder(L_base,-Tkt_Price)
  L_mnl = count(data2,vars= c("Tkt_Price"),wt_var = "MNL")
  setDT(L_mnl)[,mnl := freq/sum(freq)]
  setorder(L_mnl,-Tkt_Price)
  L_smnl = count(data2,vars= c("Tkt_Price"),wt_var = "SMNL")
  setDT(L_smnl)[,smnl := freq/sum(freq)]
  setorder(L_smnl,-Tkt_Price)
  L_rf = count(data2,vars= c("Tkt_Price"),wt_var = "RF")
  setDT(L_rf)[,rf := freq/sum(freq)]
  setorder(L_rf,-Tkt_Price)

  M_stat = merge(M_base,M_mnl,by=c("Tkt_Price"))
  M_stat = merge(M_stat,M_smnl,by=c("Tkt_Price"))
  M_stat = merge(M_stat,M_rf,by=c("Tkt_Price"))
  M_stat = M_stat[,c("Tkt_Price","base","mnl","smnl","rf")]
  setorder(M_stat,-Tkt_Price)

  N_stat = merge(N_base,N_mnl,by=c("Tkt_Price"))
  N_stat = merge(N_stat,N_smnl,by=c("Tkt_Price"))
  N_stat = merge(N_stat,N_rf,by=c("Tkt_Price"))
  N_stat = N_stat[,c("Tkt_Price","base","mnl","smnl","rf")]
  setorder(N_stat,-Tkt_Price)

  R_stat = merge(R_base,R_mnl,by=c("Tkt_Price"))
  R_stat = merge(R_stat,R_smnl,by=c("Tkt_Price"))
  R_stat = merge(R_stat,R_rf,by=c("Tkt_Price"))
  R_stat = R_stat[,c("Tkt_Price","base","mnl","smnl","rf")]
  setorder(R_stat,-Tkt_Price)

  K_stat = merge(K_base,K_mnl,by=c("Tkt_Price"))
  K_stat = merge(K_stat,K_smnl,by=c("Tkt_Price"))
  K_stat = merge(K_stat,K_rf,by=c("Tkt_Price"))
  K_stat = K_stat[,c("Tkt_Price","base","mnl","smnl","rf")]
  setorder(K_stat,-Tkt_Price)

  L_stat = merge(L_base,L_mnl,by=c("Tkt_Price"))
  L_stat = merge(L_stat,L_smnl,by=c("Tkt_Price"))
  L_stat = merge(L_stat,L_rf,by=c("Tkt_Price"))
  L_stat = L_stat[,c("Tkt_Price","base","mnl","smnl","rf")]
  setorder(L_stat,-Tkt_Price)
}


# View(R_stat)
# View(N_stat)
# View(M_stat)
# View(K_stat)
# View(L_stat)



if(enable_xlsx_output){
  stat.xlsx.name = paste0("output/",orig_selected,"-",dest_selected,
                          "_Stats_",year_selected,"_",month_selected,".xlsx")
  workbook2 <- createWorkbook()
  addWorksheet(wb = workbook2, sheetName = "L-stat")
  writeDataTable(wb = workbook2, sheet = "L-stat", x = L_stat, rowNames = FALSE)
  addWorksheet(wb = workbook2, sheetName = "N-stat")
  writeDataTable(wb = workbook2, sheet = "N-stat", x = N_stat, rowNames = FALSE)
  addWorksheet(wb = workbook2, sheetName = "R-stat")
  writeDataTable(wb = workbook2, sheet = "R-stat", x = R_stat, rowNames = FALSE)
  addWorksheet(wb = workbook2, sheetName = "M-stat")
  writeDataTable(wb = workbook2, sheet = "M-stat", x = M_stat, rowNames = FALSE)
  addWorksheet(wb = workbook2, sheetName = "K-stat")
  writeDataTable(wb = workbook2, sheet = "K-stat", x = K_stat, rowNames = FALSE)
  saveWorkbook(workbook2, file = stat.xlsx.name, overwrite = TRUE)
}

}