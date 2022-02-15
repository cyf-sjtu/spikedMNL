setDT(data.training)
setDT(data.testing)

if(enable_my_implementation){
  mnl.beta = as.numeric(mnl.model$estimate)
}else{
  mnl.beta = as.numeric(mnl.model$coefficients)
}
# mnl.x.names = names(mnl.model$coefficients)
mnl.x = data.training[, ..mnl.x.names]
data.training[, item_util := as.matrix(mnl.x) %*% mnl.beta]
data.training$item_util = data.training$item_util +
  runif(length(data.training$item_util))*1e-6 # break ties
data.training[, item_attr := exp(item_util)]
data.training[, pnr_attr := sum(item_attr), by=c("PNR") ]
data.training[, choice_prob := item_attr/pnr_attr]
data.training[, (mnl.model.name) := choice_prob]
data.training[, num_item := max(Itin_No), by=c("PNR") ]
data.training[, squared_error := 1/num_item*(Booked_Itin-choice_prob)^2]
data.training[, pick_1 := as.numeric(item_util == max(item_util)),by=c("PNR")];
data.training[, pick_3 := as.numeric(item_util %in% tail(sort(item_util),3)),by=c("PNR")];
data.training[, pick_5 := as.numeric(item_util %in% tail(sort(item_util),5)),by=c("PNR")];
# data.training[, pick_10 := as.numeric(item_util %in% tail(sort(item_util),10)),by=c("PNR")];

num_bookings = sum(data.training$Booked_Itin)
mnl.model.training.rmse = sqrt(sum(data.training$squared_error)/num_bookings)
mnl.model.training.hit_1 = sum(data.training$pick_1*data.training$Booked_Itin)/num_bookings
mnl.model.training.hit_3 = sum(data.training$pick_3*data.training$Booked_Itin)/num_bookings
mnl.model.training.hit_5 = sum(data.training$pick_5*data.training$Booked_Itin)/num_bookings
# mnl.model.training.hit_10 = sum(data.training$pick_10*data.training$Booked_Itin)/num_bookings

actual.subclass.bookings = 
  setDT(count(data.training, vars=c("Carrier","Subclass"),wt_var="Booked_Itin"))
names(actual.subclass.bookings)[names(actual.subclass.bookings) == "freq"] = "Actual"
predicted.subclass.bookings =
  setDT(count(data.training, vars=c("Carrier","Subclass"),wt_var="choice_prob"))
names(predicted.subclass.bookings)[names(predicted.subclass.bookings) == "freq"] = "Predicted"
mnl.model.training.bookings = merge(actual.subclass.bookings,
                                    predicted.subclass.bookings,
                                    by=c("Carrier","Subclass"))
setDT(mnl.model.training.bookings)[,Error:=(Predicted-Actual)/Actual]
mnl.model.training.bookings[,Error:=ifelse(is.infinite(Error),NA,Error)]

# show("In-sample accuracy: ")
# show(paste0("RMSE: ",mnl.model.training.rmse)) # hit-1
# show(paste0("Hit_1: ",mnl.model.training.hit_1)) # hit-1
# show(paste0("Hit_3: ",mnl.model.training.hit_3)) # hit-3
# show(paste0("Hit_5: ",mnl.model.training.hit_5)) # hit-5
# show(paste0("Hit_10: ",mnl.model.training.hit_10)) # hit-10

mnl.x = data.testing[, ..mnl.x.names]
data.testing[, item_util := as.matrix(mnl.x) %*% mnl.beta]
data.testing$item_util = data.testing$item_util +
  runif(length(data.testing$item_util))*1e-6 # break ties
data.testing[, item_attr := exp(item_util)]
data.testing[, pnr_attr := sum(item_attr), by=c("PNR") ]
data.testing[, choice_prob := item_attr/pnr_attr]
data.testing[, (mnl.model.name) := choice_prob]
data.testing[, num_item := max(Itin_No), by=c("PNR") ]
data.testing[, squared_error := 1/num_item*(Booked_Itin-choice_prob)^2]
data.testing[, pick_1 := as.numeric(item_util == max(item_util)),by=c("PNR")];
data.testing[, pick_3 := as.numeric(item_util %in% tail(sort(item_util),3)),by=c("PNR")];
data.testing[, pick_5 := as.numeric(item_util %in% tail(sort(item_util),5)),by=c("PNR")];
data.testing[, pick_10 := as.numeric(item_util %in% tail(sort(item_util),10)),by=c("PNR")];

num_bookings = sum(data.testing$Booked_Itin)
mnl.model.testing.loglik = sum(log(data.testing$choice_prob)*data.testing$Booked_Itin*data.testing$pnr_weight)
mnl.model.testing.rmse = sqrt(sum(data.testing$squared_error)/num_bookings)
mnl.model.testing.hit_1 = sum(data.testing$pick_1*data.testing$Booked_Itin)/num_bookings
mnl.model.testing.hit_3 = sum(data.testing$pick_3*data.testing$Booked_Itin)/num_bookings
mnl.model.testing.hit_5 = sum(data.testing$pick_5*data.testing$Booked_Itin)/num_bookings
# mnl.model.testing.hit_10 = sum(data.testing$pick_10*data.testing$Booked_Itin)/num_bookings

actual.subclass.bookings = 
  setDT(count(data.testing, vars=c("Carrier","Subclass"),wt_var="Booked_Itin"))
names(actual.subclass.bookings)[names(actual.subclass.bookings) == "freq"] = "Actual"
predicted.subclass.bookings =
  setDT(count(data.testing, vars=c("Carrier","Subclass"),wt_var="choice_prob"))
names(predicted.subclass.bookings)[names(predicted.subclass.bookings) == "freq"] = "Predicted"
mnl.model.testing.bookings = merge(actual.subclass.bookings,
                                    predicted.subclass.bookings,
                                    by=c("Carrier","Subclass"))
setDT(mnl.model.testing.bookings)[,Error:=(Predicted-Actual)/Actual]
mnl.model.testing.bookings[,Error:=ifelse(is.infinite(Error),NA,Error)]

show("Out-of-sample accuracy: ")
show(paste0("RMSE: ",mnl.model.testing.rmse)) # hit-1
show(paste0("Hit_1: ",mnl.model.testing.hit_1)) # hit-1
show(paste0("Hit_3: ",mnl.model.testing.hit_3)) # hit-3
show(paste0("Hit_5: ",mnl.model.testing.hit_5)) # hit-5
# show(paste0("Hit_10: ",mnl.model.testing.hit_10)) # hit-10

if(enable_result_output){
  # write("In-sample accuracy:", 
  #       file = result.file.name, append = TRUE)
  # write(paste0("RMSE: ",mnl.model.training.rmse),
  #       file = result.file.name, append = TRUE)
  # write(paste0("Hit_1: ",mnl.model.training.hit_1),
  #       file = result.file.name, append = TRUE)
  # write(paste0("Hit_3: ",mnl.model.training.hit_3), 
  #       file = result.file.name, append = TRUE)
  # write(paste0("Hit_5: ",mnl.model.training.hit_5), 
  #       file = result.file.name, append = TRUE)
  write("Out-of-sample accuracy: ", 
        file = result.file.name, append = TRUE)
  write(paste0("loglik: ",mnl.model.testing.loglik),
        file = result.file.name, append = TRUE)
  write(paste0("RMSE: ",mnl.model.testing.rmse), 
        file = result.file.name, append = TRUE)
  write(paste0("Hit_1: ",mnl.model.testing.hit_1), 
        file = result.file.name, append = TRUE)
  write(paste0("Hit_3: ",mnl.model.testing.hit_3), 
        file = result.file.name, append = TRUE)
  write(paste0("Hit_5: ",mnl.model.testing.hit_5), 
        file = result.file.name, append = TRUE)
  write("********************************************************************",
        file=result.file.name,append = TRUE)
}

if(enable_xlsx_output){
  if(!is.null(mnl.model)){
    temp = data.frame(test.loglik = mnl.model.testing.loglik,
                      test.rmse = mnl.model.testing.rmse,
                      test.hit1 = mnl.model.testing.hit_1,
                      test.hit3 = mnl.model.testing.hit_3,
                      test.hit5 = mnl.model.testing.hit_5
                      )
    addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-accuracy"))
    writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-accuracy"),
                   x = temp, rowNames = FALSE)
    
    # addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-subclass-training"))
    # writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-subclass-training"),
    #                x = as.data.frame(mnl.model.training.bookings),
    #                rowNames = FALSE)
    # addWorksheet(wb = workbook, sheetName = paste0(mnl.model.name,"-subclass-testing"))
    # writeDataTable(wb = workbook, sheet = paste0(mnl.model.name,"-subclass-testing"),
    #                x = as.data.frame(mnl.model.testing.bookings),
    #                rowNames = FALSE)
    
    saveWorkbook(workbook, file = result.xlsx.name, overwrite = TRUE)
  }
}
