varNames <- c ("Carrier",
               "Cabin",
               "Subclass",
               "Tkt_Price",
               "Change_Fee",
               "Cancel_Fee",
               "Mileage_Gain"
)
varFormats <- c ("character",
                 "character",
                 "character",
                 "character",
                 "character",
                 "character",
                 "character"
)

filename  <- paste0("data/",orig_selected,"-",dest_selected,
                    "_subclass_mapping_",year_selected,".csv")

subclass_mapping <- read.csv (filename,
                              stringsAsFactors = FALSE,
                              nrow = -1,
                              colClasses = varFormats,
                              col.names  = varNames)
setDT(subclass_mapping)
subclass_mapping[,Tkt_Price:= as.numeric(Tkt_Price)/1000]
subclass_mapping[,Change_Fee:= as.numeric(Change_Fee)/1000]
subclass_mapping[,Cancel_Fee:= as.numeric(Cancel_Fee)/1000]
subclass_mapping[,Mileage_Gain:= as.numeric(Mileage_Gain)/1000]
if(enable_Y_cabin_only){
  subclass_mapping = subclass_mapping[Cabin=="Y",]
}
