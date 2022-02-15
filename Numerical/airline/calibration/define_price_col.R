# features and attributes
num_book_DtD = 3 # 1-[0,6], 2-[7,14], 3-[14,inf]
num_book_hour = 2 # 1-[9,18], 2-other
num_book_weekday = 2  # 1-weekday, 2-weekend
num_channel = 3 # 1-other, 2-direct, 3-OTA
num_type_dept_hour = 2 # 1-morning, 0-afternoon

col_price_combo = c()
for (i in 1:num_book_DtD) {
  for (j in 1:num_book_hour) {
    for (k in 1:num_book_weekday) {
      for (l in 1:num_channel) {
        i = ifelse(add_DtD_to_price_combo,i,0)
        j = ifelse(add_bookhour_to_price_combo,j,0)
        k = ifelse(add_bookday_to_price_combo,k,0)
        l = ifelse(add_channel_to_price_combo,l,0)
        col_price_combo = c(col_price_combo, 
                            paste0("price_combo_",i,j,k,l))
      }
    }
  }
}
col_price_combo = unique(col_price_combo)

# Refundable Combo 
col_refundable_combo = c()
for (carrier in c("Host",OA_set)){
  for(l in 1:num_channel){
    carrier = ifelse(add_carrier_to_refundable_combo,carrier,"ALL")
    l = ifelse(add_channel_to_refundable_combo,l,0)
    col_refundable_combo = c(col_refundable_combo, 
                             paste0("refundable_",
                                    carrier,"_",l))
  }
}
col_refundable_combo = unique(col_refundable_combo)

# Cheapest Fare Spike Combo
col_spike_combo = c()
for (i in 1:dim(subclass_mapping)[1]){
  if(subclass_mapping[i]$Cabin != "Y"){next}
  carrier = ifelse(add_carrier_to_spike_combo,
                   subclass_mapping[i]$Carrier,"ALL")
  subclass = ifelse(add_subclass_to_spike_combo,
                    subclass_mapping[i]$Subclass,"ALL")
  for (l in 1:num_channel){
    l = ifelse(add_channel_to_spike_combo,l,0)
    col_spike_combo = c(col_spike_combo, 
                        paste0("spike_",carrier,"_",subclass,"_",l))
  }
}
col_spike_combo = unique(col_spike_combo)

# subinterval columns
# col_subinterval = paste0("subinterval_",0:(num_subinterval-1))
