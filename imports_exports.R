library(gdata)
library(abind)

MRIO_import_data = read.xls('~/GitHub/MRIO_BIO_SATELLITE/global_consumption_international_assessment_scale.xls')
MRIO_export_data = read.xls('~/GitHub/MRIO_BIO_SATELLITE/global_production_international_assessment_scale.xls')

sorted_export_inds = match(MRIO_import_data$Consumption_Country, MRIO_export_data$Consumption_Country)
sorted_export_threats = MRIO_export_data$Aggregated_Threats[sorted_export_inds]

sorted_import_inds = match(MRIO_export_data$Consumption_Country, MRIO_import_data$Consumption_Country)
sorted_import_threats = MRIO_import_data$Aggregated_Threats[sorted_import_inds]

concatenated_data = rbind(abind(MRIO_import_data$Aggregated_Threats, rev(sorted_import_threats)),
                          abind(-sorted_export_threats, -rev(MRIO_export_data$Aggregated_Threats)))
                          
concatenated_names = abind(MRIO_import_data$Consumption_Country, rev(MRIO_export_data$Consumption_Country))
  
last_country = 20
plot_vector = c(1:last_country, (dim(concatenated_data)[2] - last_country):dim(concatenated_data)[2])

pdf('~/GitHub/MRIO_BIO_SATELLITE/imports_exports.pdf', width = 8.3, height = 11.7)
par(mar = c(10,10,10,10)) 
barplot(concatenated_data[1, plot_vector], xlim = c(-1500, 3000), names = concatenated_names[plot_vector], col = "red", las = 2, 
        horiz = TRUE)

barplot(concatenated_data[2, plot_vector], xlim = c(-1500, 3000), names = concatenated_names[plot_vector], col = "black", las = 2, 
        horiz = TRUE, add = TRUE)


graphics.off()