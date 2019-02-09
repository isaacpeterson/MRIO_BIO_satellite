library(gdata)
library(abind)
library(rworldmap)
name_data = read.xls('~/GitHub/MRIO_BIO_SATELLITE/name_data.xls')

country_scale_data = read.xls('~/GitHub/MRIO_BIO_SATELLITE/MRIO_output_tables/country_scale_net_table.xls')
sorted_characteristics = sort(country_scale_data$net_consumption, decreasing = TRUE, index.return = TRUE)
country_scale_data = country_scale_data[sorted_characteristics$ix, ]

global_names_to_match = as.character(country_scale_data$Country)

country_num = 50

plot_names = global_names_to_match[country_num:1]

name_mapper = matrix( ncol=2, byrow=TRUE, 
                      c("USA"  , 'United States',               
                        "Taiwan"  ,   'Taiwan, Province of China',      
                        "Tanzania"  ,       'Tanzania, United Republic of'  ,  
                        "Iran"       ,   'Iran, Islamic Republic of',         
                        "Russia"    ,      'Russian Federation',       
                        "Venezuela"  ,           'Venezuela, Bolivarian Republic of',
                        "UK" , 'United Kingdom',
                        "UAE", "United Arab Emirates",
                        "Laos", "Lao People's Democratic Republic",
                        "DR Congo", "Congo, The Democratic Republic of the" ,          
                        "South Korea", "Korea, Republic of",
                        "North Korea" , "Korea, Democratic People's Republic of",
                        "Bolivia" ,  "Bolivia, Plurinational States of",
                        "Cote dIvoire"  , "Cote d'Ivoire" ,
                        "TFYR Macedonia" ,  "Macedonia, the former Yugoslav Republic of",    
                        "Syria",  "Syrian Arab Republic",              
                        "Antigua" ,  "Antigua and Barbuda" ,           
                        "British Virgin Islands", "Virgin Islands, British",
                        "Brunei", "Brunei Darussalam",
                        "Gaza Strip", "Palestinian Territory, Occupied"  ,
                        "Moldova" ,  "Republic of Moldova",
                        "Macao SAR", "Macao",
                        "Former USSR"  ,  "USSR"  
                      ))
                        
names_to_update = global_names_to_match %in% name_mapper[, 1]
global_names_to_match[names_to_update] = name_mapper[, 2]

#data_to_use = rbind(MRIO_domestic_data$Aggregated_Threats[sorted_domestic_inds], MRIO_import_data$Aggregated_Threats[sorted_import_inds])[, country_num:1]

country_scale_data$UN_codes = as.character(name_data$UN_code[match(global_names_to_match, name_data$country_name)])

mapped_data <- joinCountryData2Map(country_scale_data, joinCode = "ISO3",  nameJoinColumn = "UN_codes")

pdf('~/GitHub/MRIO_BIO_SATELLITE/MRIO_output_tables/imported_consumption_proportion.pdf', width = 8, height = 6)
par(mai=c(0,0,0.2,0), xaxs="i",yaxs="i", family="Times")
mapCountryData(mapped_data, nameColumnToPlot = "imported_consumption_proportion", colourPalette = "white2Black")
graphics.off()

pdf('~/GitHub/MRIO_BIO_SATELLITE/MRIO_output_tables/imported_proportion_from_low_income.pdf', width = 8, height = 6)
par(mai=c(0,0,0.2,0), xaxs="i",yaxs="i", family="Times")
mapCountryData(mapped_data, nameColumnToPlot = "imported_proportion_from_low_income", colourPalette = "white2Black")
graphics.off()

pdf('~/GitHub/MRIO_BIO_SATELLITE/MRIO_output_tables/proportion_exported.pdf', width = 8, height = 6)
par(mai=c(0,0,0.2,0), xaxs="i",yaxs="i", family="Times")
mapCountryData(mapped_data, nameColumnToPlot = "proportion_exported", colourPalette = "white2Black")
graphics.off()



pdf('~/GitHub/MRIO_BIO_SATELLITE/MRIO_output_tables/imports_exports.pdf', width = 8, height = 8)
par(mar = c(5,10,5,5), family="Times", font = 1)
bar_plot_lims = c(-7000, 25000)

barplot(country_scale_data$global_consumption[country_num:1], names = plot_names, xlim = bar_plot_lims, col = "white", las = 2, horiz = TRUE, width = 2, space=rep(0.3, country_num))
barplot(-country_scale_data$exported[country_num:1], names = plot_names, xlim = bar_plot_lims, col = 'gray70', las = 2, horiz = TRUE, width = 2, add = TRUE, space=rep(0.3, country_num))
barplot(country_scale_data$net_consumption[country_num:1], names = plot_names, xlim = bar_plot_lims, col = '#FF0000', las = 2, horiz = TRUE,  width = 2, add = TRUE, space=rep(0.3, country_num))
barplot(country_scale_data$imported_consumption[country_num:1], names = plot_names, xlim = bar_plot_lims, col = "#FFA500", las = 2, horiz = TRUE,  width = 2,add = TRUE, space=rep(0.3, country_num))
barplot(country_scale_data$low_income_imports[country_num:1], names = plot_names, xlim = bar_plot_lims, col = "#FFD700", las = 2, horiz = TRUE,  width = 2,add = TRUE, space=rep(0.3, country_num))
legend("bottomright", legend = c('global_footprint', 'exported_footprint', 'net_footprint', 'imported_footprint', 'low_income_imported_footprint'), 
       fill = c("white", 'gray70', '#FF0000', "#FFA500", "#FFD700"))

graphics.off()

#FF7F50	rgb(255,127,80) # coral	
#FF6347	rgb(255,99,71) # tomato	
#FF4500	rgb(255,69,0) # orangered	
#FFD700	rgb(255,215,0) # gold	
#FFA500	rgb(255,165,0) # orange	
#FF8C00	rgb(255,140,0) # darkorange	




# aggregated_tradepaths = setNames(aggregate(China_construction_data$Aggregated_Threats_Per_Tradepath~China_construction_data$Production_Country, 
#                                            sum, data = China_construction_data), 
#                                 c('Production_Country', 'Aggregated_Threats'))
# 
# aggregated_tradepaths$UN_codes = as.character(name_data$UN_code[match(aggregated_tradepaths$Production_Country, name_data$country_name)])
# 
# mapped_data <- joinCountryData2Map(aggregated_tradepaths, joinCode = "ISO3",  nameJoinColumn = "UN_codes")
# par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i", family="Times")
# 
# pdf('~/GitHub/MRIO_BIO_SATELLITE/china_construction_threats.pdf', width = 8.3, height = 5)
#   mapCountryData(mapped_data, nameColumnToPlot = "Aggregated_Threats")
# graphics.off()
# 
# 
# aggregated_tradepaths = setNames(aggregate(USA_food_data$Aggregated_Threats_Per_Tradepath~USA_food_data$Production_Country, sum, 
#                                            data = USA_food_data), 
#                                  c('Production_Country', 'Aggregated_Threats'))
# 
# aggregated_tradepaths$UN_codes = as.character(name_data$UN_code[match(aggregated_tradepaths$Production_Country, name_data$country_name)])
# 
# mapped_data <- joinCountryData2Map(aggregated_tradepaths, joinCode = "ISO3",  nameJoinColumn = "UN_codes")
# par(mai=c(0,0,0.2,0),xaxs="i",yaxs="i", family="Times")
# 
# pdf('~/GitHub/MRIO_BIO_SATELLITE/America_food_threats.pdf', width = 8.3, height = 5)
#   mapCountryData(mapped_data, nameColumnToPlot = "Aggregated_Threats")
# graphics.off()
# 
# 
# 
# 
# 
# pdf('~/GitHub/MRIO_BIO_SATELLITE/imports_exports.pdf', width = 8.3, height = 11.7)
# 
# 










# sorted_export_inds = match(MRIO_import_data$Consumption_Country, MRIO_export_data$Consumption_Country)
# sorted_export_threats = MRIO_export_data$Aggregated_Threats[sorted_export_inds]
# 
# sorted_import_inds = match(MRIO_export_data$Consumption_Country, MRIO_import_data$Consumption_Country)
# sorted_import_threats = MRIO_import_data$Aggregated_Threats[sorted_import_inds]
# 
# concatenated_data = rbind(abind(MRIO_import_data$Aggregated_Threats, rev(sorted_import_threats)),
#                           abind(-sorted_export_threats, -rev(MRIO_export_data$Aggregated_Threats)))
#                           
# concatenated_names = abind(MRIO_import_data$Consumption_Country, rev(MRIO_export_data$Consumption_Country))
#   
# last_country = 20
# plot_vector = c(1:last_country, (dim(concatenated_data)[2] - last_country):dim(concatenated_data)[2])
# 
# pdf('~/GitHub/MRIO_BIO_SATELLITE/imports_exports.pdf', width = 8.3, height = 11.7)
# par(mar = c(10,10,10,10)) 
# barplot(concatenated_data[1, plot_vector], xlim = c(-1500, 3000), names = concatenated_names[plot_vector], col = "red", las = 2, 
#         horiz = TRUE)
# 
# barplot(concatenated_data[2, plot_vector], xlim = c(-1500, 3000), names = concatenated_names[plot_vector], col = "black", las = 2, 
#         horiz = TRUE, add = TRUE)
