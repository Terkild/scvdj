library("dplyr")
library("tidyr")

sapply(list.files(path=here::here("R/"), pattern=".R$", full.names=TRUE), source)

vdj.parsed <- LoadCRVDJ(path=here::here("data/cellranger4/BCR")) %>% 
	FilterCRVDJ() %>% 
	OrderCRVDJ() %>%  
	ParseCRVDJByCell()

dim(vdj.parsed)