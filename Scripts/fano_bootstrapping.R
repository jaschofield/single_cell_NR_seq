### estimates and bootstrapping of Fano
library(tidyverse)
library(boot)

Nascent_filt <- read.csv(file = "Nascent_filt.csv")

nascent_format <- Nascent_filt %>% dplyr::select(-EF)

burst_function <- function(input_df, indices)
{
  df <- input_df[indices, ]
  c(var(df)/mean(df))
}

boot_out <- data.frame()
for(i in 1:nrow(nascent_format)){
  tryCatch({
  input_df <- nascent_format[i,] %>% t() %>% data.frame()
  bootstrap <- boot(input_df, burst_function, R = 500)
  conf <- boot.ci(boot.out = bootstrap, type ="norm")
  data_add <- data.frame(i, conf$t0, conf$normal[2], conf$normal[3])
  boot_out <- rbind(boot_out, data_add)
  print(i)
  }, error=function(e){})
}



names_df <- Nascent_filt %>% dplyr::select(EF) %>% data.frame()
names_df <- names_df %>% mutate(i = as.numeric(rownames(.)))
colnames(names_df) <- c("gene", "i")

boot_out_join <- full_join(boot_out, names_df)

write.csv(boot_out_join, file = "Fano_bootstrap.csv")
