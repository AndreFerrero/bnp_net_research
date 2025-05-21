functions_folder = here("code/functions")
funs <- list.files(functions_folder, pattern = "\\.R$", full.names = TRUE)

for(f in funs){
  source(f)
}

rm(f)
rm(funs)