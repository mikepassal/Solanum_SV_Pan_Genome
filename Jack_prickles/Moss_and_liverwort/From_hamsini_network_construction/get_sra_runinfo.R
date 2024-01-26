# help from    https://www.ncbi.nlm.nih.gov/books/NBK242621/
# list of species needed
input_file = 'query_results.csv'
if (file.info(input_file)$size>1){
  df = read.delim(input_file, sep = ',', header = T)
  df = df[df$Tumor=='no',]
  spetax = df$ScientificName[1] # df$TaxID[1]
  spetax = gsub(' ', '_', spetax)
  table1 = table(df$SRAStudy)
  
  sra = names(table1)
  sranum = as.vector(table1)
  
  len1 = sum(sranum>=10)
  
  if(len1>=25){
    newfilenm = paste0(spetax, '_sraruninfo.csv')
    file.rename(input_file, newfilenm)
  } else{
    file.remove(input_file)
  }
} else{
  file.remove(input_file)
}