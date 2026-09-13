options(stringsAsFactors=FALSE)
for (f in sort(list.files('R', pattern='[.]R$', full.names=TRUE))) {
  suppressPackageStartupMessages(source(f))
}
