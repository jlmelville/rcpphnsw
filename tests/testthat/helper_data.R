uiris <- unique(iris)
uirism <- as.matrix(uiris[, -5])

# ten iris entries where the 4 nearest neighbors are distinct
ui10 <- uirism[6:15, ]
