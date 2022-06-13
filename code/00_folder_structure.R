for (i in c("data", "output")) {
  if (!file.exists(i)) {
    dir.create(path = i, recursive = T)
  }
}
