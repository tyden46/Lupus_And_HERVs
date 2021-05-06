#### Author: Tyson Dawson
# This script uses the package formattable to make a pretty table and then
# saves it as a png with webshot. Integer values in the table are colored by
# size.
library(formattable)
library(htmltools)
library(webshot)
export_formattable <- function(f, file, background, width = "50%", height = NULL, 
                               delay = 0.2)
{
  w <- as.htmlwidget(f, width = width, height = height)
  path <- html_print(w, background = background, viewer = NULL)
  url <- paste0("file:///", gsub("\\\\", "/", normalizePath(path)))
  webshot(url,
          file = file,
          selector = ".formattable_widget",
          delay = delay)
}
setwd("C:/Users/tyson/OneDrive/Desktop/Lupus")
metadata=read.csv("StudySummaries.csv")
customGreen0 = "#007600"
customGreen = "#00ff00"
customRed="#ff0000"
customRed0="#b10000"
myTable=formattable(metadata, align =c("l","l","l","l","l", "l", "l", "l", "l"), list(
  `Indicator Name` = formatter("span", style = ~ style(color = "grey",font.weight = "bold")), 
  `SLE`= color_tile(customRed0, customRed),
  `Controls`= color_tile(customGreen0, customGreen)
))
export_formattable(myTable, "StudySummaries.png", "white !important")
