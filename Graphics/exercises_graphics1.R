## Exercise I
## ──────────

# Upload the file EconomistData.csv from your elearning platform 
##   These data consist of /Human Development Index/ and /Corruption
##   Perception Index/ scores for several countries.
#Transparency International’s (TI) Corruption Perceptions Index1 (CPI) 
#is an aggregate indicator that
#ranks countries in terms of the degree to which corruption 
#is perceived to exist among public officials and politicians. 
#It is a composite index drawing on corruption-related data from a variety of reputable
#institutions.
##   Original sources for these data are
##     [http://www.transparency.org/content/download/64476/1031428]

#TO DO
#    0. Read the data
##   1. Create a scatter plot with CPI on the x axis and HDI on the y axis.
##   2. Color the points blue.
##   3. Map the color of the the points to Region.
##   4. Make the points bigger by setting size to 2
##   5. Map the size of the points to HDI.Rank


## Exercise II
## ───────────

##   1. Re-create a scatter plot with CPI on the x axis and HDI on the y
##      axis (as you did in the previous exercise).
##   2. Overlay a smoothing line on top of the scatter plot using
##      geom_smooth.
##   3. Overlay a smoothing line on top of the scatter plot using
##      geom_smooth, but use a linear model for the predictions. Hint: see
##      `?stat_smooth'.
##   4. Overlay a smoothling line on top of the scatter plot using
##      geom_line. Hint: change the statistical transformation.
##   5. Overlay a smoothing line on top of the scatter plot using
##      the default /loess/ method, but make it less smooth. Hint: see
##      `?loess'.



#smoothing method (function) to use, eg. "lm", "glm", "gam", "loess", "rlm".
#For method = "auto" the smoothing method is chosen based on the size of the 
#largest group (across all panels). loess is used for than 1,000 observations; 
#otherwise gam is used with formula = y ~ s(x, bs = "cs"). Somewhat anecdotally, 
#loess gives a better appearance, but is O(n^2) in memory, 
#so does not work for larger datasets.

