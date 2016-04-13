# starting point: compiled df with module assignments:
head(df_moduled)

# need min and max of each Condition_Feature
library(dplyr)
by_condition_feature <- group_by(df_moduled, Condition_Feature)
stats_df <- summarise(by_condition_feature,
                      min_val = min(value),
                      max_val = max(value))
df_all <- df_moduled %>% merge(., stats_df) %>% mutate(., value_adjusted = (value - min_val)/(max_val - min_val))

# aggregated module scores:
by_module <- group_by(df_all, Donor, Module)
module_scores <- summarise(by_module,
                           module_score = mean(value_adjusted))
module_scores <- merge(dplyr::select(df_moduled, Donor, Gender), module_scores)
module_scores$Module <- as.factor(module_scores$Module)
p <- ggplot(module_scores, aes(x = Module, y = module_score, group = Donor))
p + geom_line()
p + geom_line(aes(colour = Gender))

d1 <- dplyr::filter(module_scores, Donor == "W3231")
d2 <- dplyr::filter(module_scores, Donor == "7605_2")
p + geom_line(colour = "gray") + geom_line(data = d1, colour = "red") + geom_line(data = d2, colour = "blue")

p <- ggplot(module_scores, aes(x = Module, y = module_score))
p + geom_boxplot()
p + geom_boxplot(colour = "gray") + geom_point(data = d1, colour = "red") + geom_point(data = d2, colour = "blue")
p + geom_boxplot(aes(colour = Gender))

# looking at indiv modules:
df_moduled %>% dplyr::filter(Module == 7) %>% select(Condition_Feature) %>% unique()

# test is they are different among the sexes:
 for (i in 0:11){
   x <- module_scores %>% dplyr::filter(Module == i) %>% dplyr::filter(Gender == "F") %>% dplyr::select(module_score)
   y <- module_scores %>% dplyr::filter(Module == i) %>% dplyr::filter(Gender == "M") %>% dplyr::select(module_score)
   print("Module")
   print(i)
   print(t.test(x,y)$p.value)
 }
# 2,3,4,5,6,8,9,10 are all super significant... whoa

