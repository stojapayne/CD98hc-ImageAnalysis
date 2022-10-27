install.packages("dplyr")
install.packages("readxl")
install.packages("tidyverse")
install.packages("ggplot2")
install.packages("ggpubr")
install.packages("reshape2")


library(dplyr)
library(readxl)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(reshape2)

 
####Membrane all WT####
membranesigWT<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "Membrane WT All")

na.omit.membranesigWT<- lapply(membranesigWT, na.omit)
chunk_number=100

membrane.splitallWT<-lapply(na.omit.membranesigWT, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE)))
membrane.splitall2WT<-lapply(membrane.splitallWT, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

membrane.splitall2meanWT<-sapply(membrane.splitall2WT, rowMeans) #take the average of the binned 100 pixel length sections

write.csv(membrane.splitall2meanWT, "MembraneXsectionWT.csv")

membraneWT<-sapply(seq(1, ncol(membrane.splitall2meanWT), 3), function(j) rowMeans(membrane.splitall2meanWT[, j+(0:2)])) %>% 
  as.data.frame() #take the average of each of the cross sections across corresponding pixel length (average of all 3 first value, second value, average all three, etc)

membraneWT<-membraneWT%>%
  set_names(paste0("MembraneWT",1:ncol(membraneWT))) #rename the columns to represent the mutant cell number

membraneWT$MembraneWTMean<-rowMeans(membraneWT) #calculate the mean of each row (across all 9 cells)

Binned_Pixels_membraneWT<- 1:nrow(membraneWT)

Fluorescence_IntmembraneWT<-membraneWT$MembraneWTMean

ggplot(membraneWT, aes(x=Binned_Pixels_membraneWT, y=Fluorescence_IntmembraneWT))+
  geom_line()+
  geom_point()


write.csv(membraneWT, "membraneWT.csv")

####Membrane all MUT####
membranesigMut<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "Membrane Mut All")

na.omit.membranesigMut<- lapply(membranesigMut, na.omit)
chunk_number=100

membrane.splitallMut<-lapply(na.omit.membranesigMut, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE)))
membrane.splitall2Mut<-lapply(membrane.splitallMut, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

membrane.splitall2meanMut<-sapply(membrane.splitall2Mut, rowMeans) #take the average of the binned 100 pixel length sections

write.csv(membrane.splitall2meanMut, "MembraneXsectionMut.csv")

membraneMut<-sapply(seq(1, ncol(membrane.splitall2meanMut), 3), function(j) rowMeans(membrane.splitall2meanMut[, j+(0:2)])) %>% 
  as.data.frame() #take the average of each of the cross sections across corresponding pixel length (average of all 3 first value, second value, average all three, etc)

membraneMut<-membraneMut%>%
  set_names(paste0("MembraneMut",1:ncol(membraneMut))) #rename the columns to represent the mutant cell number

membraneMut$MembraneMutMean<-rowMeans(membraneMut) #calculate the mean of each row (across all 9 cells)

Binned_Pixels_membraneMut<- 1:nrow(membraneMut)

Fluorescence_IntmembraneMut<-membraneMut$MembraneMutMean

ggplot(membraneMut, aes(x=Binned_Pixels_membraneMut, y=Fluorescence_IntmembraneMut))+
  geom_line()+
  geom_point()


write.csv(membraneMut, "membraneMut.csv")



####Mut second biological replicate####
membranetotalmut2<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "SLC3A2_mut2") #read the 3 x-sections for mutant

lmembranetotalmut2<- lapply(membranetotalmut2, na.omit)#remove all na's  #remove all na values from the file-save as an array

chunk_number=100
wtsplitallmut2<-lapply(lmembranetotalmut2, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts (so each of the 3 x-sections are now binned to 100 equal parts-you will have 3 total per cell)

wtsplitall2mut2<-lapply(wtsplitallmut2, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanmut2<-sapply(wtsplitall2mut2, rowMeans) #take the average of the binned 100 pixel length sections


mut2<-sapply(seq(1, ncol(stsplitall2meanmut2), 3), function(j) rowMeans(stsplitall2meanmut2[, j+(0:2)])) %>% 
  as.data.frame() #take the average of each of the cross sections across corresponding pixel length (average of all 3 first value, second value, average all three, etc)
  
mut2<-mut2%>%
  set_names(paste0("mut",1:ncol(mut2))) #rename the columns to represent the mutant cell number

mut2$MutMean<-rowMeans(mut2) #calculate the mean of each row (across all 9 cells)

mut2$MutSD<-apply(mut2[1:ncol(mut2)], 1, sd)

Binned_Pixelsmut2<- 1:nrow(mut2)

Fluorescence_Intmut2<-mut2$MutMean

ggplot(mut2, aes(x=Binned_Pixelsmut2, y=Fluorescence_Intmut2))+
  geom_line()+
  geom_point()
  

write.csv(mut2, "mut2.csv")
write.csv(stsplitall2meanmut2,"Cross section mutant 2.csv")



####Mut first biological replicate####
membranetotalmut1<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "SLC3A2_mut1")

lmembranetotalmut1<- lapply(membranetotalmut1, na.omit)#remove all na's 

chunk_number=100
wtsplitallmut1<-lapply(lmembranetotalmut1, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts

wtsplitall2mut1<-lapply(wtsplitallmut1, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanmut1<-sapply(wtsplitall2mut1, rowMeans)


mut1<-sapply(seq(1, ncol(stsplitall2meanmut1), 3), function(j) rowMeans(stsplitall2meanmut1[, j+(0:2)])) %>% 
  as.data.frame()

mut1<-mut1%>%
  set_names(paste0("mut",1:ncol(mut1)))

mut1$MutMean<-rowMeans(mut1)

mut1$MutSD<-apply(mut1[1:ncol(mut1)], 1, sd)

Binned_Pixelsmut1<- 1:nrow(mut1)

Fluorescence_Intmut1<-mut1$MutMean

ggplot(mut1, aes(x=Binned_Pixelsmut1, y=Fluorescence_Intmut1))+
  geom_line()+
  geom_point()

write.csv(mut1, "mut1.csv")
write.csv(stsplitall2meanmut1, "Cross section mutant 1.csv")





####Wildtype second biologycal replicate####
membranetotalwt2<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "SLC3A2_WT2")

lmembranetotalwt2<- lapply(membranetotalwt2, na.omit)#remove all na's 

chunk_number=100
wtsplitall2<-lapply(lmembranetotalwt2, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts

wtsplitall2<-lapply(wtsplitall2, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanwt2<-sapply(wtsplitall2, rowMeans)


wt2<-sapply(seq(1, ncol(stsplitall2meanwt2), 3), function(j) rowMeans(stsplitall2meanwt2[, j+(0:2)])) %>% 
  as.data.frame()

wt2<-wt2%>%
  set_names(paste0("WT",1:ncol(wt2)))

wt2$WTmean<-rowMeans(wt2)

wt2$WTsd<-apply(wt2[1:ncol(wt2)], 1, sd)

Binned_Pixelswt2=1:nrow(wt2)
Average_Fluorescence_Intensitywt2=wt2$WTmean

ggplot(wt2, aes(x=Binned_Pixelswt2, y=Average_Fluorescence_Intensitywt2))+
  geom_line(data=mut2, aes(x=Binned_Pixelsmut2, y=Fluorescence_Intmut2, color="red"))+
  geom_line()+
  geom_point()+
  theme_bw()+
  ylim(0,25)

write.csv(wt2, "wt2.csv")
write.csv(stsplitall2meanwt2, "Cross section WT 2.csv")


####Wildtype FIRST biologycal replicate####
membranetotalwt1<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "SLC3A2_WT1")

lmembranetotalwt1<- lapply(membranetotalwt1, na.omit)#remove all na's 

chunk_number=100
wtsplitall1<-lapply(lmembranetotalwt1, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts

wtsplitall1<-lapply(wtsplitall1, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanwt1<-sapply(wtsplitall1, rowMeans)


wt1<-sapply(seq(1, ncol(stsplitall2meanwt1), 3), function(j) rowMeans(stsplitall2meanwt1[, j+(0:2)])) %>% 
  as.data.frame()

wt1<-wt1%>%
  set_names(paste0("WT",1:ncol(wt1)))

wt1$WTmean<-rowMeans(wt1)

wt1$WTsd<-apply(wt1[1:ncol(wt1)], 1, sd)

Binned_Pixelswt1=1:nrow(wt1)
Average_Fluorescence_Intensitywt1=wt1$WTmean

ggplot(wt1, aes(x=Binned_Pixelswt1, y=Average_Fluorescence_Intensitywt1))+
  geom_line(data=mut1, aes(x=Binned_Pixelsmut1, y=Fluorescence_Intmut1, color="red"))+
  geom_line()+
  theme_bw()+
  ylim(0,25)



######FOR MUTANT ALL########

membranetotalmutall<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "MUT_ALL") #read the 3 x-sections for mutant

lmembranetotalmutall<- lapply(membranetotalmutall, na.omit)#remove all na's  #remove all na values from the file-save as an array

chunk_number=100
wtsplitallmutall<-lapply(lmembranetotalmutall, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts (so each of the 3 x-sections are now binned to 100 equal parts-you will have 3 total per cell)

wtsplitall2mutall<-lapply(wtsplitallmutall, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanmutall<-sapply(wtsplitall2mutall, rowMeans) #take the average of the binned 100 pixel length sections


mutall<-sapply(seq(1, ncol(stsplitall2meanmutall), 3), function(j) rowMeans(stsplitall2meanmutall[, j+(0:2)])) %>% 
  as.data.frame() #take the average of each of the cross sections across corresponding pixel length (average of all 3 first value, second value, average all three, etc)

mutall<-mutall%>%
  set_names(paste0("mut",1:ncol(mutall))) #rename the columns to represent the mutant cell number

mutall$MutMean<-rowMeans(mutall) #calculate the mean of each row (across all 9 cells)

mutall$MutSD<-apply(mutall[1:ncol(mutall)], 1, sd)

write.csv(mutall, "mutall.csv")


#####FOR WILDTYPE ALL####

membranetotalwtall<-read_excel("ImageJ slc3a2 mut vs wt analysis 2022-05-19 first 2 and last 2 folders.xlsx", sheet = "WT_ALL")

lmembranetotalwtall<- lapply(membranetotalwtall, na.omit)#remove all na's 

chunk_number=100
wtsplitall<-lapply(lmembranetotalwtall, function(x) split(x, cut(seq_along(x), chunk_number, labels=FALSE))) #split each list element to 100 equal parts

wtsplitall<-lapply(wtsplitall, function(x) do.call(rbind, x)) #order all split values to each corresponding cell as a matrix

stsplitall2meanwtall<-sapply(wtsplitall, rowMeans)


wtall<-sapply(seq(1, ncol(stsplitall2meanwtall), 3), function(j) rowMeans(stsplitall2meanwtall[, j+(0:2)])) %>% 
  as.data.frame()

wtall<-wtall%>%
  set_names(paste0("WT",1:ncol(wtall)))

wtall$WTmean<-rowMeans(wtall)

wtall$WTsd<-apply(wtall[1:ncol(wtall)], 1, sd)


write.csv(wtall, "wtall.csv")



####PLOT####

#for grid arrange : https://stackoverflow.com/questions/14743060/r-ggplot-graphs-sharing-the-same-y-axis-but-with-different-x-axis-scales
#http://www.sthda.com/english/wiki/wiki.php?id_contents=7930

#Henry's data: graph inspo - https://github.com/Bishop-Laboratory/RLoop-QC-Meta-Analysis-Miller-2022/blob/main/figures.R#L344-L437

gridplotall<-read_excel("gridplot.xlsx", sheet = "SLC3A2")

gridplotallmembrane <- read_excel("gridplot.xlsx", sheet = "Membrane")

ggp<-melt(gridplotall, id.vars="...1")

ggpmembrane <- melt(gridplotallmembrane, id.vars="...1")

colnames(ggp)[1]="Bins"

colnames(ggpmembrane)[1] = "Bins"


my_ggp_fun <- function(my_data, my_x, my_y, facetwrap, nrow, ncol, group) {    # Create user-defined function; from https://statisticsglobe.com/pass-column-names-indices-user-defined-ggplot2-function-r

  Binned_Pixels <- my_data[ , my_x]
  Average_Fluorescence_Intensity <- my_data[ , my_y]
  groupby = my_data[, group]
  
  if (missing(group) ) {

  ggplot(my_data, aes(Binned_Pixels, Average_Fluorescence_Intensity)) +
    geom_line(size = 0.2)+
    facet_wrap(vars({{facetwrap}}), nrow= nrow, ncol = ncol)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
    theme_light(base_size = 10)+
    labs(x= element_blank(), y= element_blank())+
    theme(plot.title = element_text(hjust=0.5),
          panel.border=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line (colour="black"),
          strip.background = element_blank(), 
          strip.text = element_text(colour = "black"), 
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
         panel.spacing = unit(3, "pt"),
         panel.background = element_rect(color = "black")
    ) }
  else { 
    ggplot(my_data, aes(Binned_Pixels, Average_Fluorescence_Intensity, group = groupby)) +
      geom_line(size = 0.2)+
      facet_wrap(vars({{facetwrap}}), nrow= nrow, ncol = ncol)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
      theme_light(base_size = 10)+
      labs(x= element_blank(), y= element_blank())+
      scale_y_continuous(limits = c(0, 40))+
      theme(plot.title = element_text(hjust=0.5),
            panel.border=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line (colour="black"),
            strip.background = element_blank(), 
            strip.text = element_text(colour = "black"), 
            axis.ticks.x = element_blank(),
            axis.text.x = element_blank(),
            panel.spacing = unit(3, "pt"),
            panel.background = element_rect(color = "black")
      )
    
  }#Theme from https://felixfan.github.io/ggplot2-remove-grid-background-margin/ & https://ggplot2-book.org/polishing.html

  }

#### PLOT 1 - EACH INDIVIDUAL CROSS SECTION ####
p1.1<- ggp[!grepl("WTmean|MutMean", ggp$variable),] 

p1.1mem <- ggpmembrane[!grepl("WTmean|MutMean", ggp$variable),]


p1.1$variable <- as.character(p1.1$variable)

p1.1mem$variable <- as.character(p1.1mem$variable)
p2.2<- p1.1 %>%
  mutate(condition = case_when(
    str_detect(variable, "WT")  ~ "CD98hc" ,
    str_detect(variable, "Mut") ~ "CD98hc",
    # TRUE ~ variable
  )) 

p2.2mem <-p1.1mem %>%
  mutate(condition = case_when(
    str_detect(variable, "WT")  ~ "Membrane" ,
    str_detect(variable, "Mut") ~ "Membrane",
    # TRUE ~ variable
  ))

p2.2mem$value<-p2.2mem$value / 2.6
p2.2mem <- p2.2mem %>%
  mutate(variable = case_when(
    str_detect(variable, "5.1.2_005") ~ "WT 1", 
    str_detect(variable, "5.2.2_006") ~ "WT 2",
    str_detect(variable, "5.3.2_005") ~ "WT 3",
    str_detect(variable, "4.1.2_007") ~ "WT 4",
    str_detect(variable, "4.2.2_007") ~ "WT 5",
    str_detect(variable, "4.3.2_007") ~ "WT 6",
    str_detect(variable, "6.1.2.002") ~ "Mut 1",
    str_detect(variable, "6.2.2.002") ~ "Mut 2",
    str_detect(variable, "6.3.2.002") ~ "Mut 3",
    str_detect(variable, "2.1.004") ~ "Mut 4",
    str_detect(variable, "2.2.004") ~ "Mut 5",
    str_detect(variable, "2.3.004") ~ "Mut 6",
  ))
  
p2.2 <- p2.2 %>%
  mutate(variable = case_when(
    str_detect(variable, "5.1.2_005") ~ "WT 1", 
    str_detect(variable, "5.2.2_006") ~ "WT 2",
    str_detect(variable, "5.3.2_005") ~ "WT 3",
    str_detect(variable, "4.1.2_007") ~ "WT 4",
    str_detect(variable, "4.2.2_007") ~ "WT 5",
    str_detect(variable, "4.3.2_007") ~ "WT 6",
    str_detect(variable, "6.1.2.002") ~ "Mut 1",
    str_detect(variable, "6.2.2.002") ~ "Mut 2",
    str_detect(variable, "6.3.2.002") ~ "Mut 3",
    str_detect(variable, "2.1.004") ~ "Mut 4",
    str_detect(variable, "2.2.004") ~ "Mut 5",
    str_detect(variable, "2.3.004") ~ "Mut 6",
  ))








DF <- bind_rows(p2.2mem, p2.2)

# # my_linetype <- setNames(c("solid", "solid", "dashed", "dashed"), unique(DF$condition))
# # my_color <- setNames(c("blue", "red"), unique(DF$condition))
# # 
# # my_color  
# # 
# p1<- my_ggp_fun(my_data = p2.2, my_x = "Bins", my_y = "value",
#                 facetwrap = variable, nrow = 2, ncol = 6)+
#   labs(title = "Cell Cross Sections", y = "Average Fluorescence Intensity")
# # 
# # 
# p1
# 
# p2 <- my_ggp_fun(my_data = p2.2mem, my_x = "Bins", my_y = "value",
#                  facetwrap = variable, nrow = 2, ncol = 6)+
#   labs(title = "Cell Cross Sections", y = "Average Fluorescence Intensity")
# 
# 
# p2
# 
# g1 <- ggplotGrob(p1)
# g2 <- ggplotGrob(p2)
# g <- rbind(g1, g2, size="first") # stack the two plots
# g$widths <- unit.pmax(g1$widths, g2$widths) # use the largest widths
# # center the legend vertically
# g$layout[grepl("guide", g$layout$name),c("t","b")] <- c(1,nrow(g))
# grid.newpage()
# grid.draw(g)
# 
# 
# 
# # DF %>% 
# #   ggplot(aes(x = DF[ , "Bins"], y = DF[ , "value"], 
# #              color = condition,
# #              linetype = condition))+
# #   geom_line() +
# #   facet_wrap(~ variable,  scales = "free_y") +
# #   scale_linetype_manual(values = my_linetype) +
# #   scale_color_brewer(palette = 'Dark2') +
# #   theme_bw(base_size = 14) +
# #   theme(legend.position = 'bottom') +
# #   theme(legend.key.size = unit(2.5, 'lines'))
# 
# scale = 3
# ggplot(DF, aes(DF[ , "Bins"], DF[ , "value"], linetype = condition, color = condition)) +
#   scale_color_manual(values = c("CD98hc" = "red", "Membrane" = "blue"))+
#   scale_linetype_manual(values = c("Membrane" = "dashed", "CD98hc" = "solid"))+
#   geom_line(size = 0.2)+
#   ggh4x::facet_wrap2(~variable)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
#   theme_light(base_size = 10)+
#   labs(x= element_blank(), y= element_blank())+
#   #secondary axis: https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
#    #secondary axis: https://stackoverflow.com/questions/61565304/using-secondary-y-axis-in-ggplot2-with-different-scale-factor-when-using-facet-w
#   theme(plot.title = element_text(hjust=0.5),
#         axis.text.y.right = element_text(color = "blue"),
#         axis.text.y.left = element_text(color = "red"),
#         panel.border=element_blank(),
#         panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         strip.background = element_blank(),
#         strip.text = element_text(colour = "black"),
#         axis.ticks.x = element_blank(),
#         axis.text.x = element_blank(),
#         panel.spacing = unit(3, "pt"),
#         panel.background = element_rect(color = "black")
#   )

scale = 2.6

p1<-ggplot(DF, aes(DF[ , "Bins"], DF[ , "value"], linetype = condition, color = condition)) +
  scale_color_manual(values = c("Membrane" = "blue", "CD98hc" = "red"))+
  scale_linetype_manual(values = c("Membrane" = "dashed", "CD98hc" = "solid"))+
  geom_line(size = 0.2)+
  facet_wrap(~fct_relevel(variable, "WT 1", "WT 2", "WT 3", "WT 4", "WT 5", "WT 6", "Mut 1", "Mut 2", 
                                   "Mut 3", "Mut 4", "Mut 5", "Mut 6"), nrow = 2, ncol = 6)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
  theme_light(base_size = 10)+
  labs(x= element_blank(), y= element_blank())+
  scale_y_continuous(limits = c(0, 40), sec.axis = sec_axis(~.*scale))+#secondary axis: https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
   #secondary axis: https://stackoverflow.com/questions/61565304/using-secondary-y-axis-in-ggplot2-with-different-scale-factor-when-using-facet-w
  labs(title = "Cell Cross Sections", y = "CD98hc Average F.I")+
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y.right = element_blank(),
        axis.text.y.left = element_text(color = "red"),
        axis.title.y = element_text (color = "red"),
        axis.title.y.right = element_blank(),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(3, "pt"),
        panel.background = element_rect(color = "black"), 
        legend.position = "none"
  )


p1


#### PLOT 2 & 3 - INDIVIDUAL CROSS SECTIONS OVERLAYED WT and Mut ####


# p2 <- my_ggp_fun(my_data = p2.2, my_x = "Bins", my_y = "value", 
#                  facetwrap = condition, nrow = 2, ncol = 1) +
#   theme(axis.ticks = element_blank(), 
#         axis.text = element_blank()) +
#   labs(title = "Overlayed (6)")
# 
# p2

# my_ggp_fun(my_data = DF2, my_x = "Bins", my_y = "value", facetwrap = condition3, nrow = 2, ncol = 1, group = "variable")


DF2<- DF %>%
  mutate(condition2 = case_when(
    str_detect(variable, "WT")  ~ "WT" ,
    str_detect(variable, "Mut") ~ "Mut",
    # TRUE ~ variable
  )) %>%
  mutate (condition3 = case_when(
    str_detect(variable, "WT 1") & str_detect(condition, "Membrane") ~ "MembraneWT1",
    str_detect(variable, "WT 2") & str_detect(condition, "Membrane") ~ "MembraneWT2",
    str_detect(variable, "WT 3") & str_detect(condition, "Membrane") ~ "MembraneWT3",
    str_detect(variable, "WT 4") & str_detect(condition, "Membrane") ~ "MembraneWT4",
    str_detect(variable, "WT 5") & str_detect(condition, "Membrane") ~ "MembraneWT5",
    str_detect(variable, "WT 6") & str_detect(condition, "Membrane") ~ "MembraneWT6",
    str_detect(variable, "Mut 1") & str_detect(condition, "Membrane")~ "MembraneMut1",
    str_detect(variable, "Mut 2") & str_detect(condition, "Membrane")~ "MembraneMut2",
    str_detect(variable, "Mut 3") & str_detect(condition, "Membrane")~ "MembraneMut3",
    str_detect(variable, "Mut 4") & str_detect(condition, "Membrane")~ "MembraneMut4",
    str_detect(variable, "Mut 5") & str_detect(condition, "Membrane")~ "MembraneMut5",
    str_detect(variable, "Mut 6") & str_detect(condition, "Membrane")~ "MembraneMut6",
    str_detect(variable, "WT 1") & str_detect(condition, "CD98hc") ~ "CD98hcWT1",
    str_detect(variable, "WT 2") & str_detect(condition, "CD98hc") ~ "CD98hcWT2",
    str_detect(variable, "WT 3") & str_detect(condition, "CD98hc") ~ "CD98hcWT3",
    str_detect(variable, "WT 4") & str_detect(condition, "CD98hc") ~ "CD98hcWT4",
    str_detect(variable, "WT 5") & str_detect(condition, "CD98hc") ~ "CD98hcWT5",
    str_detect(variable, "WT 6") & str_detect(condition, "CD98hc") ~ "CD98hcWT6",
    str_detect(variable, "Mut 1") & str_detect(condition, "CD98hc") ~ "CD98hcMut1",
    str_detect(variable, "Mut 2") & str_detect(condition, "CD98hc") ~ "CD98hcMut2",
    str_detect(variable, "Mut 3") & str_detect(condition, "CD98hc") ~ "CD98hcMut3",
    str_detect(variable, "Mut 4") & str_detect(condition, "CD98hc") ~ "CD98hcMut4",
    str_detect(variable, "Mut 5") & str_detect(condition, "CD98hc") ~ "CD98hcMut5",
    str_detect(variable, "Mut 6") & str_detect(condition, "CD98hc") ~ "CD98hcMut6",
  ))


scale = 2.6
p2<-ggplot(DF2, aes(DF2[ , "Bins"], DF2[ , "value"], linetype = condition, color = condition, group = condition3)) +
  scale_color_manual(values = c("Membrane" = "blue", "CD98hc" = "red"))+
  scale_linetype_manual(values = c("Membrane" = "dashed", "CD98hc" = "solid"))+
  geom_line(size = 0.2)+
  facet_wrap(~fct_rev(condition2), nrow = 2, ncol = 1)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
  theme_light(base_size = 10)+
  labs(title = "Overlayed (6)", x= element_blank(), y= element_blank())+
  scale_y_continuous(limits = c(0, 40), sec.axis = sec_axis(~.*scale))+#secondary axis: https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
  #secondary axis: https://stackoverflow.com/questions/61565304/using-secondary-y-axis-in-ggplot2-with-different-scale-factor-when-using-facet-w
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y.right = element_text(color = "blue", size = 8),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(3, "pt"),
        panel.background = element_rect(color = "black"),
        legend.position = "none",
        axis.ticks = element_blank(), 
        axis.text = element_blank()
  )


p2

#### Plot 4 and 5 - Averaged (all) ####

p3.1<- ggp[grepl("WTmean|MutMean", ggp$variable),]

p3.1mem<- ggpmembrane[grepl("WTmean|MutMean", ggp$variable),]

p3.1 <- p3.1 %>% mutate(condition = case_when(
    str_detect(variable, "WTmean")  ~ "CD98hc" ,
    str_detect(variable, "Mut") ~ "CD98hc"
)
)

p3.1mem <- p3.1mem %>% mutate(condition = case_when(
  str_detect(variable, "WTmean")  ~ "Membrane" ,
  str_detect(variable, "Mut") ~ "Membrane"
)
)

p3.1mem$value<-p3.1mem$value / 1.5
DF3 <- bind_rows(p3.1mem, p3.1)

scale2=1.5

p3<-ggplot(DF3, aes(DF3[ , "Bins"], DF3[ , "value"], linetype = condition, color = condition)) +
  scale_color_manual(values = c("Membrane" = "blue", "CD98hc" = "red"))+
  scale_linetype_manual(values = c("Membrane" = "dashed", "CD98hc" = "solid"))+
  geom_line(size = 0.2)+
  facet_wrap(~variable, nrow = 2, ncol = 1)+ #Facet wrap reference https://stackoverflow.com/questions/55290911/ggplot-facet-wrap-variable-as-an-argument-in-a-function
  theme_light(base_size = 10)+
  labs(title = "Averaged (all)", x= element_blank(), y= element_blank())+
  scale_y_continuous(limits = c(0, 40), sec.axis = sec_axis(~.*scale2, name = "Membrane Average F.I"))+#secondary axis: https://r-graph-gallery.com/line-chart-dual-Y-axis-ggplot2.html
  #secondary axis: https://stackoverflow.com/questions/61565304/using-secondary-y-axis-in-ggplot2-with-different-scale-factor-when-using-facet-w
  theme(plot.title = element_text(hjust=0.5),
        axis.text.y.right = element_text(color = "blue"),
        axis.text.y.left = element_blank(),
        axis.title.y = element_text(color = "red"),
        axis.title.y.right = element_text(color = "blue", vjust = +2),
        panel.border=element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black"),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        panel.spacing = unit(3, "pt"),
        panel.background = element_rect(color = "black"),
        legend.position = "right"
  )



p3






# p3 <- my_ggp_fun(my_data=p4.1, my_x = "Bins", my_y =  "value", 
#                  facetwrap = condition, nrow = 2, ncol = 1) +
#   theme(axis.ticks = element_blank(), 
#         axis.text = element_blank()) +
#   labs(title = "Averaged (all)")
#   
# 
# p3


ggpubr::ggarrange(p1, p2, p3, nrow = 1, widths = c(4.5, 1, 1.6)) 
                  
#save 1400 400