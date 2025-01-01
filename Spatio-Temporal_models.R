library(ggplot2)
library(tidyr)
library(dplyr)
library(gridExtra)
library(INLA)
#inla.setOption(mkl=TRUE)
library(arm)
#library(rgdal)
#library(rgeos)
#library(maptools)
#library(maps)
#library(broom)
library(spdep)
library(sf)
library(viridis)



asthma.dat <- read.csv("/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/California_Asthma/asthma-ed-visit-rates.csv", header=T)
#data with Emergency department visits for 2011-2019

head(asthma.dat)
summary(asthma.dat)
str(asthma.dat)
unique(asthma.dat$Age.Group)


asthma.dat <- asthma.dat %>% 
  filter(Age.Group=="Under 18" & Strata=="Total Population",
         Strata.Name=="Under 18")
#restrict to Under 18 ages and total population


asthma.dat <- asthma.dat %>% dplyr::select(!c("Age.Group","Strata","Strata.Name"))
#remove these columns from the data

cali.pop.dat <- read.csv("/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/California_Asthma/cali_pop_census_race.csv", header=T)
#data for population in the 2020

asthma.dat <- left_join(asthma.dat, cali.pop.dat, by=c("County_Name"="County"))
#join the two data sets


# Loading pollution data
pollution.dat <- read.csv("/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/California_Asthma/cali_pollution.csv", header=T)

# Perform the left join
asthma_pollution.dat <- asthma.dat %>%
  left_join(pollution.dat, by = "County_ID")

# Remove rows with value counts 0
asthma_pollution.dat <- asthma_pollution.dat %>% filter(Counts != 0)



# Subtract 11000 from Counts where Counts are greater than 15000
asthma_pollution.dat$Counts[asthma_pollution.dat$Counts > 15000] <- asthma_pollution.dat$Counts[asthma_pollution.dat$Counts > 15000] - 11000


# Plot a histogram of the "Counts" variable
histogram_plot <- ggplot(asthma_pollution.dat, aes(x = Counts)) + 
  geom_histogram(binwidth = 1, color = "black", fill = "blue") +
  ggtitle("Histogram of Counts") +
  xlab("Counts") +
  ylab("Frequency")

# Print the plot
print(histogram_plot)

head(asthma_pollution.dat)


summary(asthma_pollution.dat)



##### Adjacency Matrix ######

asthma.spatial <- read_sf("/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/CA_Counties/CA_Counties_TIGER2016.shp")

asthma.poly <- poly2nb(asthma.spatial)

nb2INLA("/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/CA_Counties/asthma.graph1",asthma.poly)

asthma.adj <- inla.read.graph(filename="/Users/tomasmock/LSE/Bayesian Data Analysis/Presentation Project/CA_Counties/asthma.graph1")

image(inla.graph2matrix(asthma.adj),xlab="",ylab="")



##### Adding Spatial Structure #####

library(SpatialEpi)
library(spdep)
library(broom)

asthma.sf <- st_as_sf(asthma.spatial)

asthma.sf$GEOID <- as.integer(asthma.sf$GEOID)

asthma_pollution.dat <- left_join(asthma.sf,asthma_pollution.dat, by=c("GEOID"="County_ID") )


max.count <- length(unique(asthma_pollution.dat$GEOID))

look.up <- data.frame(GEOID = unique(asthma_pollution.dat$GEOID), num_id = seq(1:max.count))

asthma_pollution.dat = left_join(asthma_pollution.dat,look.up)

head(asthma_pollution.dat)
summary(asthma_pollution.dat)

# Remove rows with value counts 0
asthma_pollution.dat <- asthma_pollution.dat %>% filter(Counts != 0)

# Check for NA values in specific columns
sum(is.na(asthma_pollution.dat$Counts))
sum(is.na(asthma_pollution.dat$Year))
sum(is.na(asthma_pollution.dat$population))  


# MODELS #

Year1 <- asthma_pollution.dat$Year

num_id1 <- asthma_pollution.dat$num_id

# The number of years
max.year <- max(asthma_pollution.dat$Year)

#Random Walk of order 2
asthma.year.rw2 <- Counts ~ 1 + white + PM2.5 +
  f(num_id, model="bym", graph=asthma.adj) +
  f(Year, model="rw2")+
  f(Year1, model="iid")

asthma.inla.rw2 <- inla(asthma.year.rw2,
                         family="nbinomial",
                         data=asthma_pollution.dat,
                         offset=log(population),
                         control.predictor = list(compute=TRUE),
                         control.compute = list(waic=T))

# Fixed effects
round(asthma.inla.rw2$summary.fixed[,1:5],5)

# WAIC
asthma.inla.rw2$waic$waic

#-----------------------------------------------------------------------

asthma.rw2.iid <- unlist(
  lapply(asthma.inla.rw2$marginals.random$Year1,
         function(x){
           inla.emarginal(exp, x)})
)

asthma.rw2.rw2 <- unlist(
  lapply(asthma.inla.rw2$marginals.random$Year,
         function(x){
           inla.emarginal(exp, x)})
)


asthma.plot <- data.frame(temp=asthma.rw2.rw2, 
                           iid=asthma.rw2.iid, 
                           
                           Year=1:9)
asthma.plot <- pivot_longer(asthma.plot, cols = !Year)

p1 <- ggplot(asthma.plot, 
             aes(x=Year, y=value, linetype=name)) + 
  geom_line()+
  ggtitle("RW2")+
  ylim(0,4.5)

print(p1)



# Subset asthma.dat to match the length of the mean vector
asthma_pollution.dat$rw2.fit <- asthma.inla.rw2$summary.fitted.values$mean





target <- c(2011,2013,2015,2017,2019) 
asthma.fitted<- filter(asthma_pollution.dat, Year %in% target)  %>%
  group_by(num_id, Year) %>%
  summarise_at(c("rw2.fit"), sum)

#we'll categorise for RW
rw.cut <- with(asthma.fitted,quantile(rw2.fit,0:10/10))

asthma.fitted <- asthma.fitted %>%
  mutate(cut.rw2.fit = cut(rw2.fit,rw.cut))

p1<- ggplot(asthma.fitted, aes(fill=rw2.fit, geometry=geometry)) + 
  geom_sf()+
  scale_fill_viridis(option="cividis") +
  theme_void() +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year)  
#categorised plot

p2<- ggplot(asthma.fitted, aes(fill=cut.rw2.fit, geometry=geometry)) + 
  geom_sf()+
  scale_fill_viridis(discrete=TRUE,option="cividis") +
  theme_void() +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year)  

grid.arrange(p1,p2,nrow=2)


#-------------------------------------------

target <- c(2011:2019)
asthma.fitted<- filter(asthma_pollution.dat, Year %in% target)  %>%
  group_by(num_id, Year) %>%
  summarise_at(c("rw2.fit"), sum)

#we'll categorise for RW
rw.cut <- with(asthma.fitted,quantile(rw2.fit,0:10/10))

asthma.fitted <- asthma.fitted %>%
  mutate(cut.rw2.fit = cut(rw2.fit,rw.cut))

p1<- ggplot(asthma.fitted, aes(fill=rw2.fit, geometry=geometry)) + 
  geom_sf()+
  scale_fill_viridis(option="cividis") +
  theme_void() +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year)  
p1

#----------------------------------------------------------------
#linear model

asthma.year.lin <- Counts ~ 1 + Year + PM2.5 + white +
  f(num_id, model="bym", graph = asthma.adj,constr=TRUE) +
  #spatial part
  f(Year1, num_id, model="iid", constr=TRUE)
#the spatial/temporal interaction term with delta coefficients


lcs <- inla.make.lincombs(Year1=diag(9), Year=rep(1,9))

asthma.inla.lin <- inla(asthma.year.lin,
                         family="nbinomial",
                         data=asthma_pollution.dat,
                         offset=log(population),
                         lincomb = lcs,
                         control.predictor = list(compute=TRUE),
                         control.compute = list(waic=T))

# Fixed effects
round(asthma.inla.lin$summary.fixed[,1:5],5)

# WAIC
asthma.inla.lin$waic$waic

# Subset asthma.dat to match the length of the mean vector
asthma_pollution.dat$lin.fit <- asthma.inla.lin$summary.fitted.values$mean

target <- c(2011:2019)
asthma.fitted<- filter(asthma_pollution.dat, Year %in% target)  %>%
  group_by(num_id, Year) %>%
  summarise_at(c("lin.fit"), sum)

#we'll categorise for RW
rw.cut <- with(asthma.fitted,quantile(lin.fit,0:10/10))

asthma.fitted <- asthma.fitted %>%
  mutate(cut.lin.fit = cut(lin.fit,rw.cut))

p1<- ggplot(asthma.fitted, aes(fill=lin.fit, geometry=geometry)) + 
  geom_sf()+
  scale_fill_viridis(option="cividis") +
  theme_void() +
  theme(legend.title = element_blank()) +
  facet_wrap(~Year)  
p1




#-------------------INTERACTIONS RW----------------------

# TYPE I RW

#another Year
asthma_pollution.dat$Year1 <- asthma_pollution.dat$Year

#interaction term
num_id.Year <- paste(asthma_pollution.dat$num_id,asthma_pollution.dat$Year,sep="_")

asthma.typeI <-Counts ~ 1 + white + PM2.5 +
  f(num_id, model="bym", graph=asthma.adj)+
  f(Year, model="rw2") +
  f(Year1, model="iid") +
  f(num_id.Year, model="iid")

asthma.typeI.inla <- inla(asthma.typeI , family="nbinomial",
                           data = asthma_pollution.dat, offset=log(population),
                           control.predictor =list(compute=TRUE),
                           control.compute = list(waic=TRUE),
                           control.inla = list(h=0.1))



# Fixed effects
round(asthma.typeI.inla$summary.fixed[,1:5],5)

# WAIC

asthma.typeI.inla$waic$waic



# TYPE II RW

#temporal 
Year.int <- asthma_pollution.dat$Year
#spatial
num_id.int <- asthma_pollution.dat$num_id




asthma.typeII <- Counts ~ 1 + white + PM2.5 +
  f(num_id, model="bym", graph=asthma.adj) +
  f(Year, model="rw2") +
  f(Year1, model="iid") +
  f(num_id.int, model="iid",
    #iid on the spatial component
    group = Year.int, 
    control.group = list(model="rw2"))
# rw2 structure on the temporal component

asthma.typeII.inla <- inla(asthma.typeII , family="nbinomial",
                            data = asthma_pollution.dat, offset=log(population),
                            control.predictor =list(compute=TRUE),
                            control.compute = list(waic=TRUE))


# Fixed effects
round(asthma.typeII.inla$summary.fixed[,1:5],5)

# WAIC

asthma.typeII.inla$waic$waic



# TYPE IV RW

asthma.typeIV <-Counts ~ 1 + white + PM2.5 +
  f(num_id, model="bym", graph=asthma.adj) +
  f(Year, model="rw2") +
  f(Year1, model="iid") +
  f(num_id.int, model="besag", graph=asthma.adj,
    #spatial component
    group=Year.int, control.group = list(model="rw2"))
#temporal component

asthma.typeIV.inla <- inla(asthma.typeIV , family="nbinomial",
                            data = asthma_pollution.dat, offset=log(population),
                            control.predictor =list(compute=TRUE),
                            control.compute = list(waic=TRUE))



# Fixed effects
round(asthma.typeIV.inla$summary.fixed[,1:5],5)

# WAIC

asthma.typeIV.inla$waic$waic


# PLOTS for RW interactions


 #function to plot the four interaction plots
 int.plot.function <- function(asthma_pollution.dat, ast.type.mean, int=int){
 
 #create the data
 delta.int <- data.frame(delta=ast.type.mean,
                       num_id=rep(unique(asthma_pollution.dat$num_id), each=6),
                        Year=rep(unique(asthma_pollution.dat$Year), 9))
 
 delta.int <- left_join(asthma_pollution.dat, delta.int)
 
 target <- c(2011,2013,2015,2017,2019)
 delta.int <- subset(delta.int, Year %in% target)
 
 ggplot(delta.int, aes(fill=delta, geometry=geometry)) +
           geom_sf()+
           scale_fill_gradient(low = "darkblue", high ="white")+
           theme_void() +
           theme(legend.title = element_blank()) +
           ggtitle(int) +
           facet_wrap(~Year)
 }
 
 p1 <- int.plot.function(ast.type.mean =
                           ast.mnth.type.inla$summary.random$num_id.Year$mean,
                   ast.mnth=ast.mnth, int=1)
 
 p2 <- int.plot.function(ast.type.mean =
                           ast.mnth.typeII.inla$summary.random$num_id.int$mean,
                   ast.mnth=ast.mnth, int=2)
 
 p4 <- int.plot.function(ast.type.mean =
                           ast.mnth.typeIV.inla$summary.random$num_id.int$mean,
                   ast.mnth=ast.mnth, int=4)
 
 grid.arrange(p1,p2,p4, nrow=2)
 
 
 
 
 
 
 # TYPE I RW plot 
 # Subset asthma.dat to match the length of the mean vector
 asthma_pollution.dat$typeI.fit <- asthma.typeI.inla$summary.fitted.values$mean
 
 target <- c(2011:2019)
 asthma.fitted<- filter(asthma_pollution.dat, Year %in% target)  %>%
   group_by(num_id, Year) %>%
   summarise_at(c("typeI.fit"), sum)
 
 #we'll categorise for RW
 rw.cut <- with(asthma.fitted,quantile(typeI.fit,0:10/10))
 
 asthma.fitted <- asthma.fitted %>%
   mutate(cut.typeI.fit = cut(typeI.fit,rw.cut))
 
 p1<- ggplot(asthma.fitted, aes(fill=typeI.fit, geometry=geometry)) + 
   geom_sf()+
   scale_fill_viridis(option="cividis") +
   theme_void() +
   theme(legend.title = element_blank()) +
   facet_wrap(~Year)  
 p1
 
 
 #----------------INTERACTIONS LINEAR----------------------
 # TYPE I Linear model + plot
 
 #another Year
 asthma_pollution.dat$Year1 <- asthma_pollution.dat$Year
 
 #interaction term
 num_id.Year <- paste(asthma_pollution.dat$num_id,asthma_pollution.dat$Year,sep="_")
 
 
 #linear
 asthma.year.lin <- Counts ~ 1 + Year + PM2.5 + white +
   f(num_id, model="bym", graph = asthma.adj,constr=TRUE) +
   #spatial part
   f(Year1, num_id, model="iid", constr=TRUE) + 
 #the spatial/temporal interaction term with delta coefficients
 f(num_id.Year, model="iid")
 #interaction term
 
 lcs <- inla.make.lincombs(Year1=diag(9), Year=rep(1,9))
 
 asthma.inla.lin <- inla(asthma.year.lin,
                         family="nbinomial",
                         data=asthma_pollution.dat,
                         offset=log(population),
                         lincomb = lcs,
                         control.predictor = list(compute=TRUE),
                         control.compute = list(waic=T), control.inla = list(h=0.1))
 
 
 # Fixed effects
 round(asthma.inla.lin$summary.fixed[,1:5],5)
 
 # WAIC
 asthma.inla.lin$waic$waic
 
 # Subset asthma.dat to match the length of the mean vector
 asthma_pollution.dat$lin.fit <- asthma.inla.lin$summary.fitted.values$mean
 
 target <- c(2011:2019)
 asthma.fitted<- filter(asthma_pollution.dat, Year %in% target)  %>%
   group_by(num_id, Year) %>%
   summarise_at(c("lin.fit"), sum)
 
 #we'll categorise for RW
 rw.cut <- with(asthma.fitted,quantile(lin.fit,0:10/10))
 
 asthma.fitted <- asthma.fitted %>%
   mutate(cut.lin.fit = cut(lin.fit,rw.cut))
 
 p1<- ggplot(asthma.fitted, aes(fill=lin.fit, geometry=geometry)) + 
   geom_sf()+
   scale_fill_viridis(option="cividis") +
   theme_void() +
   theme(legend.title = element_blank()) +
   facet_wrap(~Year)  
 p1
 
 
 
 