
# Packages ----------------------------------------------------------------

require(tidyverse)
require(magrittr)
require(cluster)
require(vegan)
require(lubridate)
require(factoextra)

# Loading data --------------------------------------------------------------------

load("Data/CapDolph_SightMatrix_2017_2020.RData")

# Manipulating data -------------------------------------------------------

# Computing the Year-Month and Year variables
data_turs <- raw_data_turs %>% mutate(YM=floor_date(Data, "month"), Y=year(Data))

# Extracting the single occasions, their order, and the days since the first one
Occs <- data_turs %>% distinct(Data) %>% mutate(Occs=1:nrow(.), Times=cumsum(c(0, diff(Data))))


# Computing the metrics ---------------------------------------------------

# Metrics based on the original occasions: Relative Span-Time, Periodicity, Resight Rate, Sight Rate
ClusVar <- data_turs %>% dplyr::select(-`Mark level`, -`Age class`, -Gender) %>% 
  left_join(Occs, by="Data") %>% mutate(maxOcc=max(Occs), maxTime=max(Times)) %>% 
  group_by(ID) %>% summarise(firstOcc = min(Occs[Avvistato==1]),
                             firstTime = min(Times[Avvistato==1]),
                             lastTime = max(Times[Avvistato==1]),
                             spanTime = (lastTime-firstTime),
                             spanTimeRel = spanTime/maxTime[1],
                             Periodicity = 1/mean(diff(Times[Avvistato==1])),
                             Resight=(sum(Avvistato)-1)/min(maxOcc-firstOcc), 
                             Sights=mean(Avvistato)) %>% 
  mutate(Periodicity = ifelse(is.nan(Periodicity), 0, Periodicity)) %>% 
  select(-firstOcc, -firstTime, -lastTime, -spanTime)

# Appending the Monthly rate (computed using the occasions grouped by month)
ClusVar$MRate <- data_turs %>% dplyr::select(-`Mark level`, -`Age class`, -Gender) %>% 
  left_join(Occs, by="Data") %>%
  group_by(ID, YM) %>% summarise(MSight=as.numeric(sum(Avvistato)>0)) %>% 
  group_by(ID) %>% summarise(MRate=mean(MSight)) %$% MRate

# Appending the Yearly rate (computed using the occasions grouped by year)
ClusVar$YRate <- data_turs %>% dplyr::select(-`Mark level`, -`Age class`, -Gender) %>% 
  left_join(Occs, by="Data") %>%
  group_by(ID, Y) %>% summarise(YSight=as.numeric(sum(Avvistato)>0)) %>% 
  group_by(ID) %>% summarise(YRate=mean(YSight)) %$% YRate


# Clustering implementation --------------------------------------------------------------

# Computing the distance metrics on the standardized numerical features
diss <- vegdist(ClusVar %>% column_to_rownames("ID")  %>% 
                   mutate(across(everything(), scale)), method = "gower")


# Selecting the optimal number of clusters K through WSS and Silhouette
fviz_nbclust(x = ClusVar[,-1], FUNcluster = hcut, diss = diss, method = "wss", k.max = 10, 
             hc_method = "ward.D2")
fviz_nbclust(x = ClusVar[,-1], FUNcluster = hcut, diss = diss, method = "silhouette", k.max = 10, 
             hc_method = "ward.D2")
# Clustering with the optimal K=3 
res <- hcut(x = diss, isdiss = T, k = 3, hc_method = "ward.D2")


# Clustering results ------------------------------------------------------

# Visualizing the results
fviz_dend(res, rect = TRUE, cex = 0.3,
          k_colors = c("firebrick", "darkorange", "darkgreen"))

# Extracting the group labels
clustcut <- res$cluster

# Joining the original data with the group labels
data_turs_clustered <- data_turs %>% left_join(tibble(ID=names(clustcut), Group=clustcut))

# Cluster groups' characterization ------------------------------------------------

# Total individual sightings colored by by group
data_turs_clustered %>% distinct(ID, TOT, Group) %>% 
  ggplot(aes(x=TOT, y = ..density.., fill=factor(Group), color = I("black"))) + 
  geom_histogram(position = "identity", alpha = 0.5) + labs(fill="Group") +
  theme_bw()

# Joining the features with the group labels
ClusVar_clustered <- ClusVar %>% left_join(tibble(ID=names(clustcut), Group=clustcut))

# Clusters' sizes
ClusVar_clustered %>% group_by(Group) %>% summarise(n=n())
# Clusters' average values
ClusVar_clustered %>% group_by(Group) %>% summarise_if(is.numeric, mean)

# Clusters interpretation on the "Principal Components Analysis" (PCA) space 
ClusVar_PCA <- ClusVar %>% select(-ID) %>% scale() %>% princomp(.)
PCA <- ClusVar_PCA$scores[, 1:2]
# Interpretation of the PCA components
ClusVar_PCA$loadings
# Visualization of the groups in the PCA space
data.frame(PCA=matrix(PCA, ncol=2), Group=clustcut) %>% ggplot() + 
  geom_point(aes(x=PCA.1, y=PCA.2, color=as.factor(Group))) + labs(color="Group") +
  theme_light() 


# Saving the output (for use in the POPAN) --------------------------------


save(data_turs_clustered, file="Data/SightMatrix_Clustered.RData")
