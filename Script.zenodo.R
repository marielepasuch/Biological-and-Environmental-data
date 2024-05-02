#ESCAPES OF Oreochromis niloticus FROM AQUACULTURE PONDS CAN ALTER THE ICHTHYOFAUNA STRUCTURE IN NEOTROPICAL STREAMS#

#packages
library(vegan)
library(BiodiversityR)
library(ggplot2)
library(ggpubr)

#Species-rank abundace - Zippin density estimates#
ambiental <- read.table("ambiental.txt", header = TRUE)
ausente <- read.table("CPUE_semZippin_ausente.txt", header = T)
aus <- ausente[,6:75]
aus1 <- log10(aus+1)

moderado <-read.table("CPUE_semZippin_moderada.txt", header = TRUE)
mod <- moderado[,6:75]
mod1 <- log10(mod+1)

intensa <- read.table("CPUE_semZippin_intensa.txt", header = TRUE)
int <- intensa[,6:75]
int1 <- log10(int+1)

ausente1 <-rankabundance(aus1, ambiental)
moderado1 <- rankabundance(mod1,ambiental)
intensa1 <- rankabundance(int1, ambiental)

rankabunplot(ausente1, scale = "logabun", scaledx = F, specnames = c(), 
             pch = 19, col = "black")

rankabunplot(moderado1, scale = "logabun", specnames = c(), pch = 19, 
             xlim = c(0,70), addit = TRUE, col = "gray")

rankabunplot(intensa1, scale = "logabun", specnames = c(), pch = 19, 
             xlim = c(0,70), addit = TRUE, col = "red")

legend(50, 20, legend = c("No propagule", "Moderate", "Intense"),
       col = c("black", "gray", "red"), lty = 1, cex = 0.8, box.lty = 0)

# PCoA and PERMANOVA - Presence and absence data, abundance data and environmental data
peixes <- read.table("CPUE_comZippin.txt", header = TRUE)
View(peixes)

#hellinger
abu <- peixes[ ,6:75]
peixes.hell <- decostand(abu, method = "hell")

#presence and absence data
peixes.pa<-ifelse(abu>0,1,0)
View(peixes.pa)

#distance matrix - Bray Curtis
peixes.bray<- vegdist(peixes.hell, method = "bray") 

#distance matrix - Jaccard
peixes.jaccard <- vegdist(peixes.pa, method = "jaccard") 

#PERMANOVA - Presence and absence data
permanovaPA <- adonis2(peixes.jaccard ~ peixes$Aquaculture)

adonis2(peixes.jaccard ~ peixes$Aquaculture * peixes$Month,
        permutations = 999)
adonis2(peixes.jaccard ~ peixes$Aquaculture)

pairwise.adonis <- function(x, factors, sim.function = 'vegdist', sim.method = 'jaccard')
{
  library(vegan)
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value)
  print("Signif. codes:  0 ????~***????T 0.001 ????~**????T 0.01 ????~*????T 0.05 ????~.????T 0.1 ????~ ????T 1")
  return(pairw.res)}

resu.permanova.PA <-pairwise.adonis(x= peixes.pa,#
                                    factors = peixes$Aquaculture, 
                                    sim.method= "jaccard") 
resu.permanova.PA

#PERMANOVA - abundance data
permanova.abu <- adonis2(peixes.bray ~ peixes$Aquaculture)

adonis2(peixes.bray ~ peixes$Aquaculture * peixes$Month,
        permutations = 999)
adonis2(peixes.bray ~ peixes$Aquaculture)

pairwise.adonis <- function(x, factors, sim.function = 'vegdist', sim.method = 'bray')
{
  library(vegan)
  co = combn(unique(as.character(factors)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    if(sim.function == 'daisy'){
      library(cluster); x1 = daisy(x[factors %in% c(co[1,elem],co[2,elem]),],metric=sim.method)
    } else{x1 = vegdist(x[factors %in% c(co[1,elem],co[2,elem]),],method=sim.method)}
    
    ad = adonis(x1 ~ factors[factors %in% c(co[1,elem],co[2,elem])] );
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
  }
  
  pairw.res = data.frame(pairs,F.Model,R2,p.value)
  print("Signif. codes:  0 ????~***????T 0.001 ????~**????T 0.01 ????~*????T 0.05 ????~.????T 0.1 ????~ ????T 1")
  return(pairw.res)}

resu.permanova.abu <-pairwise.adonis(x= peixes.hell,
                                     factors = peixes$Aquaculture, 
                                     sim.method= 'bray') 
resu.permanova.abu

#environmental data
ambiental <-read.table("ambiental.txt", header= T)
ambiental[,-c(1:4)]->abi2
decostand(abi2, method = "standardize")->abi.stand2
vegdist(abi.stand2, method = "euclidean", na.rm = TRUE)->abi.dist2

permanova.abi <- adonis2(abi.dist2 ~ ambiental$Categoria)

# ordinations - PCoA
#presence and absence
pcoa.pa <- cmdscale(peixes.jaccard, eig = T, k=2)
summary(pcoa.pa)
pcoa.pa$points
pcoa.pa$x
pcoa.pa$ac
pcoa.pa$GOF 

ordiplot(scores(pcoa.pa)[,c(1,2)], type = "t", main = "PCoA with sites")
abline(h=0, lty = 3)
abline(v = 0, lty = 3)
spe.wa<- wascores(pcoa.pa$points, peixes.pa)
site.scores <- as.data.frame(scores(pcoa.pa)[,c(1,2)])
data.scores <- cbind.data.frame(site.scores, peixes$Stream, peixes$Aquaculture) 
colnames(data.scores) <- c("PCoA1", "PCoA2", "Stream", "Aquaculture")
head(data.scores) 
species.scores <- as.data.frame(scores(pcoa.pa, "species"))  
species.scores$species <- rownames(species.scores)
head(species.scores)

plot.pa <- ggplot() + 
  geom_point(data=data.scores,aes(x= PCoA1,y= PCoA2, shape = Aquaculture, colour = Aquaculture), size = 4) + 
  scale_colour_manual(values = c("red", "gray", "black")) +
  geom_text(data=data.scores,aes(x=PCoA1,y=PCoA2,label=""), size=3, vjust= 1, hjust= 0) +  
  theme_bw() + 
  theme(legend.position = "top")+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 16))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+ 
  theme(axis.ticks = element_blank(),  
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())

#abundance
pcoa.densidade <- cmdscale(peixes.bray, eig = T, k=2)
summary(pcoa.densidade)
pcoa.densidade$points
pcoa.densidade$x
pcoa.densidade$ac
pcoa.densidade$GOF

ordiplot(scores(pcoa.densidade)[,c(1,2)], type = "t", main = "PCoA with sites")
abline(h=0, lty = 3)
abline(v = 0, lty = 3)
spe.wa<- wascores(pcoa.densidade$points, peixes.hell)
site.scores <- as.data.frame(scores(pcoa.densidade)[,c(1,2)]) 
data.scores <- cbind.data.frame(site.scores, peixes$Stream, peixes$Aquaculture) 
colnames(data.scores) <- c("PCoA1", "PCoA2", "Stream", "Aquaculture")
head(data.scores)
species.scores <- as.data.frame(scores(pcoa.densidade, "species"))
species.scores$species <- rownames(species.scores)  
head(species.scores)

plot.densidade <- ggplot() + 
  geom_point(data=data.scores,aes(x= PCoA1,y= PCoA2, shape = Aquaculture, colour = Aquaculture), size = 4) + 
  scale_colour_manual(values = c("red", "gray", "black")) +
  geom_text(data=data.scores,aes(x=PCoA1,y=PCoA2,label=""), size=3, vjust= 1, hjust= 0) +  
  theme_bw() + 
  theme(legend.position = "top")+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 16))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+ 
  theme(axis.ticks = element_blank(), 
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())

plot_ab <- ggarrange(plot.pa,
                     plot.densidade + rremove("ylab") + rremove("y.text"),#para remover o t?tulo e a informa??o textual do eixo Y
                     ncol=2, nrow=1,  #informar n?mero de linhas e colunas
                     labels= c("A", "B"), hjust = -2) #hjust para alterar a posi??o de A e B
plot_ab 

#environmental data 
pcoa.ambiental <- cmdscale(abi.dist2, eig = T, k=2)
summary(pcoa.ambiental)
pcoa.ambiental$points
pcoa.ambiental$x
pcoa.ambiental$ac
pcoa.ambiental$GOF

ordiplot(scores(pcoa.ambiental)[,c(1,2)], type = "t", main = "PCoA with sites")
abline(h=0, lty = 3)
abline(v = 0, lty = 3)
spe.wa<- wascores(pcoa.ambiental$points, abi.dist2)
site.scores <- as.data.frame(scores(pcoa.ambiental)[,c(1,2)]) 
data.scores <- cbind.data.frame(site.scores, ambiental$Riacho, ambiental$Categoria) 
colnames(data.scores) <- c("PCoA1", "PCoA2", "Stream", "Aquaculture")
head(data.scores)
species.scores <- as.data.frame(scores(pcoa.ambiental, "species"))
species.scores$species <- rownames(species.scores)  
head(species.scores)

plot.ambiental <- ggplot() + 
  geom_point(data=data.scores,aes(x= PCoA1,y= PCoA2, shape = Aquaculture, colour = Aquaculture), size = 4) + 
  scale_colour_manual(values = c("red", "gray", "black")) +
  geom_text(data=data.scores,aes(x=PCoA1,y=PCoA2,label=""), size=3, vjust= 1, hjust= 0) +  
  theme_bw() + 
  theme(legend.position = "top")+ 
  theme(legend.title = element_blank(),
        legend.text = element_text(color = "black", size = 16))+
  geom_hline(yintercept=0, linetype="dashed", color = "grey")+ 
  geom_vline(xintercept=0, linetype="dashed", color = "grey")+ 
  theme(axis.ticks = element_blank(), 
        axis.title.x = element_text(size=12), 
        axis.title.y = element_text(size=12), 
        panel.background = element_blank(), 
        panel.grid.major = element_blank(),  
        panel.grid.minor = element_blank(),  
        plot.background = element_blank())

#Influence of spatial autocorrelation on fish assemblage

read.table("Coordenadas.txt", header = TRUE)->latlong
latlong
read.table("Mantel.txt", header = TRUE)->peixes

decostand(peixes, "hellinger")->peixeshell
vegdist(latlong, "euclid")->geo.dist
vegdist(peixes, "bray")->peixes.dist
mantel(peixes.dist, geo.dist, permutations = 5000) 

