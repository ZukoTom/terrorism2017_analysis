#PRISE EN MAIN DONNEES

df <- read.csv("/home/zuko/Desktop/Maths/df_terr_small.csv",
               header = TRUE,
               stringsAsFactors = TRUE)
#on supprime colonne de l'ancien index
df <- df[,-1] 
table(df$region_txt)
quali <- which(sapply(df, is.factor))
quanti <- which(sapply(df, is.numeric))
#on supp nom de pays car redondance et lourdeur
quali <- quali[-1] #ou bien la grader en illustrative?
#comptage des modalités semble montrer des modalités fréquentes dans donnees
#barplot donne idée distribution de var aléatoires
mapply(df[,quanti],
       FUN = function(xx,name){barplot(table(xx),main = name)},
       name = names(quanti))
mapply(df[,quali],
       FUN = function(xx,name){barplot(table(xx),main = name)},
       name = names(quali))

library(stargazer)
stargazer(df[,quanti],
          summary.stat = c("n","min","p25","median","mean","p75","max","sd"),
          type = "text")

matcor.pears <-cor(df[,quanti])
#le coefficient de spearman peut être mieux adapté à données (dist assymétriques)
#mais résultas proches de ttes façons
matcor.spear <- cor(df[,quanti],method = "spearman")
#visualisation : une seule corrélation ressort!
library(corrplot)
corrplot(matcor.spear, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)
#round(matcor.spear, 2)

#ANALYSE VAR QUANTI : PCA
library(FactoMineR)
res.pca <- PCA(df[,quanti],
               graph = FALSE)
barplot(res.pca$eig[,1], las = 2)#montre que les premières PC capturent pas mal d'inertie
plot.PCA(res.pca, choix = "ind",repel=TRUE)
plot.PCA(res.pca, choix = "var")
#ANALYSE VAR QUALI : ACM
res.mca <- MCA(df[,quali],
               graph=FALSE)
eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = TRUE, ylim = c(0, 45))

var <- get_mca_var(res.mca)
head(var$coord)
# Cos2: quality on the factore map
head(var$cos2)
# Contributions to the principal components
head(var$contrib)
#graphe correlation variables/ PC
fviz_mca_var(res.mca, choice = "mca.cor", 
             repel = TRUE, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())#individus
#graphe modalités (illisible)
fviz_mca_var(res.mca, # Avoid text overlapping (slow)
             ggtheme = theme_minimal())
# Color by cos2 values: quality on the factor map
fviz_mca_var(res.mca, col.var = "cos2",
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Avoid text overlapping
             ggtheme = theme_minimal())
#etude individus illisible
ind <- get_mca_ind(res.mca)

#biplot optionnel illisible
fviz_mca_biplot(res.mca, # Avoid text overlapping (slow if many point)
                ggtheme = theme_minimal())

fviz_mca_ind(res.mca, col.ind = "cos2", 
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Avoid text overlapping (slow if many points)
             ggtheme = theme_minimal())
#ANALYSE GLOBALE : AFDM
mateta2 <- matrix(NA,length(quali),length(quanti))
rownames(mateta2) <- names(quali)
colnames(mateta2) <- names(quanti)

library(BioStatR)
for(ii in seq(nrow(mateta2))){
  for(jj in seq(ncol(mateta2))){
    mateta2[ii, jj]<-eta2(df[, colnames(mateta2)[jj]],
                          df[, rownames(mateta2)[ii]])
  }
}
corrplot(mateta2,
         args.colorlegend = list(labels = sprintf("%.1f", seq(0, 1, length = 11))),
         tl.col = "black", tl.srt = 45)
         #border = NA,
         #cols = colorRampPalette(c("white", "steelblue"), space = "rgb")(20),
         #breaks = seq(0, 1, length=21),
         #cex.lab = par("cex.lab"), cex = 0.55*par("cex"),
         #args.colorlegend = list(labels = sprintf("%.1f", seq(0, 1, length = 11)), frame = TRUE))
res.eta2<-sort(mateta2[,"nwound"])
#représentation
barplot(res.eta2,horiz = TRUE,
        las = 2,
        xlab = expression(eta^2),
        main = "link var nwound and var qual",
        cex.names  =.35)

mapply(as.data.frame(mateta2[,colnames(df)[quanti]]),
       FUN=function(xx,name){
         names(xx) = rownames(mateta2)
         res.eta2 <- sort(xx)
         barplot(res.eta2,
                 horiz = TRUE,
                 las = 2,
                 xlab = expression(eta^2),
                 main = name,
                 xlim = c(0,1))
       },
       name = colnames(df)[quanti])

res.famd <- FAMD(df,
                 ncp = Inf,
                 graph = FALSE,
                 sup.var = quanti)
round(res.famd$eig, 3)

ncp <- 55
D <- dist(res.famd$ind$coord[,1:ncp])#distance euclidienne entre observations
res.hclust  <-  hclust(D,method = "ward.D2")

barplot(sort(res.hclust$height,decreasing = TRUE)[1:15],
        names.arg = 1:15,
        xlab = "index",
        ylab = "hauteur de fusion")

#kmeans : on definit le nombre de classe sur la base de la var de l'inertie inter de CAH
nbclasse <- 8 #nb de groupes
partition <-  cutree(res.hclust, k = nbclasse) #élagage de l'arbre

#Consolidation
centres.gravite <- by(res.famd$ind$coord[,1:ncp],
                      INDICES = partition,
                      FUN = colMeans) 

centres.gravite <- do.call(rbind, centres.gravite)#donne un objet de type "matrix", nécessaire pour pouvoir utiliser ces centres comme des valeurs initiales pour la fonction kmeans

res.kmeans <- kmeans(res.famd$ind$coord[,1:ncp],
                     centers = centres.gravite)

part.finale <- as.factor(res.kmeans$cluster)
#cette table compte effectifs dans chaque groupe : les individus extrèmes semblent altérer le clustering
table(part.finale)

df_part <- cbind.data.frame(df, classe = part.finale)#on concatène le jeu de données avec la nouvelle variable classe
catdes(df_part, num = ncol(df_part))
# on voit que les couleurs matchent pas les groupes... a cause modalités rares?
res.famd <- FAMD(df_part,
                 ncp = Inf,
                 graph = FALSE,
                 sup.var =  c(ncol(df_part),quanti)) 
#write.infile(res.famd, file="resultats_afdm.csv")
library(missMDA)
res.ncp <- estim_ncpFAMD(df_part,# on retire les variables illustratves qui ne sont pas gérées par la fonction (et inutile pour déterminer le nombre d'axes
                         ncp.max = 10,
                         method.cv = "Kfold",
                         nbsim = 50
)
plot(x = as.numeric(names(res.ncp$crit)),
     y = res.ncp$crit,
     xlab = "S",
     ylab = "Erreur",
     main = "Erreur de validation croisée\n en fonction du nombre d'axes",
     type = "b")
library(factoextra)
fviz_mfa_ind(res.famd, 
             habillage = "classe")

plot(res.famd,
     choix = "quali",
     invisible = c("quali","ind")
)

plot(res.famd, choix = "ind", invisible = "quali", select = "contrib 30")

plot(res.famd,choix = "var", select = "contrib 15")

plot(res.famd, choix = "quanti", select = "coord 5", cex=.7)#ne fonctionne pas
plot(res.famd, choix = "quali", cex = .7)

res.dimdesc <- dimdesc(res.famd)
lapply(res.dimdesc$Dim.1, round, 3)

write.infile(res.dimdesc,file = "dimdesc.csv")
