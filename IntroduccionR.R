######################     PARA EJECUTAR PULSAMOS CTRL+R (WINDOWS) O CMD+R (MAC)                                    ########

# Dónde estamos
getwd()

# Poner bien el directorio de trabajo. Por consola o cambiar dir desde el menú. 
# El directorio de trabajo tiene que ser el mismo que dónde tenemos el script
setwd()

#################     INSTALAR Y CARGAR PAQUETES (no siempre necesario, ESTO ES UN EJEMPLO)                    ########
install.packages("UsingR") 
install.packages(c("devtools","pacman"))
library(UsingR)
library(devtools)
library(pacman)
library(ggplot2)

##################    MANUAL
help.start()
help(hist)

######################################################################################################################
###############################                R COMO CALCULADORA
1+1
2+3*4
3^2
exp(1)
sqrt(10)
pi
2*pi*6378 #circunferencia de la Tierra

######################################################################################################################
###############################               VARIABLES
# Numeric: contienen números con decimales
# Booleanos: TRUE or FALSE
# Strings: secuencias de caracteres
x<-1
y<-3
z<-4
x*y*z
x*y*Z #Hay distinción entre mayúsculas y minúsculas!
This.Year<-2023 #pueden incluir puntos, guiones bajos pero nunca "-". Y nunca deben empezar por números

#######################################################################################################################
###############################               VECTORES
# Series de números
# Diferentes maneras de crearlos y operar con ellos
rep(1,10)
seq(2,6)
seq(4,20,by=4)
x<-c(2,0,0,4)
y<-c(1,9,9,9)
x+y
x*4
sqrt(x)

# Navegar con ellos
x[1]
x[-1]
x[1]<-3;x
x[-1]<-5;x
y<9
y[4]=1
y<9
y[y<9]=2
y

###########################################################################################################################
############################                    DATA FRAMES
# Son agrupaciones de vectores
# Es de lo más fácil de manejar en R
bp<-read.table("bp.txt",header=T)
bp
bp["WEIGHT"]
bp[,2]
bp$WEIGHT
bp$WEIGHT[2]

bp[,-2] #aquí no mostramos parte de los datos
###########################################################################################################################
###############################             	LISTAS
# Colecciones de variables
# Se parecen mucho a las data frames

# Tenemos unos genes asociados a tres grupos
geneset1<-c(1,4)
geneset2<-c(2,5,7,8)
geneset3<-c(1,5,6)

# Almacenamos los datos en una lista
genesets<-vector("list",3)
genesets[[1]]<-geneset1
genesets[[2]]<-geneset2
genesets[[3]]<-geneset3
genesets

# Ahora le ponemos nombre a los grupos
names(genesets)=paste0("set",1:3)
genesets

# Ahora nombramos a los genes
genenames=paste0("gene",1:10)
genesets[[1]]=genenames[geneset1]
genesets[[2]]=genenames[geneset2]
genesets[[3]]=genenames[geneset3]
genesets

# Extraer elementos de la lista
genesets[1]
genesets[[1]]
genesets[[1]][1]
unlist(genesets[1])[1]
length(genesets)
length(genesets[[2]])

# Crear una lista vacía
lista_vacia<-list()
# Ahora con tres elementos
lista_vacia<-vector("list",length=3)
# Ahora la llenamos con un bucle
for(j in 1:3){
lista_vacia[[j]] <- c(1,2,3*j)
}
lista_vacia

# Añadir elementos. Para eso el data frame es lo mismo que la lista
lista_vacia[[4]] <- data.frame(x=c(8,5,3),y=c(7,9,1))
lista_vacia

# Eliminar elementos
lista_vacia[[2]]<-NULL
#lista_vacia[-2] #equivalente
lista_vacia

######################################################################################################################
#############################           IF ELSE IFELSE
if(4 > 5) {
  "Verdadero"
} else {
  "Falso"
}

num <- 1:8

ifelse(num %% 2 == 0, "Par", "Non")

######################################################################################################################
#############################            LOOPS

# For
dado <- 1:6
for(cara in dado) {
  print(cara ^ 2)
}

mi_vector<-NULL
for(cara in dado) {
  mi_vector[cara] <- cara ^ 2
}
mi_vector

# While
umbral <- 5
valor <- 0

while(valor < umbral) {
  print("Todavía no.")
  valor <- valor + 1
}

while(1 < 2) {
  print("Presiona ESC para detener")
}

######################################################################################################################
##########################             BREAK Y NEXT
for(i in 1:10) {
  if(i == 3) {
    break
  }
  print(i)
}

numero <- 20
while(numero > 5) {
  if(numero == 15) {
    break
  }
  numero <- numero - 1
}
numero

for(i in 1:4) {
  if(i == 3) {
    next
  }
  print(i)
}


#######################################################################################################################
##################################       APPLY y  LAPPLY

# APPLY apply(x,margin,fun) #en Margin, 1 es filas y 2 columnas. Fun es la funcion que como requisito siempre tiene que admitir funciones

mi_df <- data.frame(v1 = 1:3, v2 = 4:6)
mi_df
# Coerción a matriz
mi_matriz <- as.matrix(mi_df)

# Verificamos que sea matriz
is.matrix(mi_matriz)
class(mi_matriz)
mi_matriz

apply(mi_matriz,1,mean)
apply(mi_matriz,2,sum)

# LAPPLY vale para listas lapply(x,fun)
trees[1:5,] #esta por defecto en R
lapply(trees,mean)
arboles <- lapply(X = trees, FUN = mean)
class(arboles)

# Podemos usar lapply en lugar de un bucle
mi_vector <- 6:12
resultado <- NULL
posicion <- 1
for(numero in mi_vector) {
  resultado[posicion] <- sqrt(numero)
  posicion <- posicion + 1
}
resultado

resultado <- NULL
resultado <- lapply(mi_vector, sqrt)
resultado
as.numeric(resultado)

########################################################################################################################
####################                     IMPORTAR Y EXPORTAR

# Descargar datos https://archive.ics.uci.edu/ml/datasets/Breast+Cancer+Wisconsin+%28Diagnostic%29 (diagnosticos de cancer para IA)
download.file(
  url = "https://raw.githubusercontent.com/jboscomendoza/r-principiantes-bookdown/master/datos/breast-cancer-wis.data", 
  dest = "breast-cancer-wis.data"
)
bcancer <- read.table(file = "breast-cancer-wis.data")
head(bcancer)
bcancer <- read.table(file = "breast-cancer-wis.data", header = FALSE, sep = ",")
head(bcancer)
nombres <- c("id", "clump_t", "u_csize", "u_cshape", "m_adh", "spcs", "b_nuc","b_chr", "n_nuc", "mit", "class")
bcancer <- read.table(file = "breast-cancer-wis.data", header = FALSE, sep = ",",col.names = nombres)
head(bcancer)
class(bcancer)

# Exportar datos
write.table(x = iris, file = "iris.txt", sep = ",", row.names = FALSE, col.names = TRUE)
iris_txt <- read.table(file = "iris.txt", header = TRUE, sep = ",")
head(iris_txt)
write.csv(x = iris, file = "iris.csv", row.names = FALSE) 
iris_csv <- read.csv("iris.csv")
head(iris_csv)

#######################################################################################################################
##############################           GRAFICOS
#### Para los colores que valen http://www.stat.columbia.edu/~tzheng/files/Rcolor.pdf

download.file(
  url = "https://raw.githubusercontent.com/jboscomendoza/r-principiantes-bookdown/master/datos/bank.csv", 
  destfile = "bank.csv"
  )
readLines("bank.csv", n = 4) 
banco <- read.csv(file = "bank.csv", sep = ";")
head(banco)
dim(banco)
lapply(banco,class)
summary(banco)

# Histograma
hist(x=banco$age)
hist(x = banco$age, main = "Histograma de Edad", 
     xlab = "Edad", ylab = "Frecuencia")
hist(x = banco$age, main = "Histograma de Edad", 
     xlab = "Edad", ylab = "Frecuencia",
     col = "purple")
hist(x = banco$duration, main = "Histograma de Duration", 
     xlab = "Duration", ylab = "Frecuencia",
     col = "ivory")

# Plot 
plot(x = banco$education)
plot(x = banco$education, main = "Gráfica de Educacíón",
     xlab = "Nivel educativo", ylab = "Frecuencia", 
     col = c("royalblue", "seagreen", "purple", "grey"))

# Barplot
table(banco$education)
tab_banco <- table(banco$loan, banco$education)
tab_banco
barplot(tab_banco)
prop.table(tab_banco, margin = 1)
prop.table(tab_banco, margin = 2)
prop.table(tab_banco)
ptab_banco <- prop.table(tab_banco, margin = 2)
barplot(ptab_banco)
barplot(ptab_banco,  main = "Préstamos por nivel educativo",
     xlab = "Nivel educativo", ylab = "Proporción", 
     col = c("royalblue", "grey"))
unique(banco$loan)
barplot(ptab_banco,  main = "Préstamos por nivel educativo",
     xlab = "Nivel educativo", ylab = "Proporción", 
     col = c("royalblue", "grey"))
legend(x = "topright", legend = c("No", "Yes"), fill = c("royalblue", "grey"), 
       title = "Loan")

# Diagrama de dispersion
plot(x = iris$Petal.Length, y = iris$Petal.Width, col = iris$Species, 
     main = "Iris - Pétalo", xlab = "Largo", ylab = "Ancho")
legend(x = "topleft", legend = c("Setosa", "Versicolor", "Virginica"), 
       fill = c("black", "red", "green"), title = "Especie")

# Boxplot
#formula: Para esta función las fórmulas tienen el 
#formato y ~ x, donde x es el nombre de la variable continua a graficar, y la x es la variable que usaremos como agrupación.
boxplot(x = banco$age)
boxplot(formula = age ~ education, data =  banco)




