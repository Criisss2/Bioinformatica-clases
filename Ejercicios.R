######################################################################################################################

##################         UNA SOLA MUESTRA                                 ###############################################
            
# Suponemos que estamos estudiando la expresión de 10 genes

x<-c(3,6,3,6,7,1,4,7,4,4)
x
class(x)
x[1]
x[7]
x[0]
x[3:7]

# Para poner secuencias:
3:7
seq(3,7,1)

# Poner posiciones de la 1,3 y de la 5 a la 8
c(1,3,5:8)
x[c(1,3,5:8)]

# Saber qué datos cumplen una condición
x >= 4
x[x>=4] #solo nos devuelve los que dan TRUE
x[x>6]
x[x>4&x<=6]
y<-x[x>4&x<=6]
y
x[x<=4|x>6]
x[x==4]

# En qué posición tenemos ciertos valores que cumplen una condición
which(x==4)

# Ordenar un vector
sort(x) # Orden creciente
sort(x,index.return=TRUE) # Para sabe rqué posición ocupaban esos valores en el vector original
sort(x,decreasing=TRUE,index.return=TRUE) # Orden decreciente
x  # ¿Nos guarda el orden? Para eso hay que almacenarlo en una variable
z <- sort(x)
z

######################            EJERCICIO 1

# Introduce los siguientes datos:
x<-c(9,5,6,5,2,8,3,2,4,9,8,6,6,11,6,6,5,5,4,3,9,8,6,7,1,5,6,3,4,3,3,4,4,4,3,3,3,5,4,7,2,5,7,2,5,3,3,6,8,4,6,2,4,4,7,5,7,7,7,4,6,6,6,0,5,3,4,5,4,2,3,6,4,7,2,4,10,8,8,3,3,6,4,6,4,3,7,2,5,5,5,4,5,9,5,8,8,5,4,3)
x
#1.- ¿Cual es el valor medio, la varianza y la desviación estándar de estos valores?
mean(x)
var(x)
sd(x)

#2.- ¿Cuántos de los valores son 3?
length(x[x==3])
y<-x[x==3]
length(y)

#3.- ¿Cuál es el rango intercuartílico? Diferencia entre Q3 - Q1 (aprox el 50% observaciones). A más grande más dispersión
# El primer cuartil (Q1) que es el valor que separa el 25% inferior de los datos del 75% superior.
# El segundo cuartil (Q2), también conocido como mediana, que es el valor que separa el 50% inferior del 50% superior de los datos.
# El tercer cuartil (Q3) que es el valor que separa el 75% inferior del 25% superior de los datos.
IQR(x)
x
# formula Lp=(n+1)*(p/100)
y<-sort(x)#primero ordenamos los datos
y
n<-length(y)
n
q1 <- y[(n+1)*(25/100)];q1
q3 <- y[(n+1)*(75/100)];q3
q3-q1
#otra forma
q1<-quantile(x,c(0.25),type=6);q1
q3<-quantile(x,c(0.75),type=6);q3
q3-q1
#la mediana es lo mismo que el q2
mediana<-median(x)
mediana
q2<-quantile(x,c(0.50),type=6);q2

#4.- Ordena el vector y calcula las sumas acumuladas de los valores ordenados con la función base::cumsum
y<-sort(x)
suma_acumulada<-cumsum(y)
suma_acumulada

#5.- Obtén el mínimo y el máximo
min(x)
max(x)
summary(x)

##########################    EJERCICIO 2
#Generar valores aleatorios con una distribución normal de media 45 y sd de 2.3. 23489 valores a obtener.
x<-rnorm(23489,mean=45,sd=2.3)
x

#1.- ¿Cuántas valores de x son mayores o iguales de 39.4?
val<-x[x>=39.4]
length(val)

#2.- ¿Cuántos valores de x son estrictamente menores que 46?
length(x[x<46])

#3.- ¿Cuántos valores de x son estrictamente menores que 46 y mayores o iguales de 39.4?
length(x[x<46&x>=39.4])

#4.- ¿Cuántos valores son tales que su parte entera (por defecto) es igual a 40? Se puede utilizar también la función base::floor
length(x[x==40]) #¿qué pasa aquí?
help(floor)
z<-floor(x)
z
length(z[z==40])

######################################################################################################################
############################                  VARIAS MUESTRAS                              ###########################

# Suponemos que tenemos la expresión de 10 genes en 4 muestras. Distribución de Poisson. Lambda es la tasa promedio
set.seed(123)
(x=matrix(rpois(10*4,lambda=5),nrow=10))
class(x)

# Podemos añadir nombres  las columnas y filas
colnames(x)=c("m1","m2","m3","m4")
# colnames(x)=paste0("m",1:4)
rownames(x)<-paste0("gene",1:10)
x

# Buscar en la matriz
x[6,2]
x["gene6","m2"]
x[,2]
x[,"m2"]

#################################              EJERCICIO 3
# Crear una matrix con valores de distribución normal, media de 34 y sd de 3.21.
x<-matrix(rnorm(2345*122,mean=34,sd=3.21),nrow=2345)
x

#1.- ¿Cuántas columnas tiene la matriz?
ncol(x)

#2.- ¿Qué columna tiene la media más alta? Utiliza las funciones base::apply y base::which.max
colMeans(x)
y<-apply(x,2,mean) #el 2 indica columnas
which.max(y)
which.max(colMeans(x))

#3.- ¿Qué fila tiene la media más pequeña? Utiliza las funciones base::apply y base::which.min
rowMeans(x)
which.min(rowMeans(x))
y<-apply(x,1,mean)
which.min(y)

#4.- ¿Qué fila tiene el menor rango intercuartílico? Utiliza las funciones base::apply y base::which.min
which.min(apply(x,1,IQR))

#5.- ¿En cuántas ocasiones la primera columna es mayor o igual que la segunda columna de la matriz
y<-x[,1]>=x[,2]
length(y)

######################################################################################################################
#################################         DATOS GOLUB                           ######################################
# Los datos golub son datos de expresión de genes que se emplearon en "Molecular classification of cancer: class 
#discovery and class prediction by gene expression monitoring. Science 1999.
# Son 3051 genes de 38 pacientes con leucemia. 27 de esos pacientes tienen una leucemia linfoblástica aguda (ALL) y el resto (11) leucemia mielode aguda (AML).
# El tipo de tumor está indicado en el vector golub.cl donde ALL es 0 y AML es 1. Los nombres de los genes son golub.gnames
# Las filas son los genes y las columnas las muestras

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("multtest")

library("multtest")
data(golub)
summary(golub)
head(golub)
class(golub)
golub.cl
class(golub.cl)
nrow(golub);ncol(golub);dim(golub) #dim es para las dimensiones

# Vamos a crear un factor para ALL y AML. Facilita el procesado de los datos
golub.fac <- factor(golub.cl,levels=0:1,labels=c("ALL","AML"))
golub.fac
table(golub.fac) #Para ver cuántas tenemos de cada uno
summary(golub.fac)

# Creamos nuestro primer diagrama de barras
library(ggplot2) #en caso de no haberla cargado antes. Viene con el paquete pacman
primer_df <- data.frame(golub.fac)
ggplot(primer_df,aes(x=golub.fac))+geom_bar()
dev.off()
# ¿Qué pasa si sustituimos golub.fac por golub.cl?
primer_df <- data.frame(golub.cl)
ggplot(primer_df,aes(x=golub.cl))+geom_bar()
dev.off()

# Ver expresión de un gen en una muestra y en todas. Representarlo en un gráfico.
golub[2000,12]
golub[2000,]
#calculamos los ejes x e y
muestra<-1:ncol(golub)
y2000<-golub[2000,]
df<-data.frame(muestra,y2000)
ggplot(df,aes(x=muestra,y=y2000))+geom_point()
dev.off()

# Ver expresiones de los pacientes ALL y AML para el gen 2000
golub[2000,golub.fac=="ALL"]
golub[2000,golub.fac=="AML"]
df <- data.frame(muestra,y2000,tipo=golub.fac)
ggplot(df,aes(x=muestra,y=y2000,color=tipo))+geom_point()
dev.off()

# Otra opción es representar las muestras AML y ALL por separado
df<-data.frame(muestra,y2000,tipo=golub.fac)
ggplot(df,aes(x=muestra,y=y2000,color=tipo))+geom_point()+xlab("Número de muestra")+ylab("Gen en fila 2000")+facet_grid(tipo~.)
dev.off()

# Otra opción es colocar todo horizontal
df<-data.frame(muestra,y2000,tipo=golub.fac)
ggplot(df,aes(x=muestra,y=y2000,color=tipo))+geom_point()+xlab("Número de muestra")+ylab("Gen en fila 2000")+facet_grid(.~tipo)
dev.off()

# ¿Qué gen corresponde la fila 2000?
golub.gnames[2000]
class(golub.gnames)
dim(golub.gnames)
golub.gnames[2000,]
golub.gnames[2000,2]
#identificador de affymetrix
golub.gnames[2000,3]
rownames(golub)=golub.gnames[,3] #nombramos las filas con el identificador de Affy
colnames(golub)=golub.fac #de columna ponemos el tipo de leucemia
summary(golub)
head(golub)

# Ahora vamos a hacer otras gráficas. En eje x medidas de localización de ALL, en y lo mismo pero del AML
mean(golub[2000,])
median(golub[2000,])
# vemos que hay mucha asimetría porque no coinciden los valores
golub.ALL <- golub[,golub.fac == "ALL"]
golub.AML <- golub[,golub.fac == "AML"]
ALL.mean <- apply(golub.ALL,1,mean)
AML.mean <- apply(golub.AML,1,mean)
ALL.median <- apply(golub.ALL,1,median)
AML.median <- apply(golub.AML,1,median)
df<-data.frame(x0=ALL.mean,y0=AML.mean)
p=ggplot(df,aes(x=x0,y=y0))+geom_point()
p+xlab("Medias ALL")+ylab("Medias AML")
dev.off()
# vamos a poner un sombreado que ayude a distinguir dónde hay más densidad de puntos
df<-data.frame(x0=ALL.mean,y0=AML.mean)
ggplot(df,aes(x=x0,y=y0))+geom_point(alpha=.25)+xlab("Medias ALL")+ylab("Medias AML")
dev.off()
#ahora representamos la mediana
df<-data.frame(x0=ALL.median,y0=AML.median)
ggplot(df,aes(x=x0,y=y0))+geom_point(alpha=.25)+xlab("Medianas ALL")+ylab("Medianas AML")
dev.off()

################                     EJERCICIO 4
#1.- Calcular para cada gen la expresión media sin considerar el tipo de cancer
apply(golub,1,mean)

#2.- Calcular para cada gen la desviación estándar
apply(golub,1,sd)

#3.- Representar gráficamente un dibujo que, para cada gen, nos muestre en abscisas la expresión media y en ordenadas la desviación estándar
expresion_media <- apply(golub,1,mean)
desv_stdr <- apply(golub,1,sd)
df<-data.frame(x=expresion_media,y=desv_stdr)
ggplot(df,aes(x=expresion_media,y=desv_stdr))+geom_point(alpha=.25)
dev.off()

#4.- Utilizando la función grDevices::png guardar el dibujo anterior en un fichero externo
getwd()
#"C:/Users/irene/Documents/UAX/Asignaturas/Biomedicina/Bioinformatica/2023-2024/R/Golub.png"
png("C:/Users/irene/Documents/UAX/Asignaturas/Biomedicina/Bioinformatica/2023-2024/R/Golub.png")
ggplot(df,aes(x=expresion_media,y=desv_stdr))+geom_point(alpha=.25)
dev.off()

#5.- Determinar la desviación estándar de las expresiones de cada gen
desvSD<-apply(golub,1,sd)
desvSD

#6.- Representar un histograma de estas desviaciones estándar
plot(hist(desvSD))

#3.- Calcular el percentil de orden 0.9 (Q0.9) de las desviaciones calculadas en el punto anterior. 
desvSD.sort<-sort(desvSD)
n<-length(desvSD.sort)
n
q09<-desvSD.sort[(n+1)*(90/100)]
q09

#4.- Utilizando la función abline añadir al dibujo del apartado 2 una línea vertical cuya abscisa coincida con el percentil 0.9
abline(v=q09,col="red")
dev.off()
#5.- Seleccionar aquellos genes cuya desviación estándar sea mayor que el valor Q0.9
genes.sdAlta<-desvSD[desvSD>q09]
genes.sdAlta