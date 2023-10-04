#############################################################################
#
# PRACTICA 1
#
# Expresión diferencial de genes de ratón
# Microarray de Affymetrix (Affymetrix Murine Genome U74A version 2 MG_U74Av2
# Origen de los datos: GEO GSE5583 (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5583)
# Publicación: Mol Cell Biol 2006 Nov;26(21):7913-28.  16940178 (http://www.ncbi.nlm.nih.gov/pubmed/16940178)
#
# Muestras: 3 Wild Type x 3 Histone deacetylase 1 (HDAC1)
#
# R código original (credits): Ahmed Moustafa
#
## ENTREGA EL 22 OCTUBRE 23:59
## Se requiere la entrega de este script
## Adjuntar en la entrega el PDF final y el archivo con los genes
#
##############################################################################

# Instalar RCurl

if (!requireNamespace("BiocManager"))
    install.packages("BiocManager")
BiocManager::install("RCurl")

# Si esto falla, que seguro lo hace tratar de instalarlo usando el menú, Paquetes, Servidor Spain A Coruña, RCurl

# Cargamos el paquete y los datos
library(RCurl)
url = getURL ("http://bit.ly/GSE5583_data", followlocation = TRUE)
data = as.matrix(read.table (text = url, row.names = 1, header = T))

# Chequeamos las dimensiones de los datos, y vemos las primeras y las últimas filas
dim(data) #esto mide las dimensiones
head(data) #esto nos muestra las primeras filas
tail(data) #muestra últimas filas

# Hacemos un primer histograma para explorar los datos
hist(data) #histograma guardado como Primer histograma
hist(data, col = "light blue", main="GSE5583 - Histogram")#cambiamos el color con esto
 
# Transformamos los datos con un logaritmo 
data_log=log2(data) #transformar los datos y guardar los datos en esa variable

# ¿Qué pasa si hacemos una transformación logarítima de los datos? ¿Para qué sirve?
#Lo usamos para tener una distribución de los datos más normal y no tan sesgada

# Hacemos un boxplot con los datos transformados. ¿Qué significan los parámetros que hemos empleado?
boxplot(data_log) #genera el boxplot
#Los tres priemros son wildtype y los tres últimos wild type
boxplot(data_log, col=c("blue", "blue", "blue", "green", "green", "green"), main="GSE5583-boxplots",las=2) 
#col=c sirve para acceder a cada elemento, despues ponemos entre paréntesis el color que queremos que sea cada boxplot
#las=2 significa que te pone el nombre del eje x en vertical
#main sirve para poner el título de la grafica

# ¿Qué es un boxplot?
#Un diagrama de bigotes en el que vemos los cuartiles y la mediana

# Hacemos un hierarchical clustering de las muestras basándonos en un coeficiente de correlación
hc=hclust(as.dist(1-cor(data_log)))
plot(hc, main="GSE5583 - Hierarchical Clustering")

# de los valores de expresión. ¿Es correcta la separación?
#Si

#######################################
# Análisis de Expresión Diferencial 
#######################################

# Primero separamos las dos condiciones. ¿Qué tipo de datos has generado?
wt<-data[,1:3]
ko<-data[,4:6]
class(wt)
#hemos generado una tabla

# Calcula las medias de las muestras para cada condición. Usa apply
wt.mean=apply(wt, 1, mean) #apply sirve para calcular lo que yo quiera
#wt.mean es la variable donde vamos a guardar los datos
#nos aparece cada gen con su media
ko.mean=apply(ko, 1,mean)

# ¿Cuál es la media más alta?
max(wt.mean)
max(ko.mean)

# Ahora hacemos un scatter plot (gráfico de dispersión)
plot(ko.mean ~ wt.mean, xlab="WT", ylab="KO", main = "GSE5583 - Scatter")
plot(ko.mean ~ wt.mean)

# Añadir una línea diagonal con abline
abline(0,1,col="red") #ponemos 0,1 para indicar que es una recta diagonal, y = 1x + 0 (y=mx+n, donde m es la pendiente y n la ordenada en el origen)

# ¿Eres capaz de añadirle un grid?
grid() #para añadir una cuadrícula (tanto para esto como para añadir la linea diagonal hay que tener el plot abierto)

# Calculamos la diferencia entre las medias de las condiciones
diff.mean = wt.mean - ko.mean #ya tenemos calculadas las medias de cada variable, y ahora nombramos a otra variable con la diferencia de ambas
#obtenemos una tabla porque teníamos las medias de cada gen, no un único número

# Hacemos un histograma de las diferencias de medias
hist(diff.mean)
hist(diff.mean, col="light pink")

# Calculamos la significancia estadística con un t-test.
# Primero crea una lista vacía para guardar los p-values
# Segundo crea una lista vacía para guardar las estadísticas del test.

#tenemos dos matrices y por cada fila tenemos que hacer un T test
#la i significa iterar
#para cada fila desde la 1 hasta el final de filas de data coge todas las columnas de wt y todas las de ko (las 3 columnas de wt se llman x y las 3 de ko se llma y) y haz un t test con ellos. Una vez hecho el T test guardamos la información en una variable, los p values lo guardamos en un lugar de mi lista pvalue y las estadisticas en un lugar de mi lista tstat.

pvalue=NULL
tstat=NULL
for(i in 1:nrow(data)){#Para cada gen 
	x=wt[i,] #gene wt número i 
	y=ko[i,] #gene ko número i
	
	#hacemos el test
	t=t.test(x,y)
	
	#añadimos el p-value a la lista
	pvalue[i]=t$p.value
	#Añadimos las estadísticas a la lista
	tstat[i]=t$statistic
}

head(pvalue)
length(pvalue)

# OJO que aquí usamos los datos SIN TRANSFORMAR. ¿Por qué? #Porque si lo haces con un logaritmo dejan de ser datos fiables, los log solo los usamos para las gráficas.
# ¿Cuántas valores tiene cada condición? #Tenemos 6 valores por cada gen (3 WT y 3 KO), por lo que tenemos 6 muestras con 2 condiciones cada una (el WT y el KO) y 3 valores cada condición


# Ahora comprobamos que hemos hecho TODOS los cálculos


# Hacemos un histograma de los p-values.
hist(pvalue) #valores muy grandes hay que transformarlo
hist(-log10(pvalue), col="blue")
# ¿Qué pasa si le ponemos con una transformación de -log10? #nos permite visualizar mejor la gráfica, cambia la distribución de los datos, dejandonos los p values significativos en la cola, el log negativo invierte la distribución


# Hacemos un volcano plot. Aquí podemos meter la diferencia de medias y la significancia estadística
plot(diff.mean, -log10(pvalue), main="Volcano plot", col="pink") #los significativos al estar con los valores del -log10 están en la parte de arriba del volvano plot

# Queremos establecer que el mínimo para considerar una diferencia significativa, es con una diferencia de 2 y un p-value de 0.01
# ¿Puedes representarlo en el gráfico?
diff.mean_cutoff=2
pvalue_cutoff=0.01 #más restrictiva nuestra diferencia significativa
abline(v=diff.mean_cutoff, col="purple", lwd=3) #v nos hace una línea vertical y lwd es el ancho de línea
abline(h=-log10(pvalue_cutoff), col="light green", lwd=3) #h nos hace una línea horizontal, y usamos -log10 para convertir los datos, ya que como hemos hecho el volcano plot con los logaritmos no nos reconocería los datos si no pusiesemos -log10


# Ahora buscamos los genes que satisfagan estos criterios
# Primero hacemos el filtro para la diferencia de medias (fold)
filter_by_diff.mean=abs(diff.mean) >=diff.mean_cutoff #abs es valor absoluto, es decir vamos a guardar valores por encima de 2 y por debajo de -2, por eso usamos abs, para convertir la tabla a valores absolutos y sacar directamente todos los valores por encima de 2
dim(data[filter_by_diff.mean, ])

# Ahora el filtro de p-value
filter_by_pvalue=pvalue <=pvalue_cutoff #guardamos todos los valores que son menores o iguales al filtro de pvalue que habíamos puesto
dim(data[filter_by_pvalue, ]) #nos da el nº de genes que pasan el filtro
#primero creamos una veriable para el primer filtro, luego otra variable para el segundo y luego combinamos ambas

# Ahora las combinamos. ¿Cuántos genes cumplen los dos criterios? #Todos los que cumplen el filtro de pvalue, es decir 426
filter_combined=filter_by_diff.mean&filter_by_pvalue
filtered=data[filter_combined,] #nos hemos quedado con los genes comunes a los dos filtros
dim(filtered)
head(filtered)

# Ahora generamos otro volcano plot con los genes seleccionados marcados en rojo
plot(diff.mean, -log10(pvalue), main = "Volcano plot II")
points(diff.mean[filter_combined], -log10(pvalue[filter_combined]),col="red")

# Ahora vamos a marcar los que estarían sobreexpresados (rojo) y reprimidos (azul). ¿Por qué parece que están al revés? 
#Porque primero hemos calculado los WT menos los KO en la diferencia de medias y los sobreexpresados son los KO que tienen 2 veces que el WT y los reprimidos son los que tienen la mitad de WT, al ser los valores de WT mayores que los KO tenemos valores negativos, y los reprimidos aparecen en el lado positivo y los sobreexpresados en el lado negativo

plot(diff.mean, -log10(pvalue), main = "Volcano III")
points(diff.mean[filter_combined & diff.mean <0],
	-log10(pvalue[filter_combined & diff.mean <0]), col ="red")
points(diff.mean[filter_combined & diff.mean >0],
	-log10(pvalue[filter_combined & diff.mean >0]), col ="blue")

# Ahora vamos a generar un mapa. Para ello primero tenemos que hacer un cluster de las columnas y los genes 
heatmap(filtered)
rowv=as.dendrogram(hclust(as.dist(1-cor(t(filtered)))))
colv=as.dendrogram(hclust(as.dist(1-cor(filtered))))
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,labRow=FALSE) #hemos hecho una agrupación, los WT por un lado y los KO por otro lado, y nos divide cada uno por la mitad entre reprimidos (amarillos) y sobreexpresados (rojos)

# ¿Qué es cada parámetro que hemos usado dentro de la función heatmap?
#cexCol es el tamaño de la letra del eje X, Colv y Rowv son los dendogramas, el primer término es lo que vamos a representar, labRow=FALSE, es para quitar el nombre de los genes
# ¿Eres capaz de cambiar los colores del heatmap? Pista: usar el argumento col y hcl.colors
#si
heatmap(filtered, Rowv=rowv, Colv=colv, cexCol=0.7, col=hcl.colors(50))

# Ahora vamos a crear un heatmap más chulo. Para ello necesitamos dos paquetes: gplots y RcolorBrewer
#if (!requireNamespace("BiocManager"))
#    install.packages("BiocManager")
#BiocManager::install(c("gplots","RColorBrewer"))
install.packages("gplots")		
install.packages("RColorBrewer")	

library(gplots)

#la línea discontinua indica el 0. Los azules son los reprimidos y los rojos los sobreexpresados.
# Hacemos nuestro heatmap
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col=rev(redblue(256)),scale="row")

# Lo guardamos en un archivo PDF
pdf ("GSE5583_DE_Heatmap.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = rev(redblue(256)), scale = "row", labRow=FALSE)
dev.off()
pdf ("GSE5583_DE_Heatmap Green.pdf")
heatmap.2(filtered, Rowv=rowv, Colv=colv, cexCol=0.7,
	col = redgreen(75), scale ="row",labRow=FALSE)
dev.off()



# Guardamos los genes diferencialmente expresados y filtrados en un fichero
write.table(filtered, "GSE5583_DE.txt", sep="\t",
	quote=FALSE)


