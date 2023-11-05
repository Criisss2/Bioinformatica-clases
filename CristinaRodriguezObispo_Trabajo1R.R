#Primero creamos un nuevo script
#Cambiamos dir a la carpeta donde tenemos nuestros datos y el script guardado
#Comprobamos que estamos en la carpeta correcta
getwd()

#Cargamos los datos proporcionados en R
datos<-read.table("datos-trabajoR.txt", sep="\t", header=TRUE) #usamos la función read.table para que pueda leer nuestro archivo .txt, le indicamos como estan separados y le decimos que hay encabezados para que no los tome como datos

#Una vez cargados los datos empleamos las funciones:
head(datos) #esta función nos sirve para mostrar las primeras n filas del conjunto de datos.
summary(datos)#con esta función se muestra un resumen sobre las variables del archivo (mínimo, máximo, media, mediana, primer y tercer cuartil), además puede reconocer las variables categóricas, por lo que muestra la frecuencia de cada categoría.
dim(datos) #nos informa sobre las dimensiones de nuestros datos, salvo apraa vectores que debemos usar la función length().
str(datos)#esta función muestra estructuras de objetos, se utiliza principalmente para mostrar el contenido de una lista. 

#¿Cuántas variables hay? #Hay 3 variables
#¿Cuántos tratamientos? #Hay 50 tratamientos
#Ambos datos los obtenemos con la función dim() que nos proporciona las dimensiones.

#Hacemos un boxplot para nuestros datos. Uno para cada variable. Coloreando la Variable 1 y la Variable 2 de forma diferente.
boxplot(Variable1~Tratamiento, data=datos, col="light blue",
	main="Boxplot Variable 1")
boxplot(Variable2~Tratamiento, data=datos, col="light pink",
	main="Boxplot Variable 2")
#Comparamos cada variable en relación con el tratamiento (la variable en el eje y y el tratamiento en el x), después indicamos de donde se obtienen los datos para hacer el boxplot,el color con la función col y por último el título de la gráfica.

#Hacemos un gráfico de dispersión con las dos variables, cada una en un color
plot(x = datos$Variable1, y = datos$Variable2, col = datos$Tratamiento, 
     main = "Gráfico de dispersión V1-V2", xlab = "Variable 1", ylab = "Variable 2")
#en esta función primero indicamos que datos van en la x y en la y e indicamos de donde salen estos datos, luego indicamos que queremos que nos ponga colores distintos en base al tratamiento, y por último le indicamos los nombres de los ejes y el título del gráfico.

#Ahora ponemos una leyenda en el gráfico de dispersión para indicar el significado de cada color
legend(x = "bottomright", legend = c("Tto.1", "Tto.2", "Tto.3", "Tto.4", "Tto.5"), 
       fill = c("black", "red", "green", "blue", "aquamarine"), title = "Tratamientos")
#con la x indicamos la posición en la que queremos la leyenda, depués le indicamos que queremos que ponga en la leyenda, el color de cada elemento de la leyenda y por último el título de la leyenda.

#Realizamos un histograma para cada variable
hist(x = datos$Variable1, main = "Histograma de V1", 
     xlab = "Variable 1", ylab = "Frecuencia",
     col = "light blue")
hist(x = datos$Variable2, main = "Histograma de V2", 
     xlab = "Variable 2", ylab = "Frecuencia",
     col = "light pink")
#primero indicamos que datos van en el eje x y de donde salen esos datos y después le decimos que nombre poner como título del gráfico y los nombres de los ejes.

#Hacemos un factor en la columna tratamiento y lo guardamos en una variable
tratamiento_factor <- factor(datos$Tratamiento, levels=c("1","2","3","4","5"))
tratamiento_factor
#los facores se utilizan para trabajar con variables categóricas, variables que tienen un conjunto fijo y conocido de valores posibles.
#en esta función indicamos primero que variable usamos y la base de datos de donde las obtenemos y luego el número de niveles posibles,

#Calcula la media y la desviación estándar para cada tratamiento
tapply(datos$Variable1,tratamiento_factor,mean)
aggregate(Variable1~tratamiento_factor,datos,mean)

tapply(datos$Variable2,tratamiento_factor,mean)
aggregate(Variable2~tratamiento_factor,datos,mean)

#estas son dos formas para calcular la media, primero la de la variable 1 y luego la de la variable 2
#en ambas funciones usamos el facotr que hemos calculado antes.
#en ambas funciones debemos indicar la variable que queremos, el factor que hemos calculado, la base de datos de donde lo sacamos y la función que queremos realizar (media o desviación estándar).

tapply(datos$Variable1,tratamiento_factor,sd)
aggregate(Variable1~tratamiento_factor,datos,sd)

tapply(datos$Variable2,tratamiento_factor,sd)
aggregate(Variable2~tratamiento_factor,datos,sd)

#estas son dos formas para calcular la desviación típica, primero la de la variable 1 y luego la de la avriable 2

#Averiguamos cuántos elementos tiene cada tratamiento
table(tratamiento_factor) #uso el factor que he creado antes para sacar el nº de elementos más rapidamente

#Extraigo los datos para el tratamiento 1 y el tratamiento 4 y guardo cada uno en una variable diferente
Tto_1 <- subset (datos, Tratamiento ==1)
Tto_4 <- subset (datos, Tratamiento ==4)
#susbset es una función que se utiliza para obtener las filas y columnas.
#nombramos una nueva variable con la información sobre los tratamientos 1 y 4 que extraemos de nuestra tabla de datos con la función subset.
Tto_1 
Tto_4
#llamamos a cada variable apra comprobar que los datos se han guardado correctamente.

#Nuestra hipótesis nula es que las medias de tratamiento 1 y tratamiento 4 para la Variable 1 son iguales.
#Para comprobarlo primero comprobamos si los datos se distribuyen de forma normal
#Usamos la función saphiro.test() para saber si la distribución de una muestra es normal o no. Además nuestra muestra entra dentro de el rango que es de entre 3 y 5000 observaciones

Distribución_normal_Tto_1 <- shapiro.test(Tto_1$Variable1)
Distribución_normal_Tto_1

Distribución_normal_Tto_4 <- shapiro.test(Tto_4$Variable1)
Distribución_normal_Tto_4

#guardamos ambos test dentro de una nueva variable apra poder ver los resultados.
#Ambos tienen un pvalue con valores mayores de 0.05, por lo que no son valores significativamente estadisticos y por lo tanto no hay suficiente evidencia para rechazar nuestra hipótesis nula que es que los datos siguen una distribución normal. En resumen, nuestros datos siguen una distribución normal.
#En función del resultado de la prueba de normalidad, ¿qué test usarías?
#Usaría T-test, ya que es el test empleado para comparar medias de dos muestras que tienen una distribución normal

t_test <- t.test(Tto_1$Variable1, Tto_4$Variable1)
t_test

#al igual que antes guardamos los reultados del test en una nueva variable.
#Con los resultados del t test sabemos: que debido a que el p value es tan pequeño se sugiere que hay evidencia estadística de que las medias son diferentes. 
#Por otro lado la hipótesis alternativa establece que la verdadera diferencia entre las medias no es igual a 0, lo que sugiere que las medias son diferentes. 
#Y por último nos da las medias de x e y (4.0 y 50.8 respectivamente) con lo que vemos que efectivamente las medias son distintas.
#Con esto rechazamos nuestra hipótesis nula, las medias son distintas, no iguales.

#Asumimos que las muestras son independientes, pero ¿son sus varianzas iguales?
#Para comparar las varianzas usamos var.test()
varianzas_comp <- var.test (Tto_1$Variable1, Tto_4$Variable1)
varianzas_comp

#El resultado indica que hay evidencia estadística para afirmar que las varianzas de las dos variables son significativamente diferentes ya que el p value es muy pequeño. 
#Además, la estimación de la razón de varianzas es baja, lo que confirma que las varianzas son distintas.




