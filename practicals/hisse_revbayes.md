---
title: HiSSE en Revbayes
layout: home
nav_order: 8
index: true
redirect: false
parent: Tutorials
math: katex
---
Creado por Rosana Zenil-Ferguson (Mayo 2025)
Actualizado Nicolás Castillo-Rodríguez (Febrero 2026)

Creado para el taller "Métodos Comparativos para Naturalistas" dictado en la Pontificia Universidad Católica del Ecuador (PUCE) en Quito, Ecuador durante los días 3, 4 y 5 de Febrero del 2026. Los datos provienen de: di Stilio & Zenil-Ferguson. All traits considered. 2025. *Sometido*. 

Para este tutorial, crearemos un modelo de especiación y extinción dependiente del estado con caracteres escondidos (HiSSE: Hidden State-Dependent Speciation and Extinction). Evaluaremos el efecto de la polinización establecida con dos estados: Polinización por insectos, codificada como 0 en nuestra tabla de datos, y polinización por viento, codificada como 1.

## Datos

 + Lista de especies y su polinizador [here](downloads/poliniza_datos.csv)
 
 + Árbol filogenético bifurcado y ultramétrico con 103 especies [here](downloads/poliniza_arbol.tre)
 
 + Script completo: [HiSSE Rev code](downloads/hisse.Rev)
 
## Código en RevBayes
 
1. En primer lugar creamos un vector ``moves`` que guarda las propuestas de movimientos que estableceremos para cada parámetro. Así mismo, el vector ``monitors`` guarda la inferencia obtenida del modelo, incluyendo la distribución posterior para cada uno de los parámetros.
 
``` 
# Propuesta de movimientos para parámetros y monitores de la inferencia del modelo.
moves = VectorMoves()
monitors = VectorMonitors()
```

2. En este ejemplo, tenemos dos estados observados (0 y 1) y definimos dos estados escondidos (A y B). Por lo tanto, contamos con cuatro estados (0A, 1A, 0B, 1B), cada uno con su propia tasa de diversificación, es decir, estableceremos cuatro tasas de especiación y cuatro tasas de extinción.

```
NUM_STATES <- 2
NUM_HIDDEN <- 2
NUM_RATES = NUM_STATES * NUM_HIDDEN
```
 
3. Utilizando las funciones ``readTrees`` y ``readCharacterDataDelimited``, cargamos en RevBayes los archivos con la filogenia y caracteres a utilizar.

``` 
### Árbol filogenético
observed_phylogeny <- readTrees("data/poliniza_arbol.tre")[1]

## Datos
## 0 = Polinización por insectos
## 1 = Polinización por viento
data <- readCharacterDataDelimited("data/poliniza_datos.csv",
                                   stateLabels=2,
                                   type="NaturalNumbers",
                                   delimiter=",",
                                   header=TRUE)
```

4. Con la función ``expandCharacters`` extendemos nuestros datos de 0 y 1, a 0A, 0B, 1A y 1B.

```
data_exp <- data.expandCharacters(NUM_HIDDEN)
```

5. Establecemos algunos estadísticos de nuestros datos y filogenia:

```
# Qué taxones?
taxa <- observed_phylogeny.taxa()

# Qué tan antigua es nuestra filogenia?
root_age <- observed_phylogeny.rootAge()
```

6. Definir las tasas de transición:

Para las tasas de transición en nuestro modelo, utilizaremos distribuciones Gamma para nuestros *a priori*:

```
shape_pr <- 0.5
rate_pr := observed_phylogeny.treeLength()/5

#Transiciones entre estados observados:
q_0A1A ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_1A0A ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_0B1B ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_1B0B ~ dnGamma(shape=shape_pr, rate=rate_pr)

#Transitions entre estados escondidos:
q_0A0B ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_0B0A ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_1A1B ~ dnGamma(shape=shape_pr, rate=rate_pr)
q_1B1A ~ dnGamma(shape=shape_pr, rate=rate_pr)

moves.append(mvScale(q_0A1A, weight=2 ))
moves.append(mvScale(q_1A0A, weight=2 ))
moves.append(mvScale(q_0B1B, weight=2 ))
moves.append(mvScale(q_1B0B, weight=2 ))

moves.append(mvScale(q_0A0B, weight=2 ))
moves.append(mvScale(q_0B0A, weight=2 ))
moves.append(mvScale(q_1A1B, weight=2 ))
moves.append(mvScale(q_1B1A, weight=2 ))

# Crea una matriz de 4x4 porque son 4 estados
# En la Q matriz, las posiciones corresponden a: 1=0A, 2=1A, 3=0B, y 4=1B

for (i in 1:NUM_RATES) {
for (j in 1:NUM_RATES) {
q[i][j] := 0.0
}
}

q[1][2]:= q_0A1A
q[1][3]:= q_0A0B
q[2][1]:= q_1A0A
q[2][4]:= q_1A1B
q[3][4]:= q_0B1B
q[3][1]:= q_0B0A
q[4][3]:= q_1B0B
q[4][2]:= q_1B1A
```

7. Establecemos la Q matriz mediante la función ``fnFreeK``:

```
# La Q matriz es una matriz infinitesimal, lo que significa que es la derivada de la matriz de probabilidades.

rate_matrix := fnFreeK(q, rescaled=false, matrixExponentialMethod="scalingAndSquaring")
```

8. Definir las tasas de diversificación

Con el fin de establecer las tasas de diversificación, usaremos un multiplicador llamado ``speciation_alpha`` para los estados cuyo valor escondido correpsonda a A, al que sumaremos un parámetro definido como ``speciation_beta`` para los estados cuyo estado escondido corresponda a B. Definimos los valores de diversificación de esta manera ya que teóricamente los estados A y B derivan del mismo estado observado.

```
# Priors para las tasas de diversificación
total_taxa<- observed_phylogeny.ntips() # How many taxa
H = 0.5
rate_mean <- ln(ln(total_taxa/2.0) /root_age) # Magallon and Sanderson (2001)
rate_sd <- 2*H

# Especiación y extinción definidas a partir de una distribución log normal; utilizamos speciation_alpha y speciation_beta para ayudarnos a definir esa distribución.

for (i in 1:NUM_STATES) {

speciation_alpha[i] ~ dnNormal(mean=rate_mean,sd=rate_sd)
moves.append(mvSlide(speciation_alpha[i],delta=0.20,tune=true,weight=2.0))

extinction_alpha[i] ~ dnNormal(mean=rate_mean,sd=rate_sd)
moves.append(mvSlide(extinction_alpha[i],delta=0.20,tune=true,weight=2.0))
}

for (i in 1:NUM_HIDDEN) {

speciation_beta[i] ~ dnNormal(0.0,1.0)
moves.append(mvScale(speciation_beta[i],lambda=0.20,tune=true,weight=2.0))

extinction_beta[i] ~ dnNormal(0.0,1.0)
moves.append(mvSlide(extinction_beta[i],delta=0.20,tune=true,weight=2.0))

}

for (j in 1:NUM_HIDDEN) {
for (i in 1:NUM_STATES) {
if ( j == 1) {
speciation[i] := exp( speciation_alpha[i] )
extinction[i] := exp( extinction_alpha[i] )
} else {
index = i+(j*NUM_STATES)-NUM_STATES
speciation[index] := exp( speciation_alpha[i] + speciation_beta[j-1] )
extinction[index] := exp( extinction_alpha[i] + extinction_beta[j-1] )
}
}
}

net_diversification := speciation - extinction
```

9. Estimación de caracter en la raíz del árbol

Ya que desconocemos si el ancestro común más reciente de todos los taxones de nuestra filogenia se polinizaba por viento o por insectos, necesitamos estimar el valor del caracter en la raíz del árbol. En el marco de la estadística bayesiana, las frecuencias en la raíz son dos parámetros adicionales que debemos estimar. Por lo tanto, asumiremos frecuencias iniciales iguales para cada estado. Dado que hay dos estados, necesitaremos una distribución *a priori* bivariada para un vector con dos frecuencias que sumen 1. Una distribución muy útil para este propósito es la distribución de Dirichlet, una distribución de probabilidad multivariada que nos permite asignar la misma frecuencia a ambos estados.
 
``` 
root_frequencies ~ dnDirichlet(rep(1,NUM_RATES))

# Dos propuestas de movimiento para explorar los valores del caracter en la raíz del árbol

moves.append(mvDirichletSimplex(root_frequencies,tune=true,weight=2))
moves.append(mvElementSwapSimplex(root_frequencies, weight=3))
```

10. Fracción de muestreo

Para modelos de diversificación necesitamos especificar la proporción de linajes del clado que fueron muestreadas en nuestra filogenia:

```
total_clade <- 200
extant_sampling <- total_taxa/total_clade
```

11. Construir modelo HiSSE

Con el fin de conectar todos los nodos de este modelo, utilizaremos una distribución de probabilidad para árboles filogenéticos.  Esta distribución de probabilidad es la que calcula internamente la función de verosimilitud a lo largo del árbol y considera la distribución *a priori* de los parámetros en toda la estructura del árbol.

```
#  La función dnCDBDP (Character-dependent birth-death process) establece una distribución de probabilidad para un modelo de Markov que es un proceso de nacimiento y muerte dependiente de caracter a lo largo del árbol filogenético.

hisse ~ dnCDBDP(rootAge         = root_age,
                speciationRates = speciation,
                extinctionRates = extinction,
                Q               = rate_matrix,
                pi              = root_frequencies,
                rho             = extant_sampling,
                delta           = 1,
                condition       = "time")
```

12. Cálculo de la función de verosimilitud

Durante la especificación del modelo hemos ignorado los datos, sin embargo, estos son necesarios para el cálculo de verosimilitud. Para ello, podemos utilizar la función ``clamp`` que permite que la función ``dnCDBDP()`` tenga en cuenta los valores en las puntas del árbol definido para cada taxón.

``` 
# Agregar datos y filogenia al modelo
hisse.clamp(observed_phylogeny)
hisse.clampCharData(data_exp) #usamos data set expandido
```

Una vez hayamos incluido nuestros datos al modelo, obtendremos una **distribución posterior para el modelo HiSSE** y podremos realizar nuestra inferencia.

## Guardando nuestra inferencia Bayesiana

A lo largo del código hemos especificado distintos ``moves`` que corresponden a las propuestas del MCMC para cada parámetro:

1. Breve descripción de las propuestas:

+ ``mvScale(q_0A1A, weight=2)`` esta propuesta reescala el valor original de $$q_{0A1A}$$ dos veces por iteración.
+ ``mvSlide(speciation_alpha[i],delta=0.20,tune=true,weight=2.0)`` esta propuesta define una ventana que se desliza a través de los valores del parámetro y ajusta el tamaño de la ventana a medida que avanza el MCMC. 
+ ``mvDirichletSimplex(root_frequencies,tune=true,weight=2)`` esta función propone dos valores en el intervalo (0,1) cuya suma es igual a 1. Estos representan la frecuencia con la que encontraríamos el estado 0 o el estado 1 en la raíz del árbol.
+ ``mvElementSwapSimplex(root_frequencies, weight=3)`` esta función propone intercambiar las frecuencias. Por ejemplo, si tenemos (0.4, 0.6), elementswap intercambia la posición de esos valores, resultando en (0.6, 0.4).

Cada una de estas propuestas contribuye a los *momios a posteriori*, si estos mejoran, la propuesta puede ser entonces aceptada.

2. Guardar objeto gráfico en RevBayes:

Este es un paso importante en el software RevBayes. Queremos “guardar” todo el objeto gráfico para poder manipularlo posteriormente. La función `model()` nos permite hacer esto.

``` 
# Mymodel is like a box that takes care of the whole graphical model
mymodel = model(rate_matrix)
```

3. Los monitores siguen el proceso inferencial:

Los monitores almacenan toda nuestra inferencia; sin ellos no será posible recuperar las distribuciones posteriores. Existen muchos tipos de monitores, como se muestra a continuación:

```
## Este monitor guarda la distribución posterior completa para cada uno de los parámetros
monitors.append(mnModel(filename="output/hisse_pollination.log", printgen=1))
##Este monitor registra lo que ocurre con la reconstrucción de estados ancestrales en cada nodo
#monitors.append(mnJointConditionalAncestralState(tree=hisse, cdbdp=hisse, type="NaturalNumbers", printgen=500, withTips=true, withStartStates=false, filename="output/asr_hisse_pollination.log"))
##Este monitor crea mapas estocásticos de caracteres (evolución a lo largo de las ramas); es muy lento para modelos complejos.
#monitors.append(mnStochasticCharacterMap(cdbdp=hisse, printgen=500, filename="output/stochmap_hisse_pollination.log", include_simmap=true))
## Este monitor imprime en pantalla algunos parámetros, per mite conocer el avance del análisis
monitors.append(mnScreen(printgen=10, q_0A1A, q_1A0A, net_diversification))
```

+ ``mnModel`` guarda las muestras de la distribución posterior generadas por el algoritmo MCMC.
+ ``mnScreen`` imprime información en pantalla para saber que el análisis está en ejecución.
+ ``mnJointConditionalAncestralState`` guarda la reconstrucción de estados ancestrales utilizando la probabilidad posterior marginal (esto es muy diferente a lo que hacen phytools u otros programas).
+ ``mnStochasticCharacterMap`` calcula mapas estocásticos de caracteres, es decir, las transiciones que ocurren a lo largo de las ramas. Esto es muy importante, pero también muy difícil de calcular.

4. Ejecutar el MCMC

**Recuerda siempre correr dos cadenas para demostrar que el MCMC ha convergido**

En lugar de establecer un número fijo de generaciones, establecemos que el análisis se ejecute hasta que todos los parámetros en el modelo tengan un tamaño de muestra efectivo de al menos 250. Así mismo, incluimos checkpoints con el fin de retomar nuestro análisis si este es interrumpido.

```
# Crear objeto MCMC
mymcmc = mcmc(mymodel, monitors, moves, nruns=2, moveschedule="random")

# Checkpoint and stopping rule
if ( fileExists("output/hisse_pollination.state") ) {
  mymcmc.initializeFromCheckpoint("output/hisse_pollination.state")
}

stopping_rules[1] = srMinESS(250, file = "output/hisse_pollination.log", freq = 10000)

# Ejecutar MCMC
mymcmc.run(rules = stopping_rules, checkpointInterval = 1000, checkpointFile = "output/hisse_pollination.state")
```

5. Create summaries of your ancestral state reconstruction 

 + Descarga uno de los archivos log de una reconstrucción ancestral [aquí](files/anc_states_hisse_pollination_run_1.log)
 
```
# Generando una distribución posterior marginal para cada nodo para calcular su reconstrucción de estados ancestrales

#anc_states = readAncestralStateTrace("output/asr_hisse_pollination_run_1.log")
#anc_tree = ancestralStateTree(tree=observed_phylogeny, ancestral_state_trace_vector=anc_states, include_start_states=false, file="output/asr_summary_hisse_pollination_run_1.tree", burnin=0.1, summary_statistic="MAP", site = 1)

q()
```

![](images/hissegraphical.png)
