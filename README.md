# ROOT_HEP
Este reporsitorio consiste de una serie de ejercicios realizados en el desarrollo del curso de fìsica experimental de partículas. El repositorio contiene los códigos construidos en Root con el propósito de analizar los datos consignados en algunos experimentos llevados a cabo en el CERN. 

Ejercicio 1:
Se analizaron los datos del paquete Tracks_Clusters.root que contiene información del tracker y el calorímetro electromagnético.  Los resultados son:
1. El número promedio de interacciones por bunch-crossing es cercano a 39
2. Como es de esperarse, la relación entre el número promedio de interacciones y el número de vertices primarios es de proporcionalidad directa con una tendencia lineal.
3. Al comparar el número de vertices primarios y el número promedio de interacciones con el número de tracks y clusters se encuentra que la relación tambien es de proporcionalidad directa, sin embargo, cabe mencionar que mientras en el caso de los clusters se tiene una tendencia aparentemente lineal, en el caso de los tracks la tendencia no es tan evidente. 
4. Al graficar todas las variables de los tracks y clusters se encuentra que en el momento transversal existe un corte en la detecciòn el cual esta relacionado con la sensibilidad del detector, también se encuentra que la seudorapidez tiene una distribución que sugiere que hacia el centro del detector se tiene una mayor sensibilidad mientras que para ángulos cercanos a 45 grados, la detección decae drasticamente. Por último, en el caso de Phi se encuentra como es de esperarse un comportamiento periódico.
5. El PDG ID es un identificador dado a las partículas que se detectan. En el caso de los datos encontrados en el paquete se obtuvo un PDG ID de 21.93 que podrìa corresponder a un gluon o a un foton.

Ejercicio 2:
Usando el mismo paquete de datos que se usó en el ejercicio 1 se analizaron las diferentes variables asociadas con los Jets. Dentro de los datos registrados en el archivo Tracks_Clusters.root, se encontraban simulaciones hechas con Monte Carlo, las cuales permiten tener un marco de referencia para contrastar la data recopilada por los detectores. Dentro del análisis realizado se graficaron las diferentes variables cinemáticas de los jets, también se verificó el pile up existente y se utilizaron criterios basados en el JVF para reducir la influencia del pile up y mejorar las curvas resultantes. Cabe mencionar que como es de esperarse, los datos generados por el montecarlo no se muestran influenciados por el pile up, esto se debe a que al tratarse de datos simulados no existe ruido del entorno que afecten las señales.

Ejercicio 3 y 4:
En este ejercicio se se utilizó un dataset llamado Data_8TeV.root para analizar un proceso de interacción entre un par top-antitop. En esta selección se optó por el análisis del canal en el que se producen dos b-jets junto con dos bosones W, mientras que a su vez uno de los bosones decae adrónicamente y el otro leptónicamente. El propósito de este ejercicio era implementar adecuadamente los cortes necesarios sobre el dataset para seleccionar solo los eventos que coincidad con el estado final deseado. Con este objetivo se impuso las siguientes condiciones:
1. One good vertex -> Tener un buen vertice, el cual se define a aprtir de la información del branch hasGoodVertex.
2. Lepton trigger -> con los branches trigE y trigM se establece la existencia de electrones o muones en el estado final.
3. Exactly one good lepton -> Dado que en el estado final solo hay un leptón proveniente de la desintegración de uno de los bosones W, es necesario que este leptón este adecuadamente aislado y no se confunda con algun otro subproducto, para esto se imponen las condiciones Pt>25 GeV (ser suficientemente energético), lep_ptcone30/lep_pt<0.15 (aislar su recorrido), lep_etcone20/lep_pt<0.15 (aislar su presencia en el calorímetro) y encontrarse en la región del tracker que es una condición impuesta sobre el valor de seudo rapidez y que varia en función de si se trata de un muón o un electrón.
4. At least 4 good jets -> Esta condición se establece teniendo en cuenta que el estado final cuenta con dos b-jets y la desintegración hadrónica de uno de los bosones W. En este caso se exige pt > 25 GeV, |eta| < 2.5, e incluso se establece una condición para implementar el corte de JVF para reducir el pile up.
5. At least 2 b-jets -> esto surge de lo que se mencionó acerca del estado final y se traduce en la condición MV1 >= 0.7892
6. También es necesario identificar la energía perdida MET, que es energía asociada a objetos que no pudieron ser medidos en el detecto. Como en esta selección el estado final contiene una desintegración leptónica del W entonces habrá presencia de muones que no serán detectados. Para este propósito se impone la condición MET >30 GeV
7. Para finalizar, la masa invariante debe ser consecuente con la presencia de los bosones W, por lo tanto se exige que la masa transversa sea Mt > 30 GeV

Aparte de realizar el llenado de los respectivos histogramas de las diferentes variables cinématicas, se utilizaron los datos obtenidos a partir del Monte Carlo, lo cual permite establecer si la data es congruente o no con el proceso que se estima que es.

Proyecto Final
Para el proyecto final se escogió la producción WZ en colisiones de protones para el dataset de 8 TeV. Pese a que tanto para el W como para el Z la probabilidad de desintegración hadrónica esta por encima del 70%, se utilizó el canal leptónico para el análisis debido a que se trata de una señal más limpia y con un nivel de dificultad menor en su tratamiento, al menos para el caso de 8 TeV. En esta reacción se tendrá un estado final de tres leptones y una perdida de energía debido a la presencia de neitrinos. La selección se hace con los criterios:
1. One good vertex
2. El trigger de electrones y muones
3. Un buen lepton
4. Exactamente tres leptones
5. Exigir que al menos dos de los leptones sean del mismo tipo y cargas opuestas (esta condicion surge debido a la desintegración del Z en un par lepton-antilepton)
6. la masa invariante del par debe ser cercana a la masa nominal del boson Z
7. MET>30 GeV

Adicionalmente, se utilizó como background ttbar, ZZ, WW y single top que son procesos cuyo estado final puede ser similar al que deseamos. Esto se hace con el objetivo de distinguir cuantos eventos en la región de la señal pertenecen realmente al proceso que se esta estudiando.
