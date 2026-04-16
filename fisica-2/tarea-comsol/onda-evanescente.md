El **acoplamiento por ondas evanescentes** es el fenómeno físico que permite que la luz (o cualquier onda electromagnética) "salte" de una guía de onda a otra sin que haya un contacto eléctrico o físico directo entre sus núcleos.

Para entenderlo, hay que romper la idea de que la luz viaja solo *dentro* de la fibra o guía.

---

### 1. ¿Qué es una onda evanescente?
Cuando la luz viaja por una guía de onda (como los rectángulos que ves en tu modelo de COMSOL), ocurre algo llamado **Reflexión Interna Total**. Sin embargo, la física nos dice que el campo eléctrico no puede desaparecer instantáneamente en la frontera del material.

* Una pequeña parte de la energía penetra una fracción de micra hacia el material de alrededor (el revestimiento o *cladding*).
* Esta energía "sobrante" se llama **onda evanescente**.
* **Dato clave:** Su intensidad decae exponencialmente. A medida que te alejas de la guía, la energía desaparece casi por completo en una distancia muy corta.



---

### 2. El mecanismo de "salto" (Acoplamiento)
Imagina que tienes dos guías de onda paralelas (Guía A y Guía B) separadas por una distancia muy pequeña ($d$):

1.  **Guía A activa:** Envías luz por la Guía A. Su onda evanescente "flota" justo por fuera de su núcleo.
2.  **Solapamiento:** Como la Guía B está muy cerca, el campo evanescente de la Guía A alcanza a "tocar" el núcleo de la Guía B.
3.  **Transferencia:** La Guía B comienza a absorber esa energía y a propagarla. Por conservación de la energía, la potencia en la Guía A empieza a disminuir mientras que en la B aumenta.



---

### 3. Factores que determinan el acoplamiento
En tu simulación de COMSOL, el éxito del acoplamiento depende de tres variables principales que puedes modificar:

* **La Distancia (Gap):** Si separas las guías un par de micrómetros más, el acoplamiento cae a cero porque las ondas evanescentes ya no se alcanzan.
* **La Longitud de Interacción ($L$):** La energía oscila entre las guías. Existe una longitud específica donde el 100% de la luz pasa de A a B. Si la guía es el doble de larga, la luz volverá a la Guía A.
* **Sincronismo de fase:** Para que el salto sea eficiente, ambas guías deben ser (preferiblemente) idénticas en tamaño e índice de refracción.

---

### 4. ¿Por qué es importante en tu modelo?
En el **Directional Coupler** que tienes en pantalla, este fenómeno es el que define si el dispositivo funcionará como:
* **Un divisor de potencia (50/50):** La luz sale mitad por cada guía.
* **Un interruptor óptico:** Donde controlas si la luz sigue recto o cambia de carril.

¿Te gustaría que veamos cómo visualizar específicamente el decaimiento de esta onda en la sección transversal de tu modelo usando un gráfico de corte en COMSOL?
