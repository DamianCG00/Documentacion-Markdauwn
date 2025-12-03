# Tecnológico de Software
## Materia: Fundamentos de álgebra
## Alumno: COba Gongora Marco Damian
## Actividad \#20 - Documentación

---

# Resolución de sistemas de ecuaciones lineales

En este trabajo presento la resolución paso a paso de diferentes sistemas de ecuaciones lineales. Para ello apliqué métodos matriciales como la eliminación de Gauss, Gauss-Jordan, el método de la matriz inversa y la regla de Cramer. También incluí una sección donde clasifico los sistemas dependiendo de si tienen solución única, infinitas soluciones o ninguna.

---

## Ejercicio 1. Solución de un sistema (3x3) por varios métodos

Vamos a trabajar con el siguiente sistema de ecuaciones:

$$
\begin{cases}
x + y + z = 6\\
2x - y + z = 3\\
x + 2y - z = 2
\end{cases}
$$

Expresado en forma matricial, tenemos:

$$
A =
\begin{pmatrix}
1 & 1 & 1\\
2 & -1 & 1\\
1 & 2 & -1
\end{pmatrix},
\qquad
\mathbf{x} =
\begin{pmatrix}
x\\y\\z
\end{pmatrix},
\qquad
\mathbf{b} =
\begin{pmatrix}
6\\3\\2
\end{pmatrix},
\qquad
A\mathbf{x} = \mathbf{b}.
$$

### 1.1 Método de Gauss (eliminación)

Primero, armamos la matriz aumentada del sistema:

$$
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
2 & -1 & 1 & 3\\
1 & 2 & -1 & 2
\end{array}
\right]
$$

1. Empezamos eliminando la \(x\) de las filas 2 y 3 mediante operaciones entre renglones:

$$
R_2 \leftarrow R_2 - 2R_1,\quad
R_3 \leftarrow R_3 - R_1
$$

$$
\longrightarrow
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
0 & -3 & -1 & -9\\
0 & 1 & -2 & -4
\end{array}
\right]
$$

2. Intercambiamos las filas 2 y 3 para que nos quede un pivote más sencillo en la segunda fila:

$$
R_2 \leftrightarrow R_3
\Longrightarrow
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
0 & 1 & -2 & -4\\
0 & -3 & -1 & -9
\end{array}
\right]
$$

3. Ahora eliminamos la \(y\) de la fila 3:

$$
R_3 \leftarrow R_3 + 3R_2
$$

$$
\longrightarrow
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
0 & 1 & -2 & -4\\
0 & 0 & -7 & -21
\end{array}
\right]
$$

Ya que tenemos la forma escalonada, despejamos las variables de abajo hacia arriba:

$$
-7z = -21 \Rightarrow z = 3,
$$

$$
y - 2z = -4 \Rightarrow y - 6 = -4 \Rightarrow y = 2,
$$

$$
x + y + z = 6 \Rightarrow x + 2 + 3 = 6 \Rightarrow x = 1.
$$

El resultado final es:

$$
(x,y,z) = (1,2,3).
$$

---

### 1.2 Método de Gauss–Jordan

Para este método retomamos la matriz escalonada que obtuvimos en el paso anterior:

$$
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
0 & 1 & -2 & -4\\
0 & 0 & -7 & -21
\end{array}
\right]
$$

1. Convertimos el pivote de la tercera fila en 1 (normalización):

$$
R_3 \leftarrow -\frac{1}{7} R_3
\Longrightarrow
\left[
\begin{array}{ccc|c}
1 & 1 & 1 & 6\\
0 & 1 & -2 & -4\\
0 & 0 & 1 & 3
\end{array}
\right]
$$

2. Usamos ese 1 para eliminar la variable \(z\) de las filas 1 y 2 (hacemos ceros hacia arriba):

$$
R_2 \leftarrow R_2 + 2R_3,\quad
R_1 \leftarrow R_1 - R_3
$$

$$
\longrightarrow
\left[
\begin{array}{ccc|c}
1 & 1 & 0 & 3\\
0 & 1 & 0 & 2\\
0 & 0 & 1 & 3
\end{array}
\right]
$$

3. Finalmente, eliminamos la \(y\) de la primera fila:

$$
R_1 \leftarrow R_1 - R_2
$$

$$
\Longrightarrow
\left[
\begin{array}{ccc|c}
1 & 0 & 0 & 1\\
0 & 1 & 0 & 2\\
0 & 0 & 1 & 3
\end{array}
\right]
$$

La solución se lee directo de la matriz:

$$
x = 1,\quad y = 2,\quad z = 3.
$$

---

### 1.3 Método de la matriz inversa

Usamos la matriz de coeficientes \(A\):

$$
A =
\begin{pmatrix}
1 & 1 & 1\\
2 & -1 & 1\\
1 & 2 & -1
\end{pmatrix}.
$$

Calculamos el determinante usando la regla de Sarrus:

$$
\det(A) = 1(-1)(-1) + 1(1)(1) + 1(2)(2)
          - 1(-1)(1) - 1(2)(-1) - 1(1)(2)
        = 7.
$$

Como \(\det(A)\neq 0\), sabemos que la inversa existe y la calculamos:

$$
A^{-1} = \frac{1}{7}
\begin{pmatrix}
-1 & 3 & 2\\
3 & -2 & 1\\
5 & -1 & -3
\end{pmatrix}.
$$

Para encontrar la solución multiplicamos \(A^{-1}\) por el vector \(\mathbf{b}\):

$$
A^{-1}\mathbf{b}
= \frac{1}{7}
\begin{pmatrix}
-1 & 3 & 2\\
3 & -2 & 1\\
5 & -1 & -3
\end{pmatrix}
\begin{pmatrix}
6\\3\\2
\end{pmatrix}.
$$

Haciendo el producto matricial:

$$
\begin{pmatrix}
7\\14\\21
\end{pmatrix}
$$

Y dividiendo todo entre 7:

$$
\mathbf{x}=
\begin{pmatrix}
1\\2\\3
\end{pmatrix}.
$$

---

### 1.4 Regla de Cramer

Definimos las matrices y el vector de resultados:

$$
A =
\begin{pmatrix}
1 & 1 & 1\\
2 & -1 & 1\\
1 & 2 & -1
\end{pmatrix},\quad
\mathbf{b} =
\begin{pmatrix}
6\\3\\2
\end{pmatrix},
$$

El determinante general del sistema es:

$$
D = 7
$$

Ahora armamos las matrices sustituyendo la columna correspondiente por el vector \(\mathbf{b}\):

$$
A_x =
\begin{pmatrix}
6 & 1 & 1\\
3 & -1 & 1\\
2 & 2 & -1
\end{pmatrix},
$$

$$
A_y =
\begin{pmatrix}
1 & 6 & 1\\
2 & 3 & 1\\
1 & 2 & -1
\end{pmatrix},
$$

$$
A_z =
\begin{pmatrix}
1 & 1 & 6\\
2 & -1 & 3\\
1 & 2 & 2
\end{pmatrix}.
$$

Calculamos los determinantes de cada una:

$$
D_x=7,\quad D_y=14,\quad D_z=21
$$

Dividimos cada uno entre el determinante general para obtener la solución:

$$
x=1,\quad y=2,\quad z=3.
$$

---

## Ejercicio 2. Clasificación de tipos de sistemas

A continuación analizo tres casos diferentes para ver cómo se comportan las soluciones.

### 2.1 Sistema (a)

$$
\begin{cases}
x + y = 3\\
2x + 2y = 6
\end{cases}
$$

Matriz aumentada:

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
2 & 2 & 6
\end{array}
\right]
$$

Al reducirla obtenemos:

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
0 & 0 & 0
\end{array}
\right]
$$

Como se eliminó una fila completa, es un sistema compatible indeterminado (tiene **infinitas soluciones**).

---

### 2.2 Sistema (b)

$$
\begin{cases}
x + y = 3\\
2x + 2y = 7
\end{cases}
$$

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
2 & 2 & 7
\end{array}
\right]
$$

Al reducir llegamos a esto:

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
0 & 0 & 1
\end{array}
\right]
$$

La última fila implica \(0=1\), lo cual es imposible. Es un sistema incompatible (**sin solución**).

---

### 2.3 Sistema (c)

$$
\begin{cases}
x + y = 3\\
x - y = 1
\end{cases}
$$

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
1 & -1 & 1
\end{array}
\right]
$$

Paso 1 de reducción:

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
0 & -2 & -2
\end{array}
\right]
$$

Normalizamos la segunda fila:

$$
\left[
\begin{array}{cc|c}
1 & 1 & 3\\
0 & 1 & 1
\end{array}
\right]
$$

Eliminamos \(y\) de la primera:

$$
\left[
\begin{array}{cc|c}
1 & 0 & 2\\
0 & 1 & 1
\end{array}
\right]
$$

Resultado: Sistema compatible determinado con solución única:

$$
x=2,\quad y=1.
$$

---

## Ejercicio 3. Sistema más grande (4x4)

Finalmente, resolvemos un sistema de 4 variables:

$$
\begin{cases}
x + y + z + w = 10\\
2x + y - z + w = 5\\
x - y + z - w = 1\\
x + y - z + 2w = 8
\end{cases}
$$

Planteamos la matriz aumentada:

$$
\left[
\begin{array}{cccc|c}
1 & 1 & 1 & 1 & 10\\
2 & 1 & -1 & 1 & 5\\
1 & -1 & 1 & -1 & 1\\
1 & 1 & -1 & 2 & 8
\end{array}
\right]
$$

Después de aplicar la eliminación gaussiana llegamos a:

$$
\left[
\begin{array}{cccc|c}
1 & 1 & 1 & 1 & 10\\
0 & 1 & 3 & 1 & 15\\
0 & 0 & 1 & 0 & 7/2\\
0 & 0 & 0 & 1 & 5
\end{array}
\right]
$$

Hacemos la sustitución hacia atrás:

$$
w=5,\quad z=\frac{7}{2},\quad y=-\frac{1}{2},\quad x=2
$$

La solución del sistema es:

$$
(x,y,z,w)=\left(2,-\frac{1}{2},\frac{7}{2},5\right)
$$

---

Con estos ejercicios quedan demostrados los diferentes procedimientos para resolver y analizar sistemas lineales.
