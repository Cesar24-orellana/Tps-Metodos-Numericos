#---------------------------------------------------#
#    FUNCIONES PARA EL TP DE ECUACIONES LINEALES    #
#---------------------------------------------------#
#                                                   #
# [!] Nota: todas las funciones tienen arribita las #
#           librerías que necesitan para funcionar  #
#                                                   #
#---------------------------------------------------#

import numpy as np

"""   Eliminación de Gauss  """

def Eliminacion_Gauss(A,
                      b,
                      ):
    """
    A: Matriz
    B: Vector
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    M = np.concatenate((A,b.reshape(n,1)), axis=1) # matriz ampliada (A|b)
    
    for k in range(n-1):
        if M[k,k] == 0:
            raise ValueError(f"Pivote nulo en la fila: {k}. el metodo falla")
        for i in range(k+1, n):
            factor = M[i,k] / M[k,k]
            for j in range(k,n+1):
                M[i,j] -= factor*M[k,j]     # Construyendo matriz U
            M[i,k] = factor #    construyendo matriz L
    print(M)

    x = np.zeros(n)
    for i in range(n-1,-1,-1):
        suma_conocida = np.dot(M[i,i+1:n], x[i+1:n])
        x[i] = (M[i,n] - suma_conocida) / M[i,i]
    return x

"""     Eliminación de Gauss con pivote """

def gauss_pivote(A,
                 b):
    """
    A: Matriz
    B: Vector
    """
    A = np.array(A, dtype=float)
    b = np.array(b, dtype=float)
    n = len(b)
    M = np.concatenate((A,b.reshape(n,1)), axis=1) # matriz ampliada (A|b)
    p = np.arange(n-1)
    for k in range(n-1):
        valores_columna = [abs(M[p[i],k]) for i in range(k,n)]
        max_index_relativo = np.argmax(valores_columna)
        max_index_real = k + max_index_relativo
        if max_index_real != k:
            p[k], p[max_index_real] = p[max_index_real], p[k]
        fila_pivote = p[k]
        if M[fila_pivote, k] == 0:
            raise ValueError("El sistema no tiene solución unica (pivote nulo)")
        
        for i in range(k+1, n):
            fila_actual = p[i]
            factor = M[fila_actual, k] / M[fila_pivote,k]
            M[fila_actual, k:] -= factor * M[fila_actual,k:]
            M[fila_actual, k] = factor
        
        x = np.zeros(n)
        for i in range(n):
            fila_actual = p[i]
            suma_conocida = sum(M[fila_actual,j] * x[j] for j in range(i+1, n))
            x[i] = (M[fila_actual, n] - suma_conocida) / M[fila_actual, i]
        return x


"""    Metodo Gauss Seidel     """

def gauss_Seidel(A,
                b,
                x0,
                max_ite=100,
                eps=1e-03
                ):
    """
    A: Matriz
    B: Vector
    x0: Vector Inicial
    eps: Cota de error
    max_ite: Numero Maximo de iteraciones
    """
    k = 0
    n = len(A)
    x = x0.copy()
    error = 2*eps
    while(k < max_ite and error>eps):
        x_ant = x.copy()
        for i in range(n):
            suma = 0
            for j in range(n):
                if j != i:
                    suma += A[i][j]*x[j]
            x[i] = (b[i] - suma)/ A[i][i]
        k += 1
        #error = np.max(np.abs(x[i] -x_ant[i]) for i in range(n))
        error = max(abs(x[i] - x_ant[i]) for i in range(n))
    
    print(f"Converge en {k} iteraciones")
    return x, error


"""     Metodo Jacobi     """

def Jacobi(A,
                b,
                x0,
                max_ite=100,
                eps=1e-03
                ):
    """
    A: Matriz
    B: Vector
    x0: Vector Inicial
    eps: Cota de error
    max_ite: Numero Maximo de iteraciones
    """
    k = 0
    n = len(A)
    x = x0.copy()
    error = 2*eps
    while(k < max_ite and error>eps):
        x_ant = x.copy()
        for i in range(n):
            suma = 0
            for j in range(n):
                if j != i:
                    suma += A[i][j]*x_ant[j]
            x[i] = (b[i] - suma)/ A[i][i]
        k += 1
        error = max(abs(x[i] - x_ant[i]) for i in range(n))
    
    print(f"Converge en {k} iteraciones")
    return x, error


"""   Metodo SOR   """

def SOR(A,
        b,
        x0,
        omega,
        max_ite=100,
        eps=1e-03
        ):
    """
        A : Matriz
        b : vector
        x0 : valor inicial de prueba
        omega : 
        max_ite : numero maximo de iteraciones
        eps : error al que se decea llegar
    """
    n = len(A)
    x=x0.copy()
    norma_infinita=np.inf
    k=0
    while(norma_infinita >= eps and k<max_ite):
        x_ant = x.copy()
        for i in range(n):
            suma_total=0
            for j in range(i):
                if j != i:
                    suma_total = A[i][j]*x[j]
            x_gs = (b[i] - suma_total) / A[i][i]

            x[i] = (1-omega)* x_ant[i] + omega * x_gs
        norma_infinita = max(abs(x[i] - x_ant[i]) for i in range(n))
        k+=1
    print(f"Iteracion {k}: {x}")
    return x

    """
    while(norma_infinita >= eps and k<max_ite):
        x_ant = x.copy()
        for i in range(n):
            suma_total=0
            for j in range(i):
                if j != i:
                    suma_total = A[i][j]*x[j]
            x_gs = (b[i] - suma_total) / A[i][i]

            x[i] = (1-omega)* x_ant[i] + omega * x_gs
        norma_infinita = max(abs(x[i] - x_ant[i]) for i in range(n))
        k+=1
    """

    """
    while(norma_infinita >= eps and k<max_ite):
        x_ant = x.copy()
        for i in range(n):
            suma_total=0
            for j in range(i):
                suma_total += A[i][j]*x[j]
            for j in range(i+1,n):
                suma_total += A[i][j]*x_ant[j]
            x[i] = (1-omega)*x[i] + omega*(b[i]- suma_total)/A[i][i]
        norma_infinita = np.linalg.norm(x-x_ant, np.inf)
        k+=1
    """