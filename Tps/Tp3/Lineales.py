#---------------------------------------------------#
#    FUNCIONES PARA EL TP DE ECUACIONES LINEALES    #
#---------------------------------------------------#
#                                                   #
# [!] Nota: todas las funciones tienen arribita las #
#           librerías que necesitan para funcionar  #
#                                                   #
#---------------------------------------------------#

import numpy as np


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