#---------------------------------------------------#
#       FUNCIONES PARA EL TP DE NO LINEALES         #
#---------------------------------------------------#
#                                                   #
# [!] Nota: todas las funciones tienen arribita las #
#           librerías que necesitan para funcionar  #
#                                                   #
#---------------------------------------------------#

import numpy as np

""""
            Metodo de biseccion
"""
def biseccion(
        a,
        b,
        f,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i = 1
    error = 2*eps
    p = (a + b) / 2
    while i <= max_ite and error>=eps:
        i += 1
        f_a=f(a)
        f_b=f(b)
        f_p=f(p)
        if debug:
            print(f"It={i+1} | a={a} - f(a)={f_a} | b={b} - f(b)={f_b} | p={p} - f(p)={f_p}")
        if f(a)*f(p)>0:
            error = np.abs(p-a)
            a=p
        else:
            error = np.abs(b-p)
            b=p
        p = (a+b)/2
    if i > max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={p}")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

"""
            Metodo Regula Falsi
"""
def Regula_Falsi(
        a,
        b,
        f,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i=0
    error=2*eps
    while i<=max_ite and error>=eps:
        i += 1
        f_a=f(a)
        f_b=f(b)
        x=b-((f_b*(b-a))/(f_b-f_a))
        f_x = f(x)
        error = np.abs(f_x)
        if debug:
            print(f"It={i+1} | x0={a} - f(x0)={f_a} | x1={b} - f(x1)={f_b} | x2={x} - f(x2)={f_x}")
        
        if f_a*f_x>0:
            a=x
        else:
            b=x
    if i > max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={x}")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

"""
            Metodo Secante
"""

def secante(
        a,
        b,
        f,
        max_ite=100,
        eps=1e-03,
        debug=True
):
    i=0
    x0=a
    x1=b
    error =2*eps
    while i<=max_ite and error>=eps:
        x2=x1-(f(x1)*(x1-x0))/(f(x1)-f(x0))
        error=np.abs(x2-x1)
        if debug:
            print(f"It={i+1} | x0={x0} | x1={x1} | x2={x2} | error={error}")
        x0=x1
        x1=x2
        i+=1
    if i>max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={x2}")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

""""
            Metodo Newton-Raphson y Newton-Raphson Modificado
"""

def newton_raphson(
        x0,
        f,
        df,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i=0
    error=2*eps
    while i<=max_ite and error>=eps:
        x1=x0-(f(x0)/df(x0))
        error=np.abs(x1-x0)
        if debug:
            print(f"It={i+1} | x0={x0} | x1={x1} | error={error}")
        i+=1
        x0=x1
    if i>max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={x1}")

def newton_raphson_modificado(
        x0,
        f,
        df,
        df2,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i=0
    error=2*eps
    while i<=max_ite and error>=eps:
        f_x=f(x0)
        df_x=df(x0)
        df2_x=df2(x0)

        x1=x0-((f_x*df_x)/(df_x**2 - f_x*df2_x))

        error=np.abs(x1-x0)

        if debug:
            print(f"It={i+1} | x0={x0} | x1={x1} | error={error}")

        i+=1
        x0=x1
    if i>max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={x1}")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #

"""
            Metodo de Iteración de Punto Fijo
"""

def punto_fijo(
        x0,
        g,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i=0
    error=2*eps
    while i<=max_ite and error>=eps:
        x1=g(x0)
        i+=1
        error=np.abs(x1-x0)
        if debug:
            print(f"It={i+1} | x0={x0} | x1={x1} | error={error}")
        x0=x1

    if i>max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={x1}")


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - #


"""
            Metodo de Steffensen
"""

def steffensen(
        x0,
        g,
        eps=1e-03,
        max_ite=100,
        debug=True
):
    i=0
    error=2*eps
    while i<=max_ite and error>=eps:
        x1=g(x0)
        x2=g(x1)

        xs=x0-(x1-x0)**2/(x2-(2*x1)+x0)
        error=np.abs(xs-x2)

        if debug:
            print(f"It={i+1} | x0={x0} | x1={x1} | x2={x2} | xs={xs} | error={error}")
        
        i+=1
        x0=xs
    if i>max_ite:
        print(f"No converge en {i} iteraciones")
    else:
        print(f"Raíz encontrada en x={xs}")


#--------------------------------------------------------------#

"""
            Funcion Graficadora
"""

import matplotlib.pyplot as pl

def graficar(
        funciones,
        a,
        b,
        tamamio_figura=(9,5),
        puntos_intermedios=500,
        int_cerrado=True,
        y_lim=None,
        y_tope=None
):
    _validar_parametros(a,b,y_lim,y_tope)
    a = a+0.1 if not int_cerrado else a
    x = np.linspace(a,b,puntos_intermedios, endpoint=int_cerrado)

    resultados={}
    for f in funciones:
        y = _evaluar_funcion(f,x)
        if len(y[np.isfinite(y)]) == 0:
            raise RuntimeError(f"La función '{f.__name__}' no tiene valores finitos en el intervalo dado")
        resultados[f.__name__] = y
    
    _plotear(resultados,x,a,b,y_lim,y_tope,tamamio_figura)




def _validar_parametros(a,b,y_lim,y_tope):
    if a>=b:
        raise RuntimeError('El intervalo de la funcion es incorrercto (a < b)')
    if y_lim is not None:
        if y_lim[0]>y_lim[1]:
            raise RuntimeError('El límite de abscisas debe ingresarse en formato (limite_inferior, limite_superior)')
        if (y_lim[0] == y_lim[1]):
            raise RuntimeError('El límite de abscisas indicado coincide en sus dos extremos, la función no se verá correctamente')
    if y_tope is not None and y_tope < 0:
        raise RuntimeError('El tope de las abscisas no puede ser negativo, debe indicarse en unidades positivas')
    

def _evaluar_funcion(f,x):
    """Evalúa una función en x ignorando warnings matemáticos."""
    with np.errstate(all='ignore'):
        y = f(x)
    return y

def _plotear(resultados,x,a,b,y_lim,y_tope,tamamio_figura):

    y_lim_calculado = _obtener_limites(resultados)
    y_lim, y_tope= _resolver_y_lim_y_tope(y_lim,y_tope,y_lim_calculado)

    _aplicar_tope(resultados,y_tope)

    fig,ax=pl.subplots(figsize=tamamio_figura)
    for nombre,y in resultados.items():
        ax.plot(x,y,label=f'y = {nombre}(x)')
    
    _configurar_ejes(ax,a,b,y_lim)
    ax.legend()
    pl.show()



def _obtener_limites(resultados):
    todos_los_y = np.concatenate([
        y[np.isfinite(y)]
        for y in resultados.values()
    ])
    return (np.min(todos_los_y), np.max(todos_los_y))


def _resolver_y_lim_y_tope(y_lim,y_tope,y_lim_calculado):
    limites_preestablecidos=y_lim is not None
    tope_preestablecido=y_tope is not None

    if limites_preestablecidos and tope_preestablecido:
        if -y_tope > y_lim[0] or y_tope < y_lim[1]:
            print('[!] Advertencia: uno de los extremos del límite indicado supera al tope, puede que vea partes de la gráfica vacíos')
        else:
            print('[i] Información: el tope indicado es mayor al límite de abscisas indicado, puede ser ignorado')
        return y_lim, y_tope
    
    elif limites_preestablecidos:
        y_tope = max(abs(y_lim[0]), abs(y_lim[1])) * 5 
        return y_lim, y_tope
    
    elif tope_preestablecido:
        inferior = y_lim_calculado[0] if -y_tope <= y_lim_calculado[0] else -y_tope
        superior = y_lim_calculado[1] if y_tope >= y_lim_calculado[1] else y_tope
        y_tope = y_tope
        return (inferior, superior), y_tope
    
    else:
        y_tope = max(abs(y_lim_calculado[0]), abs(y_lim_calculado[1]))
        return y_lim_calculado, y_tope
    

def _aplicar_tope(resultados,y_tope):
    for y in resultados.values():
        y[np.abs(y) > y_tope] = np.nan

def _configurar_ejes(ax,a,b,y_lim):
    ax.set_ylim(y_lim[0], y_lim[1])
    ax.set_xlim(a,b)
    ax.grid(alpha=0.3)
    ax.axhline(0,color='black',linewidth=0.8,linestyle='--')
    ax.axvline(0,color='black',linewidth=0.8,linestyle='--')