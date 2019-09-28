import numpy
import time
from matplotlib import pyplot
from joblib import Parallel, delayed
import multiprocessing

#variables y vectores globales requeridos
DIM = 50
A = []
b = []
x = []
    
def generarMatriz():
    A = numpy.zeros((DIM,DIM),dtype=float)
    for i in range(0,DIM):
        if(i==0):
            A[i,i]   = 6
            A[i,i+1] = 2
        elif(i==DIM-1):
            A[i,i]   = 6
            A[i,i-1] = 2
        else:
            A[i,i]   = 6
            A[i,i+1] = 2
            A[i,i-1] = 2
    return A

def generarVectorSolucion():
    b = numpy.zeros(DIM,dtype=float)
    for i in range(0,DIM):
        if(i==0 or i==DIM-1):
            b[i] = 12  
        else:
            b[i] = 15  
    return b

def buscarError(ek):
    index = 0
    for i in range(0,DIM):
        if(ek[i] > ek[index]):
            index = i
    return index

def jacobi(x,j,A,b):
    if(j==0):           #En caso de que sea el primer elemento de la diagonal 
        return (b[j]-(0.5*x[j+1]*A[j,j+1])-(0.5*x[j+1]*A[j+1,j]))/(2*0.5*A[j,j])

    elif(j==DIM-1):     #En caso de que sea el ultimo elemento de la diagonal
        return (b[j]-(0.5*x[j-1]*A[j,j-1])-(0.5*x[j-1]*A[j-1,j]))/(2*0.5*A[j,j])

    else:               #Resto de los casos (0 < j < DIM-1)
        return (b[j]-(0.5*x[j-1]*A[j,j-1])-(0.5*x[j-1]*A[j-1,j])-(0.5*x[j+1]*A[j+1,j])-(0.5*x[j+1]*A[j,j+1]))/(2*0.5*A[j,j])

def error(x,xk,j,A,b):
    xNew = numpy.copy(x)
    xNew[j] = xk[j]
    
    return abs(numpy.dot(numpy.dot(0.5*xNew,A),numpy.transpose(xNew))-numpy.dot(b,xNew))
    

def criterioParada(A,b,x,tol,errores,itrVector,itr):
    dxk = numpy.zeros(DIM,dtype=float)
    for i in range (0,DIM):
        dxk[i] = numpy.copy(x[i]-jacobi(x,i,A,b))
    norma = numpy.linalg.norm(dxk)
    errores  = numpy.append(numpy.copy(norma),errores)
    itrVector =  numpy.append(numpy.copy(itr),itrVector)
    if(norma < tol):
        return True
    else:
        return False

def jacobi_method(j,xk,ek):
    
    xk[j] = numpy.copy(jacobi(x,j,A,b))
    ek[j] = error(x,xk,j,A,b)       #Calcular el error asociado para ese error
    return [xk[j],ek[j]]


def getResults(mat):
    xk = []
    ek = []
    for i in range(len(mat)):
        xk.append(mat[i][0])
        ek.append(mat[i][1])
    return [xk,ek]
def mmb_paralelo(tol):

    itr = 0
    errores   = numpy.array([])             #Vector de errores
    itrVector = numpy.array([])             #Vector de errores

    #num_cores = multiprocessing.cpu_count()
    while itr < 100000:                     #Maximo de iteraciones para el metodo
        xk = numpy.copy(x)
        ek = numpy.ones(DIM,dtype=float)    #Vector de errores en ek

        result = Parallel(n_jobs=1)(delayed(jacobi_method)(i,xk,ek) for i in range(0,DIM)) #se paraleliza el for loop con n_jobs

        [xk,ek] = getResults(result) #los resultados son copiados a los vectores correspondientes
        index = buscarError(ek)             #Buscar el error mas pequeno de la actualizacion y 
        x[index] = numpy.copy(xk[index])    #actualizar solo esa variable.

        if(criterioParada(A,b,x,tol,errores,itrVector,itr)):            #Comprobar el criterio de parada
            print(x)
            return x                      #Retornar el vector de soluciones
        itr = itr+1
    print(x)
    return x                      #Retornar el vector de soluciones



if __name__ == '__main__':

    DIM = 50
    #Se definen las matrices y vectores con las que se trabajara.!
    A = generarMatriz()
    b = generarVectorSolucion()
    x = numpy.ones(DIM,dtype=float)
    start_time = time.time()
    mmb_paralelo(10e-5)
    end_time = time.time()
    print("execution time: ")
    print(end_time - start_time)
