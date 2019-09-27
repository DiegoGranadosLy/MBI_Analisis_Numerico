import numpy
from matplotlib import pyplot

DIM = 50

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

def mmb_secuencial(tol):
    #Se definen las matrices y vectores con las que se trabajara.!
    A = generarMatriz()
    b = generarVectorSolucion()
    x = numpy.ones(DIM,dtype=float)
    itr = 0
    errores   = numpy.array([])             #Vector de errores
    itrVector = numpy.array([])             #Vector de errores
    while itr < 100000:                     #Maximo de iteraciones para el metodo
        xk = numpy.copy(x)
        ek = numpy.ones(DIM,dtype=float)    #Vector de errores en ek
        for j in range(0,DIM):              #Calcular Jacobi para el elemento correspondiente
            xk[j] = numpy.copy(jacobi(x,j,A,b))      
            ek[j] = error(x,xk,j,A,b)       #Calcular el error asociado para ese error
        index = buscarError(ek)             #Buscar el error mas pequeno de la actualizacion y 
        x[index] = numpy.copy(xk[index])    #actualizar solo esa variable.
        if(criterioParada(A,b,x,tol,errores,itrVector,itr)):            #Comprobar el criterio de parada
            print(x)
            return x                      #Retornar el vector de soluciones
        itr = itr+1
    print(x)
    return x                      #Retornar el vector de soluciones

mmb_secuencial(10e-10)
