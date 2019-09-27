import numpy

DIM = 100

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
    

def criterioParada(x):
    True

def mmb_secuencial():
    #Se definen las matrices y vectores con las que se trabajara.!
    A = generarMatriz()
    b = generarVectorSolucion()
    x = numpy.ones(DIM,dtype=float)
    itr = 0
    while itr < 10000:                      #Maximo de iteraciones para el metodo
        xk = numpy.copy(x)
        ek = numpy.ones(DIM,dtype=float)    #Vector de errores en ek
        for j in range(0,DIM):              #Calcular Jacobi para el elemento correspondiente
            xk[j] = numpy.copy(jacobi(x,j,A,b))      
            ek[j] = error(x,xk,j,A,b)       #Calcular el error asociado para ese error
        index = buscarError(ek)             #Buscar el error mas pequeno de la actualizacion y 
        x[index] = numpy.copy(xk[index])    #actualizar solo esa variable.
    #     if(criterioParada(x)):            #Comprobar el criterio de parada
    #         return x                      #Retornar el vector de soluciones
        itr = itr+1
    print("Sali")
    print(abs(numpy.dot(numpy.dot(0.5*x,A),numpy.transpose(x))-numpy.dot(b,x)))

mmb_secuencial()
