__precompile__()
module herramientas

export newton, newton_raices, riemann, trapecio, simpson, rk_4, euler_method, euler_imp, punto_medio, Lagrange, derivada_numerica, derivada_simetrica

 """La función newton(f,df,x0) tiene fomo argumentos la función f, su derivada (df) y la condición inicial (x0). El output de esta función es la raiz de f más cercana a x0."""
function newton(f,df,x0)
    #Comenzamos tomando nuestra condición inicial
    x=x0;
    #Hacemos el for para aplicar el método de Newton 100 veces
    for i in 1:100
        x=x-f(x)/df(x)
    end
    #La función regresa el último valor de x obtenido
    return x
end


"""newton_raices(f,df,lin) nos regresa las raices distintas de una función encontradas a partir de varias condiciones iniciales. f es la función a la que se quiere encontrar las raíces, df es su derivada y lin es la lista de condiciones iniciales."""
function newton_raices(f,df,lin)
#Definimos dos arreglos vacíos y nuestro error e
    t=[];
    s=[];
    e=(10.0)^(-8)
    
#El primer for es para ir tomando cada entrada del argumento lin
    for j in 1:length(lin)
#Tomamos nuestro valor inicial
        x=lin[j];
#El segundo for es para aplicar el método de Newton 200 veces.
        for i in 1:200
        x=x-f(x)/df(x)
        end
#Vamos guardando los valores de las raíces en el arreglo t
        push!(t,x)
    end


#Tomamos la el primer elemento de t (una de las raícez) y lo guardamos en s  
    push!(s,t[1])
    
#Hacemos un for para ir tomando elementos de t a partir del segundo
    for i in 2:length(t)
#Definimos z con un valor inicial de cero
        z=0;
#Hacemos un segundo for para ir comparando el elemento t[i] con todos los de s
       for j in 1:length(s)
#Cada que la diferencia sea mayor al error, sumamos 1 al valor de z
            if abs(t[i]-s[j])>e
                z=z+1
            end
        end
#Al final, si z es igual a la longitud de s, significa que las diferencias siempre fueron mayores que e, por lo cual el valor de t[i] que estamos comparando es una raíz diferente. Entonces, guardamos dicho valor en el arreglo s.
        if z==length(s)
            push!(s,t[i])
        end
        
    end
#La función regresa s, el arreglo obtenido, donde sólo hay raíces diferentes.
       return s
end;



"""La función riemann(f,a,b,e) regresa la integral por el método del rectángulo, donde f es la función dada, a y b son los extremos de integración y e es la longitud deseada para los subintervalos"""
function riemann(f,a,b,e)
    #Si e es la longitud deseada para los subintervalos, entonces el número de subintervalos será (b-a)/e para asegurar que sea entero usamos round() (como se muestra abajo), pero entonces el número de elementos en nuestra partición tiene que ser el número anterior más uno, de esta forma definimos h
    h=1+round((b-a)/e)
   #Definimos la lista de h elementos de "a" a "b" 
    list=linspace(a,b,h)
    #Definimos nuestro valor inicial
    int=0
    #El for es para ir sumando la contribución de cada intervalo a la integral
    for i in 2:length(list)
        #A int le vamos sumando el área de cada rectángulo.
        int=int+(list[i]-list[i-1])*f((list[i]+list[i-1])/2)
    end
    #La función regresa el valor de la integral
    return int    
end;


"""La función trapecio(f,a,b,e) regresa la integral por el método el trapecio, donde f es la función dada, a y b son los extremos de integración y e es la longitud deseada para los subintervalos"""
function trapecio(f,a,b,e)
    #Si e es la longitud deseada para los subintervalos, entonces el número de subintervalos será (b-a)/e para asegurar que sea entero usamos round() (como se muestra abajo), pero entonces el número de elementos en nuestra partición tiene que ser el número anterior más uno, de esta forma definimos h
    h=1+round((b-a)/e)
   #Definimos la lista de h elementos de "a" a "b" 
    list=linspace(a,b,h)
    #Definimos nuestro valor inicial
    int=0
    #El for es para ir sumando la contribución de cada intervalo a la integral
    for i in 2:length(list)
        #A int le vamos sumando el área de cada trapecio.
        int=int+(list[i]-list[i-1])*((f(list[i])+f(list[i-1]))/2)
    end
    #La función regresa el valor de la integral
    return int    
end;




"""La función simpson(f,a,b,e) regresa la integral por el método de Simpson, donde f es la función dada, a y b son los extremos de integración y e es la longitud deseada para los subintervalos."""
function simpson(f,a,b,e)
    #Si e es la longitud deseada para los subintervalos, entonces el número de subintervalos será (b-a)/e para asegurar que sea entero usamos round() (como se muestra abajo), pero entonces el número de elementos en nuestra partición tiene que ser el número anterior más uno, de esta forma definimos h
    h=1+round((b-a)/e)
   #Definimos la lista de h elementos de "a" a "b" 
    list=linspace(a,b,h)
    #Definimos nuestro valor inicial
    int=0
    #El for es para ir sumando la contribución de cada intervalo a la integral
    for i in 2:length(list)
        #A int le vamos sumando el área debajo de cada parábola.
        
        int=int+((list[i]-list[i-1])/6)*(f(list[i-1])+4*f((list[i]+list[i-1])/2)+f(list[i]))
        
    end
        
    #La función regresa el valor de la integral
    return int    
end;


"""La función Lagrange(X,Y,x), cuyos argumentos son: X la lista de las abscisas; Y, la lista de las correspondientes ordenadas y x, el punto a evaluar, nos regresa el polinomio de Lagrange (evaluado en x) asociado a las parejas dadas por X y Y"""
function Lagrange(X,Y,x)
    #Vamos a ir guardando cada término en el polinomio pol, que definimos como 0 al inicio
    pol=0
    #El siguiente for es para ir formando nuestros polinomios lj
    for j in 1:length(Y)
        #Defnimos l como 1 al inicio. Agregamos el z/z, para que Julia sepa que se trata de una función y no de una constante
        l=1
        
        #para m distinto de j, vamos multiplicando los términos correspondientes
        for m in 1:length(X)
            if j!=m
               l=l*(x-X[m])/(X[j]-X[m])
            end  
        end
        #A nuestro polinomio le vamos sumando la contribución de los lj, multiplicados por la respectiva yj
        pol=pol+Y[j]*l
    end
 

    #La función regresa el valor del polinomio evaluado en x
    return pol
end;


"""Método de Runge-Kutta de orden 4 para funciones de cualquier dimensión. Las entradas son: f, la función que nos da la derivada (debe ser del tipo f(x,t)); listt, el arreglo de t's; y x0, el valor inicial."""
function rk_4(f,listt,x0)
    
#Definimos h como la distancia entra las dos primeras entradas de listt.
    h=listt[2]-listt[1]
#Definimos x como nuestro valor inicial (arreglo de valores iniciales)    
    x=x0
    
#Definimos listx como un objeto de tipo indefinido. Aquí es donde iremos guardando los puntos de la solución.
    listx=[]
    
#Guardamos nuestro valor inicial en listx
    push!(listx,x)

#el siguiente for es para aplicar la relación de recurrencia
    for i in 1:length(listt)-1
        k1=f(x,listt[i])
        k2=f(x+(h/2)*k1,listt[i]+h/2)
        k3=f(x+(h/2)*k2,listt[i]+h/2)
        k4=f(x+h*k3,listt[i+1])
        
        x=x+(h/6)*(k1+2*k2+2*k3+k4)
        
 #Cada que aplicamos la relación de recurrencia, vamos guardando el punto en listx.       
        push!(listx,x)
    end
#la función regresa el arreglo de x's
    return listx
end;


"""euler_method(f,listt,x0).Método de Euler explícito para funciones de cualquier dimensión. Las entradas son: f (debe ser del tipo f(x,t)), la función que nos da la derivada;listt,el arreglo de t's; y x0, el valor inicial."""
function euler_method(f,listt,x0)
    
#Definimos h como la distancia entra las dos primeras entradas de listt.
    h=listt[2]-listt[1]
    
#Definimos x como nuestro valor inicial (arreglo de valores iniciales)
    x=x0
#Definimos listx como un objeto de tipo indefinido. Aquí es donde iremos guardando los puntos de la solución.
    listx=[]
    
#Guardamos nuestro valor inicial en listx
    push!(listx,x)

#el siguiente for es para aplicar la relación de recurrencia
    for i in 1:length(listt)-1
        
#Recordemos que el valor inicial de x es x0, por lo que el rango de i es correcto.
        x = x + f(x,listt[i])*h
        
#Cada que aplicamos la relación de recurrencia, vamos guardando el punto en listx.
        push!(listx,x) 
    end
#la función regresa el arreglo de x's
    return listx
end;


"""Método implícito de Euler para una dimensión. euler_imp(f,dxf,listt,x0). Esta función tiene como argumentos la función f (debe ser del tipo f(x,t)) , su derivada respecto a x (dxf -también debe ser función de x y t-), la lista de t's (listt) y la condición inicial (x0)."""
function euler_imp(f,dxf,listt,x0)
    
#Definimos listx como un arreglo de zeros, con el mismo numero de entradas que listt
    listx=zeros(length(listt))

#Guardamos nuestra condición inicial en la primera entrada   
    listx[1]=x0

#Definimos h como la distancia entra las dos primeras entradas de listt.
    h=listt[2]-listt[1]
    
#El siguiente for es para ir guardando las sguientes entradas de listx
    for i in 1:length(listt)-1
        
#Definimos nuestra función g, dada por el método implícito, a la cual hay que encontrarle la raíz (que será nuestro x_{i+1})
    g(z)=z-f(z,listt[i+1])*h-listx[i]
        
#Derivamos respecto a z para definir la derivada de g
    dg(z)= 1-dxf(z,listt[i+1])*h  

#Aplicamos el método de Newton a g, tomando como condición inicial nuestra x_{i} y guardando el resultado (la raíz) como nuestro x_{i+1}.
    listx[i+1]=newton(g,dg,listx[i])    
        
    end
    
#La función regresa el arreglo de valores de x
    return listx
    
end;


"""Método del punto medio. punto_medio(f,listt,x0). Las entradas son: f, la función que nos da la derivada; listt, la lista de t's; y x0, el valor inicial."""
function punto_medio(f,listt,x0)
#Definimos h como la distancia entra las dos primeras entradas de listt.
    h=listt[2]-listt[1]
    
#Definimos listx como un arreglo de zeros, con el mismo numero de entradas que listt 
    listx=zeros(length(listt))
    
#Guardamos la condición inicial en la primera entrada de nuestra lista de x's
    listx[1]=x0
    
#El siguiente for es para aplicar la relación de recurrencia
    for i in 1:length(listt)-1
       listx[i+1]= listx[i]+h*f(listx[i]+(h/2)*f(listx[i],listt[i]),listt[i]+h/2)
    end 
    
#la función regresa la lista de x's
    return listx
    
end;


"""La función derivada_numerica(f,x,h), regresa la derivada numérica de la función f evaluada en x, tomando como diferencia finita el valor h."""
function derivada_numerica(f,x,h)
    #simplemente aplicamos la definición sin la parte del límite y lo guardamos como "d"
    d=(f(x+h)-f(x))/h
    #la función regresa el valor de d
    return d
end;

"""La función derivada_simetrica(f,x,h), regresa la derivada simétrica de la función f evaluada en x, tomando como diferencia finita el valor h."""
function derivada_simetrica(f,x,h)
    d=(f(x+h)-f(x-h))/(2*h)
    return d
end;


end