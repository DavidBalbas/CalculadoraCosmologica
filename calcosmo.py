#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""

@author: David Balbás Gutiérrez, Pablo Gómez Nicolás
@date: 9-2-2019
@version: 1.3

Librería que contiene los métodos de la calculadora cosmológica.

"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp,quad

c=299792

def modeloDistancia(Cm, Cr, Cl, Ho, out):
    """
    Calcula, en función del redshift, las distancias de diámetro angular y de
    luminosidad, con los parámetros cosmológicos dados.
    
    Salida: Dibuja y guarda una representación gráfica en out'_fig.png',
    y devuelve los datos en un fichero out'.txt'.
    
    Formato de los datos: Por columnas, z, dl, da, separados por un espacio
    en blanco. Primera fila: header.
    
    """
    
    #Calculamos la distancia de luminosidad en función de E, z
    zini=0 #Empezamos a redshift 0
    zfin=15 
    npoints=20001
    
    step=(zfin-zini)/(npoints-1)
    
    z=np.linspace(zini,zfin,num=npoints)
    integral=[]
    integral.append(0) #z=0
    dl=[]
    #Se integra utilizando la regla trapezoidal, I(i+1)=I(i)+[f(i+1)+f(i)]/2
    for i in range(len(z)-1):
        integral.append(integral[i]+(1/(np.sqrt(Cl+Cm*(1+z[i+1])**3+Cr*(1+z[i+1])**4)))*0.5*step+(1/(np.sqrt(Cl+Cm*(1+z[i])**3+Cr*(1+z[i])**4))*0.5*step))
    
    #Distancias de luminosidad y diametro angular
    dl=c*(1+z)/Ho*integral
    da=dl*np.power((1+z),-2)
    
    #Plot
    plt.figure(figsize=(7,5))
    plt.plot(z,da,linewidth=2, label='$D_A$')
    plt.plot(z,dl,linewidth=2, label='$D_L$')
    plt.yscale('log')
    plt.xscale('log')
    plt.title('Distancias de luminosidad y de diámetro angular')
    plt.xlabel('Redshift')
    plt.ylabel('Distancia (Mpc)')
    plt.legend()
    
    
    #Guardado de datos
    if out!='':
        plt.savefig(out+'_fig.png')
        f = open(out+".txt", "w")
        f.write("z  dl(Mpc)  da(Mpc) \n")
        for i in range(npoints):
            f.write(str(z[i])+' '+str(dl[i])+' '+str(da[i])+"\n")
        f.close()
    
    plt.show()
    plt.close()
    

def modeloEscala(Cm, Cr, Cl, Ho, out):
    """
    Calcula, en función del tiempo, la evolución del factor de escala y de
    la constante de Hubble.
    
    Salida: Dibuja y guarda una representación gráfica en $out$_fig_x.png,
    donde x=a, h1 o h2 dependiendo de la imagen;
    y devuelve los datos en un fichero out.txt.
    
    Formato de los datos: Por columnas, t, a, H, separados por un espacio
    en blanco. Primera fila: header.
    """
    
    Ho=Ho*3.154e16/3.086e19 #in Gyr^-1
    
    def da(t,a): return a*Ho*np.sqrt(np.abs(Cl+Cm*a**-3+Cr*a**-4))
    def inta(a): return 1/(np.sqrt(np.abs(Cl+Cm*a**-3+Cr*a**-4))*a*Ho)
    
    
    universeTime=quad(inta,0,1)[0] #Tiempo del universo en Gyr
    
    #Resolvemos la ecuación diferencial utilizando un método Runge-Kutta
    #Explícito de orden 5(4).
    
    sol=solve_ivp(da, [0,-universeTime], [1], max_step=0.02)
    t=np.asarray(sol.t)+universeTime
    scaleFactor=np.asarray(sol.y[0,:])
    sol2=solve_ivp(da, [0,universeTime*3], [1], max_step=0.05)
    
    t2=np.asarray(sol2.t)+universeTime
    scaleFactor2=np.asarray(sol2.y[0,:])
    
    #Plots varios y guardado de figuras

    fig=plt.figure(figsize=(7,5))
    plt.plot(t,scaleFactor, linewidth=2, label='$a<1$')
    plt.plot(t2,scaleFactor2, linewidth=2, label='$a>1$')
    plt.title('Evolución del factor de escala')
    plt.xlabel('Tiempo desde $a=0$ (Gyr)')
    plt.ylabel('Factor de escala, $a$')
    plt.text(0.5,0.17, 'Edad actual ($a=1$): %.2f Gyr' % (universeTime), transform=fig.transFigure)
    plt.legend()
    if out!='':
        plt.savefig(out+'_fig_a.png')
    plt.show()
    plt.close()
    
    h=da(t,scaleFactor)/scaleFactor/(3.154e16/3.086e19)
    h2=da(t2,scaleFactor2)/scaleFactor2/(3.154e16/3.086e19)
    
    plt.figure(figsize=(7,5))
    plt.plot(t2,h2, linewidth=2, label='$a>1$')
    plt.title('Evolución de la Constante de Hubble ($t>t_0$)')
    plt.xlabel('Tiempo desde $a=0$ (Gyr)')
    plt.ylabel('Constante de Hubble, ((Km/s)/Mpc)')
    plt.legend()
    if out!='':
        plt.savefig(out+'_fig_h1.png')
    plt.show()
    plt.close()
    
    plt.figure(figsize=(7,5))
    plt.plot(t,h, linewidth=2, label='$a<1$')
    plt.title('Evolución de la Constante de Hubble ($t<t_0$)')
    plt.xlabel('Tiempo desde $a=0$ (Gyr)')
    plt.ylabel('Constante de Hubble, ((Km/s)/Mpc)')
    plt.yscale('log')
    plt.legend()
    
    if out!='':
        plt.savefig(out+'_fig_h2.png')
    plt.show()
    plt.close()
    
    #Guardado de datos
    if out!='':
        f = open(out+".txt", "w")
        f.write("t(Gyr)  a  H((Km/s)/Mpc)  \n")
        for i in range(len(t)):
            f.write(str(t[len(t)-i-1])+' '+str(scaleFactor[len(t)-i-1])+' '+str(h[len(t)-i-1])+"\n")
        for i in range(len(t2)):
            f.write(str(t2[i])+' '+str(scaleFactor2[i])+' '+str(h2[i])+"\n")
        f.close()
    
    
def modeloRadios(Cm, Cr, Cl, Ho, out):
    """
    Calcula, para un cierto tiempo, los valores del radio de Hubble y del 
    horizonte de sucesos, en unidades de distancia propia.
    
    Salida: Dibuja y guarda una representación gráfica en $out$_fig_1.png, y
    en $out$_fig_2.png. 
    y devuelve los datos en un fichero out.txt.
    
    Formato de los datos: Por columnas: tiempo, horizonte de partículas, radio
    de Hubble. Separados por espacios. Primera fila: Header
        
    
    """
    
    Ho=Ho*3.154e16/3.086e19 #in Gyr^-1
    
    def da(t,a): return np.sign(Cl+Cm*a**-3+Cr*a**-4)*a*Ho*np.sqrt(np.abs(Cl+Cm*a**-3+Cr*a**-4))
    def inta(a): return np.sign(Cl+Cm*a**-3+Cr*a**-4)/(np.sqrt(np.abs(Cl+Cm*a**-3+Cr*a**-4))*a*Ho)
    
    universeTime=quad(inta,0,1)[0] #Tiempo del universo en Gyr
    
    #Resolvemos la ecuación diferencial utilizando un método Runge-Kutta
    #Explícito de orden 5(4).
    
    sol=solve_ivp(da, [0,-universeTime], [1], max_step=0.02)
    t=np.asarray(sol.t)+universeTime
    scaleFactor=np.asarray(sol.y[0,:])
    sol2=solve_ivp(da, [0,universeTime*3], [1], max_step=0.05)
    
    t2=np.asarray(sol2.t)+universeTime
    scaleFactor2=np.asarray(sol2.y[0,:])
    
    h=da(t,scaleFactor)/scaleFactor/(3.154e16/3.086e19)
    h2=da(t2,scaleFactor2)/scaleFactor2/(3.154e16/3.086e19)
    
    #Reorganizamos las matrices de datos por simplicidad
    t=np.flip(t)
    scaleFactor=np.flip(scaleFactor)
    h=np.flip(h)
    t=np.concatenate((t,t2))
    scaleFactor=np.concatenate((scaleFactor,scaleFactor2))
    rH=c/h
    rH2=c/h2
    rHtot=np.concatenate((rH,rH2))/1000
    
    #Integramos numéricamente para obtener el horizonte de partículas
    dp=np.zeros(len(scaleFactor)-1)
    for i in range(len(scaleFactor)-1):
        j=0
        while j<=i:
            #Distancia en Mpc
            dp[i]=dp[i]+(c/scaleFactor[j])*(t[j+1]-t[j])*3.154e16/3.086e22
            j=j+1
        dp[i]=dp[i]*scaleFactor[i]
    
    #Plots y guardado de datos
    fig=plt.figure(figsize=(7,5))
    plt.plot(t[1:],dp, linewidth=2, label='$a<1$')
    plt.title('Horizonte de partículas')
    plt.xlabel('Tiempo desde $a=0$ (Gyr)')
    plt.ylabel('$d_P$, (Gpc)')
    plt.yscale('linear')
    plt.text(0.5,0.17, 'Horizonte actual ($a=1$): %.2f Gpc' % (dp[len(sol.t)]),transform=fig.transFigure)
    
    if out!='':
        plt.savefig(out+'_fig_1.png')
    plt.show()
    plt.close()
    
    fig=plt.figure(figsize=(7,5))
    plt.plot(t,rHtot, linewidth=2)
    plt.title('Radio de Hubble')
    plt.xlabel('Tiempo desde $a=0$ (Gyr)')
    plt.ylabel('$r_H$, Gpc')
    plt.text(0.5,0.17, 'Radio actual ($a=1$): %.2f Gpc' % (c/Ho/1000000),transform=fig.transFigure)

    if out!='':
        plt.savefig(out+'_fig_2.png')
    plt.show()
    plt.close()
    
    #Guardado de datos
    if out!='':
        f = open(out+".txt", "w")
        f.write("t(Gyr)  dP(Gpc)  rH(Gpc)  \n")
        for i in range(len(t)-1):
            f.write(str(t[i])+' '+str(dp[i])+' '+str(rHtot[i])+"\n")
        f.close()
    
    
def modeloEdad(Cm, Cr, Cl, Ho, out):
    """
    Calcula, en función del redshift, la edad del universo.
    
    Salida: Dibuja y guarda una representación gráfica en $out$_fig.png,
    y devuelve los datos en un fichero out.txt.
    Formato de los datos: por columnas, z, Edad. Separación por espacios.
    Primera fila: Header.
    """
    
    #Redshift and scale factor: 1+z=1/a(t)
    Ho=Ho*3.154e16/3.086e19 #in Gyr^-1
    
    def inta(a): return np.sign(Cl+Cm*a**-3+Cr*a**-4)/(np.sqrt(np.abs(Cl+Cm*a**-3+Cr*a**-4))*a*Ho)
    
    #Integramos numéricamente
    z=np.linspace(0,50,num=500)
    a=1/(1+z)
    universeTime=np.zeros(500)
    for i in range(500):
        universeTime[i]=quad(inta,0,a[i])[0]
        
    #Plot
    fig=plt.figure(figsize=(7,5))
    plt.plot(z,universeTime, linewidth=2)
    #plt.plot(t2,h2, linewidth=2, label='$a>1$')
    plt.title('Edad del Universo')
    plt.xlabel('Redshift, z')
    plt.ylabel('Edad del universo (Gyr)')
    plt.yscale('log')
    plt.text(0.5,0.81, 'Edad actual ($a=1$): %.2f Gyr' % (universeTime[0]),transform=fig.transFigure)
   
    #Guardado de datos
    if out!='':
        plt.savefig(out+'_fig.png')
        f = open(out+".txt", "w")
        f.write("z  Edad \n")
        for i in range(500):
            f.write(str(z[i])+' '+str(universeTime[i])+"\n")
        f.close()
    
    plt.show()
    plt.close()