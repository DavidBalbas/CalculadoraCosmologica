#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: David Balbás Gutiérrez, Pablo Gómez Nicolás
@date: 9-2-2019
@version: 1.1

*****************************
** CALCULADORA COSMOLÓGICA **
*****************************

SCRIPT DE EJECUCIÓN
Calcula, a partir de las ecuaciones de Friedmann, diversos escenarios del
Universo dados una serie de parámetros cosmológicos fundamentales.

*******************************************************************************
Parámetros necesarios:

Cm: constante de materia (bariónica+dark matter)
Cr: constante de radiación
Cl: Constante cosmológica (no confundir con los Cl)
Ho: constante de Hubble

Nota: renormaliza las constantes automáticamente para que sumen 1.

Entero indicando los resultados que desean obtenerse:
1- Distancias: diámetro angular, luminosidad en función de z (redshift)
2- Evolución de a (factor de escala) y de H (cte. Hubble) en función de t
3- Valores del radio de Hubble y horizonte de partículas, para un cierto t
4- Edad del universo en función de z

Nombre del archivo de salida. En blanco, no genera resultados.
*******************************************************************************
Salida obtenida:
    
Fichero generado .txt con los datos correspondientes.
Gráfica comparativa que representa los datos del .txt
Ver el método correspondiente (módulo calcosmo) para más información
*******************************************************************************

"""

import calcosmo

print('Calculadora Cosmológica, versión 1.00')
print('Autores: Pablo Gómez Nicolás, David Balbás Gutiérrez.')
fracaso=True
while fracaso:
    print('Introduzca los parámetros cosmológicos de su modelo.')
    try:
        Cm=float(input('Constante de materia, Cm: '))
        Cr=float(input('Constante de radiación, Cr: '))
        Cl=float(input('Constante cosmológica, Cl: '))
        Ho=float(input('Constante de Hubble, Ho [(km/s)/Mpc]: '))
        #Renormaliza las constantes
        C=Cm+Cr+Cl
        if C<=0: 
            raise TypeError
        Cm=Cm/C
        Cr=Cr/C
        Cl=Cl/C
        fracaso=False
    except (ValueError, ArithmeticError):
        print('Los datos no son correctos.')
    except (TypeError):
        print('Sus constantes deben sumar un valor mayor que 0.')
    
print('Sus constantes se han renormalizado.\n')
print('Cm={0:.4f} ; Cr={1:.4f} ; Cl={2:.4f}'.format(Cm, Cr, Cl)+'\n')
print('Seleccione el tipo de resultados deseados, introduciendo 1,2,3 o 4')
print('1- Distancias: diámetro angular, luminosidad en función de z (redshift)')
print('2- Evolución de a (factor de escala) y de H (cte. Hubble) en función de t')
print('3- Valores del radio de Hubble y horizonte de partículas, en función de t')
print('4- Edad del universo en función de z')


def desea():
    """
    Método auxiliar para realizar la pregunta de si desea realizar otro cálculo
    """
    deseado=input('¿Desea realizar otro cálculo? (s/n): ')
    if deseado=='s' or deseado=='S':
        return False
    elif deseado=='n' or deseado=='N':
        return True
    else:
        desea()
        

correcto=False
while correcto==False:
    ind=input('Introduzca 1, 2, 3 o 4: ')
    out=input('Introduzca el nombre del fichero de salida de sus resultados: ')
    if ind=='1':
        correcto=desea()
        calcosmo.modeloDistancia(Cm, Cr, Cl, Ho, out)
    elif ind=='2':
        correcto=desea()
        calcosmo.modeloEscala(Cm, Cr, Cl, Ho, out)
    elif ind=='3':
        correcto=desea()
        calcosmo.modeloRadios(Cm, Cr, Cl, Ho, out)
    elif ind=='4':
        correcto=desea()
        calcosmo.modeloEdad(Cm, Cr, Cl, Ho, out)
    else:
        print('Modelo no válido. Seleccione un número del 1 al 4')

print('Gracias por usar la Calculadora Cosmológica')


    
