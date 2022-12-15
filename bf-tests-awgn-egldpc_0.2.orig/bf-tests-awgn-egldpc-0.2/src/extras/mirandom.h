/*
 * mirandom.h
 * 
 * Copyright 2011 Fernando Pujaico Rivera <fernando.pujaico.rivera@gmail.com>
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
 * MA 02110-1301, USA.
 * 
 */

#ifndef __MI_RANDOM_H__
#define __MI_RANDOM_H__

#include <math.h>

/* Máximo valor que toman las variables aleatorias uniformemente distribuidas */
#define MY_RAND_MAX 16777216

/**\fn long int myrand(void)
 * \brief Genera un número pseudo aleatorio y uniforme entre 0 y MY_RAND_MAX-1.
 * \return Un número entre 0 y MY_RAND_MAX-1.
 */
long int myrand(void)
{
     static int z=15;
     long int a=9,c=11;
     z=(a*z+c)%MY_RAND_MAX;
     return z;
}


/**\fn long int myrand1(void)
 *  \brief Genera un número pseudo aleatorio e uniforme entre 0 y MY_RAND_MAX-1.
 *
 *  Este número es identico en estructura a myrand, pero es generado 
 *  independientemente.
 * \return Un número entre 0 y MY_RAND_MAX-1.
 */
long int myrand1(void)
{
     static int x=15;
     long int a=9,c=11;
     x=(a*x+c)%MY_RAND_MAX;
     return x;
}

/**\fn long int myrand2(void)
 * \brief Genera un número pseudo aleatorio e uniforme entre 0 y MY_RAND_MAX-1.
 *
 *  Este número es diferente en estructura a myrand1.
 * \return Un número entre 0 y MY_RAND_MAX-1.
 */
long int myrand2(void)
{
     static int y=3;
     long int a=9,c=13;
     y=(a*y+c)%MY_RAND_MAX;
     return y;
}

/**\fn float MyRandGauss(void)
 * \brief Genera un número pseudo aleatorio de distribución gausiana de media cero
 *  y varianza 1.0.
 *
 *  Este número usa en sus generación a myrand1 y myrand2.
 * \return Valores de una variable aleatoria gausiana N(0,1).
 */
float MyRandGauss(void)
{
        return sqrt(-2.0*log((1.0*myrand1())/MY_RAND_MAX))*cos((2.0*M_PI*myrand2())/MY_RAND_MAX);
}

#endif
