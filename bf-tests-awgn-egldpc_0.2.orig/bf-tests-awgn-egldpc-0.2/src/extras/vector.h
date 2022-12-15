/*
 * vector.h
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

#ifndef __VECTOR_H__
#define __VECTOR_H__

#include <stdio.h>

/** \struct Vector vector.h "extras/vector.h"
 *  \brief Vector de valores enteros.
 */
typedef struct
{
	int N;
	int *data;
}Vector;


/** \struct Linha vector.h "extras/vector.h"
 *  \brief Linha de valores enteros.
 */
typedef struct
{
	int N;
	int *data;
}Linha;

/** \struct VectorFloat vector.h "extras/vector.h"
 *  \brief Vector de valores reales.
 */
typedef struct
{
	int N;
	float *data;
}VectorFloat;


/** \fn Vector *CrearVector(unsigned short n)
 *  \brief Crea un vector de n elementos con valor cero en ellos.
 *  \param[in] n Número de elementos del vector.
 *  \return Un puntero a la estructura Vector, el vector es de longitud n y con 
 *  elementos con valor cero.
 */
Vector *CrearVector(unsigned short n)
{
	Vector *tmp=NULL;

	tmp=calloc(1, sizeof(Vector));

	tmp->N=n;
	tmp->data=calloc(n, sizeof(int));

	return tmp;
}

/** \fn void LiberaVector(Vector *tmp)
 *  \brief Destruye el vector apuntando por tmp.
 *  \param[in] tmp Puntero al vector que se desea destruir.
 */
void LiberaVector(Vector *tmp)
{
	free(tmp->data);
	free(tmp);
}

/** \fn void IniciaConZeroVector(Vector *tmp)
 *  \brief Rellena el vector con cero.
 *  \param[in] tmp Puntero al vector que se desea setear con cero.
 */
void IniciaConZeroVector(Vector *tmp)
{
	int i;
	for(i=0;i<tmp->N;i++)	tmp->data[i]=0;
}


/** \fn void IniciaConUnoZeroVector(Vector *tmp)
 *  \brief Rellena el vector con cero y uno ciclicamente.
 *  \param[in] tmp Puntero al vector que se desea setear el valor.
 */
void IniciaConUnoZeroVector(Vector *tmp)
{
	int i;
	tmp->data[0]=0;
	for(i=1;i<tmp->N;i++)	tmp->data[i]=tmp->data[i-1]^1;
}

void IniciaRandomVector(Vector *tmp)
{
	int i;

	srand(time(NULL));

	for(i=0;i<tmp->N;i++)	
	{
		if(rand()>(RAND_MAX/2))	tmp->data[i]=1;
		else		tmp->data[i]=0;
	}
}

/***********************************************/
Linha *CrearLinha(unsigned short n)
{
	Linha *tmp=NULL;

	tmp=calloc(1, sizeof(Linha));

	tmp->N=n;
	tmp->data=calloc(n, sizeof(int));

	return tmp;
}

void LiberaLinha(Linha *tmp)
{
	free(tmp->data);
	free(tmp);
}

/***********************************************/

/** \fn VectorFloat *CrearVectorFloat(unsigned short n)
 *  \brief Crea un vector de n elementos reales con valor cero(0.0) en ellos.
 *  \param[in] n Número de elementos del vector.
 *  \return Un puntero a la estructura VectorFloat, el vector es de longitud n y 
 *  con elementos con valor 0.0.
 */
VectorFloat *CrearVectorFloat(unsigned short n)
{
	VectorFloat *tmp=NULL;

	tmp=calloc(1, sizeof(VectorFloat));

	tmp->N=n;
	tmp->data=calloc(n, sizeof(float));

	return tmp;
}

/** \fn void LiberaVectorFloat(VectorFloat *tmp)
 *  \brief Destruye el vector de valores reales apuntando por tmp.
 *  \param[in] tmp Puntero al vector de valores reales que se desea destruir.
 */
void LiberaVectorFloat(VectorFloat *tmp)
{
	free(tmp->data);
	free(tmp);
}

/** \fn void IniciaConZeroVectorFloat(VectorFloat *tmp)
 *  \brief Rellena el vector de valores reales con cero.
 *  \param[in] tmp Puntero al vector de valores reales que se desea setear con cero.
 */
void IniciaConZeroVectorFloat(VectorFloat *tmp)
{
	int i;
	for(i=0;i<tmp->N;i++)	tmp->data[i]=-1.0;
}

/** \fn void IniciaConUnoZeroVectorFloat(VectorFloat *tmp)
 *  \brief Rellena el vector de valores reales con -1 y 1 ciclicamente.
 *  \param[in] tmp Puntero al vector de valores reales que se desea setear el valor.
 */
void IniciaConUnoZeroVectorFloat(VectorFloat *tmp)
{
	int i,a;
	tmp->data[0]=0.0;
	a=0;
	for(i=1;i<tmp->N;i++)	
	{
		a=a^1;
		tmp->data[i]=2.0*a-1.0;
	}
}

#endif

