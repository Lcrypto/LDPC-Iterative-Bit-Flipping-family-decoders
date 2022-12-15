/*
 * sparse.h
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

#ifndef __SPARSE_H__
#define __SPARSE_H__


#include <stdlib.h>
/**\struct SparseMatrix
 * David MacKay, Matthew Davey, and John Lafferty have all written low density 
 * parity check matrices in a format called an alist. Esta es una interpretación
 * de ese formato. para ser usado como matriz de paridad P en una LDGM sistemática.
 */
typedef struct{
	unsigned short NFil; /* Número de Filas */
	unsigned short NCol; /* Número de Columnas*/
	unsigned short N;
	unsigned short K;
	unsigned short N_K; /* N-K */
	unsigned short MaxCol; /* límite máximo de unos en todas las columnas.*/
	unsigned short MaxFil; /* límite máximo de unos en todas las filas.*/
	unsigned short *MaxUnosCol; /* máxima cantidad de unos en cada columna.*/
	unsigned short *MaxUnosFil; /* máxima cantidad de unos en cada fila.*/
	unsigned short **UnosCol;
	unsigned short **UnosFil;
}SparseMatrix;

/** \fn int CargaDatos(unsigned short int *k,unsigned short int *n,unsigned short int *n_k,const char archivo[])
 *  \brief Carga los valores n,k y n-k desde los datos de la matriz en archivo.
 *  \param[out] k Número de bits de iinformación.
 *  \param[out] n Número de bits codificados.
 *  \param[out] n_k Diferencia entre n y k.
 *  \param[in] archivo Nombre del archivo donde está guardada la matriz en el formato 
 *  <a href="http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html">Mackay</a>, 
 *  esta matriz se usará de matriz de paridad P en la matriz sistemática LDGM.
 *  \return 0 si todo fue bien o -1 en el caso contrario.
 */
int CargaDatos(unsigned short int *k,unsigned short int *n,unsigned short int *n_k,const char archivo[])
{
        FILE *fd;
	int ii;
        fd=fopen(archivo,"r");
        if(fd==NULL)          return -1;
        
        ii=fscanf(fd,"%hd %hd\n",k,n_k);
        *n=*k+*n_k;
		  
        fclose(fd);
        
        return 0;
}

/** \fn SparseMatrix *CargaMatrizG(const char archivo[])
 *  \brief Carga del archivo de la matriz P, esta matriz sigue el formato
 *  <a href="http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html">Mackay</a>.
 *  \param[in] archivo Nombre del archivo que contiene la descripción de la matriz
 *  que será usada como matriz de paridad P.
 *  \return Retorna un puntero a una estructura SparseMatrix que representa a la
 *  matriz G=[P I].
 */
SparseMatrix *CargaMatrizG(const char archivo[])
{
	int i,j,ii;
	unsigned short tmp;
	FILE *fd=NULL;
	SparseMatrix *g=NULL;

	fd=fopen(archivo,"r");
	if(fd==NULL)          return NULL;
        
	g=calloc(1,sizeof(SparseMatrix));

	////////////////////////////////////////////////////////////////////
	ii=fscanf(fd,"%hd %hd\n" ,&(g->K),&(g->N_K));
	g->N=g->K+g->N_K;
	g->NFil=g->K;
	g->NCol=g->N;
	
	g->MaxUnosCol=calloc(g->NCol,sizeof(unsigned short));
	g->MaxUnosFil=calloc(g->NFil,sizeof(unsigned short));

	////////////////////////////////////////////////////////////////////
	ii=fscanf(fd,"%hd %hd\n" ,&(g->MaxFil),&(g->MaxCol));
	g->MaxFil=g->MaxFil+1;

	////////////////////////////////////////////////////////////////////
	for(i=0;i<g->K;i++)
	{
		ii=fscanf(fd,"%hd",&(g->MaxUnosFil[i]));
		g->MaxUnosFil[i]=g->MaxUnosFil[i]+1;
	}

	////////////////////////////////////////////////////////////////////
	for(i=0;i<g->N_K;i++)		ii=fscanf(fd,"%hd",&(g->MaxUnosCol[i]));
	for(i=g->N_K;i<g->N;i++)	g->MaxUnosCol[i]=1;

	////////////////////////////////////////////////////////////////////
	g->UnosFil=calloc(g->NFil,sizeof(unsigned short*));
	for(i=0;i<g->NFil;i++)	g->UnosFil[i]=calloc(g->MaxUnosFil[i],sizeof(unsigned short));

	for(i=0;i<g->K;i++)
	{
		for(j=0;j<(g->MaxUnosFil[i]-1);j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
			g->UnosFil[i][j]=tmp-1;
		}
		for(j=(g->MaxUnosFil[i]-1);j<(g->MaxFil-1);j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
		}
		g->UnosFil[i][g->MaxUnosFil[i]-1]=i+g->N_K;
	}

	////////////////////////////////////////////////////////////////////                
	g->UnosCol=calloc(g->NCol,sizeof(unsigned short*));
	for(i=0;i<g->NCol;i++)	g->UnosCol[i]=calloc(g->MaxUnosCol[i],sizeof(unsigned short));

	for(i=0;i<g->N_K;i++)
	{
		for(j=0;j<g->MaxUnosCol[i];j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
			g->UnosCol[i][j]=tmp-1;
		}
		for(j=g->MaxUnosCol[i];j<g->MaxCol;j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
		}
	}

	for(i=g->N_K;i<g->N;i++)
	{
		g->UnosCol[i][0]=i-g->N_K;
	}

	////////////////////////////////////////////////////////////////////                

	fclose(fd);
	return g;
}


/** \fn SparseMatrix *CargaMatrizHt(const char archivo[])
 *  \brief Carga del archivo de la matriz P, esta matriz sigue el formato
 *  <a href="http://www.inference.phy.cam.ac.uk/mackay/codes/alist.html">Mackay</a>.
 *  \param[in] archivo Nombre del archivo que contiene la descripción de la matriz
 *  que será usada como matriz de paridad P.
 *  \return Retorna un puntero a una estructura SparseMatrix que representa a la
 *  matriz Htranspuesta=[I P^t]^t.
 */
SparseMatrix *CargaMatrizHt(const char archivo[])
{
	int i,j,ii;
	unsigned short tmp;
	FILE *fd=NULL;
	SparseMatrix *ht=NULL;

	fd=fopen(archivo,"r");
	if(fd==NULL)          return NULL;
        
	ht=calloc(1,sizeof(SparseMatrix));

	////////////////////////////////////////////////////////////////////
	ii=fscanf(fd,"%hd %hd\n" ,&(ht->K),&(ht->N_K));
	ht->N=ht->K+ht->N_K;
	ht->NFil=ht->N;
	ht->NCol=ht->N_K;
	
	ht->MaxUnosCol=calloc(ht->NCol,sizeof(unsigned short));
	ht->MaxUnosFil=calloc(ht->NFil,sizeof(unsigned short));

	////////////////////////////////////////////////////////////////////
	ii=fscanf(fd,"%hd %hd\n" ,&(ht->MaxFil),&(ht->MaxCol));
	ht->MaxCol=ht->MaxCol+1;

	////////////////////////////////////////////////////////////////////
	for(i=0;i<ht->N_K;i++)		ht->MaxUnosFil[i]=1;
	for(i=ht->N_K;i<ht->N;i++)	ii=fscanf(fd,"%hd",&(ht->MaxUnosFil[i]));

	////////////////////////////////////////////////////////////////////
	for(i=0;i<ht->N_K;i++)
	{
		ii=fscanf(fd,"%hd",&(ht->MaxUnosCol[i]));
		ht->MaxUnosCol[i]=ht->MaxUnosCol[i]+1;
	}

	////////////////////////////////////////////////////////////////////                
	ht->UnosFil=calloc(ht->NFil,sizeof(unsigned short*));
	for(i=0;i<ht->NFil;i++)	ht->UnosFil[i]=calloc(ht->MaxUnosFil[i],sizeof(unsigned short));

	for(i=0;i<ht->N_K;i++)
	{
		ht->UnosFil[i][0]=i;
	}

	for(i=ht->N_K;i<ht->N;i++)
	{
		for(j=0;j<ht->MaxUnosFil[i];j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
			ht->UnosFil[i][j]=tmp-1;
		}
		for(j=ht->MaxUnosFil[i];j<ht->MaxFil;j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
		}
	}

	////////////////////////////////////////////////////////////////////
	ht->UnosCol=calloc(ht->NCol,sizeof(unsigned short*));
	for(i=0;i<ht->NCol;i++)	ht->UnosCol[i]=calloc(ht->MaxUnosCol[i],sizeof(unsigned short));

	for(i=0;i<ht->N_K;i++)
	{
		ht->UnosCol[i][0]=i;

		for(j=1;j<ht->MaxUnosCol[i];j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
			ht->UnosCol[i][j]=(tmp-1)+ht->N_K;
		}
		for(j=ht->MaxUnosCol[i];j<ht->MaxCol;j++)
		{
			ii=fscanf(fd,"%hd",&tmp);
		}
	}


	////////////////////////////////////////////////////////////////////                

	fclose(fd);

	return ht;
}

#endif
