/*
 * simula.h
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

#ifndef __SIMULA_H__
#define __SIMULA_H__


#include <math.h>

#include "vector.h"
#include "sparse.h"
#include "mirandom.h"


/** \fn void multiplica(Vector *v,const Vector *u,const SparseMatrix *g)
 *  \brief Multiplica el vector u por la matriz g y lo carga en v.
 *  \param[out] v Vector codificado.
 *  \param[in] u Vector de información.
 *  \param[in] g Matriz generadora.
 */
void multiplica(Vector *v,const Vector *u,const SparseMatrix *g)
{
	int i,col,s=0;
	for(col=0;col<g->NCol;col++)
	{
		s=0;
		for(i=0;i<g->MaxUnosCol[col];i++) 
		{
			s=s^u->data[g->UnosCol[col][i]];
		}
		v->data[col]=s;
	}
}


/** \fn void canal(Vector *r,const Vector *v,const float p)
 *  \brief Pasa el vector v por un canal BSC de probabilidad p y te devuelve r.
 *  \param[out] r Vector codificado ruidoso.
 *  \param[in] v Vector codificado.
 *  \param[in] p Probabilidad de error del canal.
 */
void canal(Vector *r,const Vector *v,const float p)
{
	float disp;
	int	i;
	for(i=0;i<r->N;i++)
	{
		disp=((myrand()*1.0)/MY_RAND_MAX);
		if(disp>p)	r->data[i]=v->data[i];
		else		r->data[i]=v->data[i]^1;
	}
}

/** \fn int HardDecision(const float x,const float limear)
 *  \brief Aplica una decisión abrupta al valor x.
 *  \param[in] x valor al que se desea aplicar el limear.
 *  \param[in] limear todos los valores mayores que el limear devuelve un valor 1
 *  y cuando son menores devuelve el valor 0.
 *  \return Retorna la decisión abrupta de x.
 */
int HardDecision(const float x,const float limear)
{
    if(x>limear) return 1;
    else return 0;
}

/** \fn void canalGausiano(VectorFloat *r,const Vector *v,const float SigmaN,const float a)
 *  \brief Recibe el vector binario{0,1} v lo transforma en real{-1,1} y lo pasa 
 *  por un canal con ruido AWGN con un desvío padrón SigmaN.
 *  \param[out] r Vector a la salida del canal AWGN.
 *  \param[in] v Vector binario en la entrada del canal AWGN.
 *  \param[in] SigmaN Desvío padrón del canal AWGN.
 *  \param[in] a Amplitud de la señal de entrada en el canal AWGN.
 */
void canalGausiano(VectorFloat *r,const Vector *v,const float SigmaN,const float a)
{
	float disp;
	int	i;
	for(i=0;i<r->N;i++)
	{
		r->data[i]=a*(2.0*v->data[i]-1.0)+ SigmaN*MyRandGauss();
		
	}
}

/** \fn void imprime_vetor(const Vector *a)
 *  \brief Imprime en pantalla el vector a.
 *  \param[in] a Vector que deseo imprimir en pantalla.
 */
void imprime_vetor(const Vector *a)
{
	int	i;
	for(i=0;i<a->N;i++)
	{
		printf("%d\t",a->data[i]);
	}
}

/** \fn int encontre_error(const Vector *a)
 *  \brief busca errores en un vector, es decir valores distinto de cero
 *  \param[in] a Vector que deseo analizar.
 *  \return Retorna la cantidad de errores, es decir la cantidad de elementos
 *  distintos de cero.
 */
int encontre_error(const Vector *a)
{
	int	i,s=0;
	for(i=0;i<a->N;i++)
	{
		if(a->data[i]!=0) s=s+1;
	}
	return s;
}

/** \fn void SH_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo SHBF al vector binario r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector obtenido por decisión abrupta en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void SH_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
	
	int menor,IDmenor;
	FILE *tmp=NULL;

	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);

	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];
	
do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+(-2*S->data[ht->UnosFil[fil][j]]+1);
			}
			F->data[fil]=s;
		}


		// Encuentra el bit con mas erro

		menor=F->data[ht->N-1];
		IDmenor=ht->N-1;
		for(j=1;j<ht->N;j++)
		{
			if(menor>F->data[ht->N-1-j])	
			{
				menor=F->data[ht->N-1-j];
				IDmenor=ht->N-1-j;
			}
		}

		//Cambio un bit con menor peso
		r->data[IDmenor]=r->data[IDmenor]^1;
		
		k++;
	}
	
}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

/** \fn void PH_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo PHBF al vector binario r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector obtenido por decisión abrupta en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void PH_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
		
	int menor,IDmenor;
	FILE *tmp=NULL;
	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);
	
	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];

do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+(-2*S->data[ht->UnosFil[fil][j]]+1);
			}
			F->data[fil]=s;
		}

		// Encuentra el bit con mas erro
		menor=F->data[0];
		for(j=1;j<ht->N;j++)
		{
			if(menor>F->data[j])	
			{
				menor=F->data[j];
			}
		}

        //cambio todos los bit con menor peso
		for(j=0;j<ht->N;j++)
		{
			if(menor==F->data[j])	
			{
				r->data[j]=r->data[j]^1;
			}
		}
		k++;
	}

}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];
	

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

/** \fn void bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo BF al vector binario r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector obtenido por decisión abrupta en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int i,j,s=0,col,fil,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
			
	int mayor;
	FILE *tmp=NULL;
	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);
	
	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];

do{
	//printf("BF %7d\r",k);
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+S->data[ht->UnosFil[fil][j]];
			}
			F->data[fil]=s;
		}
		
		// Encuentra el bit con mas erro

		mayor=F->data[0];
		for(j=1;j<ht->N;j++)
		{
			if(mayor<F->data[j])	
			{
				mayor=F->data[j];
			}
		}

		//cambio todos los bit de mayor peso
		for(j=0;j<ht->N;j++)
		{
			if(mayor==F->data[j])	
			{
				r->data[j]=r->data[j]^1;
			}
		}
		k++;	
	}

}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

/** \fn void bitflipping_serial(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo BF serialmente al vector binario r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector obtenido por decisión abrupta en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void bitflipping_serial(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int IDmayor,i,j,s=0,col,fil,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
			
	int mayor;
	FILE *tmp=NULL;
	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);
	
	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];

do{
	//printf("BF %7d\r",k);
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+S->data[ht->UnosFil[fil][j]];
			}
			F->data[fil]=s;
		}
		
		// Encuentra el bit con mas erro

		mayor=F->data[0];
		for(j=1;j<ht->N;j++)
		{
			if(mayor<F->data[j])	
			{
				mayor=F->data[j];
				IDmayor=j;
			}
		}

		//Cambio un bit con mayor peso
		r->data[IDmayor]=r->data[IDmayor]^1;
		
		k++;
	}

}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

/** \fn int compara(const Vector *a,const Vector *b)
 *  \brief Compara dos vectores binarios o enteros.
 *  \param[in] a Primer vector a comparar.
 *  \param[in] b Segundo vector a comparar.
 *  \return El número de elementos distintos.
 */
int compara(const Vector *a,const Vector *b)
{
	int	i,s=0;
	for(i=0;i<a->N;i++)
	{
		if(a->data[i]!=b->data[i]) s++;
	}
	return s;
}

/** \fn int comparaSoftHard(const VectorFloat *a,const Vector *b)
 *  \brief Compara un vector real al que se le aplica una decisión abrupta con 
 *  un vector binario.
 *  \param[in] a Primer vector real a comparar.
 *  \param[in] b Segundo vector binario a comparar.
 *  \return El número de elementos distintos.
 */
int comparaSoftHard(const VectorFloat *a,const Vector *b)
{
	int	i,s=0;
	for(i=0;i<a->N;i++)
	{
		if(HardDecision(a->data[i],0)!=b->data[i]) s++;
	}
	return s;
}


/** \fn void GeneraVectorHard(Vector *b,const VectorFloat *a)
 *  \brief Genera una vector binario mediante una decisión abrupta con limear en 
 *  cero de un vector real.
 *  \param[out] b Vector binario donde se cargara la decisión abrupta.
 *  \param[in]  a Vector real que se desea convertir.
 */
void GeneraVectorHard(Vector *b,const VectorFloat *a)
{
	int	i;
	for(i=0;i<a->N;i++)
	{
		b->data[i]=HardDecision(a->data[i],0);
	}
}

/** \fn unsigned short fuente(void)
 *  \brief Genera aleatoria 1 y 0 con probabilidad 0.5.
 *  \return Un valor 1 o 0 con probabilidad 0.5.
 */
unsigned short fuente(void)
{
	if(rand()>(RAND_MAX/2))	return 1;
	else			return 0;
}

/** \fn float ProbabilidadPaquete(float p,unsigned short n)
 *  \brief Probabilidad de error de un vector a partir de la probabilidad de 
 *   error de un bit del vector.
 *  \param[in] p Probabilidad de error de bit en el vector codificado.
 *  \param[in] n Número de elementos del vector codificado.
 *  \return Retorna la probabilidad de error del paquete, osea del vector.
 */
float ProbabilidadPaquete(float p,unsigned short n)
{
	return (1.0-pow(1.0-p,n));
}

/** \fn float Qfunc(float x)
 *  \brief Evalúa la función Q(x) osea 0.5*erfc(x/sqrt(2)).
 *  \param[in] x Valor de entrada de la función Q(x).
 *  \return retorna Q(x).
 */
float Qfunc(float x)
{
      float s=0,dt,t,f,fi,fo;
      int i;
      
      if(x<0) x=0;      
      t=x;
      if(x>3) dt=0.01171875+x/(256.0);
      else
      {
          if(x>1) dt=x/(128.0);
          else    dt=(0.0078125);
      }
      
      f=exp(-t*t/2.0);
      fo=f;
      fi=fo;
      
      t=t+dt;
      for(i=0;fi>=(fo/512.0);i++)
      {
            fi=exp(-t*t/2.0);
            s=s+(fi+f);
            f=fi;
            t=t+dt;
      }
      return (s/2.0)*dt/sqrt(2.0*M_PI);
}

/** \fn void W_bitflipping(Vector *ub,const VectorFloat *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo WBF serial-mente al vector real r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector real en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void W_bitflipping(Vector *ub,const VectorFloat *r1,const SparseMatrix *ht,const float p)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	VectorFloat *E=NULL;
	VectorFloat *Ymin=NULL;
	Vector *r=NULL;
	
	int IDmayor;
	float menor,mayor,Ej;
	FILE *tmp=NULL;

	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	Ymin=CrearVectorFloat(ht->N_K);
	E=CrearVectorFloat(ht->N);
	r=CrearVector(ht->N);

	for(i=0;i<ht->N;i++) r->data[i]=HardDecision(r1->data[i],0);
	
	// Genera Ymin Vector
	for(col=0;col<ht->N_K;col++)
	{
		menor=fabs(r1->data[ht->UnosCol[col][0]]);
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			if(fabs(r1->data[ht->UnosCol[col][i]])<menor)
			menor=fabs(r1->data[ht->UnosCol[col][i]]);
		}
		Ymin->data[col]=menor;
	}
    	
do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos desde  el Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			Ej=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    Ej= Ej+(2.0*S->data[ht->UnosFil[fil][j]]-1.0)*Ymin->data[ht->UnosFil[fil][j]];
			}
			E->data[fil]=Ej;
		}


		// Encuentra el bit con mas peso
		mayor=E->data[ht->N-1];
		IDmayor=ht->N-1;
		for(j=1;j<ht->N;j++)
		{
			if(mayor<E->data[ht->N-1-j])	
			{
				mayor=E->data[ht->N-1-j];
				IDmayor=ht->N-1-j;
			}
		}

		//Cambio un bit con menor peso
		r->data[IDmayor]=r->data[IDmayor]^1;
		
		k++;
	}
	
}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=HardDecision(r1->data[j+ht->N_K],0);
	else			        for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K],0;

	LiberaVector(S);
	LiberaVectorFloat(E);
	LiberaVectorFloat(Ymin);
	LiberaVector(r);
}

/** \fn void MW_bitflipping(Vector *ub,const VectorFloat *r1,const SparseMatrix *ht,const float p,const float ALPHA)
 *  \brief Aplica el algoritmo MWBF serial-mente al vector real r1, usando la matriz h^t
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector real en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 *  \param[in] ALPHA valor de regularización.
 */
void MW_bitflipping(Vector *ub,const VectorFloat *r1,const SparseMatrix *ht,const float p,const float ALPHA)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
//	float ALPHA=0.0;
	
	Vector *S=NULL;
	VectorFloat *E=NULL;
	VectorFloat *Ymin=NULL;
	Vector *r=NULL;
	
	int IDmayor;
	float menor,mayor,Ej;
	FILE *tmp=NULL;

	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	Ymin=CrearVectorFloat(ht->N_K);
	E=CrearVectorFloat(ht->N);
	r=CrearVector(ht->N);

	for(i=0;i<ht->N;i++) r->data[i]=HardDecision(r1->data[i],0);
	
	// Genera Ymin Vector
	for(col=0;col<ht->N_K;col++)
	{
		menor=fabs(r1->data[ht->UnosCol[col][0]]);
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			if(fabs(r1->data[ht->UnosCol[col][i]])<menor)
			menor=fabs(r1->data[ht->UnosCol[col][i]]);
		}
		Ymin->data[col]=menor;
	}
    	
do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos desde  el Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			Ej=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    Ej= Ej+(2.0*S->data[ht->UnosFil[fil][j]]-1.0)*Ymin->data[ht->UnosFil[fil][j]];
			}
			E->data[fil]=Ej-ALPHA*fabs(r1->data[fil]);
		}


		// Encuentra el bit con mas peso
		mayor=E->data[ht->N-1];
		IDmayor=ht->N-1;
		for(j=1;j<ht->N;j++)
		{
			if(mayor<E->data[ht->N-1-j])	
			{
				mayor=E->data[ht->N-1-j];
				IDmayor=ht->N-1-j;
			}
		}

		//Cambio un bit con menor peso
		r->data[IDmayor]=r->data[IDmayor]^1;
		
		k++;
	}
	
}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=HardDecision(r1->data[j+ht->N_K],0);
	else			        for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K],0;

	LiberaVector(S);
	LiberaVectorFloat(E);
	LiberaVectorFloat(Ymin);
	LiberaVector(r);
}

/** \fn float MenorDiferenteDeCero(float Probabilidades[],int N)
 *  \brief Retorna la menor probabilidad en el vector de Probabilidades pero 
 *  diferente de cero.
 *  \param[in] Probabilidades Vector de probabilidades a analizar.
 *  \param[in] N Número de elementos del vector de probabilidades.
 *  \return La menor probabilidad distinta de cero del vector de probabilidades.
 */
float MenorDiferenteDeCero(float Probabilidades[],int N)
{
	float menor;
	int i;

	menor=Probabilidades[0];
      
	for(i=1;i<N;i++)
	{
		if( (Probabilidades[i]!=0)  && (Probabilidades[i]<menor) )
		menor=Probabilidades[i]; 
	}

	return menor;
}

/** \fn void P_SipSpi_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo SipSpi_BF_paralelo al vector binario r1, usando 
 *  la matriz h^t.
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector binario en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void P_SipSpi_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
	
	int menor,IDmenor;
	FILE *tmp=NULL;

	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);

	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];
	
do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+(-2*S->data[ht->UnosFil[fil][j]]+1);
			}
			F->data[fil]=s;
		}


		//Cambio los bit con mayores bits de verificacion de paridad fallida

		for(j=0;j<ht->N;j++)
		{
			if(F->data[j]<0)	
			{
				r->data[j]=r->data[j]^1;
			}
		}
		
		k++;
	}
	
}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

/** \fn void S_SipSpi_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
 *  \brief Aplica el algoritmo SipSpi_BF_serial al vector binario r1, usando 
 *  la matriz h^t.
 *  \param[out] ub Vector de información estimado.
 *  \param[in] r1 Vector binario en la salida del canal.
 *  \param[in] ht Matriz de verificación de paridad.
 *  \param[in] p Porcentaje de cantidad de iteraciones máximas con respecto a N.
 */
void S_SipSpi_bitflipping(Vector *ub,const Vector *r1,const SparseMatrix *ht,const float p)
{
	int i,j,col,fil,s=0,errores=0,MAX_INTENTOS,k=0,ErroresSindrome=0;
	
	Vector *S=NULL;
	Vector *F=NULL;
	Vector *r=NULL;
	
	int menor,IDmenor;
	FILE *tmp=NULL;

	MAX_INTENTOS=(int)(p*ht->N);

	S=CrearVector(ht->N_K);
	F=CrearVector(ht->N);
	r=CrearVector(ht->N);

	for(i=0;i<ht->N;i++) r->data[i]=r1->data[i];
	
do{
	// Genera Sindrome
	ErroresSindrome=0;
	for(col=0;col<ht->N_K;col++)
	{
		s=0;
		for(i=0;i<ht->MaxUnosCol[col];i++) 
		{
			s=s^(r->data[ht->UnosCol[col][i]]);
		}
		S->data[col]=s;
		ErroresSindrome=ErroresSindrome+s;
	}
	
	if(ErroresSindrome>0)
	{
		//printf("k:%7d S:%7d N:%7d\n",k,ErroresSindrome,S->N);
		errores=1;

		// Genera Pesos del Sindrome
		for(fil=0;fil<ht->N;fil++)
		{
			s=0;
			for(j=0;j<ht->MaxUnosFil[fil];j++) 
			{
			    s= s+(-2*S->data[ht->UnosFil[fil][j]]+1);
			}
			F->data[fil]=s;
		}


		//Cambio el primer bit con mayor cantidad de bits de verificación de paridad fallida

		for(j=0;j<ht->N;j++)
		{
			if(F->data[ht->N-1-j]<0)	
			{
				r->data[ht->N-1-j]=r->data[ht->N-1-j]^1;
				j=ht->N;
			}
		}
		
		k++;
	}
	
}while((ErroresSindrome>0)&&(k<MAX_INTENTOS));

	if(ErroresSindrome>0)	for(j=0;j<ht->K;j++)	ub->data[j]=r1->data[j+ht->N_K];
	else			for(j=0;j<ht->K;j++)	ub->data[j]=r->data[j+ht->N_K];

	LiberaVector(S);
	LiberaVector(F);
	LiberaVector(r);
}

#endif
