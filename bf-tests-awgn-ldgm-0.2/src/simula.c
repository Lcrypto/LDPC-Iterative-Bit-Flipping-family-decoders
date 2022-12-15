//gcc simula.c -o simula -lm -Os
#include <stdio.h>
#include <stdlib.h>
/*
short G[K][N]={		{1,0,1,1,0,0,  1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
		        {0,0,1,1,0,1,  0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{0,1,1,1,0,0,  0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0},
			{1,1,0,0,1,0,  0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0},
		        {1,1,1,0,0,0,  0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0},
		        {0,1,1,0,0,1,  0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0},
		        {0,1,1,0,1,0,  0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0},
	         	{1,0,0,1,1,0,  0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0},
		        {0,1,0,1,1,0,  0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0},
		        {0,0,1,1,1,0,  0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0},
		        {1,0,1,0,0,1,  0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0},
		        {1,0,0,1,0,1,  0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0},
		        {1,1,0,0,0,1,  0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0},
		        {1,0,0,0,1,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0},
		        {0,1,0,0,1,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0},
		        {0,0,0,1,1,1,  0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1} };

short Ht[N][N_K]={	{1,0,0,0,0,0},
			{0,1,0,0,0,0},
			{0,0,1,0,0,0},
			{0,0,0,1,0,0},
			{0,0,0,0,1,0},
			{0,0,0,0,0,1},

			{1,0,1,1,0,0},
			{0,0,1,1,0,1},
			{0,1,1,1,0,0},
			{1,1,0,0,1,0},
			{1,1,1,0,0,0},
			{0,1,1,0,0,1},
			{0,1,1,0,1,0},
			{1,0,0,1,1,0},
			{0,1,0,1,1,0},
			{0,0,1,1,1,0},
			{1,0,1,0,0,1},
			{1,0,0,1,0,1},
			{1,1,0,0,0,1},
			{1,0,0,0,1,1},
			{0,1,0,0,1,1},
			{0,0,0,1,1,1} };

short V[N]=	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
short U[K]=	{1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};

short R[N]=	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
short Ub2[K]=	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
short Ubf[K]=	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
short Ur[K]=	{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
*/

#include <config.h>
#include "extras/simula.h"
#include "extras/extras.h"

#define NUMERO_ALGS 9
int main(int argc,char *argv[])
{
	SparseMatrix *G=NULL;
	SparseMatrix *Ht=NULL;
	DATO *Dato=NULL;
	Vector *V=NULL;
	Vector *U=NULL;
	VectorFloat *R=NULL;
	Vector *Rhard=NULL;

	Vector *Ur=NULL;
	Vector *Ubf=NULL;
	Vector *Ubfs=NULL;
	Vector *Ushbf=NULL;
	Vector *Uphbf=NULL;
	Vector *Usssbf=NULL;
	Vector *Upssbf=NULL;
	Vector *Uwbf=NULL;
	Vector *Umwbf=NULL;

	unsigned short K,N,N_K;
	unsigned long int i,k,j;

	unsigned long int ITERACIONES;

	unsigned long int ERRORES=0;

	unsigned long int ERRORES_BIT=0;
	unsigned long int ERRORES_BIT_BF=0;
	unsigned long int ERRORES_BIT_BFS=0;
	unsigned long int ERRORES_BIT_SHBF=0;
	unsigned long int ERRORES_BIT_PHBF=0;
	unsigned long int ERRORES_BIT_SSSBF=0;
	unsigned long int ERRORES_BIT_PSSBF=0;
	unsigned long int ERRORES_BIT_WBF=0;
	unsigned long int ERRORES_BIT_MWBF=0;

	int ACTIVA_BF=1;
	int ACTIVA_BFS=1;
	int ACTIVA_SHBF=1;
	int ACTIVA_PHBF=1;
	int ACTIVA_SSSBF=1;
	int ACTIVA_PSSBF=1;
	int ACTIVA_WBF=1;
	int ACTIVA_MWBF=1;

	float BER;
	float BERbf;
	float BERbfs;
	float BERshbf;
	float BERphbf;
	float BERsssbf;
	float BERpssbf;
	float BERwbf;
	float BERmwbf;

	int m=0;
	float P=0.5,BERmin,BERminant,BERminact,sigmaN,SNR,SNRdb,dSNRdb,a=1.0; 
	float PpackTeo,PpackSim;

	float Probabilidades[NUMERO_ALGS],Rc;
	FILE *fd;
	FILE *tmp;

	Dato=Menu(argc,argv);
    
	if(Dato->Archivo==NULL) 
	{
		printf("Falto indicar el archivo de entrada en formato AList.\n");
		ImprimirAyuda();
		return EXIT_SUCCESS;
	}
	if(fpr_exist_file(Dato->Archivo)==0) 
	{
		printf("El archivo en formato AList indicado no existe.\n");
		return EXIT_SUCCESS;
	}
	if(Dato->ArchivoLog==NULL) 
	{
		Dato->ArchivoLog=(char*)calloc(1,strlen(Dato->Archivo)+8);
		sprintf(Dato->ArchivoLog,"%s.log",Dato->Archivo);
	}
	if(Dato->ArchivoOut==NULL) 
	{
		Dato->ArchivoOut=(char*)calloc(1,strlen(Dato->Archivo)+8);
		sprintf(Dato->ArchivoOut,"%s.bf",Dato->Archivo);
	}

	CargaDatos(&K,&N,&N_K,Dato->Archivo);
	G=CargaMatrizG(Dato->Archivo);
	Ht=CargaMatrizHt(Dato->Archivo);

	V=CrearVector(N);	IniciaConZeroVector(V);
	U=CrearVector(K);	IniciaConUnoZeroVector(U);
	R=CrearVectorFloat(N);	IniciaConZeroVectorFloat(R);
	Rhard=CrearVector(N);	IniciaConZeroVector(Rhard);

	Ur =CrearVector(K);	IniciaConZeroVector(Ur);
	Ubf=CrearVector(K);	IniciaConZeroVector(Ubf);
	Ubfs=CrearVector(K);	IniciaConZeroVector(Ubfs);
	Ushbf=CrearVector(K);	IniciaConZeroVector(Ushbf);
	Uphbf=CrearVector(K);	IniciaConZeroVector(Uphbf);
	Usssbf=CrearVector(K);	IniciaConZeroVector(Usssbf);
	Upssbf=CrearVector(K);	IniciaConZeroVector(Upssbf);
	Uwbf=CrearVector(K);	IniciaConZeroVector(Uwbf);
	Umwbf=CrearVector(K);	IniciaConZeroVector(Umwbf);

	tmp=fopen(Dato->ArchivoLog,"a");
	if(tmp==NULL) 
	{
		printf("Error abriendo/creando fichero: %s\n",Dato->ArchivoLog);
		return EXIT_SUCCESS;
	}
        	
	printf("\nProb_min=1/MY_RAND_MAX=%e\tMY_RAND_MAX=%e\n",1.0/MY_RAND_MAX,1.0*MY_RAND_MAX);
	fprintf(tmp,"\nProb_min=1/MY_RAND_MAX=%e\tMY_RAND_MAX=%e\n",1.0/MY_RAND_MAX,1.0*MY_RAND_MAX);


	printf("\n");
	printf("BERmin=%e\n",Dato->BERmin);
	printf("PackMax=%ld\n",Dato->PackMax);
	printf("BERminConfiable:%e\n",100.0/(N*Dato->PackMax));		
	printf("Si BERminConfiable> BERmin entonces incremente PackMax!!!!\n");
	printf("PackMax es un (long int) de %d bits,",8*sizeof(Dato->PackMax));
	printf("su valor m�ximo es de:%ld\n",(long int)pow(2.0,8.0*sizeof(Dato->PackMax)-1));
	fprintf(tmp,"\n");	
	fprintf(tmp,"BERmin=%e\n",Dato->BERmin);	
	fprintf(tmp,"PackMax=%ld\n",Dato->PackMax);
	fprintf(tmp,"BERminConfiable:%e\n",100.0/(N*Dato->PackMax));	
	fprintf(tmp,"Si BERminConfiable> BERmin entonces incremente PackMax!!!!");	
	fprintf(tmp,"PackMax es un (long int) de %d bits,",8*sizeof(Dato->PackMax));
	fprintf(tmp,"su valor m�ximo es de:%ld\n",(long int)pow(2.0,8.0*sizeof(Dato->PackMax)-1));

	printf("Presiona enter para continuar...\n");getchar();
    
	multiplica(V,U,G);

	printf("\n");fprintf(tmp,"\n");

	//SNR=Eb/No
	SNRdb=0.0;
	dSNRdb=0.3333333333333;

	SNR=pow(10,SNRdb/10);
	Rc=(1.0*Ht->K)/Ht->N;
	P=Qfunc(sqrt(SNR*2.0*Rc));

	printf("K=%d\nN=%d\nRc=%f\n\n",Ht->K,Ht->N,Rc);
	fprintf(tmp,"K=%d\nN=%d\nRc=%f\n\n",Ht->K,Ht->N,Rc);

	fclose(tmp);
    
	crear_grafico_teorico_m();
	crear_LimiteDeShannon_m();
	crear_simula1_m(Dato,Ht->K,Ht->N,Ht->MaxFil);

	BERmin=P;
	BERminact=P/0.98;
	for(k=0;SNRdb<=10;k++)
	{
		ERRORES=0;

		ERRORES_BIT=0;
		ERRORES_BIT_BF=0;
		ERRORES_BIT_BFS=0;
		ERRORES_BIT_SHBF=0;
		ERRORES_BIT_PHBF=0;
		ERRORES_BIT_SSSBF=0;
		ERRORES_BIT_PSSBF=0;
		ERRORES_BIT_WBF=0;
		ERRORES_BIT_MWBF=0;

		
		tmp=fopen(Dato->ArchivoLog,"a");
		if(tmp==NULL) 
		{
			printf("Error abriendo/creando fichero: %s\n",Dato->ArchivoLog);
			return 0;
		}
		
		fd=fopen(Dato->ArchivoOut,"a");
		if(fd==NULL) 
		{
			printf("Error creando fichero: %s\n",Dato->ArchivoOut);
			fprintf(tmp,"Error creando fichero: %s\n",Dato->ArchivoOut);
			return 0;
		}
		
		ITERACIONES=(long int)(16384.0/(BERmin*N));

		if((ITERACIONES>=Dato->PackMax)&&(BER!=0))ITERACIONES=Dato->PackMax;

		if(ITERACIONES==0) ITERACIONES=64;
		else                 
		{
			if(ITERACIONES<=512)ITERACIONES=512;
		}
				
		printf("\n\nPACOTES ENVIADOS:%e\n",ITERACIONES*1.0);
		printf("{Eb/N0}_db:%e\n",SNRdb);
		printf("P:%e\n",P);
		printf("P-pacote:%e\n",ProbabilidadPaquete(P,N));
		printf("BERminConfiable:%e\n",100.0/(N*Dato->PackMax));

		fprintf(tmp,"\n\nPACOTES ENVIADOS:%e\n",ITERACIONES*1.0);
		fprintf(tmp,"{Eb/N0}_db:%e\n",SNRdb);
		fprintf(tmp,"P:%e\n",P);
		fprintf(tmp,"P-pacote:%e\n",ProbabilidadPaquete(P,N));
		fprintf(tmp,"BERminConfiable:%e\n",100.0/(N*Dato->PackMax));

		sigmaN=a/sqrt(Rc*SNR*2.0);
		for(i=0;i<ITERACIONES;i++)
		{
			canalGausiano(R,V,sigmaN,a);
		
			if((i%(long int)(1+ITERACIONES/100))==0) printf("%6.2f %c\n",i*100.0/ITERACIONES,'%');

			GeneraVectorHard(Rhard,R);

			// UNCODED
			if(compara(Rhard,V)!=0)	ERRORES=ERRORES+1;
			/* Copia los bits de informacion de R a Ur. */
			for(j=0;j<K;j++) Ur->data[j]=Rhard->data[j+N_K];
			/* Compara Ur con el U enviado. */
			if((m=compara(U,Ur))!=0)	ERRORES_BIT=ERRORES_BIT+m;
			
            


			// Bit-Flipping
			if( (Dato->BF==1)&&(ACTIVA_BF==1) )
			{
				bitflipping(Ubf,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Ubf))!=0)
				{
					ERRORES_BIT_BF=ERRORES_BIT_BF+m;
				}
			}

			// Bit-Flipping serial
			if( (Dato->BFS==1)&&(ACTIVA_BFS==1) )
			{
				bitflipping_serial(Ubfs,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Ubfs))!=0)
				{
					ERRORES_BIT_BFS=ERRORES_BIT_BFS+m;
				}
			}

			// Serial Hard Bit-Flipping
			if( (Dato->SHBF==1)&&(ACTIVA_SHBF==1) )
			{
				SH_bitflipping(Ushbf,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Ushbf))!=0)
				{
					ERRORES_BIT_SHBF=ERRORES_BIT_SHBF+m;
				}
			}

			// Parallel Hard Bit-Flipping
			if( (Dato->PHBF==1)&&(ACTIVA_PHBF==1) )
			{
				PH_bitflipping(Uphbf,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Uphbf))!=0)
				{
					ERRORES_BIT_PHBF=ERRORES_BIT_PHBF+m;
				}
			}

			// Serial Sipser Spilman Bit-Flipping
			if( (Dato->SSSBF==1)&&(ACTIVA_SSSBF==1) )
			{
				S_SipSpi_bitflipping(Usssbf,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Usssbf))!=0) 	
				{
					ERRORES_BIT_SSSBF=ERRORES_BIT_SSSBF+m;
				}
			}

			// Parallel Sipser Spilman Bit-Flipping
			if( (Dato->PSSBF==1)&&(ACTIVA_PSSBF==1) )
			{
				P_SipSpi_bitflipping(Upssbf,Rhard,Ht,Dato->ItPor);
				if((m=compara(U,Upssbf))!=0) 	
				{
					ERRORES_BIT_PSSBF=ERRORES_BIT_PSSBF+m;
				}
			}

			// Weighted Bit-Flipping
			if( (Dato->WBF==1)&&(ACTIVA_WBF==1) )
			{
				W_bitflipping(Uwbf,R,Ht,Dato->ItPor);
				if((m=compara(U,Uwbf))!=0) 	
				{
					ERRORES_BIT_WBF=ERRORES_BIT_WBF+m;
				}
			}

			// Modified Weighted Bit-Flipping
			if( (Dato->MWBF==1)&&(ACTIVA_MWBF==1) )
			{
				MW_bitflipping(Umwbf,R,Ht,Dato->ItPor,Dato->Alpha);
				if((m=compara(U,Umwbf))!=0) 	
				{
					ERRORES_BIT_MWBF=ERRORES_BIT_MWBF+m;
				}
			}

		}
		
		printf("ERRORES-PACK: %10ld\tITERACIONES: %10ld\tProb: %e \n",ERRORES,ITERACIONES,(ERRORES*1.0)/ITERACIONES);
		fprintf(tmp,"ERRORES-PACK: %10ld\tITERACIONES: %10ld\tProb: %e \n",ERRORES,ITERACIONES,(ERRORES*1.0)/ITERACIONES);

		PpackTeo =ProbabilidadPaquete(P,N);		//prob E paquete teorico.
		PpackSim =(ERRORES*1.0)/ITERACIONES;		//prob E paquete simulado.

		BER      =ERRORES_BIT/(K*1.0*ITERACIONES);	//BER.
		BERbf    =ERRORES_BIT_BF/(K*1.0*ITERACIONES);	//BER despues de BF.
		BERbfs   =ERRORES_BIT_BFS/(K*1.0*ITERACIONES);	//BER despues de BFS.
		BERshbf  =ERRORES_BIT_SHBF/(K*1.0*ITERACIONES);	//BER despues de SHBF.
		BERphbf  =ERRORES_BIT_PHBF/(K*1.0*ITERACIONES);	//BER despues de PHBF.
		BERsssbf =ERRORES_BIT_SSSBF/(K*1.0*ITERACIONES);//BER despues de SSSBF.
		BERpssbf =ERRORES_BIT_PSSBF/(K*1.0*ITERACIONES);//BER despues de PSSBF.
		BERwbf   =ERRORES_BIT_WBF/(K*1.0*ITERACIONES);	//BER despues de WBF.
		BERmwbf  =ERRORES_BIT_MWBF/(K*1.0*ITERACIONES);	//BER despues de MWBF.
                                                        
		fprintf(fd,"%e\t%e\t%e\t%e\t"	,P		//P canal
						,BER		//BER
						,PpackTeo	//prob E paquete teorico.
						,PpackSim	//prob E paquete simulado.
						);
		if(Dato->BF==1)	  fprintf(fd,"%e\t",BERbf);	//BER después de BF
		if(Dato->BFS==1)  fprintf(fd,"%e\t",BERbfs);	//BER después de BFS
		if(Dato->SHBF==1) fprintf(fd,"%e\t",BERshbf);	//BER después de SHBF
		if(Dato->PHBF==1) fprintf(fd,"%e\t",BERphbf);	//BER después de PHBF
		if(Dato->SSSBF==1)fprintf(fd,"%e\t",BERsssbf);	//BER después de SSSBF
		if(Dato->PSSBF==1)fprintf(fd,"%e\t",BERpssbf);	//BER después de PSSBF
		if(Dato->WBF==1)  fprintf(fd,"%e\t",BERwbf);	//BER después de WBF
		if(Dato->MWBF==1) fprintf(fd,"%e\t",BERmwbf);	//BER después de MWBF
                                                        
		fprintf(fd,"\n");

		Probabilidades[0]=BER;
		Probabilidades[1]=BERbf;
		Probabilidades[2]=BERbfs;
	        Probabilidades[3]=BERshbf;
		Probabilidades[4]=BERphbf;
		Probabilidades[5]=BERsssbf;
		Probabilidades[6]=BERpssbf;
		Probabilidades[7]=BERwbf;
		Probabilidades[8]=BERmwbf;
        
		if(BERbf<=Dato->BERmin)		ACTIVA_BF=0;
		if(BERbfs<=Dato->BERmin)	ACTIVA_BFS=0;
		if(BERshbf<=Dato->BERmin)	ACTIVA_SHBF=0;
		if(BERphbf<=Dato->BERmin)	ACTIVA_PHBF=0;
		if(BERsssbf<=Dato->BERmin)	ACTIVA_SSSBF=0;
		if(BERpssbf<=Dato->BERmin)	ACTIVA_PSSBF=0;
		if(BERwbf<=Dato->BERmin)	ACTIVA_WBF=0;
		if(BERmwbf<=Dato->BERmin)	ACTIVA_MWBF=0;

		BERminant=BERminact;
		BERminact=MenorDiferenteDeCero(Probabilidades,NUMERO_ALGS);

		printf("BERmin:%e\n",BERminact);
		fprintf(tmp,"BERmin:%e\n",BERminact);

		if(BERminact<=BERminant)	BERmin=BERminact*BERminact/BERminant;
		else			BERmin=BERminact/10;
		printf("BERmin Estimado siguiente:%e\n",BERmin);
		fprintf(tmp,"BERmin Estimado siguiente:%e\n",BERmin);

		printf("BERminConfiable:%e\n",100.0/(N*Dato->PackMax));		
		fprintf(tmp,"BERminConfiable:%e\n",100.0/(N*Dato->PackMax));

		SNRdb=SNRdb+dSNRdb;
		SNR=pow(10,SNRdb/10);
		P=Qfunc(sqrt(SNR*2.0*Rc));

		fclose(fd);
		fclose(tmp);
	}

	
	return 0;
}

