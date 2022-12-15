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
	DATO *Dato=NULL;
	Vector *V=NULL;
	Vector *U=NULL;
	VectorFloat *R=NULL;
	Vector *Rhard=NULL;

	Vector *Vr=NULL;
	Vector *Vbf=NULL;
	Vector *Vbfs=NULL;
	Vector *Vshbf=NULL;
	Vector *Vphbf=NULL;
	Vector *Vsssbf=NULL;
	Vector *Vpssbf=NULL;
	Vector *Vwbf=NULL;
	Vector *Vmwbf=NULL;

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

	unsigned short NumLinhasPadrao,J;
	Linha *G=NULL;
	Linha **LinhaPadrao=NULL;

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


	CargaDatosEG(&K,&N,&J,&NumLinhasPadrao,Dato->Archivo);
	LinhaPadrao=CargaLinhaPadrao(Dato->Archivo);
	G=CargaLinhaG(Dato->Archivo);
	//getchar();

	V=CrearVector(N);	IniciaConZeroVector(V);
	U=CrearVector(K);	IniciaConUnoZeroVector(U);
	R=CrearVectorFloat(N);	IniciaConZeroVectorFloat(R);
	Rhard=CrearVector(N);	IniciaConZeroVector(Rhard);

	Vr =CrearVector(N);	IniciaConZeroVector(Vr);
	Vbf=CrearVector(N);	IniciaConZeroVector(Vbf);
	Vbfs=CrearVector(N);	IniciaConZeroVector(Vbfs);
	Vshbf=CrearVector(N);	IniciaConZeroVector(Vshbf);
	Vphbf=CrearVector(N);	IniciaConZeroVector(Vphbf);
	Vsssbf=CrearVector(N);	IniciaConZeroVector(Vsssbf);
	Vpssbf=CrearVector(N);	IniciaConZeroVector(Vpssbf);
	Vwbf=CrearVector(N);	IniciaConZeroVector(Vwbf);
	Vmwbf=CrearVector(N);	IniciaConZeroVector(Vmwbf);

	tmp=fopen(Dato->ArchivoLog,"a");
	if(tmp==NULL) 
	{
		printf("Error abriendo/creando fichero: %s\n",Dato->ArchivoLog);
		return EXIT_SUCCESS;
	}
	
	printf("N=\t%hd\n",N);
	printf("K=\t%hd\n",K);
	printf("Rc=\t%f\n",K*1.0/N);
	printf("J=\t%hd\n",J);
	printf("LINHAS-PADRAO:\t%hd\n",NumLinhasPadrao);
	printf("\n");
	fprintf(tmp,"N=\t%hd\n",N);	
	fprintf(tmp,"K=\t%hd\n",K);
	fprintf(tmp,"Rc=\t%f\n",K*1.0/N);	
	fprintf(tmp,"J=\t%hd\n",J);	
	fprintf(tmp,"LINHAS-PADRAO:\t%hd\n",NumLinhasPadrao);	
	fprintf(tmp,"\n");
        	
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
       
	GeraPalavraCodigo(V,G,U);
	//multiplica(V,U,G);

	printf("\n");fprintf(tmp,"\n");

	//SNR=Eb/No
	SNRdb=0.0;
	dSNRdb=0.3333333333333;

	SNR=pow(10,SNRdb/10);
	Rc=(1.0*K)/N;
	P=Qfunc(sqrt(SNR*2.0*Rc));

	printf("K=%d\nN=%d\nRc=%f\n\n",K,N,Rc);
	fprintf(tmp,"K=%d\nN=%d\nRc=%f\n\n",K,N,Rc);

	fclose(tmp);
    
	crear_grafico_teorico_m();
	crear_LimiteDeShannon_m();
	crear_simula1_m(Dato,K,N,LinhaPadrao[0]->N);

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
			// Copia todos los bits R a Vr.
			for(j=0;j<N;j++) Vr->data[j]=Rhard->data[j];
			// Compara Ur con el U enviado.
			if((m=compara(Vr,V))!=0)	
			{
				ERRORES=ERRORES+1;
				ERRORES_BIT=ERRORES_BIT+m;
			}
			

			// Bit-Flipping
			if( (Dato->BF==1)&&(ACTIVA_BF==1) )
			{
				bitflipping(Vbf,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vbf))!=0)
				{
					ERRORES_BIT_BF=ERRORES_BIT_BF+m;
				}
			}

			// Bit-Flipping serial
			if( (Dato->BFS==1)&&(ACTIVA_BFS==1) )
			{
				bitflipping_serial(Vbfs,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vbfs))!=0)
				{
					ERRORES_BIT_BFS=ERRORES_BIT_BFS+m;
				}
			}

			// Serial Hard Bit-Flipping
			if( (Dato->SHBF==1)&&(ACTIVA_SHBF==1) )
			{
				SH_bitflipping(Vshbf,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vshbf))!=0)
				{
					ERRORES_BIT_SHBF=ERRORES_BIT_SHBF+m;
				}
			}

			// Parallel Hard Bit-Flipping
			if( (Dato->PHBF==1)&&(ACTIVA_PHBF==1) )
			{
				PH_bitflipping(Vphbf,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vphbf))!=0)
				{
					ERRORES_BIT_PHBF=ERRORES_BIT_PHBF+m;
				}
			}

			// Serial Sipser Spilman Bit-Flipping
			if( (Dato->SSSBF==1)&&(ACTIVA_SSSBF==1) )
			{
				S_SipSpi_bitflipping(Vsssbf,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vsssbf))!=0) 	
				{
					ERRORES_BIT_SSSBF=ERRORES_BIT_SSSBF+m;
				}
			}

			// Parallel Sipser Spilman Bit-Flipping
			if( (Dato->PSSBF==1)&&(ACTIVA_PSSBF==1) )
			{
				P_SipSpi_bitflipping(Vpssbf,Rhard,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vpssbf))!=0) 	
				{
					ERRORES_BIT_PSSBF=ERRORES_BIT_PSSBF+m;
				}
			}

			// Weighted Bit-Flipping
			if( (Dato->WBF==1)&&(ACTIVA_WBF==1) )
			{
				W_bitflipping(Vwbf,R,LinhaPadrao,NumLinhasPadrao,Dato->ItPor);
				if((m=compara(V,Vwbf))!=0) 	
				{
					ERRORES_BIT_WBF=ERRORES_BIT_WBF+m;
				}
			}

			// Modified Weighted Bit-Flipping
			if( (Dato->MWBF==1)&&(ACTIVA_MWBF==1) )
			{
				MW_bitflipping(Vmwbf,R,LinhaPadrao,NumLinhasPadrao,Dato->ItPor,Dato->Alpha);
				if((m=compara(V,Vmwbf))!=0) 	
				{
					ERRORES_BIT_MWBF=ERRORES_BIT_MWBF+m;
				}
			}


		}
		
		printf("ERRORES-PACK: %10ld\tITERACIONES: %10ld\tProb: %e \n",ERRORES,ITERACIONES,(ERRORES*1.0)/ITERACIONES);
		fprintf(tmp,"ERRORES-PACK: %10ld\tITERACIONES: %10ld\tProb: %e \n",ERRORES,ITERACIONES,(ERRORES*1.0)/ITERACIONES);

		PpackTeo =ProbabilidadPaquete(P,N);		//prob E paquete teorico.
		PpackSim =(ERRORES*1.0)/ITERACIONES;		//prob E paquete simulado.

		BER      =ERRORES_BIT/(N*1.0*ITERACIONES);	//BER.
		BERbf    =ERRORES_BIT_BF/(N*1.0*ITERACIONES);	//BER despues de BF.
		BERbfs   =ERRORES_BIT_BFS/(N*1.0*ITERACIONES);	//BER despues de BFS.
		BERshbf  =ERRORES_BIT_SHBF/(N*1.0*ITERACIONES);	//BER despues de SHBF.
		BERphbf  =ERRORES_BIT_PHBF/(N*1.0*ITERACIONES);	//BER despues de PHBF.
		BERsssbf =ERRORES_BIT_SSSBF/(N*1.0*ITERACIONES);//BER despues de SSSBF.
		BERpssbf =ERRORES_BIT_PSSBF/(N*1.0*ITERACIONES);//BER despues de PSSBF.
		BERwbf   =ERRORES_BIT_WBF/(N*1.0*ITERACIONES);	//BER despues de WBF.
		BERmwbf  =ERRORES_BIT_MWBF/(N*1.0*ITERACIONES);	//BER despues de MWBF.
                                                        
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

