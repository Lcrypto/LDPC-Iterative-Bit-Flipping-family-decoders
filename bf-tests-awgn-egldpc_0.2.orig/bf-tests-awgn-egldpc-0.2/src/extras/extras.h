/*
 * extras.h
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

#ifndef _EXTRAS_H_
#define _EXTRAS_H_


#include <stdlib.h>
#include <string.h>

//#include <config.h>

#ifndef BARRA
        #define BARRA '/'
#endif

#define ITPORDEFAULT 0.01
#define ALPHADEFAULT 0.1
#define BERMINDEFAULT 1.0e-7
#define PACKMAXDEFAULT 70000
typedef struct
{
	char *Archivo;
	char *ArchivoLog;
	char *ArchivoOut;
	int BF;
	int BFS;
	int SHBF;
	int PHBF;
	int SSSBF;
	int PSSBF;
	int WBF;
	int MWBF;
	float Alpha;
	float ItPor;
	float BERmin;
	long int PackMax;
}DATO;

/**\fn int fpr_exist_file(const char ARCHIVO[])
 * \brief Verifica la existencia de un archivo.
 * \param[in] ARCHIVO Archivo a averiguar si existe.
 * \return Retorna 1 si existe o 0 si no existe.
 */
int fpr_exist_file(const char ARCHIVO[])
{
	FILE *fd=NULL;
	fd=fopen(ARCHIVO,"r");

	if(fd==NULL)	{return 0;}
	else		{fclose(fd);return 1;}
}


void ImprimirAyuda(void)
{
	printf("\n");
	printf("%s\n",PACKAGE_STRING);
	printf("\n");
	printf("%s\t[-v][--file Archivo][--log ArchivoLog]\n",PACKAGE_NAME);
	printf("\t\t\t[--out ArchivoOut][--it-por porcentaje]\n");
	printf("\t\t\t[--bf][--bfs][--shbf][--phbf][--sssbf]\n");
	printf("\t\t\t[--pssbf][--wbf][--mwbf][--alpha factor]\n");
	printf("\t\t\t[--ber-min BERminimo][--pack-max Paquetes]\n");
	printf("\n");
	printf("\t-v\n");
	printf("\tDevuelve la versión actual del programa.\n");
	printf("\t\n");
	printf("\t--file Archivo\n");
	printf("\t\tDonde \"Archivo\" es el nombre del archivo que contiene a la \n");
	printf("\t\tmatriz de paridad EG-LDPC.\n");
	printf("\t\n");
	printf("\t--log ArchivoLog\n");
	printf("\t\tDonde \"ArchivoLog\" es el nombre del archivo del log que\n");
	printf("\t\tcontiene las estadísticas en las iteraciones, por defecto si\n");
	printf("\t\tno se específica este parámetro se agrega un .log al nombre\n");
	printf("\t\tdado en --file.\n");
	printf("\t\n");
	printf("\t--out ArchivoOut\n");
	printf("\t\tDonde \"ArchivoOut\" es el nombre del archivo del salida que\n");
	printf("\t\tcontiene los datos en las iteraciones, por defecto si\n");
	printf("\t\tno se específica este parámetro se agrega un .bf al nombre\n");
	printf("\t\tdado en --file.\n");
	printf("\t\n");
	printf("\t--it-por porcentaje\n");
	printf("\t\tPorcentaje de N que será usado como cantidad de iteraciones\n");
	printf("\t\tmáximas. Si no se especifica el valor es: %e \n",ITPORDEFAULT);
	printf("\t\n");
	printf("\t--bf\n");
	printf("\t\tActivo el algoritmo bit-flipping.\n");
	printf("\t\n");
	printf("\t--bfs\n");
	printf("\t\tActivo el algoritmo bit-flipping serial.\n");
	printf("\t\n");
	printf("\t--shbf\n");
	printf("\t\tActivo el algoritmo serial hard bit-flipping.\n");
	printf("\t\n");
	printf("\t--phbf\n");
	printf("\t\tActivo el algoritmo parallel hard bit-flipping.\n");
	printf("\t\n");
	printf("\t--sssbf\n");
	printf("\t\tActivo el algoritmo secuencial Sipser Spilam bit-flipping.\n");
	printf("\t\n");
	printf("\t--pssbf\n");
	printf("\t\tActivo el algoritmo parallel  Sipser Spilam bit-flipping.\n");
	printf("\t\n");
	printf("\t--wbf\n");
	printf("\t\tActivo el algoritmo weighted bit-flipping.\n");
	printf("\t\n");
	printf("\t--mwbf\n");
	printf("\t\tActivo el algoritmo modified weighted bit-flipping.\n");
	printf("\t\n");
	printf("\t--alpha factor\n");
	printf("\t\tEl factor alpha usado en el algoritmo mwbf. Si no se especifica\n");
	printf("\t\tel valor es: %e \n",ALPHADEFAULT);
	printf("\t\n");
	printf("\t--ber-min BERminimo\n");
	printf("\t\tBERminimo es el valor mínimo que tomara cualquiera de los\n");
	printf("\t\talgoritmos BF, despues de esde valor el algoritmo se desactiva.\n");
	printf("\t\tEl valor por defecto es: %e \n",BERMINDEFAULT);
	printf("\t\n");
	printf("\t--pack-max Paquetes\n");
	printf("\t\tPaquetes es la cantidad de paquetes máximos enviados en una iteración\n");
	printf("\t\teste dato esta ligado a la confiabilidad del BER mínimo encontrado\n");
	printf("\t\ten la iteración, todo valor de BER menor que 100/(N*Paquetes) es\n");
	printf("\t\tdesconfiable dado que no tubo la suficiente cantidad de iteraciones\n");
	printf("\t\tel criterio es que para algo de probabilidad P se deben hacer 100/P\n");
	printf("\t\tpruebas para obtener una probabilidad experimental confiable, el valor\n");
	printf("\t\tpor defecto de Paquetes es: %d \n",PACKMAXDEFAULT);
	printf("\t\n");

}

void ImprimirVersion(void)
{
	printf("\n%s\n\n",PACKAGE_STRING);
}

DATO* Menu(int argc,char **argv)
{
	int i;
	DATO* Dato=NULL;

	Dato=(DATO*)calloc(sizeof(DATO),1);
	Dato->Archivo=NULL;
	Dato->ArchivoLog=NULL;
	Dato->BF=0;
	Dato->BFS=0;
	Dato->SHBF=0;
	Dato->PHBF=0;
	Dato->SSSBF=0;
	Dato->PSSBF=0;
	Dato->WBF=0;
	Dato->MWBF=0;
	Dato->ItPor=ITPORDEFAULT;
	Dato->Alpha=ALPHADEFAULT;
	Dato->PackMax=PACKMAXDEFAULT;

	if(argc==1) 
	{
		ImprimirAyuda(); 
		exit(EXIT_SUCCESS);
	}

	if( (argc==2) && (strcmp(argv[1],"-v")==0) ) 
	{
		ImprimirVersion(); 
		exit(EXIT_SUCCESS);
	}

	for(i=1;i<argc;i++)
	{
		if(strcmp(argv[i],"--file")==0) 
		{
			if(i+1<argc)	
			{
				Dato->Archivo=(char*)calloc(1,strlen(argv[i+1])+1);
				strcpy(Dato->Archivo,argv[i+1]);
			}
			else		Dato->Archivo=NULL;
			i++;
		}
		else if(strcmp(argv[i],"--log")==0)
		{
			if(i+1<argc)
			{
				Dato->ArchivoLog=(char*)calloc(1,strlen(argv[i+1])+1);
				strcpy(Dato->Archivo,argv[i+1]);
			}
			else		Dato->ArchivoLog=NULL;
			i++;
		}
		else if(strcmp(argv[i],"--out")==0)
		{
			if(i+1<argc)
			{
				Dato->ArchivoOut=(char*)calloc(1,strlen(argv[i+1])+1);
				strcpy(Dato->ArchivoOut,argv[i+1]);
			}
			else		Dato->ArchivoOut=NULL;
			i++;
		}
		else if(strcmp(argv[i],"--it-por")==0)
		{
			if(i+1<argc)
			{
				Dato->ItPor=atof(argv[i+1]);
			}
			else	Dato->ItPor=ITPORDEFAULT;
			i++;
		}
		else if(strcmp(argv[i],"--alpha")==0)
		{
			if(i+1<argc)
			{
				Dato->Alpha=atof(argv[i+1]);
			}
			else	Dato->Alpha=ALPHADEFAULT;
			i++;
		}
		else if(strcmp(argv[i],"--ber-min")==0)
		{
			if(i+1<argc)
			{
				Dato->BERmin=atof(argv[i+1]);
			}
			else	Dato->BERmin=BERMINDEFAULT;
			i++;
		}
		else if(strcmp(argv[i],"--pack-max")==0)
		{
			if(i+1<argc)
			{
				Dato->PackMax=atol(argv[i+1]);
			}
			else	Dato->PackMax=PACKMAXDEFAULT;
			i++;
		}
		else if(strcmp(argv[i],"--bf")==0)
		{
			Dato->BF=1;
		}
		else if(strcmp(argv[i],"--bfs")==0)
		{
			Dato->BFS=1;
		}
		else if(strcmp(argv[i],"--shbf")==0)
		{
			Dato->SHBF=1;
		}
		else if(strcmp(argv[i],"--phbf")==0)
		{
			Dato->PHBF=1;
		}
		else if(strcmp(argv[i],"--sssbf")==0)
		{
			Dato->SSSBF=1;
		}
		else if(strcmp(argv[i],"--pssbf")==0)
		{
			Dato->PSSBF=1;
		}
		else if(strcmp(argv[i],"--wbf")==0)
		{
			Dato->WBF=1;
		}
		else if(strcmp(argv[i],"--mwbf")==0)
		{
			Dato->MWBF=1;
		}
	}

	
	return Dato;
}

int crear_grafico_teorico_m(void)
{
	FILE *fd=NULL;
	fd=fopen("grafico_teorico.m","w");
	if(fd==NULL) 
	{printf("ERROR creando grafico_teorico.m\n");return -1;}

	fprintf(fd,"function h = grafico_teorico(p,x)\n");
	fprintf(fd,"\n");
	fprintf(fd,"n = length(p);\n");
	fprintf(fd,"h=zeros(n,1);\n");
	fprintf(fd,"\n");
	fprintf(fd,"if (mod(x,2)==0)\n");
	fprintf(fd,"     %c%cPar\n",'%','%');
	fprintf(fd,"     for II = 1:n\n");
	fprintf(fd,"\n");
	fprintf(fd,"		h(II)=1.0-binocdf(x/2.0,x,p(II))+p(II)*binopdf(x/2,x,p(II));\n");
	fprintf(fd,"     end\n");
	fprintf(fd,"     text='par'\n");
	fprintf(fd,"else\n");
	fprintf(fd,"     %c%cImpar\n",'%','%');
	fprintf(fd,"     for II = 1:n\n");
	fprintf(fd,"\n");
	fprintf(fd,"		h(II)=1.0-binocdf((x-1.0)/2.0,x,p(II));\n");
	fprintf(fd,"     end\n");
	fprintf(fd,"     text='impar'\n");
	fprintf(fd,"end\n");

	fclose(fd);
	return 0;
}

int crear_LimiteDeShannon_m(void)
{
	FILE *fd=NULL;
	fd=fopen("LimiteDeShannon.m","w");
	if(fd==NULL) 
	{printf("ERROR creando LimiteDeShannon.m\n");return -1;}

	fprintf(fd,"function y = LimiteDeShannon (P)\n");
	fprintf(fd,"	global Rc;\n");
	fprintf(fd,"\n");
	fprintf(fd,"	y=Rc-(1-P*log2(1/P)-(1-P)*log2(1/(1-P)));\n");
	fprintf(fd,"end\n");

	fclose(fd);
	return 0;
}

int crear_simula1_m(const DATO* Data,short int K,short int N,short int X)
{
	FILE *fd=NULL;
	int i;

	fd=fopen("simula1.m","w");
	if(fd==NULL) {printf("ERROR creando simula1.m\n");return -1;}

	fprintf(fd,"x=load(\'%s\');\n",Data->ArchivoOut);
	fprintf(fd,"X=x(:,1);		%c%c  Probabilidade de erro do bit (Canal)\n",'%','%');
	fprintf(fd,"Pb  =x(:,2);		%c%c  BER\n",'%','%');
	fprintf(fd,"Ppack_cal=x(:,3);	%c%c  Probabilidade de erro do pacote RX (calculado)\n",'%','%');
	fprintf(fd,"Ppack_tot=x(:,4);	%c%c  Probabilidade de erro do pacote RX (experimental)\n",'%','%');

	i=5;

	if(Data->BF==1)
	{fprintf(fd,"BERbf=x(:,%d);		%c%c  BER bit-flipping\n",i,'%','%'); i++;}
	if(Data->BFS==1)
	{fprintf(fd,"BERbfs=x(:,%d);		%c%c  BER bit-flipping sequencial\n",i,'%','%'); i++;}
	if(Data->SHBF==1)
	{fprintf(fd,"BERshbf=x(:,%d);		%c%c  BER serial hard bit-flipping\n",i,'%','%'); i++;}
	if(Data->PHBF==1)
	{fprintf(fd,"BERphbf=x(:,%d);		%c%c  BER parallel hard bit-flipping\n",i,'%','%'); i++;}
	if(Data->SSSBF==1)
	{fprintf(fd,"BERsssbf=x(:,%d);	%c%c  BER serial Sipser Spilman bit-flipping\n",i,'%','%'); i++;}
	if(Data->PSSBF==1)
	{fprintf(fd,"BERpssbf=x(:,%d);	%c%c  BER parallel Sipser Spilman bit-flipping\n",i,'%','%'); i++;}
	if(Data->WBF==1)
	{fprintf(fd,"BERwbf=x(:,%d);		%c%c  BER weighted bit-flipping\n",i,'%','%'); i++;}
	if(Data->MWBF==1)
	{fprintf(fd,"BERmwbf=x(:,%d);	%c%c  BER modified weighted bit-flipping\n",i,'%','%'); i++;}
	fprintf(fd,"\n");
	fprintf(fd,"\n");
	fprintf(fd,"global Rc;\n");
	fprintf(fd,"Rc=%d/%d;\n",K,N);
	fprintf(fd,"[P, info] = fsolve (\'LimiteDeShannon\', 0.0001)\n");
	fprintf(fd,"\n");
	fprintf(fd,"XEbNo=10*log10((qfuncinv(X).^2)/(2.0*Rc));\n");
	fprintf(fd,"PEbNo=10*log10((qfuncinv(P).^2)/(2.0*Rc));\n");
	fprintf(fd,"\n");
	fprintf(fd,"Pb=qfunc(qfuncinv(Pb)/Rc^0.5);\n");
	fprintf(fd,"AL=1.5;\n");
	fprintf(fd,"set(gcf,\'DefaultLineLineWidth\',AL);\n");
	fprintf(fd,"MZ=15;\n");
	fprintf(fd,"\n");
	fprintf(fd,"Ashannon=[PEbNo PEbNo PEbNo PEbNo PEbNo PEbNo PEbNo PEbNo PEbNo];\n");
	fprintf(fd,"Bshannon=[1e-4  3e-4  1e-3  3e-3  1e-2  3e-2  1e-1  3e-1  1    ];\n");
	fprintf(fd,"\n");
	fprintf(fd,"semilogy(XEbNo,Pb,\'k--\',\'markersize\',MZ,");

	if(Data->BF==1)
	fprintf(fd,"XEbNo,BERbf,\'k--^\',\'markersize\',MZ,");
	if(Data->BFS==1)
	fprintf(fd,"XEbNo,BERbfs,\'k--x\',\'markersize\',MZ,");
	if(Data->SHBF==1)
	fprintf(fd,"XEbNo,BERshbf,\'k--v\',\'markersize\',MZ,");
	if(Data->PHBF==1)
	fprintf(fd,"XEbNo,BERphbf,\'k--<\',\'markersize\',MZ,");
	if(Data->SSSBF==1)
	fprintf(fd,"XEbNo,BERsssbf,\'k-->\',\'markersize\',MZ,");
	if(Data->PSSBF==1)
	fprintf(fd,"XEbNo,BERpssbf,\'k--o\',\'markersize\',MZ,");
	if(Data->WBF==1)
	fprintf(fd,"XEbNo,BERwbf,\'k--s\',\'markersize\',MZ,");
	if(Data->MWBF==1)
	fprintf(fd,"XEbNo,BERmwbf,\'k--p\',\'markersize\',MZ,");

	fprintf(fd,"Ashannon,Bshannon,\'-.\',\'markersize\',MZ,");
	fprintf(fd,"XEbNo,grafico_teorico(X,%d),\'*\',\'markersize\',MZ),grid on\n",X);
	fprintf(fd,"axis([0 10 10^(-7) 1])\n");
	fprintf(fd,"\n");
	fprintf(fd,"set(gca,\'fontsize\',24,\'GridLineStyle\',\'--\'); %c sets font of numbers on axes\n",'%');
	fprintf(fd,"\n");
	fprintf(fd,"xlabel(\'Eb/No (dB)\',\'fontsize\',24);\n");
	fprintf(fd,"ylabel(\'BER\',\'fontsize\',24);\n");

	fprintf(fd,"legend(\'UNCODED\',");
	if(Data->BF==1)
	fprintf(fd,"\'BF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->BFS==1)
	fprintf(fd,"\'BFS IT=%d\',",(int)(Data->ItPor*N));
	if(Data->SHBF==1)
	fprintf(fd,"\'SHBF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->PHBF==1)
	fprintf(fd,"\'PHBF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->SSSBF==1)
	fprintf(fd,"\'SSSBF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->PSSBF==1)
	fprintf(fd,"\'PSSBF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->WBF==1)
	fprintf(fd,"\'WBF IT=%d\',",(int)(Data->ItPor*N));
	if(Data->MWBF==1)
	fprintf(fd,"\'MWBF IT=%d \\alpha=%5.3f\',",(int)(Data->ItPor*N),Data->Alpha);
	fprintf(fd,"\'Shannon L. AWGN Hard\',\'Lower Bound\',3);");

	fprintf(fd,"legend(\'boxon\');\n");
	fprintf(fd,"\n");
	fprintf(fd,"if(strcmp(computer,'PCWIN')==1)\n");
	fprintf(fd,"	print('%s.eps','-deps');\n",Data->ArchivoOut);
	fprintf(fd,"else\n");
	fprintf(fd,"	print('%s.eps','-deps','-portrait','-F:24');\n",Data->ArchivoOut);
	fprintf(fd,"end\n");

	fclose(fd);

	return 0;
}

#endif /* _EXTRAS_H_ */
