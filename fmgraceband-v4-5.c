/*BY FILIPE MATUSALEM, MARCH 2014     fmatusalem@ita.br */
//version 4.5 Sept 2022
//compilation:  g++ -o fmgraceband.x fmgraceband-v4-5.c
//use: just ./fmgraceband.x on the directory with the VASP band structure calculation.
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "filipe.h"



int main(int argc, char *argv[])
{
FILE *eigenval,*kpoints,*xmgrace,*outcar;
float x,y,z,ref,efermi,xmax,min,max,par,pi,zz,vac1,vac2;
int i,j,k,nbands,simpt,kpts,nelect,nef,spin;
char str[150],ch;


printf("\nPROGRAM TO GENERATE GRACE BAND STRUCTURES FROM VASP RESULTS (SPIN OR NON-SPIN)\n\n");
printf("-----------------------------IF ANYTHING GOES WRONG, REMEMBER----------------------------------\n");
printf("The first line of KPOINTS file must contain the symmetric points  where the calculation is done\n");
printf("The band structure calculation can be done with a line-mode or explicit KPOINTS file.\n\n"); 
printf("The reference energy and y axes interval can be entered as argument to the program.\n");
printf("To use Fermi energy as reference -> fmgraceband fermi or fmgraceband fermi -10 10\n\n");  
printf("On version 4.3: two extra horizontal lines can be added given the y positions as arguments, fmgraceband fermi -10 10 -5.62 -6.45\n\n");  
printf("To use homo or lumo energy as reference -> fmgraceband homo or fmgraceband lumo\n\n");  
printf("by Filipe Matusalem - March 2014\n");

printf("-----------------------------------------------------------------------------------------------\n\n");

pi=3.14159265359;

outcar = fopen("OUTCAR","r"); /* Arquivo ASCII, para leitura */
if(!outcar)
{
printf( "Error opening OUTCAR file\n");
exit(0);
}

eigenval = fopen("EIGENVAL","r"); /* Arquivo ASCII, para leitura */
if(!eigenval)
{
printf( "Error opening EIGENVAL file\n");
exit(0);
}

kpoints = fopen("KPOINTS","r"); /* Arquivo ASCII, para leitura */
if(!kpoints)
{
printf( "Error opening KPOINTS file\n");
exit(0);
}

xmgrace = fopen("bands.agr","w"); /* Arquivo ASCII, para escrita */
if(!xmgrace)
{
printf( "Error creating bands.agr file\n");
exit(0);
}


par=scale();
efermi=fermienergy();

if( argc > 1 ){

if(strcmp(argv[1],"fermi")==0)ref=efermi;
else {if(strcmp(argv[1],"homo")==0)ref=homo(0); else {if(strcmp(argv[1],"lumo")==0)ref=lumo(0); else{ref=atof(argv[1]);}}}}

else {ref=efermi;}


if(argc > 2){
min = atof(argv[2]);
max = atof(argv[3]);}
else {
min=(efermi-ref)-5;max=(efermi-ref)+5;
}

vac1=vac2=0;
if(argc > 4){vac1 = atof(argv[4]);}
if(argc > 5){vac2 = atof(argv[5]);}


fprintf(xmgrace,"# Grace project file\n");
fprintf(xmgrace,"#\n");
fprintf(xmgrace,"@version 50123\n");
fprintf(xmgrace,"@page size 1000, 700\n");
fprintf(xmgrace,"@map font 4 to \"Times-Roman\", \"Times-Roman\"\n");
fprintf(xmgrace,"@map font 8 to \"Symbol\", \"Symbol\"\n");
fprintf(xmgrace,"@g0 on\n");
fprintf(xmgrace,"@g0 type XY\n");
fprintf(xmgrace,"@with g0\n");
fprintf(xmgrace,"@    stack world 0, 0, 0, 0\n");
fprintf(xmgrace,"@    view 0.150000, 0.150000, 1.000000, 0.80000\n");


for(i=0;i<4;i++){
fscanf(eigenval,"%d",&spin);}


for(i=0;i<5;i++){
do
ch = getc(eigenval);              /*pula 5 linhas*/
while(ch!='\n');}

fscanf(eigenval,"%d",&kpts);
fscanf(eigenval,"%d",&kpts);
fscanf(eigenval,"%d",&nbands);




float kpoint[kpts][3],vet[kpts][3];
int aux[20],nkpsim[20];

rewind(outcar);
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra 2pi/SCALE*/
while(strcmp(str,"2pi/SCALE")!=0);

do
ch = getc(outcar);              /*chega ao fim da linha*/
while(ch!='\n');


for(i=0;i<kpts;i++){
for(j=0;j<3;j++){
fscanf(outcar,"%f",&kpoint[i][j]);
kpoint[i][j]=kpoint[i][j]*2*pi/par;}
do
ch = getc(outcar);              /*chega ao fim da linha*/
while(ch!='\n');}

simpt=2;
k=1;
for(i=0;i<kpts-1;i++){
if(fabs(kpoint[i+1][0]-kpoint[i][0])<0.00001 && fabs(kpoint[i+1][1]-kpoint[i][1])<0.00001 && fabs(kpoint[i+1][2]-kpoint[i][2])<0.00001){     /*tolerancia de 0.00001 na igualdade entre os numeros*/
kpoint[k][0]=kpoint[i][0];
kpoint[k][1]=kpoint[i][1];
kpoint[k][2]=kpoint[i][2];
nkpsim[k]=i;
k++;simpt++;}
}
nkpsim[k]=i;
nkpsim[0]=0;

kpoint[k][0]=kpoint[kpts-1][0];
kpoint[k][1]=kpoint[kpts-1][1];
kpoint[k][2]=kpoint[kpts-1][2];

float dist[simpt],abs[kpts];

dist[0]=0;
z=0;
for(i=1;i<simpt;i++){
	dist[i]=sqrt(pow(kpoint[i][0]-kpoint[i-1][0],2)+pow(kpoint[i][1]-kpoint[i-1][1],2)+pow(kpoint[i][2]-kpoint[i-1][2],2));
z=z+dist[i];
}

zz=z;

aux[0]=0;
for(i=0;i<simpt;i++) aux[i+1]=nkpsim[i+1]-nkpsim[i]-1;
aux[1]=aux[1]+1;


abs[0]=0;
k=0;
for(i=1;i<simpt;i++){
while(k<nkpsim[i])
{
abs[k+1]=abs[k]+dist[i]/aux[i]; 
k++;
}
abs[k+1]=abs[k];k++; 
}


do
ch = getc(outcar);              /*chega ao fim da linha*/
while(ch!='\n');
do
ch = getc(outcar);              /*chega ao fim da linha*/
while(ch!='\n');

for(i=0;i<kpts;i++){
for(j=0;j<3;j++) fscanf(outcar,"%f",&vet[i][j]);
do
ch = getc(outcar);              /*chega ao fim da linha*/
while(ch!='\n');
}


k=1;
for(i=0;i<kpts-1;i++){
if(fabs(vet[i+1][0]-vet[i][0])<0.00001 && fabs(vet[i+1][1]-vet[i][1])<0.00001 && fabs(vet[i+1][2]-vet[i][2])<0.00001){     /*tolerancia de 0.00001 na igualdade entre os numeros*/
vet[k][0]=vet[i][0];
vet[k][1]=vet[i][1];
vet[k][2]=vet[i][2];
k++;}
}


vet[k][0]=vet[kpts-1][0];
vet[k][1]=vet[kpts-1][1];
vet[k][2]=vet[kpts-1][2];



/*-------------------------------------------------------------------------------------------------------------------*/

if(spin==2){
float mtx[2*kpts+1][nbands];


for(j=0;j<2;j++){
do
ch = getc(eigenval);              /*chega ao fim da linha*/
while(ch!='\n');}



for(i=1;i<=kpts;i++){

fscanf(eigenval,"%s",str);

do
ch = getc(eigenval);              /*chega ao fim da linha*/
while(ch!='\n');

      for(j=0;j<nbands;j++){
              fscanf(eigenval,"%s",str);
              fscanf(eigenval,"%f",&mtx[i][j]);
              fscanf(eigenval,"%f",&mtx[i+kpts][j]);
              do
              ch = getc(eigenval);              /*chega ao fim da linha*/
              while(ch!='\n');
              mtx[i][j]=mtx[i][j]-ref;
              mtx[i+kpts][j]=mtx[i+kpts][j]-ref;}

}

fprintf(xmgrace,"@    world 0,   %f, %f,     %f\n",min,z,max);

for(i=0;i<nbands;i++){
fprintf(xmgrace,"@    s%d line linestyle 1\n",i);
fprintf(xmgrace,"@    s%d line color 1\n",i);
}
for(i=nbands;i<2*nbands;i++){
fprintf(xmgrace,"@    s%d line linestyle 2\n",i);
fprintf(xmgrace,"@    s%d line color 1\n",i);
}

for(i=2*nbands;i<2*nbands+simpt-2;i++){
fprintf(xmgrace,"@    s%d line linestyle 1\n",i);
fprintf(xmgrace,"@    s%d line color 1\n",i);
}

fprintf(xmgrace,"@    s%d line linestyle 1\n",i);
fprintf(xmgrace,"@    s%d line color 2\n",i);

if(vac1!=0){
fprintf(xmgrace,"@    s%d line linestyle 3\n",i+1);
fprintf(xmgrace,"@    s%d line color 2\n",i+1);}

if(vac2!=0){
fprintf(xmgrace,"@    s%d line linestyle 2\n",i+2);
fprintf(xmgrace,"@    s%d line color 4\n",i+2);}


fprintf(xmgrace,"@    xaxis  on\n");
fprintf(xmgrace,"@    xaxis  label \"\"\n");
fprintf(xmgrace,"@    xaxis  label place spec\n");
fprintf(xmgrace,"@    xaxis  label place 0.000000, -0.13\n");
fprintf(xmgrace,"@    xaxis  label font 4\n");
fprintf(xmgrace,"@    xaxis  label color 1\n");
fprintf(xmgrace,"@    xaxis  label char size 1.50000\n");

fprintf(xmgrace,"@    xaxis  tick on\n");
fprintf(xmgrace,"@    xaxis  tick major 0.1\n");
fprintf(xmgrace,"@    xaxis  tick minor ticks 0\n");
fprintf(xmgrace,"@    xaxis  tick default 3\n");
fprintf(xmgrace,"@    xaxis  tick place rounded true\n");
fprintf(xmgrace,"@    xaxis  tick in\n");
fprintf(xmgrace,"@    xaxis  tick major size 0.500000\n");
fprintf(xmgrace,"@    xaxis  tick major color 1\n");
fprintf(xmgrace,"@    xaxis  tick major linewidth 1.0\n");
fprintf(xmgrace,"@    xaxis  tick major linestyle 1\n");
fprintf(xmgrace,"@    xaxis  tick major grid off\n");
fprintf(xmgrace,"@    xaxis  tick minor color 1\n");
fprintf(xmgrace,"@    xaxis  tick minor linewidth 1.0\n");
fprintf(xmgrace,"@    xaxis  tick minor linestyle 1\n");
fprintf(xmgrace,"@    xaxis  tick minor grid off\n");
fprintf(xmgrace,"@    xaxis  tick minor size 0.500000\n");
fprintf(xmgrace,"@    xaxis  ticklabel on\n");
fprintf(xmgrace,"@    xaxis  ticklabel format general\n");
fprintf(xmgrace,"@    xaxis  ticklabel prec 6\n");
fprintf(xmgrace,"@    xaxis  ticklabel angle 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel skip 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel stagger 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel place normal\n");
fprintf(xmgrace,"@    xaxis  ticklabel offset auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel offset 0.000000 , 0.010000\n");
fprintf(xmgrace,"@    xaxis  ticklabel start type auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel start 0.000000\n");
fprintf(xmgrace,"@    xaxis  ticklabel stop type auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel stop 0.000000\n");
fprintf(xmgrace,"@    xaxis  ticklabel char size 1.300000\n");
fprintf(xmgrace,"@    xaxis  ticklabel font 4\n");
fprintf(xmgrace,"@    xaxis  ticklabel color 1\n");
fprintf(xmgrace,"@    xaxis  tick place normal\n");
fprintf(xmgrace,"@    xaxis  tick spec type both\n");


rewind(kpoints);
dist[0]=0;
z=0;
fprintf(xmgrace,"@    xaxis  tick spec %d\n",simpt);
for(i=0;i<simpt;i++){
z=z+dist[i]; 
fscanf(kpoints,"%s",str);

fprintf(xmgrace,"@    xaxis  tick major %d, %f\n",i,z);
fprintf(xmgrace,"@    xaxis  ticklabel %d, \"%s\"\n",i,str);

}




fprintf(xmgrace,"@    yaxis  on\n");
fprintf(xmgrace,"@    yaxis  tick on\n");
fprintf(xmgrace,"@    yaxis  tick major 5\n");
fprintf(xmgrace,"@    yaxis  tick minor ticks 1\n");
fprintf(xmgrace,"@    yaxis  label \"Energy (eV)\"\n");
fprintf(xmgrace,"@    yaxis  label font 4\n");
fprintf(xmgrace,"@    yaxis  label color 1\n");
fprintf(xmgrace,"@    yaxis  label char size 1.40000\n");
fprintf(xmgrace,"@    yaxis  ticklabel on\n");
fprintf(xmgrace,"@    yaxis  ticklabel format decimal\n");
fprintf(xmgrace,"@    yaxis  ticklabel prec 1\n");
fprintf(xmgrace,"@    yaxis  ticklabel op left\n");
fprintf(xmgrace,"@    yaxis  ticklabel font 4\n");
fprintf(xmgrace,"@    yaxis  ticklabel color 1\n");
fprintf(xmgrace,"@    yaxis  ticklabel char size 1.20000\n");
fprintf(xmgrace,"@    yaxis  tick major on\n");
fprintf(xmgrace,"@    yaxis  tick minor on\n");
fprintf(xmgrace,"@    yaxis  tick out\n");
fprintf(xmgrace,"@    yaxis  tick major color 1\n");
fprintf(xmgrace,"@    yaxis  tick major linestyle 1\n");
fprintf(xmgrace,"@    yaxis  tick minor color 1\n");
fprintf(xmgrace,"@    yaxis  tick minor linestyle 1\n");
fprintf(xmgrace,"@    yaxis  tick size 0.750000\n");
fprintf(xmgrace,"@    yaxis  tick minor size 0.500000\n");
fprintf(xmgrace,"@    yaxis  tick op left\n");
fprintf(xmgrace,"@    legend off\n");
fprintf(xmgrace,"@    frame on\n");
fprintf(xmgrace,"@    frame linestyle 1\n");
fprintf(xmgrace,"@    frame color 1\n");
fprintf(xmgrace,"@    frame background color 0\n");

fprintf(xmgrace,"@    legend on\n");
fprintf(xmgrace,"@    legend loctype view\n");
fprintf(xmgrace,"@    legend 1.1, 0.5\n");
fprintf(xmgrace,"@    legend box color 1\n");
fprintf(xmgrace,"@    legend box pattern 1\n");
fprintf(xmgrace,"@    legend box linewidth 1.0\n");
fprintf(xmgrace,"@    legend box linestyle 1\n");
fprintf(xmgrace,"@    legend box fill color 0\n");
fprintf(xmgrace,"@    legend box fill pattern 1\n");
fprintf(xmgrace,"@    legend font 4\n");
fprintf(xmgrace,"@    legend char size 1.00000\n");
fprintf(xmgrace,"@    legend color 1\n");
fprintf(xmgrace,"@    legend length 4\n");
fprintf(xmgrace,"@    legend vgap 2\n");
fprintf(xmgrace,"@    legend hgap 1\n");
fprintf(xmgrace,"@    legend invert false\n");
fprintf(xmgrace,"@    s0 legend  \"spin up\"\n");
fprintf(xmgrace,"@    s%d legend  \"spin down\"\n",nbands);
fprintf(xmgrace,"@    s%d legend  \"Fermi energy\"\n",2*nbands+simpt-2);



for(j=0;j<nbands;j++){
    fprintf(xmgrace,"@target G0.S%d\n",j);
    fprintf(xmgrace,"@type xy\n");
       for(i=1;i<=kpts;i++){
          fprintf(xmgrace,"%16.4f%16.4f\n",abs[i-1],mtx[i][j]);
          }}

for(j=0;j<nbands;j++){
    fprintf(xmgrace,"@target G0.S%d\n",j+nbands);
    fprintf(xmgrace,"@type xy\n");
       for(i=kpts+1;i<=2*kpts;i++){
          fprintf(xmgrace,"%16.4f%16.4f\n",abs[i-kpts-1],mtx[i][j]);
          }}


dist[0]=0;
z=0;
for(i=0;i<simpt;i++){
z=z+dist[i];
if(i>0 && i<simpt-1){ 
fprintf(xmgrace,"@target G0.S%d\n",nbands+i-1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",z,min);
fprintf(xmgrace,"%16.4f%16.4f\n",z,max);}
}

fprintf(xmgrace,"@target G0.S%d\n",nbands+i-2);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,efermi-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,efermi-ref);

if(vac1!=0){fprintf(xmgrace,"@target G0.S%d\n",nbands+i+1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,vac1-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,vac1-ref);}

if(vac2!=0){fprintf(xmgrace,"@target G0.S%d\n",nbands+i+1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,vac2-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,vac2-ref);}

}

/*------------------------------------------------------------------------------------------------------*/
else{
float mtx[kpts+1][nbands];


for(j=0;j<2;j++){
do
ch = getc(eigenval);              /*chega ao fim da linha*/
while(ch!='\n');}



for(i=1;i<=kpts;i++){

fscanf(eigenval,"%s",str);

do
ch = getc(eigenval);              /*chega ao fim da linha*/
while(ch!='\n');

      for(j=0;j<nbands;j++){
              fscanf(eigenval,"%s",str);
              fscanf(eigenval,"%f",&mtx[i][j]);
              do
              ch = getc(eigenval);              /*chega ao fim da linha*/
              while(ch!='\n');
              mtx[i][j]=mtx[i][j]-ref;}

}






fprintf(xmgrace,"@    world 0,   %f, %f,     %f\n",min,z,max);



for(i=0;i<nbands;i++){
fprintf(xmgrace,"@    s%d line linestyle 1\n",i);
fprintf(xmgrace,"@    s%d line color 1\n",i);
}

for(i=nbands;i<nbands+simpt-2;i++){
fprintf(xmgrace,"@    s%d line linestyle 2\n",i);
fprintf(xmgrace,"@    s%d line color 1\n",i);}





fprintf(xmgrace,"@    s%d line linestyle 1\n",i);
fprintf(xmgrace,"@    s%d line color 2\n",i);

if(vac1!=0){
fprintf(xmgrace,"@    s%d line linestyle 3\n",i+1);
fprintf(xmgrace,"@    s%d line color 2\n",i+1);}

if(vac2!=0){
fprintf(xmgrace,"@    s%d line linestyle 2\n",i+2);
fprintf(xmgrace,"@    s%d line color 4\n",i+2);}



fprintf(xmgrace,"@    xaxis  on\n");
fprintf(xmgrace,"@    xaxis  label \"\"\n");
fprintf(xmgrace,"@    xaxis  label place spec\n");
fprintf(xmgrace,"@    xaxis  label place 0.000000, -0.13\n");
fprintf(xmgrace,"@    xaxis  label font 4\n");
fprintf(xmgrace,"@    xaxis  label color 1\n");
fprintf(xmgrace,"@    xaxis  label char size 1.50000\n");

fprintf(xmgrace,"@    xaxis  tick on\n");
fprintf(xmgrace,"@    xaxis  tick major 0.1\n");
fprintf(xmgrace,"@    xaxis  tick minor ticks 0\n");
fprintf(xmgrace,"@    xaxis  tick default 3\n");
fprintf(xmgrace,"@    xaxis  tick place rounded true\n");
fprintf(xmgrace,"@    xaxis  tick in\n");
fprintf(xmgrace,"@    xaxis  tick major size 0.500000\n");
fprintf(xmgrace,"@    xaxis  tick major color 1\n");
fprintf(xmgrace,"@    xaxis  tick major linewidth 1.0\n");
fprintf(xmgrace,"@    xaxis  tick major linestyle 1\n");
fprintf(xmgrace,"@    xaxis  tick major grid off\n");
fprintf(xmgrace,"@    xaxis  tick minor color 1\n");
fprintf(xmgrace,"@    xaxis  tick minor linewidth 1.0\n");
fprintf(xmgrace,"@    xaxis  tick minor linestyle 1\n");
fprintf(xmgrace,"@    xaxis  tick minor grid off\n");
fprintf(xmgrace,"@    xaxis  tick minor size 0.500000\n");
fprintf(xmgrace,"@    xaxis  ticklabel on\n");
fprintf(xmgrace,"@    xaxis  ticklabel format general\n");
fprintf(xmgrace,"@    xaxis  ticklabel prec 6\n");
fprintf(xmgrace,"@    xaxis  ticklabel angle 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel skip 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel stagger 0\n");
fprintf(xmgrace,"@    xaxis  ticklabel place normal\n");
fprintf(xmgrace,"@    xaxis  ticklabel offset auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel offset 0.000000 , 0.010000\n");
fprintf(xmgrace,"@    xaxis  ticklabel start type auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel start 0.000000\n");
fprintf(xmgrace,"@    xaxis  ticklabel stop type auto\n");
fprintf(xmgrace,"@    xaxis  ticklabel stop 0.000000\n");
fprintf(xmgrace,"@    xaxis  ticklabel char size 1.300000\n");
fprintf(xmgrace,"@    xaxis  ticklabel font 4\n");
fprintf(xmgrace,"@    xaxis  ticklabel color 1\n");
fprintf(xmgrace,"@    xaxis  tick place normal\n");
fprintf(xmgrace,"@    xaxis  tick spec type both\n");


rewind(kpoints);
dist[0]=0;
z=0;
fprintf(xmgrace,"@    xaxis  tick spec %d\n",simpt);
for(i=0;i<simpt;i++){
z=z+dist[i]; 
fscanf(kpoints,"%s",str);

fprintf(xmgrace,"@    xaxis  tick major %d, %f\n",i,z);
fprintf(xmgrace,"@    xaxis  ticklabel %d, \"%s\"\n",i,str);

}




fprintf(xmgrace,"@    yaxis  on\n");
fprintf(xmgrace,"@    yaxis  tick on\n");
fprintf(xmgrace,"@    yaxis  tick major 5\n");
fprintf(xmgrace,"@    yaxis  tick minor ticks 1\n");
fprintf(xmgrace,"@    yaxis  label \"Energy (eV)\"\n");
fprintf(xmgrace,"@    yaxis  label font 4\n");
fprintf(xmgrace,"@    yaxis  label color 1\n");
fprintf(xmgrace,"@    yaxis  label char size 1.40000\n");
fprintf(xmgrace,"@    yaxis  ticklabel on\n");
fprintf(xmgrace,"@    yaxis  ticklabel format decimal\n");
fprintf(xmgrace,"@    yaxis  ticklabel prec 1\n");
fprintf(xmgrace,"@    yaxis  ticklabel op left\n");
fprintf(xmgrace,"@    yaxis  ticklabel font 4\n");
fprintf(xmgrace,"@    yaxis  ticklabel color 1\n");
fprintf(xmgrace,"@    yaxis  ticklabel char size 1.20000\n");
fprintf(xmgrace,"@    yaxis  tick major on\n");
fprintf(xmgrace,"@    yaxis  tick minor on\n");
fprintf(xmgrace,"@    yaxis  tick out\n");
fprintf(xmgrace,"@    yaxis  tick major color 1\n");
fprintf(xmgrace,"@    yaxis  tick major linestyle 1\n");
fprintf(xmgrace,"@    yaxis  tick minor color 1\n");
fprintf(xmgrace,"@    yaxis  tick minor linestyle 1\n");
fprintf(xmgrace,"@    yaxis  tick size 0.750000\n");
fprintf(xmgrace,"@    yaxis  tick minor size 0.500000\n");
fprintf(xmgrace,"@    yaxis  tick op left\n");
fprintf(xmgrace,"@    legend off\n");
fprintf(xmgrace,"@    frame on\n");
fprintf(xmgrace,"@    frame linestyle 1\n");
fprintf(xmgrace,"@    frame color 1\n");
fprintf(xmgrace,"@    frame background color 0\n");

fprintf(xmgrace,"@    legend on\n");
fprintf(xmgrace,"@    legend loctype view\n");
fprintf(xmgrace,"@    legend 1.1, 0.5\n");
fprintf(xmgrace,"@    legend box color 1\n");
fprintf(xmgrace,"@    legend box pattern 1\n");
fprintf(xmgrace,"@    legend box linewidth 1.0\n");
fprintf(xmgrace,"@    legend box linestyle 1\n");
fprintf(xmgrace,"@    legend box fill color 0\n");
fprintf(xmgrace,"@    legend box fill pattern 1\n");
fprintf(xmgrace,"@    legend font 4\n");
fprintf(xmgrace,"@    legend char size 1.00000\n");
fprintf(xmgrace,"@    legend color 1\n");
fprintf(xmgrace,"@    legend length 4\n");
fprintf(xmgrace,"@    legend vgap 2\n");
fprintf(xmgrace,"@    legend hgap 1\n");
fprintf(xmgrace,"@    legend invert false\n");
fprintf(xmgrace,"@    s%d legend  \"Fermi energy\"\n",nbands+simpt-2);










for(j=0;j<nbands;j++){
    fprintf(xmgrace,"@target G0.S%d\n",j);
    fprintf(xmgrace,"@type xy\n");
       for(i=1;i<=kpts;i++){
          fprintf(xmgrace,"%16.4f%16.4f\n",abs[i-1],mtx[i][j]);
          }}

dist[0]=0;
z=0;
for(i=0;i<simpt;i++){
z=z+dist[i];
if(i>0 && i<simpt-1){
fprintf(xmgrace,"@target G0.S%d\n",nbands+i-1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",z,2*min);
fprintf(xmgrace,"%16.4f%16.4f\n",z,2*max);}
}

fprintf(xmgrace,"@target G0.S%d\n",nbands+i-2);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,efermi-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,efermi-ref);

if(vac1!=0){fprintf(xmgrace,"@target G0.S%d\n",nbands+i+1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,vac1-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,vac1-ref);}

if(vac2!=0){fprintf(xmgrace,"@target G0.S%d\n",nbands+i+1);
fprintf(xmgrace,"@type xy\n");
fprintf(xmgrace,"%16.4f%16.4f\n",0.0,vac2-ref);
fprintf(xmgrace,"%16.4f%16.4f\n",z,vac2-ref);}

}/*end else*/




fclose(eigenval);
fclose(kpoints);
fclose(xmgrace);
fclose(outcar);

}






