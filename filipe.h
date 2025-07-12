#ifndef	_FILIPE_H
#define	_FILIPE_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

/*-------------------date----------------------------------*/




























/*------------------------------------homo-----------------*/
float homo(float tol)
{

FILE *outcar;
int i,j,k,l,s,nef,nbands,nkpt,band;
char str[200],ch;


outcar = fopen("OUTCAR","r"); /* Arquivo ASCII, para leitura */
if(!outcar)
{
printf( "Erro na abertura do arquivo\n");
exit(0);
}

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NKPTS*/
while(strcmp(str,"ISPIN")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&s);        

nef=0;l=0;
while (fscanf(outcar,"%s",str) != EOF){
if(strcmp(str,"E-fermi")==0)nef++;                      /*verifica se o arquivo outcar esta completo e conta qtos e-fermi*/
if(strcmp(str,"Voluntary")==0)l++;
}

if(nef == 0){printf("\n\nIncomplete OUTCAR!! bye! bye! \n\n");exit(0);} /*verifica se o calculo terminou corretamente*/


rewind(outcar);
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NKPTS*/
while(strcmp(str,"NKPTS")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&nkpt);      /*lê o numero de kpts*/

float a,b[2*nkpt],c[2*nkpt],homo,lumo;

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NBANDS=*/
while(strcmp(str,"NBANDS=")!=0);

fscanf(outcar,"%d",&nbands);      /*lê o numero de bandas*/

/*---------------------------------------------------------------------------------------------------------------------------------*/
if(s==2){
for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

for(i=0;i<7;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(k=0;k<nkpt;k++){
i=0;
do{
b[k]=c[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(k=0;k<nkpt;k++){
i=0;
do{
b[k+nkpt]=c[k+nkpt];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k+nkpt]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

homo=b[0];
lumo=c[0];

for(k=1;k<2*nkpt;k++){
if(homo<b[k])homo=b[k];
if(lumo>c[k])lumo=c[k];
}

return(homo);
}

/*---------------------------------------------------------------------------------------------------------------------------------*/
else{

for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

for(k=0;k<nkpt;k++){
i=0;
do{
b[k]=c[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

homo=b[0];
lumo=c[0];

for(k=1;k<nkpt;k++){
if(homo<b[k])homo=b[k];
if(lumo>c[k])lumo=c[k];
}


return(homo);
}

if(l == 0){
printf("\n\nWARNING!!! Incomplete OUTCAR!! Check the convergence!! \n\n");} /*verifica se o calculo terminou corretamente*/


fclose(outcar);
}

/*------------------------------------lumo-----------------*/
float lumo(float tol)
{

FILE *outcar;
int i,j,k,l,s,nef,nbands,nkpt,band;
char str[200],ch;


outcar = fopen("OUTCAR","r"); /* Arquivo ASCII, para leitura */
if(!outcar)
{
printf( "Erro na abertura do arquivo\n");
exit(0);
}

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NKPTS*/
while(strcmp(str,"ISPIN")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&s);        

nef=0;l=0;
while (fscanf(outcar,"%s",str) != EOF){
if(strcmp(str,"E-fermi")==0)nef++;                      /*verifica se o arquivo outcar esta completo e conta qtos e-fermi*/
if(strcmp(str,"Voluntary")==0)l++;
}

if(nef == 0){printf("\n\nIncomplete OUTCAR!! bye! bye! \n\n");exit(0);} /*verifica se o calculo terminou corretamente*/


rewind(outcar);
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NKPTS*/
while(strcmp(str,"NKPTS")!=0);
fscanf(outcar,"%s",str);        /*pula o =*/
fscanf(outcar,"%d",&nkpt);      /*lê o numero de kpts*/

float a,b[2*nkpt],c[2*nkpt],homo,lumo;

do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra NBANDS=*/
while(strcmp(str,"NBANDS=")!=0);

fscanf(outcar,"%d",&nbands);      /*lê o numero de bandas*/

/*---------------------------------------------------------------------------------------------------------------------------------*/
if(s==2){
for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

for(i=0;i<7;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(k=0;k<nkpt;k++){
i=0;
do{
b[k]=c[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

for(i=0;i<2;i++){
do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');
}

for(k=0;k<nkpt;k++){
i=0;
do{
b[k+nkpt]=c[k+nkpt];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k+nkpt]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

homo=b[0];
lumo=c[0];

for(k=1;k<2*nkpt;k++){
if(homo<b[k])homo=b[k];
if(lumo>c[k])lumo=c[k];
}

return(lumo);
}

/*---------------------------------------------------------------------------------------------------------------------------------*/
else{

for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

do
ch = getc(outcar);           /*vai pra proxima linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');

for(k=0;k<nkpt;k++){
i=0;
do{
b[k]=c[k];
fscanf(outcar,"%f",&a);
fscanf(outcar,"%f",&c[k]);
fscanf(outcar,"%f",&a);
i++;
}while(a>tol);


for(j=0;j<nbands-i+4;j++){
do
ch = getc(outcar);           /*pula uma linha*/
while(ch!='\n');
}
}

homo=b[0];
lumo=c[0];

for(k=1;k<nkpt;k++){
if(homo<b[k])homo=b[k];
if(lumo>c[k])lumo=c[k];
}


return(lumo);
}

if(l == 0){
printf("\n\nWARNING!!! Incomplete OUTCAR!! Check the convergence!! \n\n");} /*verifica se o calculo terminou corretamente*/


fclose(outcar);
}

/*--------------------Fermi energy---------------------*/

float fermienergy(){

FILE *outcar;
float efermi;
int i,j,l,nef;
char str[150];


outcar = fopen("OUTCAR","r"); /* Arquivo ASCII, para leitura */
if(!outcar)
{
printf( "Error opening OUTCAR file\n");
exit(0);
}


nef=0;j=0;
while (fscanf(outcar,"%s",str) != EOF){
if(strcmp(str,"E-fermi")==0)nef++;                      /*verifica se o arquivo outcar esta completo e conta qtos e-fermi*/
if(strcmp(str,"Voluntary")==0)j++;
}

if(nef == 0){printf("\n\nIncomplete OUTCAR!! bye! bye! \n\n");exit(0);} /*verifica se o calculo terminou corretamente*/


rewind(outcar);

for(i=0;i<nef;i++){
do
fscanf(outcar,"%s",str);                                      /*posiciona o ponteiro após a ultima palavra E-fermi*/
while(strcmp(str,"E-fermi")!=0);
}

fscanf(outcar,"%s",str);
fscanf(outcar,"%f",&efermi);


return(efermi);

if(l == 0){
printf("\n\nWARNING!!! Incomplete OUTCAR!! Check the convergence!! \n\n");} /*verifica se o calculo terminou corretamente*/


fclose(outcar);
}

/*-------------------------scale factor----------------------*/

float scale(){

FILE *poscar;
float par;
char ch;


poscar = fopen("POSCAR","r"); /* Arquivo ASCII, para leitura */
if(!poscar)
{
printf( "Error opening POSCAR file\n");
exit(0);
}

do
ch = getc(poscar);              /*pula linhas*/
while(ch!='\n');

fscanf(poscar,"%f",&par);
/*par=par*pow(10,-10);*/
fclose(poscar);

return(par);

}



#endif
