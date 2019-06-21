#include <stdlib.h>
#include <stdio.h>
#include <string.h>

#include "dirwalker.c"

#define nMAXLENGTH 500
#define MAXLINELENGTH 10000

//GETTD(fname,fnamelen,ttdarray,maxesPre,maxesPost,meansPre,meansPost,nMax,nElem,globalttd,messPre,microPre,protPre,messPost,microPost,protPost)
extern void htmod_mp_gettd_(unsigned char *, int *, double *,
	double *,double *,double *,double *,int *,int *,double *,double *,double *,double *,double *,double *,double *,int*);


// int main(int argc, char* argv[]){
// 	char rootpath[1024];
// 	siminfo si;
// 	siminfoarray sia = {.data = NULL, .len=0};
// 	int nRecords;
// 	if(argc==2)
// 		sprintf(rootpath,"%s",argv[1]);
// 	else
// 		sprintf(rootpath,".");
// 	nRecords = walk(rootpath,0,&sia,1);
// 	sia.data = calloc(nRecords, sizeof(si));
// 	sia.len = nRecords;
// 	nRecords = walk(rootpath,0,&sia,1);
// 	siminfo * curr = sia.data;
// 	for(;nRecords;nRecords--){
// 		printf("%s\n",printPath(curr++));
// 	}
// 	free(sia.data);
// 	sia.data = NULL;

int main(int argc,char *argv[]){
	int fnamelength,i,j;
		char rootpath[1024];
	siminfo si;
	siminfoarray sia = {.data = NULL, .len=0};
	int nRecords;
	int isProt;
	if(argc==3)
		{
			sprintf(rootpath,"%s",argv[1]);
		isProt=atoi(argv[2]);
		}
	else
		perror("must specify main [pth] [isprot]");
		//sprintf(rootpath,".");
	printf("%d\n",isProt);
	nRecords = walk(rootpath,0,&sia,1,isProt);
	sia.data = calloc(nRecords, sizeof(si));
	sia.len = nRecords;
	nRecords = walk(rootpath,0,&sia,1,isProt);
	siminfo * curr = sia.data;
		double ttarray[nMAXLENGTH];
	double maxesPre[nMAXLENGTH];
	double meansPre[nMAXLENGTH];
	double maxesPost[nMAXLENGTH];
	double meansPost[nMAXLENGTH];
	int nmax=nMAXLENGTH;
	double globalttd;
	unsigned char fname[500];
	char degreefilename[500];
	char linebuf[MAXLINELENGTH];
	int nElem;
	double messPre,microPre,protPre,messPost,microPost,protPost;
	FILE * fout;
	FILE * fin;
	unsigned char outpath[500];
	sprintf(outpath, "processedauto.%s.txt",isProt ? "prot" : "rna");
	fout = fopen(outpath, "w");
	printf("outfile=%s",outpath);
	int isProtO=1;
	fprintf(fout,"netmd\tngs\tnmicro\tstates\tsimno\tnElem\tglobalttd\tmessPre\tmicroPre\tprotPre\tmessPost\tmicroPost\tprotPost\tmaxesPre[nElem]\tmaxesPost[nElem]\tmeansPre[nElem]\tmeansPost[nElem]"
		"\tttarray\tdegrees\n");
	for(i=0;i<nRecords;i++,curr++){
		if(i%100 == 0){
			printf("%d of %d\n",i,nRecords);
		}
		sprintf(fname,"%s/%s",rootpath,printPath(curr,1,isProt));
		sprintf(degreefilename,"%s/%s",rootpath,printPath(curr, 0,isProtO));
		fnamelength = strlen(fname);
	//GETTD(fname,fnamelen,ttdarray,maxesPre,maxesPost,meansPre,meansPost,nMax,nElem,globalttd,messPre,microPre,protPre,messPost,microPost,protPost)
	htmod_mp_gettd_(fname, &fnamelength, 
		ttarray, maxesPre,maxesPost, meansPre,meansPost, &nmax, &nElem,&globalttd,
		&messPre,&microPre,&protPre,&messPost,&microPost,&protPost,&isProtO);
	if(nElem){
		fprintf(fout,"%s\t%d\t%d\t%s\t%d\t%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t",curr->md5sum,curr->ngsval,curr->mirno,curr->states,
			curr->simno,nElem,globalttd,messPre,microPre,protPre,messPost,microPost,protPost);
		
		 for(j=0;j<nElem;j++){
		 	fprintf(fout,"%f\t",maxesPre[j]);
		 }
		 for(j=0;j<nElem;j++){
		 	fprintf(fout,"%f\t",maxesPost[j]);
		 }

		for(j=0;j<nElem;j++){
			fprintf(fout,"%f\t",meansPre[j]);
		}	
		for(j=0;j<nElem;j++){
			fprintf(fout,"%f\t",meansPost[j]);
		}	
		for(j=0;j<nElem-1;j++){
			fprintf(fout,"%f\t",ttarray[j]);
		}
		fprintf(fout,"%f\t",ttarray[nElem-1]);
		fin=fopen(degreefilename,"r");
		fgets(linebuf,MAXLINELENGTH,fin);
		fgets(linebuf,MAXLINELENGTH,fin);
		char * nlpos;
		if ((nlpos=strchr(linebuf,'\n'))!=NULL)
			*nlpos = '\0';
		fprintf(fout,"%s\n",linebuf);
		fclose(fin);
	}else{
		printf("Error at %s\n",fname);
	}
	}
	
	
	fclose(fout);
	free(sia.data);
}
