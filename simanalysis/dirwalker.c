// dirwalker.c - pretty basic file for parsing directory output from genome sim
/*
    
    Copyright (C) 2018, Russell Posner

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
*/
#include <sys/types.h>
#include <unistd.h>
#include <dirent.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

typedef struct siminfo{
	char md5sum[256];
	int ngsval;
	int mirno;
	char states[256];
	int simno;
} siminfo;

char * printPath(siminfo * si,int protFile,int isProt){
	static char toRet[500];
	if(protFile && isProt)
	sprintf(toRet, "net.md.%s/ngs.%d/miRs.%d/states.%s/sim.%d.prot.txt",
		si->md5sum,si->ngsval,si->mirno,si->states,si->simno);
	else if (protFile){
	sprintf(toRet, "net.md.%s/ngs.%d/miRs.%d/states.%s/sim.%d.rna.txt",
		si->md5sum,si->ngsval,si->mirno,si->states,si->simno);

	}else{
					
			sprintf(toRet, "net.md.%s/ngs.%d/miRs.%d/states.%s/sim.%d.prot.degrees.txt",
				si->md5sum,si->ngsval,si->mirno,si->states,si->simno);
				}
	return toRet;
}


typedef struct siminfoarray{
	siminfo * data;
	unsigned int len;
} siminfoarray;

int hasSimInfoData(siminfo * sid,const char *nameIn,int isProt){
	char * pch;
	int offset;
	char name[500];
	char netmdstring[]="net.md.";
	char ngsstring[]="ngs.";
	char microstring[]="miRs.";
	char statestring[]="states.";
	char simstringP[]=".prot.txt";
	char simstringR[]=".rna.txt";
	char * simstring;
	simstring = (isProt) ? simstringP : simstringR;
	strncpy(name, nameIn, 500);
	if ((pch=strstr(name,".tar"))){
		return 0;
	}else if((pch=strstr(name,netmdstring))){
		strncpy(sid->md5sum, pch+strlen(netmdstring),256);
		return 1;
	}else if ((pch=strstr(name,ngsstring))){
		sid->ngsval=atoi(pch+strlen(ngsstring));
		return 2;
	}else if ((pch = strstr(name,microstring))){
		sid->mirno=atoi(pch+strlen(microstring));
		return 3;
	}else if ((pch = strstr(name,statestring))){
		strncpy(sid->states,pch+strlen(statestring),256);
		return 4;
	}else if ((pch = strstr(name,simstring))){
		pch=strtok(name,".");
		pch=strtok(NULL,".");
		sid->simno=atoi(pch);
		return 5;
	}else {
		return 0;
	}
}

int walk(const char *name, int level,siminfoarray * sid,int reset,int isProt){
	DIR *dir;
	struct dirent *entry;
	static int nRecords =0;
	static siminfo * curr;
	static siminfo temp = {0};
	char extensionP[]=".prot.txt";
	char extensionR[]=".rna.txt";
	char * extension = isProt ? extensionP : extensionR;
	if(reset){
		nRecords=0;
		curr=sid->data;
	}
	if(!(dir=opendir(name)))
		return 0;

	while((entry = readdir(dir)) != NULL){
		if (hasSimInfoData(&temp, entry->d_name,isProt)==0)
			continue;
		if (entry->d_type == DT_DIR) {
			char path[1024];
			if (entry->d_name[0]=='.') continue;

			snprintf(path, sizeof(path), "%s/%s",name,entry->d_name);
			//printf("%s/%s\n",path,entry->d_name);
			walk(path,level+1,sid,0,isProt);
		}else if (strstr(entry->d_name, extension)>0){
			//printf("%s\n",entry->d_name);
			if(curr){
				memmove(curr, &temp, sizeof(temp));
				curr++;
			}
			nRecords++;
		}

		
	}
	closedir(dir);
	return nRecords;
}

