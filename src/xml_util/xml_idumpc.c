/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) Per-Olof Widmark                                       *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This routine does xml dump of doubles.                                 */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/*          Lund University, Sweden                                       */
/*                                                                        */
/**************************************************************************/
#include <stdio.h>
#include "molcastype.h"
#include "xmlapi.h"
#ifdef _CAPITALS_
#define xml_idumpc XML_IDUMPC
#else
#ifndef ADD_
#define xml_idumpc xml_idumpc_
#endif
#endif
void xml_idumpc(char *name,        INT *nx_name,
                char *appear,      INT *nx_appear,
                char *units, INT *nx_units, INT *Level,
                INT *data, INT *nxx, INT *nyx) {
   FILE *f;
   char  line[256];
   char  fmt[16];
   char  mfmt[16];
   int   n_name;
   int   n_appear;
   int   n_units;
   int   nx;
   int   ny;
   int   kx,ky;
   int   k;
   int   level;
   n_name=*nx_name;
   n_appear=*nx_appear;
   n_units=*nx_units;
   nx=*nxx;
   ny=*nyx;
   level=*Level;

   sprintf(fmt," %s",INT_FORMAT);
   sprintf(mfmt,"<v> %s</v>",INT_FORMAT);
   if((f=fopen(XMLDUMP,"a"))==NULL) return;
   for(k=0; k<n_name; k++) { line[k]=name[k]; if(line[k]==' ') line[k]=0; }; line[n_name]=0;
   fprintf(f,"<%s",line);
   xml_prspec(f,"appear",appear,n_appear);
   xml_prspec(f,"units",units,n_units);
   if(level>0) fprintf(f," level=\"%i\"",level);
   fprintf(f," type=\"int\"");
   if(nx>1) fprintf(f," nx=\"%i\"",nx);
   if(ny>1) fprintf(f," ny=\"%i\"",ny);
   fprintf(f,">");
   if((ny<2)&&(nx<10)) {
      if(nx==1 && ny==1)
      {
      for(kx=0; kx<nx; kx++) for(ky=0; ky<ny; ky++) fprintf(f,fmt,data[kx*ny+ky]);
      }
      else
      {
      for(kx=0; kx<nx; kx++) for(ky=0; ky<ny; ky++) fprintf(f,mfmt,data[kx*ny+ky]);
      }
   } else {
      fprintf(f,"\n");
      for(ky=0; ky<ny; ky++) {
         for(kx=0; kx<nx; kx++) {
            if(((kx%10)==0)&&(kx!=0)) fprintf(f,"\n");
            fprintf(f,fmt,data[kx*ny+ky]);
         }
         fprintf(f,"\n");
      }
   }
   fprintf(f,"</%s>\n",line);

   fclose(f);
}
