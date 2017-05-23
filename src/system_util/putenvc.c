/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/
/*
 *  putenv
 *
 */
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#include "molcastype.h"

#ifdef _CAPITALS_
#define putenvc PUTENVC
#else
#ifndef ADD_
#define putenvc putenvc_
#endif
#endif


INT  putenvc(char *envar) {
     FILE *MYENV;

     if(envar==NULL) return(-1);

     MYENV=fopen("molcas.env","a+");
     if (MYENV==NULL) {
        fprintf(stderr,"Unable to open molcas.env file\n");
        return(-1);
     }
     fprintf(MYENV, "%s\n",envar);
     fclose(MYENV);
     return(0);
}
