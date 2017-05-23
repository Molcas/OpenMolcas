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
* Copyright (C) 2013, Victor P. Vysotskiy                              *
*               2015, Steven Vancoillie                                *
***********************************************************************/

/*
 *
 *       File:         parnell_cmd.c
 *       Author:       Victor Vysotskiy
 *       Date:         2013
 *
 *       parnell_cmd - parse a '?'-delimited
 *                     command line
 *
 *       History:
 *       Steven Vancoillie, University of Lund, September 2015
 *       Rewritten to be a bit simpler to read.
 */

#include "parnell.h"

#define PARNELL_DELIMIT '?'

parnell_status_t parnell_cmd(int argc, char *argv[]) {

        int argc_local = 0;
        char **argv_local = argv;
        parnell_status_t status;

        while ( argc-- ) {
                if (**argv == PARNELL_DELIMIT) {
                        if ((status = parnell (argc_local, argv_local)) != PARNELL_OK)
                                return status;
                        argv_local = argv;
                        argc_local = 0;
                }
                ++argc_local;
                ++argv;
        }
        status = parnell(argc_local, argv_local);
        return status;
}
