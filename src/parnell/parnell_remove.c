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
* Copyright (C) 2010, Steven Vancoillie                                *
***********************************************************************/

/* -*- mode: C -*- Time-stamp: "2010-07-02 15:38:16 stevenv"
 *
 *       File:         parnell_remove.c
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell_remove - delete files specified in a colon-separated list from
 *                       the main work directory and all subdirectories
 *
 */

#include "parnell.h"

parnell_status_t
parnell_remove (int argc, char **argv)
{
        /* go through the list, argv is guaranteed to be NULL-terminated */
        while ( argc-- ) {
                parnell_unlink(*argv++);
        }

        /* remove will always return OK, because errors on removal of files are non-critical */
        return PARNELL_OK;
}
