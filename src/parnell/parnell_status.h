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

/* -*- mode: C -*- Time-stamp: "2010-07-03 11:46:15 stevenv"
 *
 *       File:         parnell_status.h
 *       Author:       Steven Vancoillie
 *       Date:         Spring 2010
 *
 *       parnell -- basic status handling - header file
 *
 */

#ifndef PARNELL_STATUS_H
#define PARNELL_STATUS_H

enum parnell_status_e {
        PARNELL_OK = 0,
        PARNELL_START,
        PARNELL_CONTINUE,
        PARNELL_FINISHED,
        PARNELL_ERROR,
};

typedef enum parnell_status_e parnell_status_t;

#endif /* PARNELL_STATUS_H */
