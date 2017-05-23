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
* Copyright (C) 2002, Per-Olof Widmark                                 *
***********************************************************************/
/**************************************************************************/
/*                                                                        */
/* This function is to be used by fortran routines. The purpose is to     */
/* return the address of the argument, useful in f77 dymanic memory       */
/* allocation.                                                            */
/*                                                                        */
/*------------------------------------------------------------------------*/
/*                                                                        */
/* Author:  Per-Olof Widmark                                              */
/* Written: Jan. 2002                                                     */
/*                                                                        */
/**************************************************************************/
#include "molcastype.h"
#ifdef _CAPITALS_
#define illoc ILLOC
#else
#ifndef ADD_
#define illoc illoc_
#endif
#endif
INT illoc(INT x) { return x; }
