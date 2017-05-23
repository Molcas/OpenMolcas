************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 1993,2000, Markus P. Fuelscher                         *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine R8PREP(nData,Buf)
************************************************************************
*                                                                      *
*     purpose: prepare the incoming data to be packed, i.e.,           *
*              round and scale the data. The modified data will        *
*              replace the incoming data.                              *
*                                                                      *
*    Calling parameters:                                               *
*    nData : number of data elements                                   *
*    Buf   : buffer containing the data                                *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    PkCtl : Packing table                                             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher and P. O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history:                                                         *
*     replacement of packing routines                                  *
*     this subroutine is obsolete !                                    *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 2000                                 *
*                                                                      *
************************************************************************
#include "PkCtl.fh"
*
      Real*8 Buf(*)
*
*----------------------------------------------------------------------*
*
      Return
c Avoid unused argument warnings
      If (.False.) Then
         Call Unused_integer(nData)
         Call Unused_real_array(Buf)
      End If
      End
