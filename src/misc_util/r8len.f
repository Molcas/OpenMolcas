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
************************************************************************
      Subroutine R8LEN(iOpt,nData,Buf,iLen)
************************************************************************
*                                                                      *
*     Subroutine R8Len                                                 *
*                                                                      *
*     purpose: compute the packed size of a double precision           *
*              floating point number                                   *
*                                                                      *
*    Calling parameters:                                               *
*    iOpt  : Bit switch                                                *
*            bits 0-3 used for packing algorithm                       *
*    nData : number of data elements                                   *
*    Buf   : buffer containing the data                                *
*    iLen  : buffer containing the packed length                       *
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
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 2000                                 *
*                                                                      *
************************************************************************
#include "PkCtl.fh"

#include "SysDef.fh"
*
      Integer nData
      Integer iLen(*)
      Real*8  Buf(*)
*
*----------------------------------------------------------------------*
*
      If (.not.Pack) Then
         Call ICOPY(nData,[RtoB],0,iLen,1)
      Else
*        Call ERIlen(nData,PkThrs,Buf,iLen)
         Kase=iAnd(iOpt,15)
         If(Kase.eq.0) Then
            Call tcl_r8(Buf,nData, iLen, PkThrs,Init_do_setup_l)
            Init_do_setup_l=0
         Else
            iZero=8
            call iCopy(nData,[8],0,iLen,1)
            Do i=1,nData
               If(Abs(Buf(i)).lt.PkThrs) Then
                  iLen(i)=iZero
                  iZero=0
               Else
c                  iLen(i)=8
                  iZero=8
               End If
            End Do
         End If
      End If
*
*----------------------------------------------------------------------*
*
      Return
      End
