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
      Subroutine UPKR8(iOpt,nData,nByte,InBuf,OutBuf)
************************************************************************
*                                                                      *
*     purpose: decode pack double precision floating point numbers     *
*                                                                      *
*    Calling parameters:                                               *
*    iOpt  : Bit switch                                                *
*            bits 0-3 used for packing algorithm                       *
*    nData : number of unpacked double precision floating point numbers*
*    nByte : length of pack array in bytes                             *
*    OutBuf: contains on output unpacked numbers                       *
*    InBuf : contains packed integral on input                         *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    PkCtl : Packing table                                             *
*                                                                      *
*    Local data declarations: none                                     *
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

      Integer iOpt,nData,nByte
      Real*8  InBuf(*)
      Real*8  OutBuf(*)
*
*----------------------------------------------------------------------*
*     If data have not been packed copy them from                      *
*     the input buffer to the output buffer                            *
*----------------------------------------------------------------------*
*
      If (.not.Pack) Then
        call dcopy_(nData,InBuf,1,OutBuf,1)
        nByte=8*nData
        Return
      End If
*
*----------------------------------------------------------------------*
*     Unpack the data and transfer the result from                     *
*     the input buffer to the output buffer                            *
*----------------------------------------------------------------------*
*
*     Call upkERI(nData,PkThrs,nByte,InBuf,OutBuf)
      Kase=iAnd(iOpt,15)
      If(Kase.eq.0) Then
         Call tcd_r8(InBuf,nComp, OutBuf,nData, PkThrs,Init_do_setup_d)
         Init_do_setup_d=0
         nByte=nComp
      Else
         Call rld_r8(InBuf,nComp, OutBuf,nData, PkThrs)
         nByte=8*nComp
      End If

*
*----------------------------------------------------------------------*
*
      Return
      End
