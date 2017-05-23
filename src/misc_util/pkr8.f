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
      Subroutine PKR8(iOpt,nData,nByte,InBuf,OutBuf)
************************************************************************
*                                                                      *
*     purpose: pack the incoming buffer of double prcision floating    *
*              point number into the outgoing buffer.                  *
*                                                                      *
*    Calling parameters:                                               *
*    iOpt  : Bit switch                                                *
*            bits 0-3 used for packing algorithm                       *
*    nData : number of unpacked integrals                              *
*    nByte : length of pack array in bytes                             *
*    InBuf : on input unpacked double precision floating point numbers *
*            (the content of this buffer is destroyed on output!)      *
*    OutBuf: contains packed numbers  on output                        *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    PkCtl : Packing table                                             *
*                                                                      *
*    Local data declarations:                                          *
*    Round : Table of constant used for rounding                       *
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
*
      Integer iOpt,nData,nByte
      Real*8 InBuf(*)
      Real*8 OutBuf(*)

*----------------------------------------------------------------------*
*     If no packing is desired copy the data from                      *
*     the input buffer to the output buffer                            *
*----------------------------------------------------------------------*

      If (.not.Pack) Then
        call dcopy_(nData,InBuf,1,OutBuf,1)
        nByte=8*nData
        Return
      End If

*----------------------------------------------------------------------*
*     Pack the data and transfer the result from                       *
*     the input buffer to the output buffer                            *
*----------------------------------------------------------------------*

*     Call pkERI(nData,PkThrs,nByte,InBuf,OutBuf)
      Kase=iAnd(iOpt,15)
      If(Kase.eq.0) Then
         Call tce_r8(InBuf,nData, OutBuf,nComp, PkThrs,Init_do_setup_e)
         Init_do_setup_e=0
         nByte=nComp
      Else
         Call rle_r8(InBuf,nData, OutBuf,nComp, PkThrs)
         nByte=8*nComp
      End If

*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*

      Return
      End
