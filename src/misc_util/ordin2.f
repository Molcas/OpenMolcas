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
* Copyright (C) 1996, Markus P. Fuelscher                              *
************************************************************************
      Subroutine OrdIn2(iOpt,Buf,lBuf,iBatch)
************************************************************************
*                                                                      *
*    Purpose: Read a buffer of ordered two electron integrals          *
*                                                                      *
*    Calling arguments:                                                *
*    Buf    : contains on output the integrals                         *
*    lBuf   : number of integrals to be transfered                     *
*    iOpt   : option code (iOpt=1:start reading at first integral)     *
*                         (iOpt=2:continue reading)                    *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat : table of contents and auxiliary information              *
*                                                                      *
*    Local data declarations: none                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher                                                  *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
*
#include "Molcas.fh"
#include "TwoDat.fh"
*
      Real*8 Buf(*)
*---------------------------------------------------------------------*
*     If this is the first block of a symmetry batch                  *
*     get the disk disk start adress and load the buffer              *
*---------------------------------------------------------------------*
      If ( iOpt.eq.1 ) then
        iOff = RAMD_adr(iBatch)
        RAMD_next = iOff
        Call dCopy_(lBuf,RAMD_ints(iOff),1,Buf,1)
        RAMD_next = RAMD_next + lBuf
*---------------------------------------------------------------------*
*     If the number of requested integrals is larger than             *
*     the current buffer first drain the current buffer and           *
*     read as many subsequent buffers as needed                       *
*---------------------------------------------------------------------*
      Else
        iOff = RAMD_next
        Call dCopy_(lBuf,RAMD_ints(iOff),1,Buf,1)
        RAMD_next = RAMD_next + lBuf
      End If
*---------------------------------------------------------------------*
*     exit                                                            *
*---------------------------------------------------------------------*
      Return
      End
