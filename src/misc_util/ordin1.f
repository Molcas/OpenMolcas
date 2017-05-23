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
* Copyright (C) 1993, Markus P. Fuelscher                              *
*               1993, Per-Olof Widmark                                 *
************************************************************************
      Subroutine OrdIn1(iOpt,Buf0,lBuf0,iBatch)
************************************************************************
*                                                                      *
*     Purpose: Read a buffer of ordered two electron integrals         *
*                                                                      *
*     Note:    This subroutine has internal buffers.                   *
*                                                                      *
*    Calling parameters:                                               *
*    Buf0   : contains on output the integrals                         *
*    lBuf0  : number of integrals to be transfered                     *
*    iOpt   : option code (iOpt=1:start reading at first integral)     *
*                         (iOpt=2:continue reading)                    *
*    rc     : return code                                              *
*                                                                      *
*    Global data declarations (Include files) :                        *
*    TwoDat : table of contents and auxiliary information              *
*    TowRc  : Table of return code                                     *
*    TwoDef : definitions of record structure                          *
*    TwoBuf : save area for buffering of two electron integrals        *
*                                                                      *
*    Local data declarations: none                                     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.O. Widmark                                 *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Integer (A-Z)
      External C2R8
      Real*8 C2R8
*
#include "Molcas.fh"
#include "TwoDat.fh"
#include "TwoDef.fh"
*---------------------------------------------------------------------*
*                                                                     *
*     This is a buffer exculsively used for I/O buffering             *
*     of the ordered 2-el. integrals file                             *
*                                                                     *
*     !!!     The current version uses one buffer only     !!!        *
*     !!! Double buffering and asynchronous I/O is preseen !!!        *
*                                                                     *
*---------------------------------------------------------------------*
      Character*1 Buf1
      Common /TWOBUF/ Buf1(8*lStRec)

#include "SysDef.fh"
*
      Real*8 Buf0(*)
      Save kOpt
*     Call qEnter('OrdIn1')
*---------------------------------------------------------------------*
*     Fetch the unit number, disk start adress and pointers           *
*---------------------------------------------------------------------*
      LuTwo=AuxTwo(isUnit)
      iDisk1=AuxTwo(isDaDa)
      isBuf1=AuxTwo(isUpk8)
      lBuf1=AuxTwo(islBf1)
*---------------------------------------------------------------------*
*     If this is the first block of a symmetry batch                  *
*     get the disk disk start adress and load the buffer              *
*---------------------------------------------------------------------*
      If ( iOpt.eq.1 ) then
        iDisk1=TocTwo(isDAdr+iBatch-1)
        jOpt=2
        Call dDAFILE(LuTwo,jOpt,Buf1,lStRec,iDisk1)
*---------------------------------------------------------------------*
*                                                                     *
*       Note:                                                         *
*       If the records are organized in sequential order              *
*       (SORT3 in SEWARD is activated)                                *
*       deactivate the update of the disk adress                      *
*                                                                     *
*       iDisk1=NINT(C2R8(Buf1(1)))                                    *
*---------------------------------------------------------------------*
        lBuf1=nint(C2R8(Buf1(17)))
        kOpt=nint(C2R8(Buf1(25)))
        isBuf1=lTop*RtoB+1
      End If
*---------------------------------------------------------------------*
*     If the number of requested integrals is smaller than            *
*     the current buffer transfer data                                *
*---------------------------------------------------------------------*
      If ( lBuf0.le.lBuf1 ) then
        Call UPKR8(kOpt,lBuf0,nByte,Buf1(isBuf1),Buf0)
        isBuf1=isBuf1+nByte
        lBuf1=lBuf1-lBuf0
*---------------------------------------------------------------------*
*     If the number of requested integrals is larger than             *
*     the current buffer first drain the current buffer and           *
*     read as many subsequent buffers as needed                       *
*---------------------------------------------------------------------*
      Else
        Call UPKR8(kOpt,lBuf1,nByte,Buf1(isBuf1),Buf0)
        isBuf0=lBuf1+1
        nleft=lBuf0-lBuf1
        Do while ( nleft.gt.0 )
          jOpt=2
          Call dDAFILE(LuTwo,jOpt,Buf1,lStRec,iDisk1)
*---------------------------------------------------------------------*
*                                                                     *
*         Note:                                                       *
*         If the records are organized in sequential order            *
*         (SORT3 in SEWARD is activated)                              *
*         deactivate the update of the disk adress                    *
*                                                                     *
*         iDisk1=NINT(C2R8(Buf1(1)))                                  *
*---------------------------------------------------------------------*
          lBuf1=nint(C2R8(Buf1(17)))
          ncopy=min(nleft,lBuf1)
          kOpt=nint(C2R8(Buf1(25)))
          isBuf1=lTop*RtoB+1
          Call UPKR8(kOpt,ncopy,nByte,Buf1(isBuf1),Buf0(isBuf0))
          isBuf0=isBuf0+ncopy
          isBuf1=lTop*RtoB+1+nByte
          lBuf1=lBuf1-ncopy
          nleft=nleft-ncopy
        End Do
      End If
*---------------------------------------------------------------------*
*     Update pointer to next disk adress and integral to unpack       *
*---------------------------------------------------------------------*
      AuxTwo(isDaDa)=iDisk1
      AuxTwo(isUpk8)=isBuf1
      AuxTwo(islBf1)=lBuf1
*---------------------------------------------------------------------*
*     exit                                                            *
*---------------------------------------------------------------------*
*     Call qExit('OrdIn1')
      Return
      End
