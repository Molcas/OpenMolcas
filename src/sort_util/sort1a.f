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
* Copyright (C) 1993,1996, Markus P. Fuelscher                         *
*               1993, Per Ake Malmqvist                                *
************************************************************************
      Subroutine SORT1A(nUt,vInt,nSqNum,nSyBlk)
************************************************************************
*                                                                      *
*     Purpose: Phase 1 of the bin sorting algorithm                    *
*             Distribute the integrals into the bins where the indices *
*             and the integral values are stored in different buffers. *
*             When a bin is completed append the buffer containing the *
*             integral values to the file LuTwo and the buffer         *
*             containing the indice to LuTmp. Before draining the      *
*             buffers they are sorted in ascending order using the     *
*             ESSL routine ISORTX.                                     *
*                                                                      *
*     Called from: PLF2,INDSFT2                                        *
*                                                                      *
*     Calls to : PKI4,PKR8,SetVec,ErrTra,ISORTX,I4Len,R8Len            *
*                                                                      *
*     Calling Parameters:                                              *
*     nUt    : number of 2el integrals in the buffers vInt,nSqNum      *
*              and nSyBlk                                              *
*     vInt   : Buffer of 2el integral values                           *
*     nSqNum : sequence number of the integral relative to             *
*              the first adress of the symmetry block                  *
*     nSyBlk : symmetry block number of an integral                    *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     TwoDat : definitions of sorting flags and address tables         *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M. P. Fuelscher and P.-AA. Malmqvist                             *
*     University of Lund, Sweden, 1993                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*     - modified to use a virtual disk                                 *
*       M. P. Fuelscher, University of Lund, Sweden, 1996              *
*                                                                      *
************************************************************************
*
      use srt2
      Implicit Real*8 (A-H,O-Z)
*
#include "Molcas.fh"
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "stdalloc.fh"
#include "print.fh"
*
      Real*8 vInt(nUt),nSqNum(nUt),nSyBlk(nUt)
*
*----------------------------------------------------------------------*
*     pick up print level                                              *
*----------------------------------------------------------------------*
*
      iRout = 81
      iPrint = nPrint(iRout)
      If ( iPrint.ge.99 ) then
        Write(6,*) ' >>> Enter SORT1A <<<'
        Call dVcPrt('nSqNum',' ',nSqNum,nUt)
        Call dVcPrt('nSyBlk',' ',nSyBlk,nUt)
        Call dVcPrt('vInt',' ',vInt,nUt)
      End If
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
C     Call QEnter('Sort1A')
      If ( RAMD ) then
        Call SORT1C(nUt,vInt,nSqNum,nSyBlk)
C       Call QExit('Sort1A')
        Return
      End If
*
      iOpt=0 ! Always tight!
*
*----------------------------------------------------------------------*
*     put the 2el integrals into the bins                              *
*----------------------------------------------------------------------*
*
      Do iUt=1,nUt
         iBin=INT(nSyBlk(iUt))
         next=nInt(iBin)+1
         lwVBin(nOffV(iBin)+next)=vInt(iUt)
         lwIBin(nOffI(iBin)+next)=INT(nSqNum(iUt))
         nInt(iBin)=next
         mInt(1,iBin)=mInt(1,iBin)+1
*
*----------------------------------------------------------------------*
*     If a bin is completed pack the buffer and put it to disc         *
*     First, sort the the integrals such that the indices are given in *
*     strictly ascending order.                                        *
*----------------------------------------------------------------------*
*
         If ( next.ge.(lBin-1) ) Then
            Call SaveBin(iBin,iOpt)
         End If
      End Do
*
*----------------------------------------------------------------------*
*     Turn timing OFF and exit                                         *
*----------------------------------------------------------------------*
*
C     Call QExit('Sort1A')
      Return
      End
