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
*               1998, Roland Lindh                                     *
************************************************************************
      Subroutine SORT2
************************************************************************
*                                                                      *
*     Purpose: Control phase 2 of bin sorting algorithm. First,        *
*              reload all integrals belonging to the same slice        *
*              and sort them ( c.f. SORT2A ). While reading the        *
*              keep a trace of the records which have been read.       *
*              They are used to generate forward chaining indices.     *
*              Also integrals which are of zero value are placed       *
*              back into the list. In the second step put the sorted   *
*              integrals back onto disk. If new records are needed     *
*              append them to the end of the list.                     *
*                                                                      *
*     Called from: Seward_main                                         *
*                                                                      *
*     Calls to : Sort2A,Sort2B,Sort3                                   *
*                                                                      *
*     Calling parameters: none                                         *
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
*     Srt3    : dynamic stack to control inout and output of           *
*               integral buffers                                       *
*                                                                      *
*     local data declarations: none                                    *
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
*     - modified to possibly use a virtual disk                        *
*       M. P. Fuelscher,University of Lund, Sweden, 1996               *
*     - modified to compute MxOrd                                      *
*       R. Lindh,University of Lund, Sweden, 1998                      *
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
*---------------------------------------------------------------------*
*                                                                     *
*     Stack information pertinent to phase2 of sorting                *
*     (to keep track of records that have been read in and/or         *
*      written out                                                    *
*                                                                     *
*     lStk   : size of the stack                                      *
*                                                                     *
*---------------------------------------------------------------------*
C     Parameter ( lStk = 64*1024 )
C     Integer IOStk(lStk)
C     Common /SRT3/ nStk,IOStk
C     Common /SRT3/ nStk,ip_IOStk,lStk
      Integer nStk,lStk
*
#include "stdalloc.fh"
#include "SysDef.fh"
#include "print.fh"
      Integer, Allocatable:: IOStk(:)
      Real*8, Allocatable:: SrtA(:), Scr(:)
*
*----------------------------------------------------------------------*
*     pick up the print level                                          *
*----------------------------------------------------------------------*
*
      iRout = 84
      iPrint = nPrint(iRout)
      If ( iPrint.ge.10) Write (6,*) ' >>> Enter SORT2 <<<'
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
      Call qEnter('Sort2')
*
*----------------------------------------------------------------------*
*     Initialize the IO-stack                                          *
*----------------------------------------------------------------------*
*
      Call mma_MaxINT(lStk_Max)
      lStk=Max(64*1024,lStk_Max/5)
      Call mma_allocate(IOStk,lStk,Label='IOStk')
      IOStk(:)=0
      nStk=0
*
*----------------------------------------------------------------------*
*     loop over all submatrices of 2el integrals (slices)              *
*----------------------------------------------------------------------*
*
      iOrd=0
      iBin=0
      nSym=nSyOp
      Do iSymi=1,nSym
        ib=nBs(iSymi)
        iSkip=nSkip(iSymi)
        Do jSymj=1,iSymi
          iSymj=1+ieor(iSymi-1,jSymj-1)
          iSyblj=jSymj+iSymi*(iSymi-1)/2
          jb=nBs(jSymj)
          ibj=ib*jb
          If( iSymi.eq.jSymj ) ibj=ib*(ib+1)/2
          jSkip=nSkip(jSymj)
          kSymMx=iSymi
          If( Square ) kSymMx=nSym
          Do kSymk=1,kSymMx
            kb=nBs(kSymk)
            kSkip=nSkip(kSymk)
            lSymMx=jSymj
            If( kSymk.ne.iSymi .or. Square ) lSymMx=kSymk
            Do lSyml=1,lSymMx
              kSyml=1+ieor(kSymk-1,lSyml-1)
              kSybll=lSyml+kSymk*(kSymk-1)/2
              If( ieor(iSymj-1,kSyml-1).eq.0 ) then
                lb=nBs(lSyml)
                kbl=kb*lb
                If( kSymk.eq.lSyml ) kbl=kb*(kb+1)/2
                lSkip=nSkip(lSyml)
                iSyBlk=kSybll+mxSyP*(iSyblj-1)
*
                If ( (iSkip+jSkip+kSkip+lSkip).eq.0 .and.
     &               (ibj*kbl).ne.0 ) Then
*                                                                      *
************************************************************************
*                                                                      *
                   If (RAMD) Then
*
*                     RAMD option
*
                      lSrtA=ibj*kbl
                      iBin = iBin + 1
                      iBatch=nBatch(iSyBlk)
                      iOff = RAMD_Adr(iBatch)+1
                      Call SORT2B(iBin,lSrtA,iOrd,lSrtA,RAMD_Ints(iOff),
     &                            IOStk,lStk,nStk)
*                                                                      *
************************************************************************
*                                                                      *
                   Else
*                                                                      *
************************************************************************
*                                                                      *
*                    Sorting option
*
                      nSlice=nSln(iSyBlk)
                      lSlice=lSll(iSyBlk)
                      nij=lSlice/kbl
                      mxij=ibj
*
                      lSrtA_=min(mxij,nij)*kbl
                      Call mma_allocate(SrtA,lSrtA_,Label='SrtA')
*
*----------------------------------------------------------------------*
*     Allocate space to keep all integrals of a slice in memory,       *
*     fetch them from disk and sort them.                              *
*----------------------------------------------------------------------*
*
                      Do iSlice=1,nSlice
                         iBin=iBin+1
                         lSrtA=min(mxij,nij)*kbl
                         SrtA(1:lSrtA)=0.0D0
                         Call SORT2A(iBin,lSrtA,SrtA,IOStk,lStk,nStk)
*
*----------------------------------------------------------------------*
*     Sort the IO-stack in ascending order, i.e., preference is        *
*     always given to the lowest available disk adresses.              *
*----------------------------------------------------------------------*
*
                         Call ILASRT('D',nStk,IOStk,iErr)
*
*----------------------------------------------------------------------*
*     Restore integrals on LuTwo                                       *
*----------------------------------------------------------------------*
*
                         Call SORT2B(iBin,lSrtA,iOrd,lSrtA,SrtA,
     &                               IOStk,lStk,nStk)
                         mxij=mxij-nij
                      End Do
*
                      Call mma_deallocate(SrtA)
*                                                                      *
************************************************************************
*                                                                      *
                   End If
*                                                                      *
************************************************************************
*                                                                      *
                End If
              End If
            End Do
          End Do
        End Do
      End Do
*
*----------------------------------------------------------------------*
*     Clean the remaining records in the I/O Stack                     *
*----------------------------------------------------------------------*
*
      Call mma_allocate(Scr,lStRec,Label='Scr')
      Scr(:)=0.0D0
      Do iStk=1,nStk
        iOrd=iOrd+1
        iDisk=IOStk(iStk)
        Scr(2)=DBLE(iOrd)
        iOpt=1
        Call dDAFILE(LuTwo,iOpt,Scr,lStRec,iDisk)
      End Do
      MxOrd=iOrd
      Call mma_deallocate(Scr)
      Call mma_deallocate(IOStk)
*
      If (.Not.RAMD) Then
         Call mma_deallocate(ValBin)
         Call mma_deallocate(IndBin)
      End If
*
*----------------------------------------------------------------------*
*     Turn timing OFF and exit                                         *
*----------------------------------------------------------------------*
*
      Call qExit('Sort2')
      Return
      End
