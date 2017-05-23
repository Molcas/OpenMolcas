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
* Copyright (C) 1991, Markus P. Fuelscher                              *
*               1991, Per Ake Malmqvist                                *
************************************************************************
      Subroutine SORT2A(iBin,lSrtA,SrtArr,ValBin,IndBin,l_Bin,
     &                  IOStk,lStk,nStk)
************************************************************************
*                                                                      *
*     Purpose: Reload all integral belonging to a given slice          *
*              and sort them.                                          *
*                                                                      *
*     Called from: Sort2                                               *
*                                                                      *
*     Calls to : UPkI4,UPkR8,DaFile                                    *
*                                                                      *
*     Calling parameters:                                              *
*     iBin   : Slice number                                            *
*     SrtArr : Work space to keep the 2el integrals                    *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*     Srt3    : dynamic stack to control inout and output of           *
*               integral buffers                                       *
*                                                                      *
*     local data declarations:                                         *
*     PkVBin : I/O buffer contains packed integral values              *
*     PkIBin : I/O buffer contains packed ordering numbers             *
*     ValBin : contains unpacked integral values                       *
*     IndBin : contains unpacked ordering numbers                      *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      Implicit Real*8 (A-H,O-Z)
*
#include "TwoDat.fh"
#include "TwoDef.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "srt2.fh"
#include "SysDef.fh"
#include "print.fh"
#include "PkCtl.fh"
#include "warnings.fh"
*
      Dimension SrtArr(lSrtA)
      Dimension PkVBin(lStRec)
      Integer   PkIBin(lStRec)
      Dimension ValBin(l_Bin),IndBin(l_Bin)
      Integer IOStk(lStk)
*
*----------------------------------------------------------------------*
*     as the packed integral labels add an extra 1-2 Byte              *
*     disk space per integral we have to adjust the record             *
*     length of LuTmp to the different machines.                       *
*----------------------------------------------------------------------*
*
      idiv = ItoB/2
      If ( Pack ) idiv = idiv/2
      mStRec=(lStRec/idiv)
      mDaRec=(lDaRec/idiv)
*
*----------------------------------------------------------------------*
*     pick up the print level                                          *
*----------------------------------------------------------------------*
*
      iRout = 85
      iPrint = nPrint(iRout)
      If ( iPrint.ge.10) then
        Write (6,*) ' >>> Enter SORT2A <<<'
        Write (6,*) ' iBin  ',iBin
        Write (6,*) ' lSrtA ',lSrtA
      End If
*
      iZero=lSrtA-mInt(1,iBin)
      iInt=mInt(2,iBin)*RtoB
      iInd=mInt(3,iBin)
      iP_Storage=(iZero+iInt+RtoB)/RtoB
      iI_Storage=( iInd+iInt+RtoB)/RtoB
      If (iP_Storage .le. iI_Storage) Then
         iDVBin(4,iBin) = 0  ! Dense mode.
C        Write (*,*) 'Mode: Dense'
      Else
         iDVBin(4,iBin) = 1  ! Sparse mode.
C        Write (*,*) 'Mode: Sparse'
      End If
*
#ifdef _DEBUG_
      Write (6,*)
      Write (6,*) 'Processing slice                   :',iBin
      Write (6,*) 'Actual number of non-zero integrals:', mInt(1,iBin)
      Write (6,*) 'Effective number of integrals      :', mInt(2,iBin)
      Write (6,*) 'Effective number of indicies       :', mInt(3,iBin)
      Write (6,*) 'Total number of integrals          :',lSrtA
      Write (6,*) 'Packed storage                     :',iP_Storage
      Write (6,*) 'Indexed storage                    :',iI_Storage
#endif
*
*----------------------------------------------------------------------*
*     Start reading packed buffers                                     *
*----------------------------------------------------------------------*
*
      iDaTmp=iDIBin(2,iBin)
      iDaTwo=iDVBin(2,iBin)
      Do while ( iDaTmp.ge.0 )
        nStk=nStk+1
        If ( nStk.gt.lStk ) then
          iRC=001
          Write(6,*)
          Write(6,'(2X,A,I3.3,A)')
     &    '*** Error in SORT2A ***'
          Write(6,'(2X,A)') 'nStk exceeds limits (nStk>lStk)'
          Write(6,'(2X,A,I8)') 'nStk =',nStk
          Write(6,'(2X,A,I8)') 'lStk =',lStk
          Write(6,'(2X,A,I8)') 'iBin =',iBin
          Write(6,*)
          Write(6,*) 'Action: rerun with a larger MOLCAS_MEM'
          Call qTrace
          Call Quit(_RC_MEMORY_ERROR_)
        End If
        IOStk(nStk)=iDaTwo
        iOpt=2
        If ( iPrint.ge.10 ) then
          Write (6,*) ' read records: iDaTmp,iDaTwo ',iDaTmp,iDaTwo
        End If
        Call iDAFILE(LuTmp,iOpt,PkIBin,mStRec,iDaTmp)
        Call dDAFILE(LuTwo,iOpt,PkVBin,lStRec,iDaTwo)
        ist1=lTop
        ist2=lTop
        Do iSec=1,nSect
          nInts1=PkIBin(ist1-1)
          nInts2=Int(PkVBin(ist2-1))
          If ( nInts1.ne.nInts2 ) then
            iRC=002
            Write(6,*)
            Write(6,'(2X,A,I3.3,A)')
     &      '*** Error in SORT2A ***'
            Write(6,'(2X,A)') 'An inconsistency has been deteced'
            Write(6,'(2X,A)') 'nInts1#nInts2'
            Write(6,*)
            Call qTrace
            Call xFlush(6)
            Call Abend
          End If
          nInts=nInts1
          If ( nInts.gt.l_Bin ) then
            iRC=003
            Write(6,*)
            Write(6,'(2X,A,I3.3,A)')
     &      '*** Error in SORT2A ***'
            Write(6,'(2X,A)') 'An inconsistency has been deteced'
            Write(6,'(2X,A)') 'nInts>l_Bin'
            Write(6,*)
            Call qTrace
            Call xFlush(6)
            Call Abend
          End If
          If ( nInts.gt.0 ) then
*
*----------------------------------------------------------------------*
*     Unpack buffers                                                   *
*----------------------------------------------------------------------*
*
            Call UPKI4(nInts,lIBin,PkIBin(ist1+1),IndBin)
            iOpt=0 ! Always tight mode
            Call UPKR8(iOpt,nInts,lVBin,PkVBin(ist2+1),ValBin)
*
*----------------------------------------------------------------------*
*     Sort 2el integrals                                               *
*----------------------------------------------------------------------*
*
            Do indx=1,nInts
              SrtArr(IndBin(indx))=ValBin(indx)
            End Do
            ist1=ist1+mDaRec
            ist2=ist2+lDaRec
          End If
        End Do
*
*----------------------------------------------------------------------*
*     Get the disk adress of the next record                           *
*----------------------------------------------------------------------*
*
        iDaTmp=PkIBin(1)
        iDaTwo=Int(PkVBin(1))
      End Do
      If ( iPrint.ge.99 ) Call dVcPrt('sorted ERIs',' ',SrtArr,lSrtA)
*     Call GetMem('Sort2a','Check',' ',junk,junk)
*
      Return
      End
