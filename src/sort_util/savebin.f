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
*               1998, Roland Lindh                                     *
************************************************************************
      Subroutine SaveBin(iBin,ValBin,IndBin,lIndx,lInts,l_Bin,iOpt)
************************************************************************
*                                                                      *
*     Purpose: Phase 1 of the bin sorting algorithm                    *
*             Once a bin is filled it has to be dumped to disk.        *
*             First, however, we like to sort and pack the buffers.    *
*             A further complication arises due to the fact that a     *
*             bin is shorter than the records as used in the wave-     *
*             function codes.                                          *
*                                                                      *
*     Called from: SORT1A and SORT1B                                   *
*                                                                      *
*     Calls to : PKI4,PKR8,DaFile,SetVec,ErrTra,ISORTX,I4Len,R8Len     *
*                                                                      *
*     Calling Parameters:                                              *
*     iBin   : Bin numberto be saved                                   *
*                                                                      *
*     Global data declarations (Include files) :                       *
*     TwoDef  : definitions of the record structure                    *
*     Srt0    : common block containing information pertinent to       *
*               the calculation of 2el integral sequence numbers       *
*     Srt1    : common block containing information the number of      *
*               bins and partitioning of symmetry blocks               *
*     Srt2    : common block containing information pertinent to       *
*               the bin sorting algorithm                              *
*     WSpc    : dynamic work space                                     *
*                                                                      *
*     local data declarations:                                         *
*     IndBin : temporary array which contains the packed labels        *
*     IndBin : temporary array which contains the packed 2el integrals *
*     lIndx  : temporary array used to store the length                *
*              of the packed indices                                   *
*     lInts  : temporary array used to store the length                *
*              of the packed integral                                  *
*     Scr    : temporary array used for initializing records           *
*                                                                      *
*     Modified to remove sorting step, R. Lindh, March '98             *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      Implicit Real*8 (A-H,O-Z)
*
#include "TwoDef.fh"
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"
#include "srt2.fh"
#include "WrkSpc.fh"

#include "SysDef.fh"
#include "PkCtl.fh"

      Dimension ValBin(l_Bin),IndBin(l_Bin)
      Dimension lIndx(l_Bin),lInts(l_Bin)
      Integer iScr(lStRec)
      Real*8   Scr(lStRec)
      Integer rc
*
*----------------------------------------------------------------------*
*     Turn timing ON                                                   *
*----------------------------------------------------------------------*
*
C     Call QEnter('Savebin')
*
*----------------------------------------------------------------------*
*         as the packed integral labels add an extra 1-2 Byte          *
*         disk space per integral we have to adjust the record         *
*         length of LuTmp to the different machines.                   *
*----------------------------------------------------------------------*
      idiv = ItoB/2
      If ( Pack ) idiv = idiv/2
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      nInts=nInt(iBin)
      iOff1=nOff1(iBin)
      iOff2=nOff2(iBin)
*----------------------------------------------------------------------*
*         precompute the packed length of a record                     *
*----------------------------------------------------------------------*
      Call I4LEN(nInts,iwork(iOff2+1),lIndx)
      Call R8LEN(iOpt,nInts,work(iOff1+1),lInts)
      mxIRec=((lDaRec/idiv)-lTop-1)*ItoB
      mxVRec=(lDaRec-lTop-1)*RtoB
      nSave=0
      lIRec=0
      lVRec=0
      Do iSave=1,nInts
         lIRec=lIRec+lIndx(iSave)
         lVRec=lVRec+lInts(iSave)
         If ( lIRec.lt.mxIRec .and. lVRec.lt.mxVRec ) nSave=iSave
      End Do
      lIRec=0
      lVRec=0
      Do iSave=1,nSave
         lIRec=lIRec+lIndx(iSave)
         lVRec=lVRec+lInts(iSave)
      End Do
*----------------------------------------------------------------------*
*         now pack and check the packed record length again            *
*----------------------------------------------------------------------*
      Call PKI4(nSave,lIBin,iwork(iOff2+1),IndBin(lTop+1))
      mInds=(lIBin+ItoB-1)/ItoB
      mInt(3,iBin)=mInt(3,iBin)+mInds
      If ( lIBin.gt.mxIRec ) then
         rc=001
         Write(6,*)
         Write(6,'(2X,A,I3.3,A)')
     &   '*** Error in SAVEBIN ***'
         Write(6,'(2X,A)') 'An inconsistency has been deteced'
         Write(6,'(2X,A)') 'lIRec > mxIRec '
         Write(6,*)
         Call qTrace
         Call xFlush(6)
         Call Abend
      End If
      If ( lIBin.ne.lIRec ) then
         rc=002
         Write(6,*)
         Write(6,'(2X,A,I3.3,A)')
     &   '*** Error in SAVEBIN ***'
         Write(6,'(2X,A)') 'An inconsistency has been deteced'
         Write(6,'(2X,A)') 'lIBin # lIRec'
         Write(6,*)
         Call qTrace
         call xFlush(6)
         Call Abend
      End If
      Call PKR8(iOpt,nSave,lVBin,Work(iOff1+1),ValBin(lTop+1))
      mInts=(lVBin+RtoB-1)/RtoB
      mInt(2,iBin)=mInt(2,iBin)+mInts
      If ( lVBin.gt.mxVRec ) then
         rc=003
         Write(6,*)
         Write(6,'(2X,A,I3.3,A)')
     &   '*** Error in SAVEBIN ***'
         Write(6,'(2X,A)') 'An inconsistency has been deteced'
         Write(6,'(2X,A)') 'lVRec > mxVRec '
         Write(6,*)
         Call qTrace
         Call xFlush(6)
         Call Abend
      End If
      If ( lVBin.ne.lVRec ) then
         rc=004
         Write(6,*)
         Write(6,'(2X,A,I3.3,A)')
     &   '*** Error in SAVEBIN ***'
         Write(6,'(2X,A)') 'An inconsistency has been deteced'
         Write(6,'(2X,A)') 'lVBin # lVRec'
         Write(6,*)
         Call qTrace
         call xFlush(6)
         Call Abend()
      End If
*----------------------------------------------------------------------*
*     write the record header                                          *
*----------------------------------------------------------------------*
      IndBin(1)=iDIBin(2,iBin)
      IndBin(2)=lIBin+lTop*ItoB
      IndBin(3)=nSave
      ValBin(1)=iDVBin(2,iBin)
      ValBin(2)=lVBin+lTop*RtoB
      ValBin(3)=nSave
*----------------------------------------------------------------------*
*     save the record on disk                                          *
*----------------------------------------------------------------------*
      iDaTmp=iDIBin(1,iBin)
      iDaTwo=iDVBin(1,iBin)
      iOptIO=1
*
*     If a new sector mark it on the disk to make sure that it is
*     continuous.
*
      If ( mod(nRec(iBin),nSect).eq.0 ) then
         iDaTmp=mDaTmp
         iDaTwo=mDaTwo
         Zero=0.0d0
*
         Call ICopy(lStRec,0,0,iScr,1)
         Call iDAFILE(LuTmp,iOptIO,iScr,(lStRec/idiv),mDaTmp)
*
         call dcopy_(lStRec,Zero,0,Scr,1)
         Call dDAFILE(LuTwo,iOptIO,Scr,lStRec,mDaTwo)
*
         iDIBin(2,iBin)=iDaTmp
         iDVBin(2,iBin)=iDaTwo
         If ( iDVBin(3,iBin).lt.0 ) iDVBin(3,iBin)=iDaTwo
      End If
*
      Call iDAFILE(LuTmp,iOptIO,IndBin,(lDaRec/idiv),iDaTmp)
      iDIBin(1,iBin)=iDaTmp
      Call dDAFILE(LuTwo,iOptIO,ValBin,lDaRec,iDaTwo)
      iDVBin(1,iBin)=iDaTwo
      nRec(iBin)=nRec(iBin)+1
*----------------------------------------------------------------------*
*         cleanup                                                      *
*----------------------------------------------------------------------*
      nKeep=nInts-nSave
      If ( nKeep.gt.0 ) Then
         Do i=1,nKeep
            iwork(iOff2+i)=iwork(iOff2+nSave+i)
            work(iOff1+i)=work(iOff1+nSave+i)
         End Do
      End If
      nInt(iBin)=nKeep
*
*----------------------------------------------------------------------*
*     Turn timing OFF and exit                                         *
*----------------------------------------------------------------------*
*
C     Call QExit('SaveBin')
      Return
      End
