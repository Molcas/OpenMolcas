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
      Subroutine SORT2B(iBin,nInts,iOrd,lSrtA,SrtArr,IOStk,lStk,nStk)
************************************************************************
*                                                                      *
*     Purpose: Store the sorted integrals on disk.                     *
*              Every record that is written contains a header of       *
*              length lTop which containd the number of integrals      *
*              on the record and its ordering number.                  *
*                                                                      *
*     Called from: Sort2                                               *
*                                                                      *
*     Calls to : PkI4,PkR8,DaFile                                      *
*                                                                      *
*     Calling parameters:                                              *
*     iBin   : Slice number                                            *
*     nInts  : no. of 2el integrals in slice                           *
*     iOrd   : ordering number of record                               *
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
*     PkVal  : I/O buffer contains packed integral values              *
*     IntLen : this buffer contains the length of the integrals        *
*              to be pack into the next buffer                         *
*                                                                      *
**** M. Fuelscher and P.-Aa. Malmqvist, Univ. of Lund, Sweden, 1991 ****
*
      use srt2
      Implicit Real*8 (A-H,O-Z)
*
#include "TwoDat.fh"
#include "srt0.fh"
#include "srt1.fh"

#include "SysDef.fh"
#include "print.fh"
*
      Dimension PkVal(lStRec)
      Dimension IntLen(4*lStRec)
      Dimension SrtArr(lSrtA)
      Integer IOStk(lStk)
*
*----------------------------------------------------------------------*
*     pick up the print level                                          *
*----------------------------------------------------------------------*
*
#ifdef _DEBUGPRINT_
      iRout = 86
      iPrint = nPrint(iRout)
      If ( iPrint.gt.5 ) then
        Write (6,*) ' >>> Enter SORT2B <<<'
        Write (6,*) ' iBin  ',iBin
        Write (6,*) ' LSrtA ',lSrtA
        Write (6,*) ' nInts ',nInts
      End If
      If ( iPrint.ge.10 ) then
        Call iVcPrt('stack of free records',' ',IOStk,nStk)
      End If
#endif
*
*----------------------------------------------------------------------*
*     Transfer integrals from the sorting area to I/O buffer           *
*----------------------------------------------------------------------*
      nRec(iBin)=0
      kStk=0
      nSaved=0
      iOpt=iDVBin(4,iBin)
      Do while ( nSaved.lt.nInts )
         iStart=nSaved+1
         iEnd=Min(nSaved+4*lStRec,nInts)
         Call R8Len(iOpt,1+iEnd-iStart,SrtArr(iStart),IntLen)
         mxVRec=(lStRec-lTop-1)*RtoB
         nSave=0
         lVRec=0
         Do iSave=1,1+iEnd-iStart
            lVRec=lVRec+IntLen(iSave)
            If ( lVRec.lt.mxVRec ) nSave=iSave
         End Do
         If ( nSave.eq.0 ) then
            iRC=001
            Write(6,*)
            Write(6,'(2X,A,I3.3,A)')
     &      '*** Error in SORT2B ***'
            Write(6,'(2X,A)') 'nSave = 0'
            Write(6,*)
            Call xFlush(6)
            Call Abend()
         End If
         lVRec=0
         Do iSave=1,nSave
            lVRec=lVRec+IntLen(iSave)
         End Do
         If ( lVRec.gt.mxVRec ) then
            iRC=002
            Write(6,*)
            Write(6,'(2X,A,I3.3,A)')
     &      '*** Error in SORT2B ***'
            Write(6,'(2X,A)') 'An inconsistency has been deteced'
            Write(6,'(2X,A)') 'lVRec > mxVRec '
            Write(6,*)
            Call xFlush(6)
            Call Abend()
         End If
*----------------------------------------------------------------------*
*     Pack integrals                                                   *
*----------------------------------------------------------------------*
         Call PKR8(iOpt,nSave,llVRec,SrtArr(iStart),PkVal(lTop+1))
         If ( llVRec.ne.lVRec ) then
            iRC=003
            Write(6,*)
            Write(6,'(2X,A,I3.3,A)')
     &      '*** Error in SORT2B ***'
            Write(6,'(2X,A)') 'An inconsistency has been deteced'
            Write(6,'(2X,A)') 'llVBin # lVRec'
            Write(6,*)
            call xFlush(6)
            Call Abend()
         End If
         iOrd=iOrd+1
         PkVal(1)=0
         PkVal(2)=iOrd
         PkVal(3)=nSave
         PkVal(4)=iOpt
*----------------------------------------------------------------------*
*     Write out the new record                                         *
*----------------------------------------------------------------------*
         If ( kStk.lt.nStk) then
*---------- use an old record, get address from IOStk
            kStk=kStk+1
            iDaTwo=IOStk(kStk)
         Else
            iDaTwo=mDaTwo
            iOptIO=0
            Call dDAFILE(LuTwo,iOptIO,[0.0d0],lStRec,mDaTwo)
         End If
#ifdef _DEBUGPRINT_
         If ( iPrint.ge.10 ) then
           Write (6,*) ' write record: iOrd,iDaTwo ',iOrd,iDaTwo
         End If
#endif
*
*------- Write out the buffer
*
         iOptIO=1
         Call dDAFILE(LuTwo,iOptIO,PkVal,lStRec,iDaTwo)
*
         nRec(iBin)=nRec(iBin)+1
         nSaved=nSaved+nSave
*
      End Do
*----------------------------------------------------------------------*
*     cleanup the stack of Record addresses                            *
*----------------------------------------------------------------------*
      jStk=0
*
*     Do this only if the number of records was reduced
*
      If ( kStk.lt.nStk ) then
*
*        Move unused record addresses to the start of IOStk!
*
         Do iStk=kStk+1,nStk
            jStk=jStk+1
            IOStk(jStk)=IOStk(iStk)
         End Do
      End If
      nStk=jStk
*
*----------------------------------------------------------------------*
*     exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
