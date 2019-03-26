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
* Copyright (C) 2003, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
* This routine writes a record into the runfile.                       *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine gxWrRun(iRc,Label, Data,nData, iOpt, RecTyp)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Character*(*) Label
      Character     Data(*)
      Integer       nData
      Integer       iOpt
      Integer       RecTyp
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Integer      Lu
      Integer      iDisk
      Integer      DataAdr
      Logical      ok
      Logical      remove
      Integer      item
      Integer      i
      Integer      NewLen
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      DataAdr=-99999
      ok=.false.
      If(RecTyp.eq.TypInt) ok=.true.
      If(RecTyp.eq.TypDbl) ok=.true.
      If(RecTyp.eq.TypStr) ok=.true.
      If(RecTyp.eq.TypLgl) ok=.true.
      If(.not.ok) Then
         Call SysAbendMsg('gxWrRun',
     &                    'Argument RecTyp is of wrong type',
     &                    'Aborting')
      End If
      If(nData.lt.0) Then
         Call SysAbendMsg('gxWrRun',
     &                    'Number of data items less than zero',
     &                    'Aborting')
      End If
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('gxWrRun',ErrMsg,' ')
      End If
      iRc=0
*----------------------------------------------------------------------*
* Does the runfile exist? If not create it.                            *
*----------------------------------------------------------------------*
      Call f_inquire(RunName,ok)
      If(.not.ok) Call MkRun(iRc, iOpt)
*----------------------------------------------------------------------*
* Open runfile.                                                        *
*----------------------------------------------------------------------*
      Call OpnRun(iRc,Lu,iOpt)
*----------------------------------------------------------------------*
* Do we have space left on file?                                       *
*----------------------------------------------------------------------*
      If(RunHdr(ipItems).ge.nToc) Then
         Call DaClos(Lu)
         Call SysFilemsg('gxWrRun',
     &                   'Ran out of ToC record in RunFile',
     &                   Lu,
     &                   ' ')
         Call Abend()
      End If
*----------------------------------------------------------------------*
* Read the ToC                                                         *
*----------------------------------------------------------------------*
      iDisk=RunHdr(ipDaLab)
      Call cDaFile(Lu,icRd,TocLab,16*nToc,iDisk)
      iDisk=RunHdr(ipDaPtr)
      Call iDaFile(Lu,icRd,TocPtr,nToc,iDisk)
      iDisk=RunHdr(ipDaLen)
      Call iDaFile(Lu,icRd,TocLen,nToc,iDisk)
      iDisk=RunHdr(ipDaMaxLen)
      Call iDaFile(Lu,icRd,TocMaxLen,nToc,iDisk)
      iDisk=RunHdr(ipDaTyp)
      Call iDaFile(Lu,icRd,TocTyp,nToc,iDisk)
*----------------------------------------------------------------------*
* Reuse old field?                                                     *
*----------------------------------------------------------------------*
      item=-1
      Do i=1,nToc
         If(TocLab(i).eq.Label) item=i
      End Do
      NewLen=0
      If(item.ne.-1) Then
         remove=.false.
         If(TocTyp(item).ne.RecTyp) remove=.true.
         If(TocMaxLen(item).lt.nData)  remove=.true.
         If(remove) Then
c            write (6,*) '*******************************************'
c            write (6,'(a,a,a,i10)') 'Label=',Label,
c     &      ' expands in RUNFILE with size=',nData
c            write (6,*) '*******************************************'
c            call Abend
            TocLab(item)='Empty   '
            TocPtr(item)=NulPtr
            TocLen(item)=0
            TocTyp(item)=TypUnk
            item=-1
         Else
            DataAdr=TocPtr(item)
            NewLen=TocLen(item)
         End If
         RunHdr(ipItems)=RunHdr(ipItems)-1
      End If
*----------------------------------------------------------------------*
* Use new field?                                                       *
*----------------------------------------------------------------------*
      If(item.eq.-1) Then
         Do i=nToc,1,-1
            If(TocPtr(i).eq.NulPtr) item=i
         End Do
         If(item.eq.-1) Then
            Call DaClos(Lu)
            Call SysFilemsg('gxWrRun',
     &                      'Internal inconsistency handling RunFile',
     &                      Lu,
     &                      ' ')
            Call Abend()
         End If
         DataAdr=RunHdr(ipNext)
      End If
*----------------------------------------------------------------------*
* Write data to runfile and update header.                             *
*----------------------------------------------------------------------*
      RunHdr(ipItems)=RunHdr(ipItems)+1
      TocLab(item)=Label
      TocPtr(item)=DataAdr
      TocLen(item)=nData
      TocMaxLen(item)=max(NewLen,nData)
      TocTyp(item)=RecTyp

      iDisk=DataAdr
      Call gzRWRun(Lu,icWr,Data,nData,iDisk,RecTyp)

      if(iDisk.gt.RunHdr(ipNext)) RunHdr(ipNext)=iDisk
      iDisk=0
      Call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)
      iDisk=RunHdr(ipDaLab)
      Call cDaFile(Lu,icWr,TocLab,16*nToc,iDisk)
      iDisk=RunHdr(ipDaPtr)
      Call iDaFile(Lu,icWr,TocPtr,nToc,iDisk)
      iDisk=RunHdr(ipDaLen)
      Call iDaFile(Lu,icWr,TocLen,nToc,iDisk)
      iDisk=RunHdr(ipDaMaxLen)
      Call iDaFile(Lu,icWr,TocMaxLen,nToc,iDisk)
      iDisk=RunHdr(ipDaTyp)
      Call iDaFile(Lu,icWr,TocTyp,nToc,iDisk)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call DaClos(Lu)
      Return
      End
