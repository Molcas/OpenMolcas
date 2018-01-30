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
* This routine reads a record from the runfile.                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine gxRdRun(iRc,Label, Data,nData, iOpt, RecTyp)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Character*(*) Label
      Real*8        Data(*)
      Integer       nData
      Integer       iOpt
      Integer       RecTyp
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer      Lu
      Integer      iDisk
      Integer      DataAdr
      Logical      ok
      Integer      item
      Integer      i
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      ok=.false.
      If(RecTyp.eq.TypInt) ok=.true.
      If(RecTyp.eq.TypDbl) ok=.true.
      If(RecTyp.eq.TypStr) ok=.true.
      If(RecTyp.eq.TypLgl) ok=.true.
      If(.not.ok) Then
         Call SysAbendMsg('gxRdRun',
     &                    'Argument RecTyp is of wrong type',
     &                    'Aborting')
      End If
      If(nData.lt.0) Then
         Call SysAbendMsg('gxRdRun',
     &                    'Number of data items less than zero',
     &                    'Aborting')
      End If
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('gxRdRun',ErrMsg,' ')
      End If
      iRc=0
*----------------------------------------------------------------------*
* Does the runfile exist? If not abort.                                *
*----------------------------------------------------------------------*
      Call f_inquire(RunName,ok)
      If(.not.ok) Then
         Call SysFilemsg('gxRdRun',
     &                   'RunFile does not exist',
     &                   Lu,
     &                   ' ')
      End If
*----------------------------------------------------------------------*
* Open runfile.                                                        *
*----------------------------------------------------------------------*
      Call OpnRun(iRc,Lu,iOpt)
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
* Find field.                                                          *
*----------------------------------------------------------------------*
      item=-1
      Do i=1,nToc
         CmpLab1=TocLab(i)
         CmpLab2=Label
c         Call Upcase(CmpLab1)
c         Call Upcase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
      If(item.eq.-1) Then
         Call DaClos(Lu)
         Write(ErrMsg,'(a,a)') 'Record not found in runfile: ',
     &                         Label
         Call SysFilemsg('gxRdRun',
     &                   ErrMsg,
     &                   Lu,
     &                   ' ')
      End If
      DataAdr=TocPtr(item)
*----------------------------------------------------------------------*
* Read data from runfile.                                              *
*----------------------------------------------------------------------*
      iDisk=DataAdr
      Call gzRWRun(Lu,icRd,Data,nData,iDisk,RecTyp)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call DaClos(Lu)
      Return
      End
