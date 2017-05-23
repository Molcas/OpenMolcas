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
* This routine locates a field in the runfile.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine ffxRun(iRc,Label, nData,RecTyp, iOpt)
#include "runinfo.fh"
#include "runtypes.fh"
#include "runrc.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Character*(*) Label
      Integer       nData
      Integer       RecTyp
      Integer       iOpt
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Character*16 CmpLab1
      Character*16 CmpLab2
      Integer Lu
      Integer iDisk
      Logical ok
      Integer item
      Integer i
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('ffxRun',ErrMsg,' ')
      End If
      iRc=0
*----------------------------------------------------------------------*
* Does the runfile exist? If not abend.                                *
*----------------------------------------------------------------------*
      Call f_inquire(RunName,ok)
*
* Do not return error for querying a runfile that does not exist,
* but rather return "not found", patch 6.7.263
*
*      If(.not.ok) Then
*         Call SysAbendMsg('ffxRun','RunFile does not exist',' ')
*      End If
      If(.not.ok) Then
c         write(6,*) ' Warning! In ffxRun: runfile does not exist!'
         iRc=rcNotFound
         nData=0
         RecTyp=TypUnk
         Return
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
* Locate record.                                                       *
*----------------------------------------------------------------------*
      item=-1
      Do i=1,nToc
         CmpLab1=TocLab(i)
         CmpLab2=Label
         Call UpCase(CmpLab1)
         Call UpCase(CmpLab2)
         If(CmpLab1.eq.CmpLab2) item=i
      End Do
      If(item.eq.-1) Then
         iRc=rcNotFound
         nData=0
         RecTyp=TypUnk
         Call DaClos(Lu)
         Return
      End If
      nData=TocLen(item)
      RecTyp=TocTyp(item)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call DaClos(Lu)
      Return
      End
