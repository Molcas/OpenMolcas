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
* This routine creates a runfile.                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* iOpt: Bitswitch                                                      *
*    1 -- Do not create if exist.                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine MkRun(iRc, iOpt)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare dummy arguments                                              *
*----------------------------------------------------------------------*
      Integer      iRc
      Integer      iOpt
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Integer      Lu
      Integer      iDisk
      Integer      iAllow
      Logical      ok
      Integer      i
*----------------------------------------------------------------------*
* Declare external entry points                                        *
*----------------------------------------------------------------------*
      Integer      isFreeUnit
      External     isFreeUnit
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      iAllow=-1
      iAllow=iEor(iAllow,1)
      If(iAnd(iOpt,iAllow).ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('MkRun',ErrMsg,' ')
      End If
      iRc=0
*----------------------------------------------------------------------*
* Optionally do not create.                                            *
*----------------------------------------------------------------------*
      If(iAnd(iOpt,1).ne.0) Then
         Call f_inquire(RunName,ok)
         If(ok) Then
c            Write(6,*) '*** NOT creating runfile ',RunName
            Return
         End If
      End If
*     Write(6,*) '*** Creating runfile ',RunName
*----------------------------------------------------------------------*
* Create it.                                                           *
*----------------------------------------------------------------------*
      Lu=11
      Lu=isFreeUnit(Lu)

      RunHdr(ipID)=IDrun
      RunHdr(ipVer)=VNrun
      RunHdr(ipNext)=0
      RunHdr(ipItems)=0
      Call DaName(Lu,RunName)
      iDisk=0
      Call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)
      RunHdr(ipNext)=iDisk
      iDisk=0
      Call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)

      iDisk=RunHdr(ipNext)
      Do i=1,nToc
         TocLab(i)='Empty   '
         TocPtr(i)=NulPtr
         TocLen(i)=0
         TocMaxLen(i)=0
         TocTyp(i)=TypUnk
      End Do
      RunHdr(ipDaLab)=iDisk
      Call cDaFile(Lu,icWr,TocLab,16*nToc,iDisk)
      RunHdr(ipDaPtr)=iDisk
      Call iDaFile(Lu,icWr,TocPtr,nToc,iDisk)
      RunHdr(ipDaLen)=iDisk
      Call iDaFile(Lu,icWr,TocLen,nToc,iDisk)
      RunHdr(ipDaMaxLen)=iDisk
      Call iDaFile(Lu,icWr,TocMaxLen,nToc,iDisk)
      RunHdr(ipDaTyp)=iDisk
      Call iDaFile(Lu,icWr,TocTyp,nToc,iDisk)
      RunHdr(ipNext)=iDisk
      iDisk=0
      Call iDaFile(Lu,icWr,RunHdr,nHdrSz,iDisk)

      Call DaClos(Lu)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
