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
      Subroutine DumpRun(iRc, iOpt)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Integer       iOpt
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Integer      Lu
      Integer      iDisk
      Integer      i
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('DumpRun',ErrMsg,' ')
      End If
      iRc=0
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
* Print record information.                                            *
*----------------------------------------------------------------------*
      Write(6,*)
      Write(6,'(2a)') '------------------------------------',
     &                '------------------'
      Write(6,'(a)')  'Contents in RunFile'
      Write(6,'(2a)') '------------------------------------',
     &                '------------------'
      Write(6,'(2a)') '  Slot        Label       Disk loc. ',
     &                '  Field len.  Type'
      Write(6,'(2a)') '  ----  ----------------  ----------',
     &                '  ----------  ----'
      Do i=1,nToc
         If(TocPtr(i).ne.NulPtr) Then
            Write(6,'(i6,2x,a16,i12,2i12,i6)')
     &            i,TocLab(i),TocPtr(i),TocLen(i),TocMaxLen(i),TocTyp(i)
         End If
      End Do
      Write(6,'(2a)') '------------------------------------',
     &                '------------------'
      Write(6,*)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call DaClos(Lu)
      Return
      End
