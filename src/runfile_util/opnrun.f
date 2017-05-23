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
* This routine opens the runfile.                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine OpnRun(iRc,Lu, iOpt)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Integer       Lu
      Integer       iOpt
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64 ErrMsg
      Integer      iDisk
      Logical      ok
*----------------------------------------------------------------------*
* Declare external entry points                                        *
*----------------------------------------------------------------------*
      Integer isFreeUnit
      External isFreeUnit
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('OpnRun',ErrMsg,' ')
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
* Open runfile and check that file is ok.                              *
*----------------------------------------------------------------------*
      Lu=11
      Lu=isFreeUnit(Lu)

      RunHdr(ipID)=-77
      RunHdr(ipVer)=-77
      Call DaName(Lu,RunName)
      iDisk=0
      Call iDaFile(Lu,icRd,RunHdr,nHdrSz,iDisk)
      If(RunHdr(ipID).ne.IDrun) Then
         Call DaClos(Lu)
         Call SysFilemsg('gxWrRun',
     &                   'Wrong file type, not a RunFile',
     &                   Lu,
     &                   ' ')
         Call Abend()
      End If
      If(RunHdr(ipVer).ne.VNrun) Then
         Call DaClos(Lu)
         Call SysFilemsg('gxWrRun',
     &                   'Wrong version of RunFile',
     &                   Lu,
     &                   ' ')
         Call Abend()
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
