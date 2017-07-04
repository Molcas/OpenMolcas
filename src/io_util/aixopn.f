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
* Copyright (C) 1990,1991, Per-Olof Widmark                            *
************************************************************************
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixOpn(Handle,FileName)                                           *
*                                                                      *
* A file is opened for read/write operations. If the file does not     *
* exist it is automatically created.                                   *
*                                                                      *
* Input:  FileName - A character string specifying a complete path     *
*                    name or relative path name for the file. The      *
*                    name must be shorter than 128 characters.         *
*                                                                      *
* Output: Handle   - When a file is sucessfully opened, a unique file  *
*                    id is supplied by the routine. This is used for   *
*                    specifying the file to other routines.            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          S&TC, ACIS, IBM Sweden                                      *
* Written: November 1990                                               *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
* 911021 - If return code is 13 when opening file, try to open it as a *
*          read only file. Per-Olof Widmark.                           *
* 931217 - Flags obtained from system calls.                           *
* 010710 - Change using iRand to incremental number                    *
* 120304 - Use handle for the passing of 'Lu' value needed by FiM      *
*                                                                      *
************************************************************************
      Integer Function AixOpn(handle,name,translate)
      Implicit Integer (a-z)
#include "switch.fh"
#include "ctl.fh"
      Character*(*) name
      Character*256 tmp
      Logical Translate
c      Integer Length,Prog_Length
      Integer StrnLn
      External StrnLn
      Character*256 tmp1
      Character*80 ErrTxt
      External Get_Progname
      Character*100 Get_Progname
      save NVV
      data NVV /666/
*----------------------------------------------------------------------*
* Entry to AixOpn                                                      *
*----------------------------------------------------------------------*
      AixOpn=0
*----------------------------------------------------------------------*
* Check if slot in table is available                                  *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pStat,n).ne.0) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixOpn=eTmF
            Call SysWarnMsg('Aixopn','Too many opened files\n',
     *      'try to increase MxFile')
            Return
         End If
         Go To 100
      End If
      nFile=n
*----------------------------------------------------------------------*
* Strip file name and append string terminator                         *
*----------------------------------------------------------------------*
      n=Len(name)
200   If(name(n:n).eq.' ') Then
         n=n-1
         If(n.le.0) Then
            AixOpn=eBlNme
            Return
         End If
         Go To 200
      End If
      n=n+1
      If(n.ge.Len(tmp)) Then
         AixOpn=eTlFn
         Return
      End If
      tmp=name
      tmp(n:n)=Char(0)
*----------------------------------------------------------------------*
* Attempt to open file.                                                *
*----------------------------------------------------------------------*
      rc = 0
      tmp1=tmp
      ltmp=StrnLn(tmp1)
      if(translate) then
        call PrgmTranslate(tmp1,tmp,ltmp)
      endif
c       write(*,*) 'DEBUG AIXOPN: what would happen ',tmp(1:ltmp),'<'
       tmp=tmp(1:ltmp)
c
c        print *,'in=',tmp1
c        print *,'len=',ltmp
c        print *,'res=',tmp
      tmp(ltmp+1:ltmp+1)=Char(0)
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      If(handle.eq.0) Then
         rc=c_open(tmp)
      Else
         rc=c_FimOpen(tmp,handle)
         If(rc.lt.0) Then
            rc=-rc
            AixOpn=eFiMFo
         End If
      End If
      If(handle.eq.0) Then
#else
      rc=c_open(tmp)
#endif
      If(rc.lt.0) Then
         rc=AixErr(ErrTxt)
        Call SysWarnFileMsg('AixOpn',name,
     *            'MSG: open',ErrTxt)
      call SysPutsEnd()
      Call Abend()
      End If
#if defined (_HAVE_EXTRA_) && ! defined (_GA_)
      End If
#endif
      desc=rc
*----------------------------------------------------------------------*
* Attempt sucessful, update control blocks.                            *
*----------------------------------------------------------------------*
c      handle=iRand()
       NVV=NVV+100
       handle=NVV
      CtlBlk(pHndle,nFile)=handle
      CtlBlk(pDesc ,nFile)=desc
      CtlBlk(pStat ,nFile)=1
      CtlBlk(pWhere,nFile)=0
      FCtlBlk(nFile)=name
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
