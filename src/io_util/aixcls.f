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
* Copyright (C) 1990, Per-Olof Widmark                                 *
************************************************************************
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixCls(Handle)                                                    *
*                                                                      *
* A file is closed.                                                    *
*                                                                      *
* Input:  Handle   - This is the unique file identifier associated     *
*                    with the file. It is created by AixOpn, and must  *
*                    be used on subsequent references to the file.     *
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
*                                                                      *
************************************************************************
      Integer Function AixCls(handle)
      Implicit Integer (a-z)
#include "switch.fh"
#include "ctl.fh"
      Character*80 ErrTxt
*----------------------------------------------------------------------*
* Entry to AixCls                                                      *
*----------------------------------------------------------------------*
      AixCls=0
*----------------------------------------------------------------------*
* Check if file is opened.                                             *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pHndle,n).ne.handle) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixCls=eNtOpn
            Return
         End If
         Go To 100
      End If
      nFile=n
      desc=CtlBlk(pDesc,nFile)
*----------------------------------------------------------------------*
* Close file                                                           *
*----------------------------------------------------------------------*
      rc=c_close(desc)
      If(rc.lt.0) Then
         AixCls=AixErr(ErrTxt)
      Call SysAbendFileMsg('AixCls',FCtlBlk(nFile),
     &                    'MSG: close', ErrTxt)
      End If
*----------------------------------------------------------------------*
* Update control block                                                 *
*----------------------------------------------------------------------*
      CtlBlk(pHndle,nFile) = 0
      CtlBlk(pDesc ,nFile) = 0
      CtlBlk(pWhere,nFile) = 0
      CtlBlk(pStat,nFile)  = 0
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
