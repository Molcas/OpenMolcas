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
* Copyright (C) 2012,2013, Victor P. Vysotskiy                         *
************************************************************************
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixFsz(Handle)                                                    *
*                                                                      *
* Return size of file.                                                 *
*                                                                      *
* Input:  Handle   - This is the unique file identifier associated     *
*                    with the file. It is created by AixOpn, and must  *
*                    be used on subsequent references to the file.     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Victor P. Vysotskiy                                         *
*          Lund University, Sweden                                     *
* Written: 2012-2013                                                   *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*                                                                      *
************************************************************************
      Integer  Function AixFsz(handle)
      Implicit Integer (a-z)
#include "switch.fh"
#include "ctl.fh"
      Character*80 ErrTxt

*----------------------------------------------------------------------*
* Entry to AixCls                                                      *
*----------------------------------------------------------------------*
      AixFsz=0
*----------------------------------------------------------------------*
* Check if file is opened.                                             *
*----------------------------------------------------------------------*
      n=1
100   If(CtlBlk(pHndle,n).ne.handle) Then
         n=n+1
         If(n.gt.MxFile) Then
            AixFsz=eNtOpn
            Return
         End If
         Go To 100
      End If
      nFile=n
      desc=CtlBlk(pDesc,nFile)
*----------------------------------------------------------------------*
* Get file size                                                        *
*----------------------------------------------------------------------*
      rc=c_stat(desc)
      If(rc.lt.0) Then
         AixFsz=AixErr(ErrTxt)
      Call SysAbendFileMsg('AixFsz',FCtlBlk(nFile),
     &                    'MSG: close', ErrTxt)
      End If
      AixFsz=rc
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
