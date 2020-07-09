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
* rc=AixRm(FileName)                                                   *
*                                                                      *
* Erase file with specified file name.                                 *
*                                                                      *
* Input:  FileName - This file will be erased. Given as a character    *
*                    string.                                           *
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
      Integer Function AixRm(name)
      Implicit Integer (a-z)
#include "switch.fh"
#include "ctl.fh"
      External Get_Progname
      Character*100 Get_Progname
      Integer StrnLn
      External StrnLn

      Character*(*) name
      Character*256 tmp, out
      Character*80 ErrTxt
*----------------------------------------------------------------------*
* Entry to AixRm                                                       *
*----------------------------------------------------------------------*
      AixRm=0
*----------------------------------------------------------------------*
* Strip file name and append string terminator                         *
*----------------------------------------------------------------------*
      n=Len(name)
100   If(name(n:n).eq.' ') Then
         n=n-1
         If(n.le.0) Then
            AixRm=eBlNme
            Return
         End If
         Go To 100
      End If
      n=n+1
      If(n.ge.Len(tmp)) Then
         AixRm=eTlFn
         Return
      End If
      tmp=name
      tmp(n:n)=Char(0)
*----------------------------------------------------------------------*
* erase file                                                           *
*----------------------------------------------------------------------*
      out=' '

      Call PrgmTranslate(Name,out,ltmp)
      out(ltmp+1:ltmp+1)=Char(0)
      rc=c_remove(out)
      If(rc.ne.0) Then
         AixRm=AixErr(ErrTxt)
      Call SysAbendMsg('AixRm','MSG: delete', ErrTxt)
         Return
      End If
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
