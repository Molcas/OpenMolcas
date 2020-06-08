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
*               2020, Ignacio Fdez. Galvan                             *
************************************************************************
************************************************************************
*                                                                      *
*                             A I X - I / O                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* rc=AixMv(FileName,NewName)                                           *
*                                                                      *
* Rename (move) file.                                                  *
*                                                                      *
* Input:  FileName - This file will be renamed. Given as a character   *
*                    string.                                           *
* Input:  NewName  - New name for the file. Given as a character       *
*                    string.                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Ignacio Fdez. Galvan                                        *
* Written: April 2020                                                  *
*          (based on AixRm)                                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* History:                                                             *
*                                                                      *
************************************************************************
      Integer Function AixMv(FileName,NewName)
      Implicit Integer (a-z)
      External Get_Progname
      Character*100 Get_Progname
      Integer StrnLn
      External StrnLn

      Character*(*) FileName, NewName
      Character*256 out1, out2
      Character*80 ErrTxt
*----------------------------------------------------------------------*
* Entry to AixMv                                                       *
*----------------------------------------------------------------------*
      AixMv=0
*----------------------------------------------------------------------*
* rename file                                                          *
*----------------------------------------------------------------------*
      out1=' '
      out2=' '

      Call PrgmTranslate(FileName,out1,ltmp)
      out1(ltmp+1:ltmp+1)=Char(0)
      Call PrgmTranslate(NewName,out2,ltmp)
      out2(ltmp+1:ltmp+1)=Char(0)
      rc=c_rename(out1,out2)
      If(rc.ne.0) Then
         AixMv=AixErr(ErrTxt)
         Call SysAbendMsg('AixMv','MSG: rename', ErrTxt)
         Return
      End If
*----------------------------------------------------------------------*
* Finished so return to caller                                         *
*----------------------------------------------------------------------*
      Return
      End
