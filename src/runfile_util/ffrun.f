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
      Subroutine ffRun(Label,nData,RecTyp)
#include "runinfo.fh"
#include "runtypes.fh"
#include "runrc.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Character*(*) Label
      Integer       nData
      Integer       RecTyp
*----------------------------------------------------------------------*
* Declare local data                                                   *
*----------------------------------------------------------------------*
      Character*64  ErrMsg
      Integer       iRc
      Integer       iOpt
*----------------------------------------------------------------------*
* Call extended version routine.                                       *
*----------------------------------------------------------------------*
      iRc=0
      iOpt=0
      Call ffxRun(iRc,Label,nData,RecTyp,iOpt)
      if(iRc.eq.rcNotFound) Then
         nData=0
         RecTyp=TypUnk
      Else If(iRc.ne.0) Then
         Write(ErrMsg,'(3a)') 'Error locating field "',
     &                        Label,
     &                        '" in runfile'
         Call SysAbendMsg('ffRun',ErrMsg,' ')
      End If
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
