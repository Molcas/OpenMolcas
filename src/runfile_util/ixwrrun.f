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
* This routine write a record into the runfile.                        *
* Data type is Integer.                                                *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund university, Sweden                                     *
* Written: August 2003                                                 *
*                                                                      *
************************************************************************
      Subroutine ixWrRun(iRc,Label, Data,nData, iOpt)
#include "runinfo.fh"
#include "runtypes.fh"
*----------------------------------------------------------------------*
* Declare arguments                                                    *
*----------------------------------------------------------------------*
      Integer       iRc
      Character*(*) Label
      Integer       Data(*)
      Integer       nData
      Integer       iOpt
*----------------------------------------------------------------------*
* Local variables.                                                     *
*----------------------------------------------------------------------*
      Character*64 Errmsg
*----------------------------------------------------------------------*
* Check that arguments are ok.                                         *
*----------------------------------------------------------------------*
      If(iOpt.ne.0) Then
         Write(ErrMsg,*) 'Illegal option flag:',iOpt
         Call SysAbendMsg('ixWrRun',ErrMsg,' ')
      End If
      iRc=0
*----------------------------------------------------------------------*
* Call generic writing routine.                                        *
*----------------------------------------------------------------------*
      Call gxWrRun(iRc,Label, Data,nData, iOpt, TypInt)
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Return
      End
