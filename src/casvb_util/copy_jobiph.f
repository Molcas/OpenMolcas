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
* Copyright (C) 1996-2006, T. Thorsteinsson and D. L. Cooper           *
************************************************************************
      SUBROUTINE Copy_JobIph(InFile,OutFile)
c This is very nasty routine, actually a temporary hack!!!
c should be replaced to a proper copy of JobIph files!
      Implicit Real*8 (a-h,o-z)
#include "WrkSpc.fh"
      Character*(*) InFile,OutFile

      call fcopy(InFile,OutFile,ierr)
        if(iErr.ne.0) call abend
c      Junk='cp -p '//File1(1:lFile1)//' '//File2(1:lFile2)
      Return
      End
