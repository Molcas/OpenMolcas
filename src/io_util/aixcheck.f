************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine AixCheck()
      Implicit Integer (a-z)
#include "ctl.fh"
      Logical Opened
      character name*256
*----------------------------------------------------------------------*
* Check if slot in table is available, it should NOT!                  *
*----------------------------------------------------------------------*
      Do n = 1, MxFile
         If (CtlBlk(pStat,n).ne.0) Then
#ifndef _DEVEL_
            Call SysAbendFileMsg('AixCheck',FCtlBlk(n),'Active unit.',
     &                           'Should have been closed!')
#else
            Call SysWarnMsg('AixCheck','Active unit: '//FCtlBlk(n),
     &                           ', should have been closed!')
#endif
         End If
         Inquire(unit=n,Opened=Opened)
         If (Opened ) Then
         if(.not. (n.eq.6.or.n.eq.5))then
          inquire(unit=n,name=name)
            Write (6,*) 'Fortran file:', n, '(',name(1:index(name,' ')),
     *        ')  is still open!'
#ifndef _DEVEL_
            Call Abend()
#endif
         End If
         endif
      End Do
      Return
      End
