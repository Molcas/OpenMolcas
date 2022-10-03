!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine Get_dExcdRa(ipdExcdRa,ndExcdRa)
      Implicit Real*8 (A-H,O-Z)
#include "WrkSpc.fh"

      Character(LEN=24) Label
      Logical      Found

      Label='dExcdRa'
      Call qpg_dArray(Label,Found,ndExcdRa)
      If(.not.Found .or. ndExcdRa.eq.0) Then
         Call SysAbendmsg('Get_dExcdRa','Did not find:',Label)
      End If
      Call GetMem('dExcdRa','Allo','Real',ipdExcdRa,ndExcdRa)
      Call Get_dArray(Label,Work(ipdExcdRa),ndExcdRa)

      Return
      End
      Subroutine Get_dExcdRa_X(dExcdRa,ndExcdRa)
      Implicit None
      Character(LEN=24) Label
      Logical      Found
      Integer :: mdExcdRa=-1
      Integer, Intent(In) :: ndExcdRa
      Real*8,  Intent(Out) :: dExcdRa(ndExcdRa)

      Label='dExcdRa'
      Call qpg_dArray(Label,Found,mdExcdRa)
      If(.not.Found .or. mdExcdRa.eq.0) Then
         Call SysAbendmsg('Get_dExcdRa','Did not find:',Label)
      End If
      If (mdExcdRa/=ndExcdRa) Then
         Write (6,*) 'mdExcdRa/=ndExcdRa'
         Write (6,*)  mdExcdRa,'/=',ndExcdRa
         Call AbEnd()
      End If
      Call Get_dArray(Label,dExcdRa,ndExcdRa)

      Return
      End
