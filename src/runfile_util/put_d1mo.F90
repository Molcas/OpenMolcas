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

subroutine Put_D1MO(D1MO,nD1MO)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
real*8 D1MO(nD1MO)
character*24 Label

Label = 'D1mo'
call Put_dArray(Label,D1MO,nD1MO)

return

end subroutine Put_D1MO
