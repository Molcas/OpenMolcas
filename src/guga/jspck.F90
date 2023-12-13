!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine JSPCK(INTSYM,L,ISYM)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: INTSYM(*)
integer(kind=iwp), intent(in) :: L, ISYM
integer(kind=iwp) :: IPOW, MY

!INTW = INTSYM((L+9)/10)
!IPOW = 2**(27-3*mod(L-1,10))
!INTW = INTW+(ISYM-1)*IPOW
!INTSYM((L+9)/10) = INTW
MY = (L+9)/10
IPOW = 2**(27-3*mod(L-1,10))
INTSYM(MY) = INTSYM(MY)+(ISYM-1)*IPOW

return

end subroutine JSPCK
