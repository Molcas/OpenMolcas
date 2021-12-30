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

subroutine ICPCK(ICSPCK,L,ICASE)

use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(inout) :: ICSPCK(*)
integer(kind=iwp), intent(in) :: L, ICASE
integer(kind=iwp) :: IPOW, MY

!INTW = ICSPCK((L+14)/15)
!IPOW = 2**(28-2*mod(L-1,15))
!INTW = INTW+ICASE*IPOW
!ICSPCK((L+14)/15) = INTW
MY = (L+14)/15
IPOW = 2**(28-2*mod(L-1,15))
ICSPCK(MY) = ICSPCK(MY)+ICASE*IPOW

return

end subroutine ICPCK
