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

subroutine WrCheck(Label,Arr,n)

<<<<<<< HEAD
use Constants, only: One
=======
>>>>>>> upstream-openmolcas/master
use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(in) :: Label
integer(kind=iwp), intent(in) :: n
real(kind=wp), intent(in) :: Arr(n)
<<<<<<< HEAD
real(kind=wp), external :: DDot_

write(u6,*) Label,DDot_(n,Arr,1,Arr,1),DDot_(n,Arr,1,[One],0),n
=======

write(u6,*) Label,sum(Arr(:)**2),sum(Arr(:))
>>>>>>> upstream-openmolcas/master

end subroutine WrCheck
