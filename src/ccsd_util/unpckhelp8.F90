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

subroutine unpckhelp8(a,b,dimp,dimef,eadd,noe,bb,dimb)
! this routine does:
! b(ef,_Bb) = a(pe,qf)-a(qf,pe) for symp=symq

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dimp, dimef, eadd, noe, bb, dimb
real(kind=wp), intent(in) :: a(dimp,dimp)
real(kind=wp), intent(inout) :: b(dimef,dimb)
integer(kind=iwp) :: e, ef

ef = 0
do e=2,noe
  b(ef+1:ef+e-1,bb) = a(eadd+e,eadd+1:eadd+e-1)-a(eadd+1:eadd+e-1,eadd+e)
  ef = ef+e-1
end do

return

end subroutine unpckhelp8
