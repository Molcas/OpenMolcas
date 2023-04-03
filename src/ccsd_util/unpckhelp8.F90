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
integer(kind=iwp) :: dimp, dimef, eadd, noe, bb, dimb
real(kind=wp) :: a(dimp,dimp), b(dimef,dimb)
integer(kind=iwp) :: ef, pe, qf

ef = 0
do pe=eadd+2,eadd+noe
  do qf=eadd+1,pe-1
    ef = ef+1
    b(ef,bb) = a(pe,qf)-a(qf,pe)
  end do
end do

return

end subroutine unpckhelp8
