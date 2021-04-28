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

subroutine facab(binom,na1,nb1,crda,crdb,xab)

use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: binom(*), crda(*), crdb(*)
integer(kind=iwp), intent(in) :: na1, nb1
real(kind=wp), intent(out) :: xab(na1+nb1-1)
integer(kind=iwp) :: ia1, ib1, naind, nbind

xab(:) = Zero
naind = (na1*(na1-1))/2
nbind = (nb1*(nb1-1))/2
do ia1=1,na1
  do ib1=1,nb1
    xab((ia1-1)+ib1) = xab((ia1-1)+ib1)+(binom(naind+ia1)*crda((na1+1)-ia1))*binom(nbind+ib1)*crdb((nb1+1)-ib1)
  end do
end do

return

end subroutine facab
