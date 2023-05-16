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

subroutine MakeWwdHlp1(Ww,W1,dima,dimbe,dimbega)
! this routine does:
! Make Ww(+)((aa)",(bega)") from W1(a",be",b",ga")
! for the case beSGrp=gaSGrp
! N.B. algoritmus nieje prilis vymakany
!
! parameter description
! Ww   - array for Ww+(-) (O)
! W1   - array for W1(a",be",b",ga") (I)
! dimx - dimension of a",ga",bega" (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, dimbega
real(kind=wp), intent(out) :: Ww(dima,dimbega)
real(kind=wp), intent(in) :: W1(dima,dimbe,dima,dimbe)
integer(kind=iwp) :: a, be, bega, ga

bega = 0
do be=1,dimbe
  do ga=1,be
    bega = bega+1
    do a=1,dima
      Ww(a,bega) = W1(a,be,a,ga)
    end do
  end do
end do

return

end subroutine MakeWwdHlp1
