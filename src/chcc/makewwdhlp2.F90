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

subroutine MakeWwdHlp2(Ww,W1,dima,dimbe,dimga)
! this routine does:
! Make Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
! for the case a"=b" , beSGrp/=gaSGrp
! N.B. algoritmus nieje prilis vymakany
!
! parameter description
! Ww   - array for Ww+(-) (O)
! W1   - array for W1(a",be",b",ga") (I)
! dimx - dimension of a",be",ga" (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimbe, dimga
real(kind=wp), intent(out) :: Ww(dima,dimbe,dimga)
real(kind=wp), intent(in) :: W1(dima,dimbe,dima,dimga)
integer(kind=iwp) :: a

do a=1,dima
  !DEC$ VECTOR UNALIGNED
  Ww(a,:,:) = W1(a,:,a,:)
end do

return

end subroutine MakeWwdHlp2
