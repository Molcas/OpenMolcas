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

subroutine MakeWwHlp1(Ww,W1,dima,dimb,dimab,dimbe,dimga,dimbega,key)
! this routine does:
! Make Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
! for the case a"=b" , be"=ga"
! N.B. algoritmus nieje prilis vymakany
!
! parameter description
! Ww   - array for Ww+(-) (O)
! W1   - array for W1(a",be",b",ga") (I)
! dimx - dimension of a",b",ab",be",ga",bega" (I)
! key  - 1 - calc Ww+, 2 - calc Ww- (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: dima, dimb, dimab, dimbe, dimga, dimbega, key
real(kind=wp), intent(out) :: Ww(dimab,dimbega)
real(kind=wp), intent(in) :: W1(dima,dimbe,dimb,dimga)
integer(kind=iwp) :: a, ab, be, bega

if (key == 1) then
  bega = 0
  do be=1,dimbe
    ab = 0
    do a=2,dima
      Ww(ab+1:ab+a-1,bega+1:bega+be) = W1(a,be,1:a-1,1:be)+W1(1:a-1,be,a,1:be)
      ab = ab+a-1
    end do
    bega = bega+be
  end do
else
  bega = 0
  do be=2,dimbe
    ab = 0
    do a=2,dima
      Ww(ab+1:ab+a-1,bega+1:bega+be-1) = W1(a,be,1:a-1,1:be-1)-W1(1:a-1,be,a,1:be-1)
      ab = ab+a-1
    end do
    bega = bega+be-1
  end do
end if

! Cely clen ma Faktor 2, tu teda nevydelim 2
!Ww(:,:) = Half*Ww(:,:)

return

end subroutine MakeWwHlp1
