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

subroutine MakeWwHlp3(Ww,W1,W2,dima,dimb,dimbe,dimga,dimbega,key)
! this routine does:
! Make Ww(+,-)((ab)",(bega)") from W1(a",be",b",ga")
!                              and W2(b",be",a",ga")
! for the case a"/=b" , be"=ga"
! N.B. algoritmus nieje prilis vymakany, prva cast sa da
!      urobit pomodou dcopy
!
! parameter description
! Ww   - array for Ww+(-) (O)
! W1   - array for W1(a",be",b",ga") (I)
! W2   - array for W1(b",be",a",ga") (I)
! dimx - dimension of a",b",be",ga",be>(>=)ga" (I)
! key  - 1 - calc Ww+, 2 - calc Ww- (I)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: dima, dimb, dimbe, dimga, dimbega, key
real(kind=wp) :: Ww(dima,dimb,dimbega), W1(dima,dimbe,dimb,dimga), W2(dimb,dimbe,dima,dimga)
integer(kind=iwp) :: a, b, be, bega, ga

if (key == 1) then
  bega = 0
  do be=1,dimbe
    do ga=1,be
      bega = bega+1
      do b=1,dimb
        do a=1,dima
          Ww(a,b,bega) = W1(a,be,b,ga)+W2(b,be,a,ga)
        end do
      end do
    end do
  end do
else
  bega = 0
  do be=2,dimbe
    do ga=1,be-1
      bega = bega+1
      do b=1,dimb
        do a=1,dima
          Ww(a,b,bega) = W1(a,be,b,ga)-W2(b,be,a,ga)
        end do
      end do
    end do
  end do
end if

! Cely clen ma Faktor 2, tu teda nevydelim 2
!call mv0sv(dima*dimb*dimbega,dima*dimb*dimbega,Ww(1,1,1),Half)

return

end subroutine MakeWwHlp3
