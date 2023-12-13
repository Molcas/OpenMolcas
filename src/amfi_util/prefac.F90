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

subroutine prefac(Lmax,preroots,clebsch)

use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: Lmax
real(kind=wp), intent(out) :: preroots(2,0:Lmax), clebsch(3,2,-Lmax:Lmax,0:Lmax)
integer(kind=iwp) :: L, L2, M, M2
real(kind=wp) :: fact
real(kind=wp), external :: getCG

!bs the roots appearing in front of all
!bs the contributions
!write(u6,*) 'begin of prefac'
do L=0,Lmax
  fact = One/sqrt(real(L+L+1,kind=wp))
  preroots(1,L) = sqrt(real(L,kind=wp))*fact
  preroots(2,L) = sqrt(real(L+1,kind=wp))*fact
end do
!bs there are Clebsch-Gordan-Coefficients
!bs which always appear:
!bs
!bs   -----                       ------
!bs  |                                 |
!bs  |  l +/- 1     1        |      l  |
!bs  |                       |         |
!bs  |                       |         |
!bs  |  m+/-1,0   -1,1,0     |      m  |
!bs  |                       |         |
!bs  |                                 |
!bs   -----                       -----
!bs
!bs
!bs array clebsch (3,2,-Lmax:Lmax,0:Lmax)
!bs first index    1:  m-1
!bs                2:  m
!bs                3:  m+1
!bs second index   1:  l-1
!bs                2:  l+1
!bs third index        m
!bs fourth index       l

!write(u6,*),'start to generate CGs'
do L=0,Lmax
  L2 = 2*L
  do M=-L,L
    !write(u6,*) 'L,M: ',L,M
    M2 = 2*M
    !bs getCG calculates CG-coeffecients. In order to avoid fractions,
    !bs e.g. for spins, arguments are doubled values...
    clebsch(1,1,M,L) = getCG(L2-2,2,L2,M2-2,2,M2)
    clebsch(2,1,M,L) = getCG(L2-2,2,L2,M2,0,M2)
    clebsch(3,1,M,L) = getCG(L2-2,2,L2,M2+2,-2,M2)
    clebsch(1,2,M,L) = getCG(L2+2,2,L2,M2-2,2,M2)
    clebsch(2,2,M,L) = getCG(L2+2,2,L2,M2,0,M2)
    clebsch(3,2,M,L) = getCG(L2+2,2,L2,M2+2,-2,M2)
  end do
end do

return

end subroutine prefac
