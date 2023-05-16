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

subroutine saamphlp1(t24a,t24b,t22b,noa,nob,nva,nvb,key)
! adaptation routine for T2 amplitudes for symi=symj, syma=symb
!
! t24a - amplitudes T2(ab,ij)aaaa (I/O)
! t24b - amplitudes T2(ab,ij)bbbb (I/O)
! t22b - amplitudes T2(a,b,i,j)abab (I/O)
! noa  - number of alpha occupied orbitals in symi (I)
! nob  - number of beta occupied orbitals in symi (I)
! nva  - number of alpha virtual orbitals in syma (I)
! nvb  - number of beta virtual orbitals in syma (I)
! key  - type of adaptation (I)
!        0 - no adaptation
!        1 - T2 DDVV adaptation
!        2 - T2 DDVV + T1 DV adaptation
!        3 - full T1 and T2 adaptation (only for doublets)
!        4 - full T2 without SDVS (only for doublets)

use Index_Functions, only: nTri_Elem
use Constants, only: Zero, Two, Six, Half, Quart
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: noa, nob, nva, nvb, key
real(kind=wp), intent(inout) :: t24a(nTri_Elem(nva),nTri_Elem(noa)), t24b(nTri_Elem(nvb),nTri_Elem(nob)), t22b(nva,nvb,noa,nob)
integer(kind=iwp) :: a, ab, ab1, b, i, ij, j, nd, nsa, nsi, nv
real(kind=wp) :: t1, t2, taaaa, tabab, tabba, tbaab, tbaba, tbbbb

if (key == 0) return

nd = nob
nv = nva
nsi = noa-nob
nsa = nvb-nva

! T2 adatption
!
! case I) DDVV i>=j,a>=b
! turn od in any case of type of adaptation
!
! aaaa=t24a(ijab) DDVV
! bbbb=t24b(ijab) DDVV
! abab=t22b(ijab) DDVV
! abba=-t22b(ijba) DDVV
! baab=-t22b(jiab) DDVV
! baba=t22b(ji,ba) DDVV
!
! direct
! t1 = (abab+baba-abba-baab)/4
! t2 = (2aaaa+2bbbb+abab+baba+abba+baab)/12
!
! reverse
!
! aaaa = 2t2
! bbbb = 2t2
! abab = t2+t1
! abba = t2-t1
! baab = t2-t1
! baba = t2+t1

do a=2,nv
  do b=1,a-1
    ab = (a-1)*(a-2)/2+b
    ab1 = (a+nsa-1)*(a+nsa-2)/2+b+nsa
    do i=2,nd
      do j=1,i-1
        ij = (i-1)*(i-2)/2+j

        taaaa = t24a(ab,ij)
        tbbbb = t24b(ab1,ij)
        tabab = t22b(a,b+nsa,i,j)
        tabba = -t22b(b,a+nsa,i,j)
        tbaab = -t22b(a,b+nsa,j,i)
        tbaba = t22b(b,a+nsa,j,i)

        t1 = Quart*(tabab+tbaba-tabba-tbaab)
        t2 = (Two*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/12.0_wp

        t24a(ab,ij) = Two*t2
        t24b(ab1,ij) = Two*t2
        t22b(a,b+nsa,i,j) = t2+t1
        t22b(b,a+nsa,i,j) = t1-t2
        t22b(a,b+nsa,j,i) = t1-t2
        t22b(b,a+nsa,j,i) = t2+t1

      end do
    end do
  end do
end do

do a=2,nv
  do b=1,a-1
    ab = (a-1)*(a-2)/2+b
    ab1 = (a+nsa-1)*(a+nsa-2)/2+b+nsa
    do i=1,nd
      j = i
      ij = (i-1)*(i-2)/2+j

      !taaaa = t24a(ab,ij)
      !tbbbb = t24b(ab1,ij)
      tabab = t22b(a,b+nsa,i,j)
      tabba = -t22b(b,a+nsa,i,j)
      tbaab = -t22b(a,b+nsa,j,i)
      tbaba = t22b(b,a+nsa,j,i)

      t1 = Quart*(tabab+tbaba-tabba-tbaab)
      !t2 = (Two*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/12.0_wp
      t2 = Zero

      !t24a(ab,ij) = Two*t2
      !t24b(ab1,ij) = Two*t2
      t22b(a,b+nsa,i,j) = t2+t1
      t22b(b,a+nsa,i,j) = t1-t2
      t22b(a,b+nsa,j,i) = t1-t2
      t22b(b,a+nsa,j,i) = t2+t1

    end do
  end do
end do

do a=1,nv
  b = a
  ab = (a-1)*(a-2)/2+b
  ab1 = (a+nsa-1)*(a+nsa-2)/2+b+nsa
  do i=2,nd
    do j=1,i-1
      ij = (i-1)*(i-2)/2+j

      !taaaa = t24a(ab,ij)
      !tbbbb = t24b(ab1,ij)
      tabab = t22b(a,b+nsa,i,j)
      tabba = -t22b(b,a+nsa,i,j)
      tbaab = -t22b(a,b+nsa,j,i)
      tbaba = t22b(b,a+nsa,j,i)

      t1 = Quart*(tabab+tbaba-tabba-tbaab)
      !t2 = (Two*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/12.0_wp
      t2 = Zero

      !t24a(ab,ij) = Two*t2
      !t24b(ab1,ij) = Two*t2
      t22b(a,b+nsa,i,j) = t2+t1
      t22b(b,a+nsa,i,j) = t1-t2
      t22b(a,b+nsa,j,i) = t1-t2
      t22b(b,a+nsa,j,i) = t2+t1

    end do
  end do
end do

do a=1,nv
  b = a
  ab = (a-1)*(a-2)/2+b
  ab1 = (a+nsa-1)*(a+nsa-2)/2+b+nsa
  do i=1,nd
    j = i
    ij = (i-1)*(i-2)/2+j

    !taaaa = t24a(ab,ij)
    !tbbbb = t24b(ab1,ij)
    tabab = t22b(a,b+nsa,i,j)
    tabba = -t22b(b,a+nsa,i,j)
    tbaab = -t22b(a,b+nsa,j,i)
    tbaba = t22b(b,a+nsa,j,i)

    t1 = Quart*(tabab+tbaba-tabba-tbaab)
    !t2 = (Two*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/12.0_wp
    t2 = Zero

    !t24a(ab,ij) = Two*t2
    !t24b(ab1,ij) = Two*t2
    t22b(a,b+nsa,i,j) = t2+t1
    t22b(b,a+nsa,i,j) = t1-t2
    t22b(a,b+nsa,j,i) = t1-t2
    t22b(b,a+nsa,j,i) = t2+t1

  end do
end do

if (((key == 3) .or. (key == 4)) .and. (nsa > 0)) then

  ! case II) DDVS i>=j,a,b
  !
  ! bbbb = t24b(ij,ab)
  ! abab = t22b(i,j,a,b)
  ! baab =-t22b(j,i,a,b)
  !
  ! direct :
  !
  ! t1=(abab-baab)/2
  ! t2=(2bbbb+abab+baab)/6
  !
  ! reverse:
  !
  ! bbbb = 2t2
  ! abab = t2+t1
  ! baab = t2-t1

  b = nsa
  do a=1,nv
    ab = a*(a-1)/2+b
    do i=2,nd
      do j=1,i-1
        ij = (i-1)*(i-2)/2+j

        tbbbb = t24b(ab,ij)
        tabab = t22b(a,b,i,j)
        tbaab = -t22b(a,b,j,i)

        t1 = Half*(tabab-tbaab)
        t2 = (Two*tbbbb+tabab+tbaab)/Six

        t24b(ab,ij) = Two*t2
        t22b(a,b,i,j) = t2+t1
        t22b(a,b,j,i) = t1-t2

      end do
    end do
  end do

  b = nsa
  do a=1,nv
    ab = a*(a-1)/2+b
    do i=1,nd
      j = i
      ij = (i-1)*(i-2)/2+j

      !tbbbb = t24b(ab,ij)
      tabab = t22b(a,b,i,j)
      tbaab = -t22b(a,b,j,i)

      t1 = Half*(tabab-tbaab)
      !t2 = (Two*tbbbb+tabab+tbaab)/Six
      t2 = Zero

      !t24b(ab,ij) = Two*t2
      t22b(a,b,i,j) = t2+t1
      t22b(a,b,j,i) = t1-t2

    end do
  end do

end if

if (((key == 3) .or. (key == 4)) .and. (nsi > 0)) then

  ! case III) SDVV i,j,a>=b
  !
  ! aaaa=t24a(ij,ab)
  ! abab=t22b(i,j,a,b)
  ! abba=-t22b(i,j,b,a)
  !
  ! direct
  !
  ! t1=(abab-abba)/2
  ! t2=(2aaaa+abab+abba)/6
  !
  ! reverse
  !
  ! aaaa=2t2
  ! abab=t2+t1
  ! abba=t2-t1

  i = nd+nsi
  do a=2,nv
    do b=1,a-1
      ab = (a-1)*(a-2)/2+b
      do j=1,nd
        ij = (i-1)*(i-2)/2+j

        taaaa = t24a(ab,ij)
        tabab = t22b(a,b+nsa,i,j)
        tabba = -t22b(b,a+nsa,i,j)

        t1 = Half*(tabab-tabba)
        t2 = (Two*taaaa+tabab+tabba)/Six

        t24a(ab,ij) = Two*t2
        t22b(a,b+nsa,i,j) = t2+t1
        t22b(b,a+nsa,i,j) = t1-t2

      end do
    end do
  end do

  i = nd+nsi
  do a=1,nv
    b = a
    ab = (a-1)*(a-2)/2+b
    do j=1,nd
      ij = (i-1)*(i-2)/2+j

      !taaaa = t24a(ab,ij)
      tabab = t22b(a,b+nsa,i,j)
      tabba = -t22b(b,a+nsa,i,j)

      t1 = Half*(tabab-tabba)
      !t2 = (Two*taaaa+tabab+tabba)/Six
      t2 = Zero

      !t24a(ab,ij) = Two*t2
      t22b(a,b+nsa,i,j) = t2+t1
      t22b(b,a+nsa,i,j) = t1-t2

    end do
  end do

end if

return

end subroutine saamphlp1
