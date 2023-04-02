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

subroutine saamphlp2(t2aaaa,t2bbbb,t2abab,t2baba,t2abba,t2baab,noai,noaj,nobi,nobj,nvaa,nvab,nvba,nvbb,key)
! adaptation routine for T2 amplitudes for symi>symj, syma>symb
!
! t2aaaa - amplitudes T2(a,b,i,j)aaaa (I/O)
! t2bbbb - amplitudes T2(a,b,i,j)bbbb (I/O)
! t2abab - amplitudes T2(a,b,i,j)abab (I/O)
! t2baba - amplitudes T2(b,a,j,i)abab (I/O)
! t2abba - amplitudes T2(b,a,i,j)abab (I/O)
! t2baab - amplitudes T2(a,b,j,i)abab (I/O)
! noai   - number of alpha occupied orbitals in symi (I)
! noaj   - number of alpha occupied orbitals in symj (I)
! nobi   - number of beta occupied orbitals in symi (I)
! nobj   - number of beta occupied orbitals in symj (I)
! nvaa   - number of alpha virtual orbitals in syma (I)
! nvab   - number of alpha virtual orbitals in symb (I)
! nvba   - number of beta virtual orbitals in syma (I)
! nvbb   - number of beta virtual orbitals in symb (I)
! key    - type of adaptation (I)
!          0 - no adaptation
!          1 - T2 DDVV adaptation
!          2 - T2 DDVV + T1 DV adaptation
!          3 - full T1 and T2 adaptation (only for doublets)
!          4 - full T2 without SDVS (only for doublets)

integer noai, noaj, nobi, nobj, nvaa, nvab, nvba, nvbb
!T2
real*8 t2aaaa(1:nvaa,1:nvab,1:noai,1:noaj)
real*8 t2bbbb(1:nvba,1:nvbb,1:nobi,1:nobj)
real*8 t2abab(1:nvaa,1:nvbb,1:noai,1:nobj)
real*8 t2baba(1:nvab,1:nvba,1:noaj,1:nobi)
real*8 t2abba(1:nvab,1:nvba,1:noai,1:nobj)
real*8 t2baab(1:nvaa,1:nvbb,1:noaj,1:nobi)
integer key
! help variables
integer ndi, ndj, nva, nvb, nsi, nsj, nsa, nsb
integer i, j, a, b
real*8 taaaa, tbbbb, tabab, tabba, tbaab, tbaba, t1, t2

if (key == 0) return

ndi = nobi
ndj = nobj
nva = nvaa
nvb = nvab
nsi = noai-nobi
nsj = noaj-nobj
nsa = nvba-nvaa
nsb = nvbb-nvab

! T2 adatption
!
! case I) DDVV i>=j,a>=b
! turn od in any case of type of adaption
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

do j=1,ndj
  do i=1,ndi
    do b=1,nvb
      do a=1,nva

        taaaa = t2aaaa(a,b,i,j)
        tbbbb = t2bbbb(a+nsa,b+nsb,i,j)
        tabab = t2abab(a,b+nsb,i,j)
        tbaba = t2baba(b,a+nsa,j,i)
        tabba = -t2abba(b,a+nsa,i,j)
        tbaab = -t2baab(a,b+nsb,j,i)

        t1 = (tabab+tbaba-tabba-tbaab)/4.0d0
        t2 = (2.0d0*(taaaa+tbbbb)+tabab+tbaba+tabba+tbaab)/1.2d1

        taaaa = 2.0d0*t2
        tbbbb = 2.0d0*t2
        tabab = t1+t2
        tbaba = t1+t2
        tabba = t2-t1
        tbaab = t2-t1

        t2aaaa(a,b,i,j) = taaaa
        t2bbbb(a+nsa,b+nsb,i,j) = tbbbb
        t2abab(a,b+nsb,i,j) = tabab
        t2baba(b,a+nsa,j,i) = tbaba
        t2abba(b,a+nsa,i,j) = -tabba
        t2baab(a,b+nsb,j,i) = -tbaab

      end do
    end do
  end do
end do

if (((key == 3) .or. (key == 4)) .and. (nsb > 0)) then

  ! case II) DDVS
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

  b = nsb
  do j=1,ndj
    do i=1,ndi
      do a=1,nva

        tbbbb = t2bbbb(a+nsa,b,i,j)
        tabab = t2abab(a,b,i,j)
        tbaab = -t2baab(a,b,j,i)

        t1 = (tabab-tbaab)/2.0d0
        t2 = (2.0d0*tbbbb+tabab+tbaab)/6.0d0

        tbbbb = 2.0d0*t2
        tabab = t2+t1
        tbaab = t2-t1

        t2bbbb(a+nsa,b,i,j) = tbbbb
        t2abab(a,b,i,j) = tabab
        t2baab(a,b,j,i) = -tbaab

      end do
    end do
  end do

end if

if (((key == 3) .or. (key == 4)) .and. (nsa > 0)) then

  ! case III) DDSV
  !
  ! bbbb = - t24b(ij,ab)
  ! abab = - -t22b(i,j,b,a)
  ! baab = - t22b(j,i,b,a)
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

  a = nsa
  do j=1,ndj
    do i=1,ndi
      do b=1,nvb

        tbbbb = -t2bbbb(a,b+nsb,i,j)
        tabab = t2abba(b,a,i,j)
        tbaab = -t2baba(b,a,j,i)

        t1 = (tabab-tbaab)/2.0d0
        t2 = (2.0d0*tbbbb+tabab+tbaab)/6.0d0

        tbbbb = 2.0d0*t2
        tabab = t2+t1
        tbaab = t2-t1

        t2bbbb(a,b+nsb,i,j) = -tbbbb
        t2abba(b,a,i,j) = tabab
        t2baba(b,a,j,i) = -tbaab

      end do
    end do
  end do

end if

if (((key == 3) .or. (key == 4)) .and. (nsi > 0)) then

  ! case IV) SDVV
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

  i = ndi+nsi
  do j=1,ndj
    do b=1,nvb
      do a=1,nva

        taaaa = t2aaaa(a,b,i,j)
        tabab = t2abab(a,b+nsb,i,j)
        tabba = -t2abba(b,a+nsa,i,j)

        t1 = (tabab-tabba)/2.0d0
        t2 = (2.0d0*taaaa+tabab+tabba)/6.0d0

        taaaa = 2.0d0*t2
        tabab = t2+t1
        tabba = t2-t1

        t2aaaa(a,b,i,j) = taaaa
        t2abab(a,b+nsb,i,j) = tabab
        t2abba(b,a+nsa,i,j) = -tabba

      end do
    end do
  end do

end if

if (((key == 3) .or. (key == 4)) .and. (nsj > 0)) then

  ! case V) DSVV
  !
  ! aaaa= - t24a(ij,ab)
  ! abab= - -t22b(j,i,a,b)
  ! abba= - t22b(j,i,b,a)
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

  j = ndj+nsj
  do i=1,ndi
    do b=1,nvb
      do a=1,nva

        taaaa = -t2aaaa(a,b,i,j)
        tabab = t2baab(a,b+nsb,j,i)
        tabba = -t2baba(b,a+nsa,j,i)

        t1 = (tabab-tabba)/2.0d0
        t2 = (2.0d0*taaaa+tabab+tabba)/6.0d0

        taaaa = 2.0d0*t2
        tabab = t2+t1
        tabba = t2-t1

        t2aaaa(a,b,i,j) = -taaaa
        t2baab(a,b+nsb,j,i) = tabab
        t2baba(b,a+nsa,j,i) = -tabba

      end do
    end do
  end do

end if

return

end subroutine saamphlp2
