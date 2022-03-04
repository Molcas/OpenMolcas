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

subroutine tosigY(m1,m2,m3,m4,angint,mcombina,ncontl1,ncontl2,ncontl3,ncontl4,carteY,preY,interxyz,isgnprod,cleaner)
!bs this subroutine combines the angular integrals
!bs to the integrals for the real-valued linear
!bs combinations for the sigma_X part
!bs definition of the real-valued linear combinations:
!bs
!bs M=0  is the same as   Y(L,0)
!bs
!bs M > 0
!bs
!bs | L,M,+> = 2**(-0.5) ( (-1)**M Y(L,M) + Y(L,-M))
!bs
!bs | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M)) ($$$$)
!bs
!bs due to symmetry, there can be only integrals
!bs with one or three (sigma_+ and sigma_-)  - combinations

implicit real*8(a-h,o-z)
#include "para.fh"
parameter(fine=7.29735308D-03) !TO_BE_CHECKED
!bs at least it's identical with Odd's valuE
parameter(speed=1d0/fine)
parameter(speed2=speed*speed)
!bs !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
dimension mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), angint(ncontl1,ncontl2,ncontl3,ncontl4,*), &
          carteY(ncontl1,ncontl3,ncontl2,ncontl4), preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), interxyz(*), &
          isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), isgnM(-1:1,-1:1,-1:1,-1:1)
logical cleaner

!write(6,*) 'begin tosigy '
!bs cleaning up the integral-array
irun = ncontl4*ncontl2*ncontl3*ncontl1
call dzero(carteY,irun)
!bs set some signs
!bs isgnM will give an additonal minus-sign if both m-values
!bs (cartesian and angular) are negative  see $$$$
do irun4=-1,1
  do irun3=-1,1
    do irun2=-1,1
      do irun1=-1,1
        isgnM(irun1,irun2,irun3,irun4) = 1
      end do
    end do
  end do
end do
if (m1 < 0) then
  do irun4=-1,1
    do irun3=-1,1
      do irun2=-1,1
        isgnM(-1,irun2,irun3,irun4) = -isgnM(-1,irun2,irun3,irun4)
      end do
    end do
  end do
end if
if (m2 < 0) then
  do irun4=-1,1
    do irun3=-1,1
      do irun1=-1,1
        isgnM(irun1,-1,irun3,irun4) = -isgnM(irun1,-1,irun3,irun4)
      end do
    end do
  end do
end if
if (m3 < 0) then
  do irun4=-1,1
    do irun2=-1,1
      do irun1=-1,1
        isgnM(irun1,irun2,-1,irun4) = -isgnM(irun1,irun2,-1,irun4)
      end do
    end do
  end do
end if
if (m4 < 0) then
  do irun3=-1,1
    do irun2=-1,1
      do irun1=-1,1
        isgnM(irun1,irun2,irun3,-1) = -isgnM(irun1,irun2,irun3,-1)
      end do
    end do
  end do
end if
!bs define absolute m-values
Mabs1 = abs(m1)
Mabs2 = abs(m2)
Mabs3 = abs(m3)
Mabs4 = abs(m4)
irun = 0
if (interxyz(1) == 0) then
  write(6,*) 'tosigy: no interaction: ',m1,m2,m3,m4
  call Abend()
end if
prey1234 = preY(m1,m2,m3,m4)
!write(6,*) 'prey ',prey1234
do while (interxyz(irun+1) > 0)
  irun = irun+1
  !write(6,*) 'tosigy: ',irun,interxyz(irun)

  select case (interxyz(irun))

    case (1)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 1',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
      factor = isgnM(1,1,1,1)*prey1234*dble(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (2)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 2',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(-1,-1,-1,-1)*prey1234*dble(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (3)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 3',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = isgnM(1,1,1,-1)*prey1234*dble(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (4)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 4',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = isgnM(-1,-1,-1,1)*prey1234*dble(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (5)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 5',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = isgnM(1,1,-1,1)*prey1234*dble(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (6)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 6',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = isgnM(-1,-1,1,-1)*prey1234*dble(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (7)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 7',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = isgnM(1,-1,1,1)*prey1234*dble(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (8)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 8',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(-1,1,-1,-1)*prey1234*dble(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (9)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 9',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
      factor = isgnM(-1,1,1,1)*prey1234*dble(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (10)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 10',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(1,-1,-1,-1)*prey1234*dble(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (11)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 11',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(1,1,-1,-1)*prey1234*dble(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (12)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 12',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = isgnM(-1,-1,1,1)*prey1234*dble(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (13)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 13',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = isgnM(1,-1,1,-1)*prey1234*dble(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (14)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 14',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = isgnM(-1,1,-1,1)*prey1234*dble(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (15)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 15',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = isgnM(1,-1,-1,1)*prey1234*dble(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (16)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 16',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = isgnM(-1,1,1,-1)*prey1234*dble(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4))
      if (ityp == 3) factor = -factor
      call daxpint(angint(1,1,1,1,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

  end select
end do
if (cleaner) then
  do irun4=1,ncontl4
    do irun2=1,ncontl2
      do irun1=1,ncontl1
        cartey(irun1,irun1,irun2,irun4) = 0d0
      end do
    end do
  end do
end if

return

end subroutine tosigY
