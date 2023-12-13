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

use AMFI_global, only: Lmax
use Constants, only: Zero, One, speed => c_in_au
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: m1, m2, m3, m4, mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), ncontl1, ncontl2, &
                                 ncontl3, ncontl4, interxyz(*), isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
real(kind=wp), intent(in) :: angint(ncontl1,ncontl2,ncontl3,ncontl4,*), preY(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
!bs !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
real(kind=wp), intent(out) :: carteY(ncontl1,ncontl3,ncontl2,ncontl4)
logical(kind=iwp), intent(in) :: cleaner
integer(kind=iwp) :: iblock, irun, isgnM(-1:1,-1:1,-1:1,-1:1), ityp, Mabs1, Mabs2, Mabs3, Mabs4
real(kind=wp) :: factor, prey1234
real(kind=wp), parameter :: fine = One/speed, speed2 = speed*speed

!write(u6,*) 'begin tosigy '
!bs cleaning up the integral-array
carteY(:,:,:,:) = Zero
!bs set some signs
!bs isgnM will give an additonal minus-sign if both m-values
!bs (cartesian and angular) are negative  see $$$$
isgnM(:,:,:,:) = 1
if (m1 < 0) isgnM(-1,:,:,:) = -isgnM(-1,:,:,:)
if (m2 < 0) isgnM(:,-1,:,:) = -isgnM(:,-1,:,:)
if (m3 < 0) isgnM(:,:,-1,:) = -isgnM(:,:,-1,:)
if (m4 < 0) isgnM(:,:,:,-1) = -isgnM(:,:,:,-1)
!bs define absolute m-values
Mabs1 = abs(m1)
Mabs2 = abs(m2)
Mabs3 = abs(m3)
Mabs4 = abs(m4)
irun = 0
if (interxyz(1) == 0) then
  write(u6,*) 'tosigy: no interaction: ',m1,m2,m3,m4
  call Abend()
end if
prey1234 = preY(m1,m2,m3,m4)
!write(u6,*) 'prey ',prey1234
do while (interxyz(irun+1) > 0)
  irun = irun+1
  !write(u6,*) 'tosigy: ',irun,interxyz(irun)

  select case (interxyz(irun))

    case (1)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 1',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
      factor = isgnM(1,1,1,1)*prey1234*real(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (2)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 2',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(-1,-1,-1,-1)*prey1234*real(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (3)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 3',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = isgnM(1,1,1,-1)*prey1234*real(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (4)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 4',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = isgnM(-1,-1,-1,1)*prey1234*real(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (5)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 5',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = isgnM(1,1,-1,1)*prey1234*real(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (6)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 6',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = isgnM(-1,-1,1,-1)*prey1234*real(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (7)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 7',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = isgnM(1,-1,1,1)*prey1234*real(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (8)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 8',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(-1,1,-1,-1)*prey1234*real(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (9)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 9',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
      factor = isgnM(-1,1,1,1)*prey1234*real(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (10)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 10',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(1,-1,-1,-1)*prey1234*real(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (11)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 11',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = isgnM(1,1,-1,-1)*prey1234*real(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (12)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 12',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = isgnM(-1,-1,1,1)*prey1234*real(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (13)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 13',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = isgnM(1,-1,1,-1)*prey1234*real(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (14)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 14',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = isgnM(-1,1,-1,1)*prey1234*real(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (15)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 15',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = isgnM(1,-1,-1,1)*prey1234*real(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (16)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigy','wrong ityp in tosigY 16',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = isgnM(-1,1,1,-1)*prey1234*real(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4),kind=wp)
      if (ityp == 3) factor = -factor
      call daxpint(angint(:,:,:,:,iblock),carteY,factor,ncontl1,ncontl2,ncontl3,ncontl4)

  end select
end do
if (cleaner) then
  do irun=1,ncontl1
    carteY(irun,irun,:,:) = Zero
  end do
end if

return

end subroutine tosigY
