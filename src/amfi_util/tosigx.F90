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

subroutine tosigX(m1,m2,m3,m4,angint,mcombina,ncontl1,ncontl2,ncontl3,ncontl4,carteX,preXZ,interxyz,isgnprod,cleaner)
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
!bs | L,M,-> = -i 2**(-0.5) ( (-1)**M Y(L,M) - Y(L,-M))  ($$$$)
!bs
!bs due to symmetry, there can be only integrals
!bs with indices one or three  (sigma_+ and sigma_-)- combinations

use AMFI_global, only: Lmax
use Constants, only: Zero, One, speed => c_in_au
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: m1, m2, m3, m4, mcombina(2,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax), ncontl1, ncontl2, &
                                 ncontl3, ncontl4, interxyz(*), isgnprod(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
real(kind=wp), intent(in) :: angint(ncontl1,ncontl2,ncontl3,ncontl4,*), preXZ(-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax,-Lmax:Lmax)
!bs !!!!!!!!!!!changing now to the order of chemists notation!!!!!!!!!!
real(kind=wp), intent(out) :: carteX(ncontl1,ncontl3,ncontl2,ncontl4)
logical(kind=iwp), intent(in) :: cleaner
integer(kind=iwp) :: iblock, irun, isgnM(-1:1,-1:1,-1:1,-1:1), ityp, Mabs1, Mabs2, Mabs3, Mabs4
real(kind=wp) :: factor, prexz1234
real(kind=wp), parameter :: fine = One/speed, speed2 = speed**2

!write(u6,*) ' begin tosigx'
!bs cleaning up the integral-array
carteX(:,:,:,:) = Zero
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
  write(u6,*) 'tosigx: no interaction: ',m1,m2,m3,m4
  call Abend()
end if
prexz1234 = preXZ(m1,m2,m3,m4)
do while (interxyz(irun+1) > 0)
  irun = irun+1
  !write(u6,*) 'tosigx: ',irun,interxyz(irun)

  select case (interxyz(irun))

    case (1)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 1',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,Mabs4)
      factor = real(isgnM(1,1,1,1),kind=wp)*prexz1234*real(isgnprod(Mabs1,Mabs2,Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (2)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 2',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = real(isgnM(-1,-1,-1,-1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,-Mabs2,-Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (3)
      ityp = mcombina(1,Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 3',' ')
      iblock = mcombina(2,Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = real(isgnM(1,1,1,-1),kind=wp)*prexz1234*real(isgnprod(Mabs1,Mabs2,Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (4)
      ityp = mcombina(1,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 4',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = real(isgnM(-1,-1,-1,1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,-Mabs2,-Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (5)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 5',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = real(isgnM(1,1,-1,1),kind=wp)*prexz1234*real(isgnprod(Mabs1,Mabs2,-Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (6)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 6',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = real(isgnM(-1,-1,1,-1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,-Mabs2,Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (7)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 7',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = real(isgnM(1,-1,1,1),kind=wp)*prexz1234*real(isgnprod(Mabs1,-Mabs2,Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (8)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 8',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = real(isgnM(-1,1,-1,-1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,Mabs2,-Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (9)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 9',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,Mabs4)
      factor = real(isgnM(-1,1,1,1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,Mabs2,Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (10)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 10',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,-Mabs4)
      factor = real(isgnM(1,-1,-1,-1),kind=wp)*prexz1234*real(isgnprod(Mabs1,-Mabs2,-Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (11)
      ityp = mcombina(1,Mabs1,Mabs2,-Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 11',' ')
      iblock = mcombina(2,Mabs1,Mabs2,-Mabs3,-Mabs4)
      factor = real(isgnM(1,1,-1,-1),kind=wp)*prexz1234*real(isgnprod(Mabs1,Mabs2,-Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (12)
      ityp = mcombina(1,-Mabs1,-Mabs2,Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 12',' ')
      iblock = mcombina(2,-Mabs1,-Mabs2,Mabs3,Mabs4)
      factor = real(isgnM(-1,-1,1,1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,-Mabs2,Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (13)
      ityp = mcombina(1,Mabs1,-Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 13',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,Mabs3,-Mabs4)
      factor = real(isgnM(1,-1,1,-1),kind=wp)*prexz1234*real(isgnprod(Mabs1,-Mabs2,Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (14)
      ityp = mcombina(1,-Mabs1,Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 14',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,-Mabs3,Mabs4)
      factor = real(isgnM(-1,1,-1,1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,Mabs2,-Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (15)
      ityp = mcombina(1,Mabs1,-Mabs2,-Mabs3,Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 15',' ')
      iblock = mcombina(2,Mabs1,-Mabs2,-Mabs3,Mabs4)
      factor = real(isgnM(1,-1,-1,1),kind=wp)*prexz1234*real(isgnprod(Mabs1,-Mabs2,-Mabs3,Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

    case (16)
      ityp = mcombina(1,-Mabs1,Mabs2,Mabs3,-Mabs4)
      if ((ityp /= 1) .and. (ityp /= 3)) call SysAbendMsg('tosigx','wrong ityp in tosigX 16',' ')
      iblock = mcombina(2,-Mabs1,Mabs2,Mabs3,-Mabs4)
      factor = real(isgnM(-1,1,1,-1),kind=wp)*prexz1234*real(isgnprod(-Mabs1,Mabs2,Mabs3,-Mabs4),kind=wp)
      call daxpint(angint(:,:,:,:,iblock),carteX,factor,ncontl1,ncontl2,ncontl3,ncontl4)

  end select
end do
if (cleaner) then
  do irun=1,ncontl1
    carteX(irun,irun,:,:) = Zero
  end do
end if

return

end subroutine tosigX
