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

subroutine contone(L,oneoverR3,onecontr,Lmax,contcoeff,nprim,ncont,MxcontL,dummy,onecartX,onecartY,onecartZ,charge,oneonly)
!bs contracts one-electron integrals and multiplies with l,m-dependent
!bs factors for L-,L0,L+

use index_functions, only: iTri
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: L, Lmax, nprim, ncont, MxcontL
real(kind=wp), intent(in) :: oneoverR3(*), contcoeff(nprim,ncont), charge
real(kind=wp), intent(out) :: onecontr(MxcontL,MxcontL,-Lmax:Lmax,3), dummy(ncont,ncont)
real(kind=wp), intent(inout) :: onecartX(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1)), &
                                onecartY(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1)), onecartZ(MxcontL,MxcontL,(Lmax+Lmax+1)*(Lmax+1))
logical(kind=iwp), intent(in) :: oneonly
integer(kind=iwp) :: icont1, icont2, iprim1, iprim2, M
real(kind=wp) :: factor0, factormin, factorplus

!bs first of all cleaning dummy and onecontr
dummy(:,:) = Zero
if (oneonly) then
  onecartx(:,:,:) = Zero
  onecarty(:,:,:) = Zero
  onecartz(:,:,:) = Zero
end if
onecontr(:,:,:,:) = Zero
!bs contract onto dummy
do icont1=1,ncont
  do icont2=1,ncont
    do iprim2=1,nprim
      do iprim1=1,nprim
        dummy(icont2,icont1) = dummy(icont2,icont1)+contcoeff(iprim1,icont1)*contcoeff(iprim2,icont2)*oneoverR3(iTri(iprim1,iprim2))
      end do
    end do
  end do
end do
!bs start to add l,m dependent factors
do M=-L,L
  factormin = charge*sqrt(real(L*L-M*M+L+M,kind=wp))
  factor0 = charge*real(M,kind=wp)
  factorplus = charge*sqrt(real(L*L-M*M+L-M,kind=wp))
  onecontr(1:ncont,1:ncont,M,1) = dummy*factormin  ! L-minus
  onecontr(1:ncont,1:ncont,M,2) = dummy*factor0    ! L-0
  onecontr(1:ncont,1:ncont,M,3) = dummy*factorplus ! L-plus
end do
!bs make the final cartesian integrals
call cartoneX(L,Lmax,onecontr,ncont,MxcontL,onecartX)
call cartoneY(L,Lmax,onecontr,ncont,MxcontL,onecartY)
call cartoneZ(L,Lmax,onecontr,ncont,MxcontL,onecartZ)

return

end subroutine contone
