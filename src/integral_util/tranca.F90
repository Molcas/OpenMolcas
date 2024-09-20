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

subroutine Tranca(Cavxyz,Cavsph,lMax,CarSph)

use Index_Functions, only: nTri_Elem1, nTri3_Elem1
use Real_Spherical, only: ipSph, RSph
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: lMax
real(kind=wp), intent(inout) :: Cavxyz(nTri3_Elem1(lMax)), Cavsph((lMax+1)**2)
logical(kind=iwp), intent(in) :: CarSph
integer(kind=iwp) :: iOff1, iOff2, m, nElem

iOff1 = 1
iOff2 = 1
do m=0,lMax
  nElem = nTri_Elem1(m)
  !call RecPrt('Car-->Sph',' ',RSph(ipSph(m)),nElem,nElem)
  if (CarSph) then
    Cavsph(iOff2:iOff2+2*m) = Zero
    !call RecPrt('Cartesian',' ',Cavxyz(iOff1),1,nElem)
    call dGeMV_('T',nElem,2*m+1,One,RSph(ipSph(m)),nElem,Cavxyz(iOff1),1,Zero,CavSph(iOff2),1)
    !call RecPrt('Spherical',' ',Cavsph(iOff2),1,2*m+1)
  else
    Cavxyz(iOff1:iOff1+nElem-1) = Zero
    !call RecPrt('Spherical',' ',Cavsph(iOff2),1,2*m+1)
    call dGeMV_('N',nElem,2*m+1,One,RSph(ipSph(m)),nElem,Cavsph(iOff2),1,Zero,Cavxyz(iOff1),1)
    !call RecPrt('Cartesian',' ',Cavxyz(iOff1),1,nElem)
  end if
  iOff1 = iOff1+nElem
  iOff2 = iOff2+2*m+1
end do

end subroutine Tranca
