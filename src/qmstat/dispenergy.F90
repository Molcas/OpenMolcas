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

subroutine DispEnergy(EEDisp,BoMaH,BoMaO,dAtO1,dAtH1,dAtH2,Rab13i,Rab23i,Rab33i,indQAt)

use qmstat_global, only: DispDamp, uDisp
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: EEDisp
real(kind=wp), intent(in) :: BoMaH, BoMaO, dAtO1, dAtH1, dAtH2, Rab13i, Rab23i, Rab33i
integer(kind=iwp), intent(in) :: indQAt
integer(kind=iwp) :: k, kFac
real(kind=wp) :: BM(3), DampBarn(3), EfromH1, EfromH2, EfromO1

! If not damping, set factors to 1.0

DampBarn(:) = One

if (DispDamp) then

  ! If Damping, do it, otherwise easy stuff

  ! Get the damping, for the relevant QM-atom with all solvent atoms.
  BM(1) = BoMaH*dAtH1
  BM(2) = BoMaH*dAtH2
  BM(3) = BoMaO*dAtO1

  kFac = 1
  do k=1,6
    kFac = kFac*k
    DampBarn(:) = DampBarn(:)+BM(:)**k/real(kFac,kind=wp)
  end do
  DampBarn(:) = One-DampBarn(:)*exp(-BM(:))

end if

! Now evaluate the Dispersion energy.

EfromH1 = Rab23i**2*DampBarn(1)*uDisp(2,indQAt)
EfromH2 = Rab33i**2*DampBarn(2)*uDisp(2,indQAt)
EfromO1 = Rab13i**2*DampBarn(3)*uDisp(1,indQAt)

EEDisp = EEdisp+EfromO1+EfromH1+EfromH2

return

end subroutine DispEnergy
