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

subroutine GeoDer(nAt,Cond,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerDM,Grd,DerTes,DerPunt,DerRad,DerCentr)

use Constants, only: Zero, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAt, nTs, nS, ISphe(nTs), NOrd(nS)
logical(kind=iwp), intent(in) :: Cond
real(kind=wp), intent(in) :: Eps, Sphere(4,nS), Tessera(4,nTs), Q(2,nTs), DerTes(nTs,nAt,3), DerPunt(nTs,nAt,3,3), &
                             DerRad(nS,nAt,3), DerCentr(nS,nAt,3,3)
real(kind=wp), intent(out) :: DerDM(nTs,nTs), Grd(3,nAt)
integer(kind=iwp) :: IAtom, iTs, IXYZ, jTs
real(kind=wp) :: GeoGrd, Qi, Qj

! Compute the PCM geometric contribution to gradients

Grd(:,:) = Zero
DerDM(:,:) = Zero
do IAtom=1,nAt
  do IXYZ=1,3
    ! Dielectric model
    if (.not. Cond) then
      call Over(IAtom,IXYZ,GeoGrd,nAt,nTs,nS,Eps,Sphere,ISphe,NOrd,Tessera,Q,DerRad,DerCentr)
    ! Conductor model
    else if (Cond) then
      GeoGrd = Zero
      call DerD(IAtom,IXYZ,Tessera,ISphe,DerDM,DerTes,DerPunt,DerCentr,nTs,nAt,nS)
      do iTs=1,nTs
        Qi = Q(1,iTs)+Q(2,iTs)
        do jTs=1,nTs
          Qj = Q(1,jTs)+Q(2,jTs)
          GeoGrd = GeoGrd+Qi*DerDM(iTs,jTs)*Qj
        end do
      end do
    end if
    Grd(IXYZ,IAtom) = GeoGrd*Half
  end do
end do

return

end subroutine GeoDer
