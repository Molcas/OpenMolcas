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

subroutine MatPCM(NTs,Eps,Conductor,ISphe,Coor_Sph,Tessera,DMat,SMat,SDMat,TMat,RMat)
! Compute PCM matrix with the formalism in
! M. Cossi, N. Rega, G. Scalmani, V. Barone JCP in press;
! D.M. Chipman JCP 112, 5558 (2000).
! Solvation charges are defined through
! Tq = RV
! where V is the solute electrostatic potential. Here T^-1*R is computed
! and finally returned in DMat.

use PCM_Arrays, only: DiagScale
use Constants, only: Zero, One, Two, Four, Pi
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: NTs, ISPhe(NTs)
real(kind=wp), intent(in) :: Eps, Coor_Sph(4,*), Tessera(4,NTs)
real(kind=wp), intent(out) :: DMat(NTs,NTs), SMat(NTs,NTs), SDMat(NTs,NTs), TMat(NTs,NTs), RMat(NTs,NTs)
logical(kind=iwp), intent(in) :: Conductor
integer(kind=iwp) :: ITs, JTs, KTs, LI
real(kind=wp) :: EpsFac, Fac, Prod, RIJ, XI, XJ, XNI, YI, YJ, YNI, ZI, ZJ, ZNI
real(kind=wp), parameter :: FPI = Four*Pi, TPI = Two*Pi

if (Conductor) then

  ! Conductor model

  ! S matrix
  EpsFac = Eps/(Eps-One)
  SMat(:,:) = Zero
  do ITs=1,NTs
    XI = Tessera(1,iTs)
    YI = Tessera(2,iTs)
    ZI = Tessera(3,iTs)
    SMat(ITs,ITs) = -DiagScale*EpsFac*sqrt(FPI/Tessera(4,ITs))
    do JTs=1,ITs-1
      XJ = Tessera(1,jTs)
      YJ = Tessera(2,jTs)
      ZJ = Tessera(3,jTs)
      RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      SMat(ITs,JTs) = -EpsFac*One/RIJ
      SMat(JTs,ITs) = SMat(ITs,JTs)
    end do
  end do

  ! Invert S matrix and store it in D

  if (Eps > One) then
    call MatInvert(SMat,nTs)
    DMat(:,:) = SMat(:,:)
  else
    DMat(:,:) = Zero
  end if

else

  ! Dielectric model:

  ! S and D* matrices
  DMat(:,:) = Zero
  do ITs=1,NTs
    XI = Tessera(1,iTs)
    YI = Tessera(2,iTs)
    ZI = Tessera(3,iTs)
    LI = ISphe(ITs)
    XNI = (XI-Coor_Sph(1,LI))/Coor_Sph(4,LI)
    YNI = (YI-Coor_Sph(2,LI))/Coor_Sph(4,LI)
    ZNI = (ZI-Coor_Sph(3,LI))/Coor_Sph(4,LI)
    SMat(ITs,ITs) = DiagScale*sqrt(FPI/Tessera(4,ITs))
    DMat(ITs,ITs) = DMat(ITs,ITs)-TPI/Tessera(4,ITs)
    do JTs=1,NTs
      if (JTs == ITs) cycle
      XJ = Tessera(1,jTs)
      YJ = Tessera(2,jTs)
      ZJ = Tessera(3,jTs)
      RIJ = sqrt((XI-XJ)**2+(YI-YJ)**2+(ZI-ZJ)**2)
      SMat(ITs,JTs) = One/RIJ
      Prod = (XI-XJ)*XNI+(YI-YJ)*YNI+(ZI-ZJ)*ZNI
      DMat(ITs,JTs) = -Prod/RIJ**3
      DMat(JTs,JTs) = DMat(JTs,JTs)-DMat(ITs,JTs)*Tessera(4,ITs)/Tessera(4,JTs)
    end do
  end do

  ! S*A*D matrix
  SDMat(:,:) = Zero
  do ITs=1,NTs
    do JTs=1,NTs
      do KTs=1,NTs
        SDMat(ITs,JTs) = SDMat(ITs,JTs)+SMat(ITs,KTs)*Tessera(4,KTs)*DMat(KTs,JTs)
      end do
    end do
  end do

  ! The charges are defined as
  ! q = T-1 R V,         T = f(e)*S - SAD / 2p

  ! T and R matrices
  Fac = (Eps+One)/(Eps-One)
  TMat(:,:) = Fac*SMat(:,:)-SDMat(:,:)/TPI
  do ITs=1,NTs
    do JTs=1,NTs
      RMat(ITs,JTs) = DMat(JTs,ITs)*Tessera(4,JTs)/TPI
    end do
    RMat(ITs,ITs) = RMat(ITs,ITs)-One
  end do

  ! Invert T matrix

  if (Eps > One) then
    call MatInvert(TMat,nTs)
  else
    TMat(:,:) = Zero
  end if

  ! Form T^-1 * R and store it in D

  call DGEMM_('N','N',nTs,nTs,nTs,One,TMat,nTs,RMat,nTs,Zero,DMat,nTs)

end if

return

end subroutine MatPCM
