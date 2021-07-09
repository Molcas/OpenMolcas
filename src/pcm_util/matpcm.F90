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
use Constants, only: Zero, One, Two, Pi
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: NTs, ISPhe(*)
real(kind=wp), intent(in) :: Eps, Coor_Sph(4,*), Tessera(4,NTs)
real(kind=wp), intent(_OUT_) :: DMat(NTs,*), SMat(NTs,*), SDMat(NTs,*), TMat(NTs,*), RMat(NTs,*)
logical(kind=iwp), intent(in) :: Conductor
integer(kind=iwp) :: ITs, JTs, KTs, LI
real(kind=wp) :: EpsFac, Fac, FPI, Prod, RIJ, TPI, XI, XJ, XNI, YI, YJ, YNI, ZI, ZJ, ZNI

TPI = Two*PI
FPI = Two*TPI
if (Conductor) then

  ! Conductor model

  ! S matrix
  EpsFac = Eps/(Eps-One)
  call dcopy_(nTs*nTs,[Zero],0,SMat,1)
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
    call dcopy_(nTs*nTs,SMat,1,DMat,1)
  else
    call FZero(DMat,nTs**2)
  end if

else

  ! Dielectric model:

  ! S and D* matrices
  call dcopy_(nTs*nTs,[Zero],0,DMat,1)
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
  call dcopy_(nTs*nTs,[Zero],0,SDMat,1)
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
  do ITs=1,NTs
    TMat(ITs,ITs) = Fac*SMat(ITs,ITs)-SDMat(ITs,ITs)/TPI
    RMat(ITs,ITs) = -One+DMat(ITs,ITs)*Tessera(4,ITs)/TPI
    do JTs=1,NTs
      if (JTs == ITs) cycle
      TMat(ITs,JTs) = Fac*SMat(ITs,JTs)-SDMat(ITs,JTs)/TPI
      RMat(ITs,JTs) = Tessera(4,JTs)*DMat(JTs,ITs)/TPI
    end do
  end do

  ! Invert T matrix

  if (Eps > One) then
    call MatInvert(TMat,nTs)
  else
    call FZero(TMat,nTs**2)
  end if

  ! Form T^-1 * R and store it in D

  call DGEMM_('N','N',nTs,nTs,nTs,One,TMat,nTs,RMat,nTs,Zero,DMat,nTs)

end if

return

end subroutine MatPCM
