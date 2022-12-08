!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************

! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 22, 2021, created this file.               *
! ****************************************************************
subroutine Calc_Pot1(Pot1,TabMO,mAO,mGrid,nMOs,P2_ontop,nP2_ontop,MOs)

use nq_Grid, only: GradRho, vRho, vSigma, Weights
use nq_pdft, only: d2RdRho2, d2RdRhodPi, d2ZdR2, dEdRho, dEdRhox, dEdRhoy, dEdRhoz, dF_dRhoapb, dF_dRhoamb, dF_dRhoxamb, &
                   dF_dRhoxapb, dF_dRhoyamb, dF_dRhoyapb, dF_dRhozamb, dF_dRhozapb, dRdPi, dRdRho, dRhodx, dRhody, dRhodz, dZdR, &
                   dZdRho, fta, ftb, ftc, GradPidFdRho, GradRdFdRho, GradRhodFdRho, lft, lGGA, Pass1, Pass2, Pass3, RatioA, RhoAB, &
                   ThrsNT, ZetaA
use nq_Info, only: mIrrep, mOrb, nOrbt, nPot1, OffBasFro, OffOrb, OffOrb2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Six, Twelve, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(inout) :: Pot1(nPot1)
integer(kind=iwp), intent(in) :: mAO, mGrid, nMOs, nP2_ontop
real(kind=wp), intent(in) :: TabMO(mAO,mGrid,nMOs), MOs(mGrid,nOrbt), P2_ontop(nP2_ontop,mGrid)
integer(kind=iwp) :: iGrid, iIrrep, iMO, iOff1, IOff2, iOrb
real(kind=wp) :: dEdRhop2, dF_dRhoax, dF_dRhoay, dF_dRhoaz, dF_dRhobx, dF_dRhoby, dF_dRhobz, Diff1, dRdx, dRdy, dRdz
! PreMO is MO multiplied with things needed for potential calculation
real(kind=wp), allocatable :: PreMO(:,:)

call mma_allocate(PreMO,mGrid,nOrbt)
PreMO(:,:) = MOs

! black terms in the notes
do iGrid=1,mGrid
  if (Pass1(iGrid)) then
    dF_dRhoapb(iGrid) = vRho(1,iGrid)+vRho(2,iGrid)
    dF_dRhoamb(iGrid) = vRho(1,iGrid)-vRho(2,iGrid)
    dRdRho(iGrid) = RatioA(iGrid)*(-Two/RhoAB(iGrid))
    dZdRho(iGrid) = dZdR(iGrid)*dRdRho(iGrid)
    dEdRho(iGrid) = dF_dRhoapb(iGrid)+dF_dRhoamb(iGrid)*(ZetaA(iGrid)+RhoAB(iGrid)*dZdRho(iGrid))
    dRdPi(iGrid) = Four/RhoAB(iGrid)**2
  else
    dRdPi(iGrid) = Zero
    dF_dRhoapb(iGrid) = Zero
    dF_dRhoamb(iGrid) = Zero
    dEdRho(iGrid) = Zero
  end if
end do

! red terms in the notes
if (lGGA) then
  do iGrid=1,mGrid
    if (Pass1(iGrid)) then
      dF_dRhoax = Two*vSigma(1,iGrid)*GradRho(1,iGrid)+vSigma(2,iGrid)*GradRho(4,iGrid)
      dF_dRhobx = Two*vSigma(3,iGrid)*GradRho(4,iGrid)+vSigma(2,iGrid)*GradRho(1,iGrid)
      dF_dRhoay = Two*vSigma(1,iGrid)*GradRho(2,iGrid)+vSigma(2,iGrid)*GradRho(5,iGrid)
      dF_dRhoby = Two*vSigma(3,iGrid)*GradRho(5,iGrid)+vSigma(2,iGrid)*GradRho(2,iGrid)
      dF_dRhoaz = Two*vSigma(1,iGrid)*GradRho(3,iGrid)+vSigma(2,iGrid)*GradRho(6,iGrid)
      dF_dRhobz = Two*vSigma(3,iGrid)*GradRho(6,iGrid)+vSigma(2,iGrid)*GradRho(3,iGrid)
      dF_dRhoxapb(iGrid) = dF_dRhoax+dF_dRhobx
      dF_dRhoxamb(iGrid) = dF_dRhoax-dF_dRhobx
      dF_dRhoyapb(iGrid) = dF_dRhoay+dF_dRhoby
      dF_dRhoyamb(iGrid) = dF_dRhoay-dF_dRhoby
      dF_dRhozapb(iGrid) = dF_dRhoaz+dF_dRhobz
      dF_dRhozamb(iGrid) = dF_dRhoaz-dF_dRhobz

      dRhodx(iGrid) = GradRho(1,iGrid)+GradRho(4,iGrid)
      dRhody(iGrid) = GradRho(2,iGrid)+GradRho(5,iGrid)
      dRhodz(iGrid) = GradRho(3,iGrid)+GradRho(6,iGrid)

      GradRhodFdRho(iGrid) = (dF_dRhoxamb(iGrid)*dRhodx(iGrid)+dF_dRhoyamb(iGrid)*dRhody(iGrid)+dF_dRhozamb(iGrid)*dRhodz(iGrid))
      dEdRho(iGrid) = dEdRho(iGrid)+dZdRho(iGrid)*GradRhodFdRho(iGrid)

      dEdRhox(iGrid) = dF_dRhoxapb(iGrid)+ZetaA(iGrid)*dF_dRhoxamb(iGrid)
      dEdRhoy(iGrid) = dF_dRhoyapb(iGrid)+ZetaA(iGrid)*dF_dRhoyamb(iGrid)
      dEdRhoz(iGrid) = dF_dRhozapb(iGrid)+ZetaA(iGrid)*dF_dRhozamb(iGrid)

    else
      dF_dRhoxapb(iGrid) = Zero
      dF_dRhoxamb(iGrid) = Zero
      dF_dRhoyapb(iGrid) = Zero
      dF_dRhoyamb(iGrid) = Zero
      dF_dRhozapb(iGrid) = Zero
      dF_dRhozamb(iGrid) = Zero
      GradRhodFdRho(iGrid) = Zero
      dEdRhox(iGrid) = Zero
      dEdRhoy(iGrid) = Zero
      dEdRhoz(iGrid) = Zero
    end if
  end do
  ! green and blue terms in the notes
  if (lft) then
    do iGrid=1,mGrid
      if (Pass1(iGrid)) then
        dRdX = dRdRho(iGrid)*dRhodX(iGrid)+dRdPi(iGrid)*P2_ontop(2,iGrid)
        dRdY = dRdRho(iGrid)*dRhodY(iGrid)+dRdPi(iGrid)*P2_ontop(3,iGrid)
        dRdZ = dRdRho(iGrid)*dRhodZ(iGrid)+dRdPi(iGrid)*P2_ontop(4,iGrid)
        GradRdFdRho(iGrid) = dRdX*dF_dRhoxamb(iGrid)+dRdY*dF_dRhoyamb(iGrid)+dRdZ*dF_dRhozamb(iGrid)
        GradPidFdRho(iGrid) = P2_ontop(2,iGrid)*dF_dRhoxamb(iGrid)+P2_ontop(3,iGrid)*dF_dRhoyamb(iGrid)+ &
                              P2_ontop(4,iGrid)*dF_dRhozamb(iGrid)
        d2RdRho2(iGrid) = Six*RatioA(iGrid)/RhoAB(iGrid)**2
        d2RdRhodPi(iGrid) = -Two*dRdPi(iGrid)/RhoAB(iGrid)
        if (Pass2(iGrid)) then
          d2ZdR2(iGrid) = Two*dZdR(iGrid)**3
        else if (Pass3(iGrid)) then
          Diff1 = RatioA(iGrid)-ThrsNT
          d2ZdR2(iGrid) = (20.0_wp*fta*Diff1**2+Twelve*ftb*Diff1+Six*ftc)*Diff1
        else
          d2ZdR2(iGrid) = Zero
        end if
        dEdRho(iGrid) = dEdRho(iGrid)+(dZdR(iGrid)+RhoAB(iGrid)*d2ZdR2(iGrid)*dRdRho(iGrid))*GradRdFdRho(iGrid)+ &
                        RhoAB(iGrid)*dZdR(iGrid)*d2RdRho2(iGrid)*GradRhodFdRho(iGrid)+ &
                        RhoAB(iGrid)*dZdR(iGrid)*d2RdRhodPi(iGrid)*GradPidFdRho(iGrid)
        dEdRhop2 = RhoAB(iGrid)*dZdRho(iGrid)
        ! end of dEdRho term
        ! now dEdRhoprime terms
        dEdRhox(iGrid) = dEdRhox(iGrid)+dEdRhop2*dF_dRhoxamb(iGrid)
        dEdRhoy(iGrid) = dEdRhoy(iGrid)+dEdRhop2*dF_dRhoyamb(iGrid)
        dEdRhoz(iGrid) = dEdRhoz(iGrid)+dEdRhop2*dF_dRhozamb(iGrid)
      else
        GradRdFdRho(iGrid) = Zero
        GradPidFdRho(iGrid) = Zero
        d2RdRho2(iGrid) = Zero
        d2RdRhodPi(iGrid) = Zero
        d2ZdR2(iGrid) = Zero
      end if
    end do
  end if
end if

dEdRho(:) = Half*dEdRho

do iGrid=1,mGrid
  PreMO(iGrid,:) = PreMO(iGrid,:)*dEdRho(iGrid)
end do

if (lGGA) then
  do iIrrep=0,mIrrep-1
    do iOrb=1,mOrb(iIrrep)
      IOff1 = iOrb+OffOrb(iIrrep)
      iMO = iOrb+OffBasFro(iIrrep)
      do iGrid=1,mGrid
        PreMO(iGrid,IOff1) = PreMO(iGrid,IOff1)+TabMO(2,iGrid,iMO)*dEdRhox(iGrid)+TabMO(3,iGrid,iMO)*dEdRhoy(iGrid)+ &
                             TabMO(4,iGrid,iMO)*dEdRhoz(iGrid)
      end do
    end do
  end do
end if

do iGrid=1,mGrid
  PreMO(iGrid,:) = PreMO(iGrid,:)*Weights(iGrid)
end do

do iIrrep=0,mIrrep-1
  IOff1 = OffOrb(iIrrep)+1
  IOff2 = OffOrb2(iIrrep)+1
  call DGEMM_('T','N',mOrb(iIrrep),mOrb(iIrrep),mGrid,One,PreMO(:,IOff1:),mGrid,MOs(:,IOff1:),mGrid,One,Pot1(iOff2:),mOrb(iIrrep))
end do

call mma_deallocate(PreMO)

return

end subroutine Calc_Pot1
