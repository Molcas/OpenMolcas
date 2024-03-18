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

subroutine dWdR(R,iNQ,Weights,list_p,nlist_p,invlist,dW_dR,nGrad_Eff,iTab,dW_Temp,dPB,nGrid)

use NQ_Structure, only: NQ_data
use nq_Grid, only: Pax
use Constants, only: Zero, One, Two, Three, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: iNQ, nlist_p, list_p(nlist_p), invlist(*), nGrad_Eff, iTab(4,nGrad_Eff), nGrid
real(kind=wp), intent(in) :: R(3,nGrid), Weights(nGrid)
real(kind=wp), intent(out) :: dW_dR(nGrad_Eff,nGrid), dW_Temp(3,nlist_p), dPB(3,nlist_p,nlist_p)
integer(kind=iwp) :: iA, iB, iC, iCar, iD, iGrad, iGrid, iiB, jNQ, kNQ, lNQ
real(kind=wp) :: dmu_BC_dA(3), dmu_BC_dB(3), dmu_BC_dC(3), dOdxs(3), dZ_dB(3), Fact, Osxyz(3), P_A, P_B, r_B, R_BC, R_BCxyz(3), &
                 r_Bxyz(3), r_C, r_Cxyz(3), rMU_BC, s_MU_BC, sxyz(3), temp, tMU_BC, xdiff0, xdiff1, xdiff2, xdiff3, Z
real(kind=wp), parameter :: Thrs = 1.0e-20_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! See J. Chem. Phys. 98 (1993) 5612 doi:10.1063/1.464906
!                                                                      *
!***********************************************************************
!                                                                      *
! iNQ is the index of the current atomic grid to which these grid points belong.

iA = invlist(iNQ)
!                                                                      *
!***********************************************************************
!                                                                      *
dW_dR(:,1:nGrid) = Zero
do_grid: do iGrid=1,nGrid
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The current grid point is associated with center A and the
  ! "atomic" displacement vector relative to center A is computed.

  Osxyz(:) = R(:,iGrid)-NQ_data(iNQ)%Coor(:)

  sxyz(1) = sum(Pax(:,1)*Osxyz(:))
  sxyz(2) = sum(Pax(:,2)*Osxyz(:))
  sxyz(3) = sum(Pax(:,3)*Osxyz(:))

  Z = Zero
  dPB(:,:,:) = Zero

  ! Compute all P_B and corresponding derivatives.

  !P_A = Zero
  P_A = One ! Dummy set
  do iiB=1,nlist_p
    if (iiB == 1) then
      iB = iA
    else if (iiB == iA) then
      iB = 1
    else
      iB = iiB
    end if

    kNQ = list_p(iB)
    r_Bxyz(:) = R(:,iGrid)-NQ_Data(kNQ)%Coor(:)
    r_B = sqrt(sum(r_Bxyz(:)**2))

    ! loop over C=/=B for all s(mu_BC), see Eq. B3

    P_B = One
    do iC=1,nlist_p

      if (iC /= iB) then
        lNQ = list_p(iC)
        r_Cxyz(:) = R(:,iGrid)-NQ_Data(lNQ)%Coor(:)
        r_C = sqrt(sum(r_Cxyz(:)**2))
        R_BCxyz(:) = NQ_Data(kNQ)%Coor(:)-NQ_Data(lNQ)%Coor(:)
        R_BC = sqrt(sum(R_BCxyz(:)**2))

        ! Eq. B6

        rMU_BC = (r_B-r_C)/R_BC
        if (rMU_BC <= Half) then
          xdiff0 = rMU_BC
          xdiff1 = (xdiff0*Half)*(Three-xdiff0**2)
          xdiff2 = (xdiff1*Half)*(Three-xdiff1**2)
          xdiff3 = (xdiff2*Half)*(Three-xdiff2**2)

          ! Eq. B4

          s_MU_BC = Half*(One-xdiff3)

          P_B = P_B*s_MU_BC
          if (P_B <= Thrs) exit
          tMU_BC = -27.0_wp*(One-xdiff2**2)*(One-xdiff1**2)*(One-rMU_BC**2)/(16.0_wp*max(s_MU_BC,1.0e-99_wp))
        else
          xdiff0 = rMU_BC-One
          xdiff1 = (-OneHalf-Half*xdiff0)*xdiff0**2
          xdiff2 = (-OneHalf-Half*xdiff1)*xdiff1**2
          xdiff3 = (-OneHalf-Half*xdiff2)*xdiff2**2
          s_MU_BC = -Half*xdiff3

          P_B = P_B*s_MU_BC
          if (P_B <= Thrs) exit
          tMU_BC = 27.0_wp*(Two+xdiff2)*xdiff2*(Two+xdiff1)*xdiff1*(Two+xdiff0)*xdiff0/(16.0_wp*max(s_MU_BC,1.0e-99_wp))
        end if

        ! Differentiate mu_BC with respect to the center, D.

        do iD=1,nlist_p
          !jNQ = list_p(iD)

          if (iD == iB) then

            ! d mu_BC(r_A) / dB, Eq. B10

            if (r_B == Zero) then
              dmu_BC_dB(:) = r_C*R_BCxyz(:)/R_BC**3
            else
              dmu_BC_dB(:) = -r_Bxyz(:)/(r_B*R_BC)-(r_B-r_C)*R_BCxyz(:)/R_BC**3
            end if

            dPB(:,iB,iB) = dPB(:,iB,iB)+tMU_BC*dmu_BC_dB(:)

          else if (iD == iC) then

            ! d mu_BC(r_A) / dC, Eq, B10

            if (r_C == Zero) then
              dmu_BC_dC(:) = r_B*R_BCxyz(:)/R_BC**3
            else
              dmu_BC_dC(:) = r_Cxyz(:)/(r_C*R_BC)+(r_B-r_C)*R_BCxyz(:)/R_BC**3
            end if
            dPB(:,iC,iB) = dPB(:,iC,iB)+tMU_BC*dmu_BC_dC(:)

          end if

          ! d mu_BC(r_A) / dr_A

          if (r_B == Zero) then
            dmu_BC_dA(:) = (-r_Cxyz(:)/r_C)/R_BC
          else if (r_C == Zero) then
            dmu_BC_dA(:) = (r_Bxyz(:)/r_B)/R_BC
          else
            dmu_BC_dA(:) = (r_Bxyz(:)/r_B-r_Cxyz(:)/r_C)/R_BC
          end if

          ! The direct term

          if (iD == iA) dPB(:,iA,iB) = dPB(:,iA,iB)+tMU_BC*dmu_BC_dA(:)

          jNQ = list_p(iD)
          do iCar=1,3
            dOdxs(1) = sum(NQ_Data(jNQ)%dOdx(1,:,iCar)*sxyz(:))
            dOdxs(2) = sum(NQ_Data(jNQ)%dOdx(2,:,iCar)*sxyz(:))
            dOdxs(3) = sum(NQ_Data(jNQ)%dOdx(3,:,iCar)*sxyz(:))
            temp = tMU_BC*sum(dmu_BC_dA(:)*dOdxs(:))

            dPB(iCar,iD,iB) = dPB(iCar,iD,iB)-temp
          end do

        end do ! iD

      end if
    end do     ! iC

    ! Multiply derivatives with P_B as in Eq. B8

    dPB(:,:,iB) = P_B*dPB(:,:,iB)

    if (iB == iA) P_A = P_B
    if (P_A <= Thrs) cycle do_grid

    ! Denominator Eq. B2
    Z = Z+P_B
  end do ! iB

  if (P_A == Zero) then
    Fact = Zero
  else
    Fact = Weights(iGrid)*Z/P_A
  end if
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Assemble the gradient

  do iB=1,nlist_p

    dZ_dB(1) = sum(dPB(1,iB,:))
    dZ_dB(2) = sum(dPB(2,iB,:))
    dZ_dB(3) = sum(dPB(3,iB,:))

    ! Eq. B7

    dW_Temp(:,iB) = dPB(:,iB,iA)/Z-(P_A*dZ_dB(:))/Z**2
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pick up the relevant gradients

  do iGrad=1,nGrad_Eff
    iCar = iTab(1,iGrad)
    kNQ = iTab(3,iGrad)
    dW_dR(iGrad,iGrid) = Fact*dW_Temp(iCar,invlist(kNQ))
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end do do_grid ! iGrid
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine dWdR
