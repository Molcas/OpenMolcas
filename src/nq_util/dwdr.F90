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

subroutine dWdR(R,ilist_p,Weights,list_p,nlist_p,dW_dR,nGrad_Eff,iTab,dW_Temp,dPB,nGrid)

use NQ_Structure, only: NQ_data
use nq_Grid, only: Pax
use Constants, only: Zero, One, Two, Three, Half, OneHalf
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ilist_p, nlist_p, list_p(nlist_p), nGrad_Eff, iTab(4,nGrad_Eff), nGrid
real(kind=wp), intent(in) :: R(3,nGrid), Weights(nGrid)
real(kind=wp), intent(out) :: dW_dR(nGrad_Eff,nGrid), dW_Temp(3,nlist_p), dPB(3,nlist_p,nlist_p)
integer(kind=iwp) :: iA, iB, iC, iCar, iD, iGrad, iGrid, iiB, iNQ, jNQ, kNQ, lNQ
real(kind=wp) :: dmu_BC_dAx, dmu_BC_dAy, dmu_BC_dAz, dmu_BC_dBx, dmu_BC_dBy, dmu_BC_dBz, dmu_BC_dCx, dmu_BC_dCy, dmu_BC_dCz, &
                 dOdx_11, dOdx_12, dOdx_13, dOdx_21, dOdx_22, dOdx_23, dOdx_31, dOdx_32, dOdx_33, dOdxs(3), dZ_dBx, dZ_dBy, &
                 dZ_dBz, Fact, Fact0, O11, O12, O13, O21, O22, O23, O31, O32, O33, Osxyz(3), p1, p2, p3, P_A, P_B, r_B, R_BC, &
                 R_BCx, R_BCy, R_BCz, r_Bx, r_By, r_Bz, r_C, r_Cx, r_Cy, r_Cz, rMU_BC, s_MU_BC, sxyz(3), temp, tMU_BC, xdiff0, &
                 xdiff1, xdiff2, Z
real(kind=wp), parameter :: Thrs = 1.0e-20_wp

!                                                                      *
!***********************************************************************
!                                                                      *
! iNQ is the index of the current atomic grid to which these grid
! points belong.

iNQ = list_p(ilist_p)
iA = ilist_p
O11 = Pax(1,1)
O12 = Pax(2,1)
O13 = Pax(3,1)
O21 = Pax(1,2)
O22 = Pax(2,2)
O23 = Pax(3,2)
O31 = Pax(1,3)
O32 = Pax(2,3)
O33 = Pax(3,3)
!                                                                      *
!***********************************************************************
!                                                                      *
do_grid: do iGrid=1,nGrid
  dW_dR(:,iGrid) = Zero
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! The current grid point is associated with center A and the
  ! "atomic" displacement vector relative to center A is computed.

  Osxyz(:) = R(:,iGrid)-NQ_data(iNQ)%Coor(:)

  sxyz(1) = O11*Osxyz(1)+O12*Osxyz(2)+O13*Osxyz(3)
  sxyz(2) = O21*Osxyz(1)+O22*Osxyz(2)+O23*Osxyz(3)
  sxyz(3) = O31*Osxyz(1)+O32*Osxyz(2)+O33*Osxyz(3)

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
    r_Bx = R(1,iGrid)-NQ_Data(kNQ)%Coor(1)
    r_By = R(2,iGrid)-NQ_Data(kNQ)%Coor(2)
    r_Bz = R(3,iGrid)-NQ_Data(kNQ)%Coor(3)
    r_B = sqrt(r_Bx**2+r_By**2+r_Bz**2)

    ! loop over C=/=B for all s(mu_BC), see Eq. B3

    P_B = One
    do iC=1,nlist_p

      if (iC /= iB) then
        lNQ = list_p(iC)
        r_Cx = R(1,iGrid)-NQ_Data(lNQ)%Coor(1)
        r_Cy = R(2,iGrid)-NQ_Data(lNQ)%Coor(2)
        r_Cz = R(3,iGrid)-NQ_Data(lNQ)%Coor(3)
        r_C = sqrt(r_Cx**2+r_Cy**2+r_Cz**2)
        R_BCx = NQ_Data(kNQ)%Coor(1)-NQ_Data(lNQ)%Coor(1)
        R_BCy = NQ_Data(kNQ)%Coor(2)-NQ_Data(lNQ)%Coor(2)
        R_BCz = NQ_Data(kNQ)%Coor(3)-NQ_Data(lNQ)%Coor(3)
        R_BC = sqrt(R_BCx**2+R_BCy**2+R_BCz**2)

        ! Eq. B6

        rMU_BC = (r_B-r_C)/R_BC
        if (rMU_BC <= Half) then
          p1 = (rMU_BC*Half)*(Three-rMU_BC**2)
          p2 = (p1*Half)*(Three-p1**2)
          p3 = (p2*Half)*(Three-p2**2)

          ! Eq. B4

          s_MU_BC = Half*(One-p3)

          P_B = P_B*s_MU_BC
          if (P_B <= Thrs) exit
          tMU_BC = -27.0_wp*(One-p2**2)*(One-p1**2)*(One-rMU_BC**2)/(16.0_wp*max(s_MU_BC,1.0e-99_wp))
        else
          xdiff0 = rMU_BC-One
          xdiff1 = (-OneHalf-Half*xdiff0)*xdiff0**2
          xdiff2 = (-OneHalf-Half*xdiff1)*xdiff1**2
          p3 = (OneHalf+Half*xdiff2)*xdiff2**2
          s_MU_BC = Half*p3

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
              dmu_BC_dBx = -(r_B-r_C)*R_BCx/R_BC**3
              dmu_BC_dBy = -(r_B-r_C)*R_BCy/R_BC**3
              dmu_BC_dBz = -(r_B-r_C)*R_BCz/R_BC**3
            else
              dmu_BC_dBx = -r_Bx/(r_B*R_BC)-(r_B-r_C)*R_BCx/R_BC**3
              dmu_BC_dBy = -r_By/(r_B*R_BC)-(r_B-r_C)*R_BCy/R_BC**3
              dmu_BC_dBz = -r_Bz/(r_B*R_BC)-(r_B-r_C)*R_BCz/R_BC**3
            end if

            dPB(1,iB,iB) = dPB(1,iB,iB)+tMU_BC*dmu_BC_dBx
            dPB(2,iB,iB) = dPB(2,iB,iB)+tMU_BC*dmu_BC_dBy
            dPB(3,iB,iB) = dPB(3,iB,iB)+tMU_BC*dmu_BC_dBz

          else if (iD == iC) then

            ! d mu_BC(r_A) / dC, Eq, B10

            if (r_C == Zero) then
              dmu_BC_dCx = +(r_B-r_C)*R_BCx/R_BC**3
              dmu_BC_dCy = +(r_B-r_C)*R_BCy/R_BC**3
              dmu_BC_dCz = +(r_B-r_C)*R_BCz/R_BC**3
            else
              dmu_BC_dCx = r_Cx/(r_C*R_BC)+(r_B-r_C)*R_BCx/R_BC**3
              dmu_BC_dCy = r_Cy/(r_C*R_BC)+(r_B-r_C)*R_BCy/R_BC**3
              dmu_BC_dCz = r_Cz/(r_C*R_BC)+(r_B-r_C)*R_BCz/R_BC**3
            end if
            dPB(1,iC,iB) = dPB(1,iC,iB)+tMU_BC*dmu_BC_dCx
            dPB(2,iC,iB) = dPB(2,iC,iB)+tMU_BC*dmu_BC_dCy
            dPB(3,iC,iB) = dPB(3,iC,iB)+tMU_BC*dmu_BC_dCz

          end if

          ! d mu_BC(r_A) / dr_A

          if (r_B == Zero) then
            dmu_BC_dAx = (-r_Cx/r_C)/R_BC
            dmu_BC_dAy = (-r_Cy/r_C)/R_BC
            dmu_BC_dAz = (-r_Cz/r_C)/R_BC
          else if (r_C == Zero) then
            dmu_BC_dAx = (r_Bx/r_B)/R_BC
            dmu_BC_dAy = (r_By/r_B)/R_BC
            dmu_BC_dAz = (r_Bz/r_B)/R_BC
          else
            dmu_BC_dAx = (r_Bx/r_B-r_Cx/r_C)/R_BC
            dmu_BC_dAy = (r_By/r_B-r_Cy/r_C)/R_BC
            dmu_BC_dAz = (r_Bz/r_B-r_Cz/r_C)/R_BC
          end if

          if (iD == iA) then

            ! The direct term

            dPB(1,iA,iB) = dPB(1,iA,iB)+tMU_BC*dmu_BC_dAx
            dPB(2,iA,iB) = dPB(2,iA,iB)+tMU_BC*dmu_BC_dAy
            dPB(3,iA,iB) = dPB(3,iA,iB)+tMU_BC*dmu_BC_dAz

          end if

          jNQ = list_p(iD)
          do iCar=1,3
            dOdx_11 = NQ_Data(jNQ)%dOdx(1,1,iCar)
            dOdx_21 = NQ_Data(jNQ)%dOdx(2,1,iCar)
            dOdx_31 = NQ_Data(jNQ)%dOdx(3,1,iCar)
            dOdx_12 = NQ_Data(jNQ)%dOdx(1,2,iCar)
            dOdx_22 = NQ_Data(jNQ)%dOdx(2,2,iCar)
            dOdx_32 = NQ_Data(jNQ)%dOdx(3,2,iCar)
            dOdx_13 = NQ_Data(jNQ)%dOdx(1,3,iCar)
            dOdx_23 = NQ_Data(jNQ)%dOdx(2,3,iCar)
            dOdx_33 = NQ_Data(jNQ)%dOdx(3,3,iCar)
            dOdxs(1) = dOdx_11*sxyz(1)+dOdx_12*sxyz(2)+dOdx_13*sxyz(3)
            dOdxs(2) = dOdx_21*sxyz(1)+dOdx_22*sxyz(2)+dOdx_23*sxyz(3)
            dOdxs(3) = dOdx_31*sxyz(1)+dOdx_32*sxyz(2)+dOdx_33*sxyz(3)
            temp = tMU_BC*(dmu_BC_dAx*dOdxs(1)+dmu_BC_dAy*dOdxs(2)+dmu_BC_dAz*dOdxs(3))

            dPB(iCar,iD,iB) = dPB(iCar,iD,iB)-temp
          end do

        end do ! iD

      end if
    end do     ! iC

    ! Multiply derivatives with P_B as in Eq. B8

    do iD=1,nlist_p
      dPB(1,iD,iB) = P_B*dPB(1,iD,iB)
      dPB(2,iD,iB) = P_B*dPB(2,iD,iB)
      dPB(3,iD,iB) = P_B*dPB(3,iD,iB)
    end do

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

    dZ_dBx = Zero
    dZ_dBy = Zero
    dZ_dBz = Zero
    do iC=1,nlist_p
      dZ_dBx = dZ_dBx+dPB(1,iB,iC)
      dZ_dBy = dZ_dBy+dPB(2,iB,iC)
      dZ_dBz = dZ_dBz+dPB(3,iB,iC)
    end do

    ! Eq. B7

    dW_Temp(1,iB) = dPB(1,iB,iA)/Z-(P_A*dZ_dBx)/Z**2
    dW_Temp(2,iB) = dPB(2,iB,iA)/Z-(P_A*dZ_dBy)/Z**2
    dW_Temp(3,iB) = dPB(3,iB,iA)/Z-(P_A*dZ_dBz)/Z**2
  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Pick up the relevant gradients

  do iGrad=1,nGrad_Eff
    iCar = iTab(1,iGrad)
    kNQ = iTab(3,iGrad)
    Fact0 = Fact
    do iB=1,nlist_p
      lNQ = list_p(iB)
      if (kNQ == lNQ) dW_dR(iGrad,iGrid) = Fact0*dW_Temp(iCar,iB)
    end do
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
