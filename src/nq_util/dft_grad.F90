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
! Copyright (C) 2002, Roland Lindh                                     *
!***********************************************************************

subroutine DFT_Grad(Grad,nGrad,nD,Grid,mGrid,dRho_dR,ndRho_dR,nGrad_Eff,Weights,iNQ)
!***********************************************************************
!                                                                      *
!     Object: to trace the correct parts to get the contributions to   *
!             the gradient due to the DFT energy.                      *
!                                                                      *
!     Author: Roland Lindh, Dept. of Chemical Physics, University of   *
!             Lund, Sweden.  May 2002 in Bologna, Italy.               *
!***********************************************************************

use nq_Grid, only: dW_dR, F_xc, GradRho, IndGrd, iTab, Pax, Temp, vLapl, vRho, vSigma, vTau
use nq_Structure, only: NQ_data
use nq_Info, only: Functional_type, GGA_Type, Grid_Type, LDA_Type, meta_GGA_type1, meta_GGA_type2, Moving_Grid, Off
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half, Quart
use Definitions, only: wp, iwp
!define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: nGrad, nD, mGrid, ndRho_dR, nGrad_Eff, iNQ
real(kind=wp), intent(inout) :: Grad(nGrad)
real(kind=wp), intent(in) :: Grid(3,mGrid), dRho_dR(ndRho_dR,mGrid,nGrad_Eff), Weights(mGrid)
integer(kind=iwp) :: i, i_Eff, iCar, iGrad, ixyz, j, jGrad, jNQ
real(kind=wp) :: dF_dr, Fact, gxa, gxb, gya, gyb, gza, gzb, OV(3,3), OVT(3), R_Grid(3), tmp, V(3,3)
real(kind=wp), allocatable :: Aux(:,:)
real(kind=wp), external :: DDot_

!                                                                      *
!***********************************************************************
!                                                                      *
R_Grid(:) = NQ_Data(iNQ)%Coor(:)
#ifdef _DEBUGPRINT_
call RecPrt('R_Grid',' ',R_Grid,1,3)
call RecPrt('Grid',' ',Grid,3,mGrid)
call RecPrt('Weights',' ',Weights,1,mGrid)
call RecPrt('dW_dR',' ',dW_dR,nGrad_Eff,mGrid)
call RecPrt('dRho_dR(1)',' ',dRho_dR,ndRho_dR,mGrid)
call RecPrt('dF_dRho',' ',dF_dRho,ndF_dRho,mGrid)
do iEff=1,nGrad_Eff
  write(u6,*) 'iTab=',iTab(1,iEff),iTab(2,iEff),iTab(3,iEff),iTab(3,iEff)
  write(u6,*) 'IndGrd=',IndGrd(iEff)
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! We have that the DFT energy is expressed as
!
! E_DFT = Sum_Gg  w(r_g(G))  f(G,r_g(G))
!
! r_g = R_G + O s_g
!
! The first derivative is computed as
!
! E_DFT^x = Sum w^x f  + w f^x
!
! where
!
! f^x = f^(x) + <nabla_r  f * r^x >
!
! where
!
! nabla_r f is the functional differentiated with respect to a
! displacement of a grid point.
!
! and
!
! r^x = delta_AG e_i + O^x s
!                                                                      *
!***********************************************************************
!                                                                      *
! Add the contributions
!
! w * f^x to centers other than the origin of the current
!         set of grid points.
!
! Note that for x being one of the cartesian components of the
! center of the atomic grid we do not have f^x but rather f^(x).
! The correct contribution will be added below.
!
! Here we also accumulate contributions for the rotational
! invariance.
!                                                                      *
!***********************************************************************
!                                                                      *
select case (Functional_type)
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  case (LDA_type)
    !                                                                  *
    !*******************************************************************
    !                                                                  *

    call mma_Allocate(Aux,nD,mGrid,Label='Aux')
    if (nD == 1) then
      do j=1,mGrid
        Aux(1,j) = vRho(1,j)
      end do
    else
      do j=1,mGrid
        Aux(1,j) = vRho(1,j)
        Aux(2,j) = vRho(2,j)
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (GGA_type)
    !                                                                  *
    !*******************************************************************
    !
    call mma_Allocate(Aux,4*nD,mGrid,Label='Aux')
    if (nD == 1) then
      do j=1,mGrid
        Aux(1,j) = vRho(1,j)
        Aux(2,j) = Two*vSigma(1,j)*Gradrho(1,j)
        Aux(3,j) = Two*vSigma(1,j)*Gradrho(2,j)
        Aux(4,j) = Two*vSigma(1,j)*Gradrho(3,j)
      end do
    else
      do j=1,mGrid
        gxa = Gradrho(1,j)
        gya = Gradrho(2,j)
        gza = Gradrho(3,j)
        gxb = Gradrho(4,j)
        gyb = Gradrho(5,j)
        gzb = Gradrho(6,j)

        Aux(1,j) = vRho(1,j)
        Aux(2,j) = vRho(2,j)
        Aux(3,j) = Two*vSigma(1,j)*gxa+vSigma(2,j)*gxb
        Aux(4,j) = Two*vSigma(1,j)*gya+vSigma(2,j)*gyb
        Aux(5,j) = Two*vSigma(1,j)*gza+vSigma(2,j)*gzb
        Aux(6,j) = Two*vSigma(3,j)*gxb+vSigma(2,j)*gxa
        Aux(7,j) = Two*vSigma(3,j)*gyb+vSigma(2,j)*gya
        Aux(8,j) = Two*vSigma(3,j)*gzb+vSigma(2,j)*gza
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type1)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call mma_Allocate(Aux,5*nD,mGrid,Label='Aux')
    if (nD == 1) then
      do j=1,mGrid
        Aux(1,j) = vRho(1,j)
        Aux(2,j) = Two*vSigma(1,j)*Gradrho(1,j)
        Aux(3,j) = Two*vSigma(1,j)*Gradrho(2,j)
        Aux(4,j) = Two*vSigma(1,j)*Gradrho(3,j)
        Aux(5,j) = Quart*vTau(1,j)
      end do
    else
      do j=1,mGrid
        gxa = Gradrho(1,j)
        gya = Gradrho(2,j)
        gza = Gradrho(3,j)
        gxb = Gradrho(4,j)
        gyb = Gradrho(5,j)
        gzb = Gradrho(6,j)

        Aux(1,j) = vRho(1,j)
        Aux(2,j) = vRho(2,j)
        Aux(3,j) = Two*vSigma(1,j)*gxa+vSigma(2,j)*gxb
        Aux(4,j) = Two*vSigma(1,j)*gya+vSigma(2,j)*gyb
        Aux(5,j) = Two*vSigma(1,j)*gza+vSigma(2,j)*gzb
        Aux(6,j) = Two*vSigma(3,j)*gxb+vSigma(2,j)*gxa
        Aux(7,j) = Two*vSigma(3,j)*gyb+vSigma(2,j)*gya
        Aux(8,j) = Two*vSigma(3,j)*gzb+vSigma(2,j)*gza
        Aux(9,j) = Half*vTau(1,j)
        Aux(10,j) = Half*vTau(2,j)
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case (meta_GGA_type2)
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call mma_Allocate(Aux,6*nD,mGrid,Label='Aux')
    if (nD == 1) then
      do j=1,mGrid
        Aux(1,j) = vRho(1,j)
        Aux(2,j) = Two*vSigma(1,j)*Gradrho(1,j)
        Aux(3,j) = Two*vSigma(1,j)*Gradrho(2,j)
        Aux(4,j) = Two*vSigma(1,j)*Gradrho(3,j)
        Aux(5,j) = Quart*vTau(1,j)
        Aux(6,j) = vLapl(1,j)
      end do
    else
      do j=1,mGrid
        gxa = Gradrho(1,j)
        gya = Gradrho(2,j)
        gza = Gradrho(3,j)
        gxb = Gradrho(4,j)
        gyb = Gradrho(5,j)
        gzb = Gradrho(6,j)

        Aux(1,j) = vRho(1,j)
        Aux(2,j) = vRho(2,j)
        Aux(3,j) = Two*vSigma(1,j)*gxa+vSigma(2,j)*gxb
        Aux(4,j) = Two*vSigma(1,j)*gya+vSigma(2,j)*gyb
        Aux(5,j) = Two*vSigma(1,j)*gza+vSigma(2,j)*gzb
        Aux(6,j) = Two*vSigma(3,j)*gxb+vSigma(2,j)*gxa
        Aux(7,j) = Two*vSigma(3,j)*gyb+vSigma(2,j)*gya
        Aux(8,j) = Two*vSigma(3,j)*gzb+vSigma(2,j)*gza
        Aux(9,j) = Half*vTau(1,j)
        Aux(10,j) = Half*vTau(2,j)
        Aux(11,j) = vLapl(1,j)
        Aux(12,j) = vLapl(2,j)
      end do
    end if
    !                                                                  *
    !*******************************************************************
    !                                                                  *
  case default
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    call WarningMessage(2,'Do_Grad: wrong functional type!')
    call Abend()
    !                                                                  *
    !*******************************************************************
    !                                                                  *
end select
!                                                                      *
!***********************************************************************
!                                                                      *
OV(:,:) = Zero
do i_Eff=1,nGrad_Eff
  tmp = Zero
  OVT(:) = Zero
  do j=1,mGrid
    dF_dr = Weights(j)*DDot_(ndRho_dR,Aux(:,j),1,dRho_dR(:,j,i_Eff),1)
    tmp = tmp+dF_dr

    ! Accumulate stuff for rotational invariance

    OVT(:) = OVT(:)+dF_dr*Grid(:,j)
  end do
  ixyz = iTab(1,i_Eff)
  OV(ixyz,:) = OV(ixyz,:)+OVT(:)-tmp*R_Grid(:)
  Temp(i_Eff) = -tmp
end do

call mma_deAllocate(Aux)

do i_Eff=1,nGrad_Eff
  if (iTab(2,i_Eff) == Off) Temp(i_Eff) = Zero
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
call RecPrt('w * f^x before translational contributions',' ',Temp,1,nGrad_Eff)
call RecPrt('OV',' ',OV,3,3)
#endif
!                                                                      *
!***********************************************************************
!***********************************************************************
!                                                                      *
! Here we compute the term
!
! w * f^x
!
! for x being a cartesian component of the center of the atomic
! grid. This is done using the translational invariance condition.

if (Grid_Type == Moving_Grid) then
  do ixyz=1,3
    iGrad = 0
    do jGrad=1,nGrad_Eff
      if ((iTab(1,jGrad) == ixyz) .and. (iTab(2,jGrad) == Off) .and. (IndGrd(jGrad) > 0)) iGrad = jGrad
    end do
#   ifdef _DEBUGPRINT_
    write(u6,*) 'iGrad=',iGrad
#   endif
    if (iGrad /= 0) then

      ! Evaluate indirectly via the translational invariance
      ! the sum of the direct and indirect term

      do jGrad=1,nGrad_Eff
        if ((jGrad /= iGrad) .and. (iTab(1,jGrad) == ixyz)) then

          Temp(iGrad) = Temp(iGrad)-Temp(jGrad)

#         ifdef _DEBUGPRINT_
          write(u6,*) 'jGrad,Temp(jGrad)=',jGrad,Temp(jGrad)
#         endif
        end if
      end do

    end if
  end do
# ifdef _DEBUGPRINT_
  call RecPrt('w * f^x',' ',Temp,1,nGrad_Eff)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! For a "moving" grid add contributions due to the derivative with
  ! respect to the partitioning.
  !
  ! w^x * f

  call DGEMM_('N','N',nGrad_Eff,1,mGrid,One,dW_dR,nGrad_Eff,F_xc,mGrid,One,Temp,nGrad_Eff)
# ifdef _DEBUGPRINT_
  call RecPrt('w * f^x + w^x * f',' ',Temp,1,nGrad_Eff)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Add the rotational invariance term
  !
  ! First transform back to the cartesian coordinates system.

  Fact = real(2-(nD/2),kind=wp)
  call DGEMM_('N','N',3,3,3,Fact,OV,3,Pax,3,Zero,V,3)
# ifdef _DEBUGPRINT_
  call RecPrt('V',' ',V,3,3)
# endif

  do i_Eff=1,nGrad_Eff
    iCar = iTab(1,i_Eff)
    jNQ = iTab(3,i_Eff)

    ! Compute < nabla_r f * r^x > as Tr (O^x V)

    Tmp = DDot_(9,NQ_Data(jNQ)%dOdx(:,:,iCar),1,V,1)*Half
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'iCar,jNQ=',iCar,jNQ
    call RecPrt('dOdx',' ',NQ_Data(jNQ)%dOdx(:,:,iCar),3,3)
    write(u6,*) 'Tmp=',Tmp
#   endif
    Temp(i_Eff) = Temp(i_Eff)-Tmp
  end do

end if !moving grid
#ifdef _DEBUGPRINT_
call RecPrt('Gradient contribution from this block',' ',Temp,1,nGrad_Eff)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Accumulate and symmetry adapt.

do i_Eff=1,nGrad_Eff
  i = IndGrd(i_Eff)
  if (i >= 1) then
    Fact = real(iTab(4,i_Eff),kind=wp)
    Grad(i) = Grad(i)+Fact*Temp(i_Eff)
  end if
end do
#ifdef _DEBUGPRINT_
call RecPrt('Gradient accumulated so far',' ',Grad,1,nGrad)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *

return

end subroutine DFT_Grad
