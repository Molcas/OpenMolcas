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

subroutine NDSD_Ts(mGrid,nDmat)
!***********************************************************************
!                                                                      *
! Object:  compute Func for Thomas-Fermi KE functional                 *
!          compute non-TF part (rho_B dependent) of NDSD potential     *
!                                                                      *
!          (see J.-M. Garcia Lastra, J. W. Kaminski, T. A. Wesolowski, *
!                                   J. Chem. Phys.  129 (2008) 074107.)*
!                                                                      *
!          Note: for a spin-polarized rho_B (environment density), the *
!                NDSD potential is computed using the alpha+beta       *
!                density, gradient and laplacian.                      *
!                                                                      *
!***********************************************************************

use nq_Grid, only: F_xc, GradRho, Lapl, Rho, vRho
use Constants, only: Zero, One, Two, Three, Five, Ten, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: mGrid, nDmat
integer(kind=iwp) :: iGrid, k
real(kind=wp) :: Cf, d_sys, da_sys, db_sys, dfunc_NDSD, dfunc_NDSD_alpha, dfunc_NDSD_beta, DTot, functional, Rho_min, &
                 wGradRho(1:3), wLaplRho
real(kind=wp), parameter :: Coeff = One, Five3 = Five/Three, T_X = 1.0e-20_wp, Two3 = Two/Three
real(kind=wp), external :: Fexp, Vt_lim

!                                                                      *
!***********************************************************************
!                                                                      *
vRho(:,:) = Zero
Cf = (Three/Ten)*(three*Pi**Two)**Two3
Rho_min = T_X*1.0e-2_wp
!                                                                      *
!***********************************************************************
!                                                                      *
! Compute value of energy and integrand on the grid
!                                                                      *
!***********************************************************************
!                                                                      *
if (nDmat == 1) then
  do iGrid=1,mGrid
    d_sys = Two*Rho(1,iGrid)
    if (d_sys < T_X) cycle

    ! Kinetic energy contributions

    functional = Cf*d_sys**Five3
    F_xc(iGrid) = F_xc(iGrid)+Coeff*functional

    ! Contributions to the potential

    do k=1,3
      wGradRho(k) = Two*GradRho(k,iGrid)
    end do
    wLaplRho = Two*Lapl(1,iGrid)

    dfunc_NDSD = Fexp(d_sys,wGradRho(1))*Vt_lim(d_sys,wGradRho(1),wLaplRho)
    vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_NDSD

  end do

else if (nDmat == 2) then

  Cf = Cf*(Two**Two3)

  do iGrid=1,mGrid
    da_sys = max(Rho_Min,Rho(1,iGrid))
    db_sys = max(Rho_Min,Rho(2,iGrid))
    DTot = da_sys+db_sys
    if (DTot < T_X) cycle

    ! Kinetic energy contributions

    functional = Cf*(da_sys**Five3+db_sys**Five3)
    F_xc(iGrid) = F_xc(iGrid)+Coeff*functional

    ! Contributions to the potential

    do k=1,3
      wGradRho(k) = Rho(k,iGrid)+Rho(k+3,iGrid)
    end do
    wLaplRho = Lapl(1,iGrid)+Lapl(2,iGrid)

    dfunc_NDSD_alpha = Fexp(DTot,wGradRho(1))*Vt_lim(DTot,wGradRho(1),wLaplRho)
    dfunc_NDSD_beta = dfunc_NDSD_alpha

    vRho(1,iGrid) = vRho(1,iGrid)+Coeff*dfunc_NDSD_alpha
    vRho(2,iGrid) = vRho(2,iGrid)+Coeff*dfunc_NDSD_beta

  end do

else
  write(u6,*) 'In NDSD_Ts: invalid # of densities. nDmat=  ',nDmat
  call Abend()
end if

return

end subroutine NDSD_Ts
