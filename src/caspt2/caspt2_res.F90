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
! Copyright (C) 2021, Yoshio Nishimoto                                 *
!***********************************************************************

subroutine CASPT2_Res(VECROT,nVECROT)

use caspt2_global, only: imag_shift, iVecL, jStLag, real_shift, sigma_p_epsilon
use caspt2_module, only: HZERO, IFMSCOUP, NSTATE
use EQSOLV, only: IRHS, IVECC, IVECC2, IVECR, IVECW, IVECX
use SC_NEVPT2, only: SC_prop, SC_NEVPT2_res
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nVECROT
real(kind=wp), intent(in) :: VECROT(nVECROT)
integer(kind=iwp) :: ICONV, iRHSbk, iStLag, iVecRbk, iVecXbk
real(kind=wp) :: SAV, SAVI, savreg, Scal
integer(kind=iwp), parameter :: iVecG = 8

!! Some subroutines (e.g., PSCAVEC) use NINDEP, so they cannot be
!! directly used for SC-NEVPT2. It should be easier to write from scratch
if (SC_prop) then
  call SC_NEVPT2_res(VECROT)
  return
end if

if (HZERO == 'DYALL' .and. real_shift == Zero .and. &
    imag_shift == Zero .and. sigma_p_epsilon == Zero) then
  !! The correlation energy of the target state for the plane NEVPT2 is stationary
  !! Note that IVECL for SC-NEVPT2 is in an MO (contravariant) basis
  call RHS_ZERO(iVecL)
else
  !! 1) Calculate the derivative of the CASPT2 energy with respect
  !!    to the amplitude.
  !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
  !!    diagonal CASPT2, compute the lambda directly.
  !!
  !! L_S = U_{TS}*H_{TU}*U_{US}
  !!     + U_{SS}*(H_{SS} + <\Psi_S|H0-E0|\Psi_S>)*U_{SS}
  !!     + <lambda|H|\Psi0> + <lambda|H0-E0+Eshift|\Psi_S>

  !write(u6,*) 'in CASPT2_res'
  call PSCAVEC(One,IRHS,iVecL)

  !! Construct the partial derivative of the target state
  !! The derivative is constructed in iVecL
  !! The shift parameters are set to zero, because the actual energy
  !! is computed without them. The reference state has to be
  !! multiplied by two, from the above equation for L_S.
  !! For MS-CASPT2, the rotation is mutiplied later.
  SAV = real_shift
  SAVI = imag_shift
  savreg = sigma_p_epsilon
  real_shift = Zero
  imag_shift = Zero
  sigma_p_epsilon = Zero
  call SIGMA_CASPT2(Two,Two,IVECX,iVecL)
  real_shift = SAV
  imag_shift = SAVI
  sigma_p_epsilon = savreg
end if

!! Add the partial derivative contribution for MS-CASPT2
!! (off-diagonal elements). The derivative is taken with IVECW
!! and put in IVECC.
!write(u6,*) 'Ifmscoup = ',ifmscoup,nstate
if (IFMSCOUP) then
  call RHS_ZERO(IVECC)
  call PSCAVEC(VECROT(jStLag),iVecL,iVecL)
  do iStLag=1,nState
    Scal = VECROT(iStLag)
    if (iStLag == jStLag) Scal = Zero
    if (abs(VECROT(iStLag)) <= 1.0e-12_wp) cycle
    call MS_Res(1,iStLag,jStLag,Scal)
  end do
  !! Transform to SR representatin (IRHS).
  call PTRTOSR(0,IVECC,IRHS)
  !! Add to iVecL
  call PLCVEC(One,One,IRHS,iVecL)
end if

!! Finally, solve the lambda equation.
!! The following is just a copy-and-paste of eqctl2 and pcg,
!! but some unnecessary lines (comuptation of energy etc.) are
!! omitted.
!! The lambda equation is solved with the shift parameters,
!! as is the case for the T-amplitude.

iVecXbk = iVecX
iVecRbk = iVecR
iRHSbk = iRHS
iVecX = iVecR
iRHS = iVecL !! = 7
iVecR = iVecG !! = 8

call PCG_RES(ICONV)
if (ICONV /= 0) then
  write(u6,*) 'Lambda equation did not converge...'
  write(u6,*) 'Continue anyway?'
end if

iVecX = iVecXbk
iVecR = iVecRbk
iRHS = iRHSbk

!! For implicit derivative of S
if (IFMSCOUP) then
  call PTRTOSR(1,IVECW,IRHS)
  call RHS_ZERO(IVECC)
  do iStLag=1,nState
    Scal = VECROT(iStLag)*Half
    if (iStLag == jStLag) Scal = Scal*Two
    if (abs(VECROT(iStLag)) <= 1.0e-12_wp) cycle
    call MS_Res(1,iStLag,jStLag,Scal)
  end do
  call PTRTOSR(0,IVECC,iVecL)
end if

!! Restore contravariant and covariant representations of the
!! non-variational T-amplitude
call PTRTOC(0,IVECX,IVECC)
call PTRTOC(1,IVECX,IVECC2)

return

end subroutine CASPT2_Res
