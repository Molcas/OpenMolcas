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
      Subroutine CASPT2_Res(VECROT,nVECROT)

      use caspt2_global, only: real_shift, imag_shift, sigma_p_epsilon
      use caspt2_global, only: jStLag,iVecL
      use EQSOLV, only: IRHS, IVECX, IVECR, IVECC, IVECC2, IVECW
      use caspt2_module, only: IFMSCOUP, NSTATE
      use Constants, only: Zero, One, Half, Two
      use definitions, only: wp, iwp, u6

      implicit none

      integer(kind=iwp), intent(in) :: nVECROT
      real(kind=wp), intent(in) :: VECROT(nVECROT)

      real(kind=wp) :: SAV, SAVI, savreg, Scal
      integer(kind=iwp) :: iVecXbk, iVecRbk, iRHSbk, iStLag, ICONV
      integer(kind=iwp), parameter :: iVecG = 8

      !! 1) Calculate the derivative of the CASPT2 energy with respect
      !!    to the amplitude.
      !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
      !!    diagonal CASPT2, compute the lambda directly.
      !!
      !! L_S = U_{TS}*H_{TU}*U_{US}
      !!     + U_{SS}*(H_{SS} + <\Psi_S|H0-E0|\Psi_S>)*U_{SS}
      !!     + <lambda|H|\Psi0> + <lambda|H0-E0+Eshift|\Psi_S>

!     write(u6,*) 'in CASPT2_res'
      CALL PSCAVEC(One,IRHS,iVecL)

      !! Construct the partial derivative of the target state
      !! The derivative is constructed in iVecL
      !! The shift parameters are set to zero, because the actual energy
      !! is computed without them. The reference state has to be
      !! multiplied by two, from the above equation for L_S.
      !! For MS-CASPT2, the rotation is mutiplied later.
      SAV=real_shift
      SAVI=imag_shift
      savreg=sigma_p_epsilon
      real_shift=Zero
      imag_shift=Zero
      sigma_p_epsilon=Zero
      CALL SIGMA_CASPT2(Two,Two,IVECX,iVecL)
      real_shift=SAV
      imag_shift=SAVI
      sigma_p_epsilon=savreg

      !! Add the partial derivative contribution for MS-CASPT2
      !! (off-diagonal elements). The derivative is taken with IVECW
      !! and put in IVECC.
!     write (u6,*) 'Ifmscoup = ', ifmscoup, nstate
      IF (IFMSCOUP) Then
        Call RHS_ZERO(IVECC)
        Call PSCAVEC(VECROT(jStLag),iVecL,iVecL)
        Do iStLag = 1, nState
          Scal = VECROT(iStLag)
          If (iStLag == jStLag) Scal = Zero
          If (ABS(VECROT(iStLag)) <= 1.0e-12_wp) Cycle
          Call MS_Res(1,iStLag,jStLag,Scal)
        End Do
        !! Transform to SR representatin (IRHS).
        CALL PTRTOSR(0,IVECC,IRHS)
        !! Add to iVecL
        Call PLCVEC(One,One,IRHS,iVecL)
      End If

      !! Finally, solve the lambda equation.
      !! The following is just a copy-and-paste of eqctl2.f and pcg.f,
      !! but some unnecessary lines (comuptation of energy etc.) are
      !! omitted.
      !! The lambda equation is solved with the shift parameters,
      !! as is the case for the T-amplitude.

      iVecXbk = iVecX
      iVecRbk = iVecR
      iRHSbk  = iRHS
      iVecX   = iVecR
      iRHS    = iVecL !! = 7
      iVecR   = iVecG !! = 8

      Call PCG_RES(ICONV)
      IF (ICONV /= 0) THEN
        WRITE (u6,'(" Lambda equation did not converge...")')
        WRITE (u6,'(" Continue anyway?")')
      END IF

      iVecX   = iVecXbk
      iVecR   = iVecRbk
      iRHS    = iRHSbk

      !! For implicit derivative of S
      IF (IFMSCOUP) THEN
        CALL PTRTOSR(1,IVECW,IRHS)
        Call RHS_ZERO(IVECC)
        Do iStLag = 1, nState
          Scal = VECROT(iStLag)*Half
          If (iStLag == jStLag) Scal = Scal*Two
          If (ABS(VECROT(iStLag)) <= 1.0e-12_wp) Cycle
          Call MS_Res(1,iStLag,jStLag,Scal)
        End Do
        CALL PTRTOSR(0,IVECC,iVecL)
      END IF

      !! Restore contravariant and covariant representations of the
      !! non-variational T-amplitude
      CALL PTRTOC(0,IVECX,IVECC)
      CALL PTRTOC(1,IVECX,IVECC2)

      RETURN

      end subroutine CASPT2_Res
