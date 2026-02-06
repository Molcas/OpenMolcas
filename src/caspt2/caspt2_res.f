************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Yoshio Nishimoto                                 *
************************************************************************
      Subroutine CASPT2_Res(VECROT)

      use caspt2_global, only: real_shift, imag_shift, sigma_p_epsilon
      use caspt2_global, only: jStLag,iVecL,iVecG
      use EQSOLV, only: IRHS, IVECX, IVECR, IVECC, IVECC2, IVECW
      use caspt2_module, only: IFMSCOUP, NSTATE
      use Constants, only: Zero, One, Half, Two
      use definitions, only: wp, iwp, u6

      implicit none

      real(kind=wp), intent(in) :: VECROT(*)

      real(kind=wp) :: SAV, SAVI, savreg, Scal
      integer(kind=iwp) :: iVecXbk, iVecRbk, iRHSbk, iStLag, ICONV

      !! 1) Calculate the derivative of the CASPT2 energy with respect
      !!    to the amplitude.
      !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
      !!    diagonal CASPT2, compute the lambda directly.
      !!
      !! L_S = U_{TS}*H_{TU}*U_{US}
      !!     + U_{SS}*(H_{SS} + <\Psi_S|H0-E0|\Psi_S>)*U_{SS}
      !!     + <lambda|H|\Psi0> + <lambda|H0-E0+Eshift|\Psi_S>

!     write(u6,*) "in CASPT2_res"
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
!     write (*,*) "Ifmscoup = ", ifmscoup, nstate
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
!
!-----------------------------------------------------------------------
!
      !! RHS_SGMDIA
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W1,lg_W2,DIN,DIS)
#ifdef _MOLCAS_MPP_
      USE Para_Info, ONLY: Is_Real_Par
#endif
      use fake_GA, only: GA_Arrays
      use definitions, only: wp, iwp

      implicit none

      integer(kind=iwp), intent(in) :: Mode, NIN, NIS, lg_W1, lg_W2
      real(kind=wp), intent(in) :: DIN(*), DIS(*)

#ifdef _MOLCAS_MPP_
      integer(kind=iwp) :: myRank, iLo1, iHi1, jLo1, jHi1, iLo2, iHi2,
     &                     jLo2, jHi2, NROW, NCOL, mW1, LDW1, mW2, LDW2
#endif

! Apply the resolvent of the diagonal part of H0 to an RHS array

#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
!-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W1,myRank,iLo1,iHi1,jLo1,jHi1)
        CALL GA_Distribution (lg_W2,myRank,iLo2,iHi2,jLo2,jHi2)
        !! Well, assume the same dimension
        IF (iLo1 > 0 .AND. jLo1 > 0 .AND. iLo2 > 0 .AND. jLo2 > 0) THEN
          NROW=iHi1-iLo1+1
          NCOL=jHi1-jLo1+1
          CALL GA_Access (lg_W1,iLo1,iHi1,jLo1,jHi1,mW1,LDW1)
          CALL GA_Access (lg_W2,iLo2,iHi2,jLo2,jHi2,mW2,LDW2)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW1),DBL_MB(mW2),
     &                      LDW1,DIN(iLo1),DIS(jLo1))
          CALL GA_Release_Update (lg_W1,iLo1,iHi1,jLo1,jHi1)
          CALL GA_Release_Update (lg_W2,iLo2,iHi2,jLo2,jHi2)
        END IF
        CALL GA_Sync()
      ELSE
#endif
        CALL CASPT2_ResD2(MODE,NIN,NIS,GA_Arrays(lg_W1)%A,
     &                                 GA_Arrays(lg_W2)%A,
     &                                 NIN,DIN,DIS)
#ifdef _MOLCAS_MPP_
      END IF
#endif

      END SUBROUTINE CASPT2_ResD
!
!-----------------------------------------------------------------------
!
      !! RESDIA
      subroutine CASPT2_ResD2(Mode,nRow,nCol,W1,W2,LDW,dIn,dIs)

      use Constants, only: Zero, One
      use definitions, only: wp, iwp
      use caspt2_global, only: real_shift, imag_shift,
     &                         sigma_p_epsilon, sigma_p_exponent

      implicit none

      integer(kind=iwp), intent(in)    :: Mode, nRow, nCol, LDW
      real(kind=wp),     intent(inout) :: W1(LDW,*), W2(LDW,*)
      real(kind=wp),     intent(in)    :: dIn(*), dIs(*)

      integer(kind=iwp)                :: i, j, p
      real(kind=wp)                    :: scal, delta, delta_inv,
     &                                    sigma, epsilon, expscal,
     &                                    delta_ps

      epsilon = Zero
      if (sigma_p_epsilon /= Zero) then
        epsilon = sigma_p_epsilon
      end if

      do j = 1, nCol
        do i = 1, nRow
          scal = Zero
          select case (mode)
          case (1)
            ! energy denominator plus real shift
            delta = dIn(i) + dIs(j) + real_shift
            ! inverse denominator plus imaginary shift
            delta_inv = delta/(delta**2 + imag_shift**2)
            ! multiply by (inverse) sigma-p regularizer
            epsilon = sigma_p_epsilon
            p = sigma_p_exponent
            if (epsilon > Zero) then
              sigma = One/epsilon**p
              delta_inv
     &          = delta_inv * (One - exp(-sigma*abs(delta)**p))
            end if
            !! The following SCAL is the actual residual
            scal = One - (dIn(i)+dIs(j))*delta_inv
            !! Another scaling is required for lambda
            scal = -scal*delta_inv

            W1(i,j) = scal*W1(i,j)
            W2(i,j) = scal*W2(i,j)
          case (2)
            if (imag_shift /= Zero) then
              scal = imag_shift/(dIn(i)+dIs(j))
              W1(i,j) = scal*W1(i,j)
              W2(i,j) = scal*W2(i,j)
            else if (epsilon /= Zero) then
              ! derivative of the denominator of sigma-p CASPT2
              ! always real_shift = imag_shift = 0
              delta   = dIn(i)+dIs(j)
              ! multiply by (inverse) sigma-p regularizer
              p = sigma_p_exponent
              sigma   = One/epsilon**p
              delta_ps= (delta**p)*sigma
              expscal = exp(-abs(delta_ps))
              delta_inv = One/(One - expscal)

              W1(i,j) = delta_inv*p*(delta_ps)*expscal*W1(i,j)
     &                  *(SIGN(One,delta)**p)
              W2(i,j) = delta_inv*W2(i,j)
!             W2(i,j) = delta_inv*p*(delta_ps)*expscal*W2(i,j)
!    *                  *(SIGN(One,delta)**p)
!             W1(i,j) = delta_inv*W1(i,j)
            end if
          case (3)
            ! derivative of the numerator of sigma-p CASPT2
            ! always real_shift = imag_shift = 0
            delta     = dIn(i)+dIs(j)
            ! multiply by (inverse) sigma-p regularizer
            p = sigma_p_exponent
            sigma = One/epsilon**p
            expscal = exp(-sigma*abs(delta)**p)
            delta_inv = One/(One - expscal)
            W1(i,j) = delta_inv*W1(i,j)
          end select
        end do
      end do

      end subroutine CASPT2_ResD2
!
!-----------------------------------------------------------------------
!
      SUBROUTINE PCG_RES(ICONV)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: terse, usual
      use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecX
      use caspt2_module, only: MxCase, MaxIt, rNorm, ThrConv
      use Constants, only: Zero, One
      use definitions, only: wp, iwp, u6

      IMPLICIT NONE

      integer(kind=iwp), intent(out) :: ICONV

      integer(kind=iwp) :: I, ITER, IVECP, IVECT, IVECU
      real(kind=wp) :: ALPHA,BETA,PR,PT,UR,ECORR(0:8,0:MXCASE),EAIVX,
     &                 EATVX,EBJAI,EBJAT,EBVAT,EVJAI,EVJTI,EVJTU,E2NONV,
     &                 OVLAPS(0:8,0:MXCASE),DSCALE
      logical(kind=iwp) :: converged

! Flag to tell wether convergence was obtained
      ICONV = 0
      converged = .false.


! Mnemonic names for vectors stored on LUSOLV, see EQCTL.
! Here, we use the local names IVECP, IVECT, IVECU which are thus
! to be seen as overlayed areas. The true vectors IVECC and IVECC2
! are computed on return from this routine, so for a while we use them
! for temporaries.
      IVECP=IVECC
      IVECT=IVECC2
      IVECU=IVECC2


      ITER=0
      RNORM=Zero

! Solve equations for the diagonal case, in eigenbasis:
! Current solution vector X, Current residual vector R
      CALL PSCAVEC(-One,IRHS,IVECR)
      CALL PRESDIA(IVECR,IVECX,OVLAPS)
      IF(MAXIT == 0) THEN
        IF(IPRGLB >= TERSE) THEN
          WRITE(u6,*)
          WRITE(u6,'(23A5)')('-----',i=1,23)
          WRITE(u6,*)' DIAGONAL CASPT2 APPROXIMATION:'
        END IF
        converged = .true.
      else

! Pre-conditioned conjugate gradient:
! R <- R - (H0-E0)*X
        CALL SIGMA_CASPT2(-One,One,IVECX,IVECR)
        CALL POVLVEC(IVECR,IVECR,OVLAPS)
        RNORM=SQRT(OVLAPS(0,0))
        IF(RNORM < THRCONV) then
          converged = .true.
        ELSE

          IF(IPRGLB >= USUAL) THEN
            WRITE(u6,*)
        WRITE(u6,*) "Solving the Lambda equation for analytic gradients"
            WRITE(u6,'(25A5)')('-----',I=1,25)
            WRITE(u6,'(2X,A,A)')
     & 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI ',
     & '       BVAT        BJAT        BJAI        TOTAL       RNORM  '
            WRITE(u6,'(25A5)')('-----',I=1,25)
          END IF
          CALL PRESDIA(IVECR,IVECP,OVLAPS)
! PCG iteration loops:
!---------------------
          do
            CALL POVLVEC(IVECP,IVECP,OVLAPS)
            DSCALE=One/SQRT(OVLAPS(0,0))
            CALL PSCAVEC(DSCALE,IVECP,IVECP)
            CALL POVLVEC(IVECP,IVECR,OVLAPS)
            PR=OVLAPS(0,0)
            CALL SIGMA_CASPT2(One,Zero,IVECP,IVECT)
            CALL POVLVEC(IVECP,IVECT,OVLAPS)
            PT=OVLAPS(0,0)
            ALPHA=PR/PT
            CALL PLCVEC(ALPHA,One,IVECP,IVECX)
            CALL PLCVEC(-ALPHA,One,IVECT,IVECR)
            CALL POVLVEC(IVECR,IVECR,OVLAPS)
            RNORM=SQRT(OVLAPS(0,0))
            IF(RNORM < THRCONV) then
              converged = .true.
              exit
            end if
            ITER=ITER+1
            CALL POVLVEC(IRHS,IVECX,ECORR)
            EVJTU=ECORR(0,1)
            EVJTI=ECORR(0,2)+ECORR(0,3)
            EATVX=ECORR(0,4)
            EAIVX=ECORR(0,5)
            EVJAI=ECORR(0,6)+ECORR(0,7)
            EBVAT=ECORR(0,8)+ECORR(0,9)
            EBJAT=ECORR(0,10)+ECORR(0,11)
            EBJAI=ECORR(0,12)+ECORR(0,13)
            E2NONV=ECORR(0,0)
            IF(IPRGLB >= USUAL) THEN
            WRITE(u6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,
     &                            EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM
              CALL XFLUSH(u6)
            END IF
            IF(ITER >= MAXIT) exit
            CALL PRESDIA(IVECR,IVECU,OVLAPS)
            UR=OVLAPS(0,0)
            BETA=PR/UR
            CALL PLCVEC(BETA,One,IVECU,IVECP)
          end do
        end if
!---------------------

        if (.not.converged) then
          IF(IPRGLB >= TERSE) THEN
            WRITE(u6,*)
            WRITE(u6,*)' NOT CONVERGED AFTER MAX ITERATIONS.'
          END IF
          ICONV = 16
        end if
      end if

      IF(IPRGLB >= TERSE) THEN
       WRITE(u6,'(25A5)')('-----',I=1,25)
       WRITE(u6,*)
      END IF

      END SUBROUTINE PCG_RES
