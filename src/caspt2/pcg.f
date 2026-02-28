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
* Copyright (C) 1996,1999, Per Ake Malmqvist                           *
************************************************************************
*--------------------------------------------*
* 1996  PER-AAKE MALMQUIST                   *
* DEPARTMENT OF THEORETICAL CHEMISTRY        *
* UNIVERSITY OF LUND                         *
* SWEDEN                                     *
* 1999: GEMINAL-R12 ENABLED                  *
*--------------------------------------------*
      SUBROUTINE PCG(ICONV)
      USE INPUTDATA, ONLY: INPUT
      use caspt2_global, only: EMP2
      use caspt2_global, only: iPrGlb
      use caspt2_global, only: sigma_p_epsilon,imag_shift,real_shift
      use caspt2_global, only: do_grad, nStpGrd
      use caspt2_global, only: LISTS
      use PrintLevel, only: TERSE, USUAL
      use stdalloc, only: mma_allocate, mma_deallocate
      use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecX, NLSTOT
      use caspt2_module, only: DeNorm, E2Corr, E2Tot, MxCase, ERef,
     &                         IfChol, MaxIt, nSym, rNorm,
     &                         ThrConv, Cases
      use constants, only: Zero, One, Two
      use definitions, only: iwp, wp, u6
      IMPLICIT NONE

      integer(kind=iwp), intent(out):: ICONV

      integer(kind=iwp) I,IC,IS,ITER
      integer(kind=iwp) IVECP,IVECT,IVECU
      integer(kind=iwp) LAXITY
      integer(kind=iwp), EXTERNAL :: Cho_X_GetTol
      real(kind=wp) ALPHA,BETA,PR,PT,UR
      real(kind=wp) ECORR(0:8,0:MXCASE)
      real(kind=wp) EAIVX,EATVX,EBJAI,EBJAT,EBVAT,EVJAI,EVJTI,EVJTU
      real(kind=wp) E2NONV,ESHIFT
      real(kind=wp) OVLAPS(0:8,0:MXCASE)
      real(kind=wp) SAV,SAVI,savreg,DSCALE,REFWGT

C Flag to tell whether convergence was obtained
      ICONV = 0

C Lists of coupling coefficients, used for sigma vector
C generation from non-diagonal blocks of H0.
      CALL mma_allocate(LISTS,NLSTOT,LABEL='LISTS')
      CALL MKLIST(LISTS,NLSTOT)


C Mnemonic names for vectors stored on LUSOLV, see EQCTL.
C Here, we use the local names IVECP, IVECT, IVECU which are thus
C to be seen as overlayed areas. The true vectors IVECC and IVECC2
C are computed on return from this routine, so for a while we use them
C for temporaries.
      IVECP=IVECC
      IVECT=IVECC2
      IVECU=IVECC2

C Solve equations for the diagonal case, in eigenbasis:
C Current solution vector X, Current residual vector R
      CALL PSCAVEC(-One,IRHS,IVECR)
      CALL PRESDIA(IVECR,IVECX,OVLAPS)

      IF(MAXIT.EQ.0) THEN
       IF(IPRGLB.GE.TERSE) THEN
        WRITE(u6,*)
        WRITE(u6,'(23A5)')('-----',i=1,23)
        WRITE(u6,*)' DIAGONAL CASPT2 APPROXIMATION:'
        GOTO 900
       END IF
      END IF

      RNORM=Zero
C Pre-conditioned conjugate gradient:
C R <- R - (H0-E0)*X
      CALL SIGMA_CASPT2(-One,One,IVECX,IVECR)
      CALL POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM=SQRT(OVLAPS(0,0))

      IF(RNORM>=THRCONV) THEN

      IF(IPRGLB.GE.USUAL) THEN
       WRITE(u6,*)
       WRITE(u6,*)'The contributions to the second order'//
     &     ' correlation energy in atomic units.'
       WRITE(u6,'(25A5)')('-----',I=1,25)
       WRITE(u6,'(2X,A,A)')
     & 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI ',
     & '       BVAT        BJAT        BJAI        TOTAL       RNORM  '
       WRITE(u6,'(25A5)')('-----',I=1,25)
      END IF

      CALL PRESDIA(IVECR,IVECP,OVLAPS)
C PCG iteration loops:
C---------------------
      ITER=0
      Do
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

         IF(RNORM.LT.THRCONV) EXIT

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
         IF(IPRGLB.GE.USUAL) THEN
          WRITE(u6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,
     &                        EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM
         END IF

         IF(ITER>=MAXIT) THEN
            IF(IPRGLB.GE.TERSE) THEN
             WRITE(6,*)
             WRITE(6,*)' NOT CONVERGED AFTER MAX ITERATIONS.'
            END IF
            ICONV = 16
            EXIT
         END IF

         CALL PRESDIA(IVECR,IVECU,OVLAPS)
         UR=OVLAPS(0,0)
         BETA=PR/UR
         CALL PLCVEC(BETA,One,IVECU,IVECP)
      END Do
C---------------------
      END IF

 900  CONTINUE

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
      CALL POVLVEC(IVECX,IVECX,OVLAPS)
      DENORM=One+OVLAPS(0,0)
      REFWGT=One/DENORM
CPAM Insert: Compute the variational second-order energy.
CPAM Use unshifted H0. Save any shifts, then restore them.
      SAV=real_shift
      SAVI=imag_shift
      savreg=sigma_p_epsilon
      real_shift=Zero
      imag_shift=Zero
      sigma_p_epsilon=Zero
      CALL SIGMA_CASPT2(One,Zero,IVECX,IVECT)
      real_shift=SAV
      imag_shift=SAVI
      sigma_p_epsilon=savreg
      CALL POVLVEC(IVECX,IVECT,OVLAPS)
      E2CORR=Two*E2NONV+OVLAPS(0,0)
CPAM End of insert.
      ESHIFT=E2CORR-E2NONV
      E2TOT=EREF+E2CORR

      IF(IPRGLB.GT.USUAL) THEN
        WRITE(u6,*)
        WRITE(u6,*)' Correlation energy /Case, /Symm, and sums:'
        DO IC=1,13
         WRITE(u6,'(1X,A8,9F12.8)')
     &      CASES(IC),(ECORR(IS,IC),IS=1,NSYM),ECORR(0,IC)
        END DO
        WRITE(u6,'(1X,A8,9F12.8)')
     &    'Summed: ', (ECORR(IS,0),IS=1,NSYM),ECORR(0,0)
      ENDIF

      if (nStpGrd == 1 .or. (nStpGrd == 2 .and. .not.do_grad)) then
      IF (IPRGLB.GE.TERSE) THEN
         WRITE(u6,'(25A5)')('-----',I=1,25)
         WRITE(u6,*)
         WRITE(u6,*)' FINAL CASPT2 RESULT:'
         WRITE(u6,*)

         If (.not.Input % LovCASPT2) Then
            WRITE(u6,'(6x,a,f18.10)')'Reference energy:     ',EREF
            WRITE(u6,'(6x,a,f18.10)')'E2 (Non-variational): ',E2NONV
            IF(real_shift.NE.0.0d0.or.imag_shift.ne.0.0d0
     &       .or.sigma_p_epsilon.ne.0.0d0) THEN
              WRITE(u6,'(6x,a,f18.10)')'Shift correction:     ',ESHIFT
            END IF
            WRITE(u6,'(6x,a,f18.10)')'E2 (Variational):     ',E2CORR
            If (.not.Input % FnoCASPT2) Then
               WRITE(u6,'(6x,a,f18.10)')'Total energy:         ',E2TOT
            Else
               WRITE(u6,'(6x,a,f18.10,a)')'FNO correction:       ',EMP2,
     &              '   (estimate)   '
               WRITE(u6,'(6x,a,f13.5)')
               E2TOT=E2TOT+EMP2
               WRITE(u6,'(6x,a,f18.10,a)')'Total energy:         ',
     &                   E2TOT, '   (FNO-CASPT2) '
            EndIf
            WRITE(u6,'(6x,a,f18.10)')'Residual norm:        ',RNORM
            WRITE(u6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT
         Else
            WRITE(u6,'(6x,a,f18.10)')
     &              'Reference energy:                 ',EREF
            WRITE(u6,'(6x,a,f18.10)')
     &              'Active-Site E2 (Non-variational): ',E2NONV
            IF(real_shift.NE.0.0d0.or.imag_shift.ne.0.0d0) THEN
              WRITE(u6,'(6x,a,f18.10)')
     &              'Shift correction:                 ',ESHIFT
            END IF
            WRITE(u6,'(6x,a,f18.10)')
     &              'Active-Site E2 (Variational):     ',E2CORR
            WRITE(u6,'(6x,a,f18.10)')
     &              'Frozen region E2 :                ',EMP2
            WRITE(u6,'(6x,a,f18.10)')
     &              'Residual norm:                    ',RNORM
            WRITE(u6,'(6x,a,f13.5)')
     &              'Reference weight:                 ',REFWGT
            WRITE(u6,'(6x,a,f13.5)')
            E2TOT=E2TOT+EMP2
            WRITE(u6,'(6x,a,f18.10)')
     &              'Total energy (LovCASPT2):         ',E2TOT
         EndIf
      END IF
      end if

* In automatic verification calculations, the precision is lower
* in case of Cholesky calculation.
      LAXITY=8
      IF(IfChol) LAXITY=Cho_X_GetTol(LAXITY)
      Call Add_Info('E_CASPT2',[E2TOT],1,LAXITY)

      if (nStpGrd == 1 .or. (nStpGrd == 2 .and. .not.do_grad)) then
      IF(IPRGLB.GE.USUAL) THEN
       WRITE(u6,*)
       WRITE(u6,'(6x,a)')
     &  'Contributions to the CASPT2 correlation energy'
       WRITE(u6,'(6x,a,F18.10)')
     &  'Active & Virtual Only:    ',EATVX+EBVAT
       WRITE(u6,'(6x,a,F18.10)')
     &  'One Inactive Excited:     ',EVJTU+EAIVX+EBJAT
       WRITE(u6,'(6x,a,F18.10)')
     &  'Two Inactive Excited:     ',EVJTI+EVJAI+EBJAI
       WRITE(u6,*)
      END IF
      end if
      CALL mma_deallocate(LISTS)

      END SUBROUTINE PCG
