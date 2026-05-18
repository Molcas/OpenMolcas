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
      SUBROUTINE PCG_RES(ICONV)
      use caspt2_global, only:iPrGlb
      use PrintLevel, only: TERSE, USUAL
      use EQSOLV, only: iRHS, iVecc, iVecc2, iVecR, iVecX
      use caspt2_module, only: MxCase, MaxIt, rNorm, ThrConv
      use Constants, only: Zero, One
      use definitions, only: wp, iwp, u6

      IMPLICIT NONE

      integer(kind=iwp), intent(out) :: ICONV

      integer(kind=iwp) :: I, ITER, IVECP, IVECT, IVECU
      real(kind=wp) :: ALPHA,BETA,PR,PT,UR,ECORR(0:8,0:MXCASE),EAIVX,   &
     &                 EATVX,EBJAI,EBJAT,EBVAT,EVJAI,EVJTI,EVJTU,E2NONV,&
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
        WRITE(u6,*) 'Solving the Lambda equation for analytic gradients'
            WRITE(u6,'(25A5)')('-----',I=1,25)
            WRITE(u6,'(2X,A,A)')                                        &
     & 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI ',&
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
            WRITE(u6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,&
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
