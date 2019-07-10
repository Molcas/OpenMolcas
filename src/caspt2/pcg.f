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
      USE INPUTDATA
      IMPLICIT NONE

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

      INTEGER ICONV

      INTEGER I,IC,IS,ITER
      INTEGER IVECP,IVECT,IVECU
      INTEGER LAXITY
      INTEGER Cho_X_GetTol
      EXTERNAL Cho_X_GetTol
      REAL*8 ALPHA,BETA,PR,PT,UR
      REAL*8 ECORR(0:8,0:MXCASE)
      REAL*8 EAIVX,EATVX,EBJAI,EBJAT,EBVAT,EVJAI,EVJTI,EVJTU
      REAL*8 E2NONV,ESHIFT
      REAL*8 OVLAPS(0:8,0:MXCASE)
      REAL*8 SAV,SAVI,DSCALE

      CALL QENTER('PCG')
C Flag to tell wether convergence was obtained
      ICONV = 0

C Lists of coupling coefficients, used for sigma vector
C generation from non-diagonal blocks of H0.
      CALL GETMEM('LISTS','ALLO','INTE',LLISTS,NLSTOT)
      CALL MKLIST(iWORK(LLISTS))


C Mnemonic names for vectors stored on LUSOLV, see EQCTL.
C Here, we use the local names IVECP, IVECT, IVECU which are thus
C to be seen as overlayed areas. The true vectors IVECC and IVECC2
C are computed on return from this routine, so for a while we use them
C for temporaries.
      IVECP=IVECC
      IVECT=IVECC2
      IVECU=IVECC2


      ITER=0
      RNORM=0.0d0

C Solve equations for the diagonal case, in eigenbasis:
C Current solution vector X, Current residual vector R
      CALL PSCAVEC(-1.0D00,IRHS,IVECR)
      CALL PRESDIA(IVECR,IVECX,OVLAPS)
      IF(MAXIT.EQ.0) THEN
       IF(IPRGLB.GE.TERSE) THEN
        WRITE(6,*)
        WRITE(6,'(23A5)')('-----',i=1,23)
        WRITE(6,*)' DIAGONAL CASPT2 APPROXIMATION:'
        GOTO 900
       END IF
      END IF

C Pre-conditioned conjugate gradient:
C R <- R - (H0-E0)*X
      CALL SIGMA_CASPT2(-1.0D00,1.0D00,IVECX,IVECR)
      CALL POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM=SQRT(OVLAPS(0,0))
      IF(RNORM.LT.THRCONV) GOTO 900
      IF(IPRGLB.GE.USUAL) THEN
       WRITE(6,*)
       WRITE(6,*)'The contributions to the second order'//
     &     ' correlation energy in atomic units.'
       WRITE(6,'(25A5)')('-----',I=1,25)
       WRITE(6,'(2X,A,A)')
     & 'IT.      VJTU        VJTI        ATVX        AIVX        VJAI ',
     & '       BVAT        BJAT        BJAI        TOTAL       RNORM  '
       WRITE(6,'(25A5)')('-----',I=1,25)
      END IF
      CALL PRESDIA(IVECR,IVECP,OVLAPS)
C PCG iteration loops:
C---------------------
 100  CONTINUE
      CALL POVLVEC(IVECP,IVECP,OVLAPS)
      DSCALE=1.0D00/SQRT(OVLAPS(0,0))
      CALL PSCAVEC(DSCALE,IVECP,IVECP)
      CALL POVLVEC(IVECP,IVECR,OVLAPS)
      PR=OVLAPS(0,0)
      CALL SIGMA_CASPT2(1.0D00,0.0D00,IVECP,IVECT)
      CALL POVLVEC(IVECP,IVECT,OVLAPS)
      PT=OVLAPS(0,0)
      ALPHA=PR/PT
      CALL PLCVEC(ALPHA,1.0D00,IVECP,IVECX)
      CALL PLCVEC(-ALPHA,1.0D00,IVECT,IVECR)
      CALL POVLVEC(IVECR,IVECR,OVLAPS)
      RNORM=SQRT(OVLAPS(0,0))
      IF(RNORM.LT.THRCONV) GOTO 900
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
       WRITE(6,'(1X,I3,1X,10F12.6)') ITER,EVJTU,EVJTI,EATVX,EAIVX,
     &                     EVJAI,EBVAT,EBJAT,EBJAI,E2NONV,RNORM
       CALL XFLUSH(6)
      END IF
      IF(ITER.GE.MAXIT) GOTO 800
      CALL PRESDIA(IVECR,IVECU,OVLAPS)
      UR=OVLAPS(0,0)
      BETA=PR/UR
      CALL PLCVEC(BETA,1.0D00,IVECU,IVECP)
      GOTO 100
C---------------------

 800  CONTINUE
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,*)
       WRITE(6,*)' NOT CONVERGED AFTER MAX ITERATIONS.'
      END IF
      ICONV = 16
 900  CONTINUE
      IF(IPRGLB.GE.TERSE) THEN
       WRITE(6,'(23A5)')('-----',I=1,23)
       WRITE(6,*)
       WRITE(6,*)' FINAL CASPT2 RESULT:'
      END IF
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
      DENORM=1.0D0+OVLAPS(0,0)
      REFWGT=1.0D00/DENORM
CPAM Insert: Compute the variational second-order energy.
CPAM Use unshifted H0. Save any shifts, then restore them.
      SAV=SHIFT
      SAVI=SHIFTI
      SHIFT=0.0d0
      SHIFTI=0.0d0
      CALL SIGMA_CASPT2(1.0d0,0.0d0,IVECX,IVECT)
      SHIFT=SAV
      SHIFTI=SAVI
      CALL POVLVEC(IVECX,IVECT,OVLAPS)
      E2CORR=2.0D0*E2NONV+OVLAPS(0,0)
CPAM End of insert.
      ESHIFT=E2CORR-E2NONV
      E2TOT=EREF+E2CORR

      IF(IPRGLB.GT.USUAL) THEN
        WRITE(6,*)
        WRITE(6,*)' Correlation energy /Case, /Symm, and sums:'
        DO IC=1,13
         WRITE(6,'(1X,A8,9F12.8)')
     &      CASES(IC),(ECORR(IS,IC),IS=1,NSYM),ECORR(0,IC)
        END DO
        WRITE(6,'(1X,A8,9F12.8)')
     &    'Summed: ', (ECORR(IS,0),IS=1,NSYM),ECORR(0,0)
      ENDIF

      IF (IPRGLB.GE.TERSE) THEN
         WRITE(6,*)

         If (.not.Input % LovCASPT2) Then
            WRITE(6,'(6x,a,f18.10)')'Reference energy:     ',EREF
            WRITE(6,'(6x,a,f18.10)')'E2 (Non-variational): ',E2NONV
            IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
              WRITE(6,'(6x,a,f18.10)')'Shift correction:     ',ESHIFT
            END IF
            WRITE(6,'(6x,a,f18.10)')'E2 (Variational):     ',E2CORR
            If (.not.Input % FnoCASPT2) Then
               WRITE(6,'(6x,a,f18.10)')'Total energy:         ',E2TOT
            Else
               WRITE(6,'(6x,a,f18.10,a)')'FNO correction:       ',EMP2,
     &              '   (estimate)   '
               WRITE(6,'(6x,a,f13.5)')
               E2TOT=E2TOT+EMP2
               WRITE(6,'(6x,a,f18.10,a)')'Total energy:         ',E2TOT,
     &              '   (FNO-CASPT2) '
            EndIf
            WRITE(6,'(6x,a,f18.10)')'Residual norm:        ',RNORM
            WRITE(6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT
         Else
            WRITE(6,'(6x,a,f18.10)')
     &              'Reference energy:                 ',EREF
            WRITE(6,'(6x,a,f18.10)')
     &              'Active-Site E2 (Non-variational): ',E2NONV
            IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
              WRITE(6,'(6x,a,f18.10)')
     &              'Shift correction:                 ',ESHIFT
            END IF
            WRITE(6,'(6x,a,f18.10)')
     &              'Active-Site E2 (Variational):     ',E2CORR
            WRITE(6,'(6x,a,f18.10)')
     &              'Frozen region E2 :                ',EMP2
            WRITE(6,'(6x,a,f18.10)')
     &              'Residual norm:                    ',RNORM
            WRITE(6,'(6x,a,f13.5)')
     &              'Reference weight:                 ',REFWGT
            WRITE(6,'(6x,a,f13.5)')
            E2TOT=E2TOT+EMP2
            WRITE(6,'(6x,a,f18.10)')
     &              'Total energy (LovCASPT2):         ',E2TOT
         EndIf
      END IF

* In automatic verification calculations, the precision is lower
* in case of Cholesky calculation.
      LAXITY=8
      IF(IfChol) LAXITY=Cho_X_GetTol(LAXITY)
      Call Add_Info('E_CASPT2',[E2TOT],1,LAXITY)

      IF(IPRGLB.GE.USUAL) THEN
       WRITE(6,*)
       WRITE(6,'(6x,a)')
     &  'Contributions to the CASPT2 correlation energy'
       WRITE(6,'(6x,a,F18.10)')
     &  'Active & Virtual Only:    ',EATVX+EBVAT
       WRITE(6,'(6x,a,F18.10)')
     &  'One Inactive Excited:     ',EVJTU+EAIVX+EBJAT
       WRITE(6,'(6x,a,F18.10)')
     &  'Two Inactive Excited:     ',EVJTI+EVJAI+EBJAI
       WRITE(6,*)
      END IF
      CALL GETMEM('LISTS','FREE','INTE',LLISTS,NLSTOT)
      CALL QEXIT('PCG')
      RETURN
      END
