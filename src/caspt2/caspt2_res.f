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
      Subroutine CASPT2_Res
C
      Implicit Real*8 (A-H,O-Z)
C
#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"

C#include "SysDef.fh"
C
      !! 1) Calculate the derivative of the CASPT2 energy with respect
      !!    to the amplitude.
      !! 2) In the standard CASPT2, solve the CASPT2 equation. In the
      !!    diagonal CASPT2, compute the lambda directly.
C
C     write(6,*) "in CASPT2_res"
      If (MAXIT.ne.0) THEN !  .and. (SHIFT.NE.0.0D+00.or.SHIFTI.ne.0.0D+00)) Then
        iRHS2  = 7
      End If
C
      !! Copy the solution vector to the residual space
      Do iCase = 1, 13
C       write(6,*) "icase=",icase
C       if (icase.ne.12.and.icase.ne.13) cycle
C       if (icase.ne.10.and.icase.ne.11) cycle
C       if (icase.ne. 8.and.icase.ne. 9) cycle
        Do iSym = 1, nSym
          nIN = nINDEP(iSym,iCase)
          IF(NIN.EQ.0) Cycle
          nAS = nASUP(iSym,iCase)
          nIS = nISUP(iSym,iCase)
C Remember: NIN values in BDIAG, but must read NAS for correct
C positioning.
          Call GETMEM('LBD','ALLO','REAL',LBD,nAS)
          Call GETMEM('LID','ALLO','REAL',LID,nIS)
          iD = iDBMat(iSym,iCase)
          Call dDaFile(LUSBT,2,Work(LBD),nAS,iD)
          Call dDaFile(LUSBT,2,Work(LID),nIS,iD)
C         if (icase.eq.4) then
C           do i = 1, nis
C             write(6,*) "ir = ",i
C             do j = 1, nin
C               write(6,'(i3,3f20.10)') i,work(lbd+j-1),work(lid+i-1),
C    *          1.0d+00/(work(lbd+j-1)+work(lid+i-1))
C             end do
C           end do
C         end if

          Call RHS_ALLO(nIN,nIS,lg_V1)
          Call RHS_ALLO(nIN,nIS,lg_V2)
          !! Read the solution vector
          Call RHS_Read(nIN,nIS,lg_V2,iCase,iSym,iVecX)
          !! Save it in the residual vector space immediately
C         Call RHS_Save(nIN,nIS,lg_V2,iCase,iSym,iVecR)
          !! Read the RHS vector
          Call RHS_Read(nIN,nIS,lg_V1,iCase,iSym,iRHS)
          !! Scale the RHS vector appropriately (compute lambda)
          Call CASPT2_ResD(1,nIN,nIS,lg_V1,Work(LBD),Work(LID))
          !! T <- T + lambda
          Call DScal_(nIN*nIS,2.0D+00,Work(lg_V1),1)
          If (MaxIt.eq.0) Then
C           call dcopy_(nin*nis,0.0d+00,0,work(lg_v2),1)
C           Call DaXpY_(nIN*nIS,2.0D+00,Work(lg_V1),1,Work(lg_V2),1)
            !! Save the modified T in the original T
C           Call RHS_Save(nIN,nIS,lg_V2,iCase,iSym,iVecX)
C           write(6,*) "lambda"
C       do i = 1, nin*nis
C       write(6,'(i3,f20.10)') i,work(lg_v1+i-1)
C       end do
C       do i = 1, 10
C       write(6,*) i,work(lg_v1)
C       end do
            Call RHS_Save(nIN,nIS,lg_V1,iCase,iSym,iVecR)
          Else
C           Call RHS_Save(nIN,nIS,lg_V1,iCase,iSym,iRHS2)
              Call RHS_ALLO(NAS,NIS,lg_V3)
              CALL RHS_READ(NAS,NIS,lg_V3,ICASE,ISYM,iRHS)
              CALL RHS_SAVE(NAS,NIS,lg_V3,ICASE,ISYM,iRHS2)
              Call RHS_FREE(NAS,NIS,lg_V3)
          End If
          Call RHS_Free(nIN,nIS,lg_V1)
          Call RHS_Free(nIN,nIS,lg_V2)

          Call GETMEM('LBD','FREE','REAL',LBD,nAS)
          Call GETMEM('LID','FREE','REAL',LID,nIS)
        End Do
      End Do
C
C     Now, going to solve the Lambda for non-variational CASPT2, i.e.
C     with real/imaginary shift, without the diagonal approximation.
C     The following is just a copy-and-paste of eqctl2.f and pcg.f,
C     but some unnecessary lines (comuptation of energy etc.) are
C     omitted.
C
C Transform RHS of CASPT2 equations to eigenbasis for H0:
C     CALL PTRTOSR(1,IVECW,IRHS)
C
      !! We need IRHS,IVECR,IVECX,IVECC,IVECC2
      !! The original IRHS is no longer needed (?),
      !! but has to be modified for (X)MS
      !! IVECR is also not needed
      !! IVECX is needed, so use a different array
      !! IVECC and IVECC2 are later transformed
      If (MAXIT.ne.0) THEN !  .and. SHIFT.NE.0.0D+00.or.SHIFTI.ne.0.0D+00) Then
      SAV=SHIFT
      SAVI=SHIFTI
      SHIFT=0.0d0
      SHIFTI=0.0d0
      CALL SIGMA_CASPT2(2.0d+00,2.0d+00,IVECX,iRHS2)
      SHIFT=SAV
      SHIFTI=SAVI
C
      iVecXbk = iVecX
      iVecRbk = iVecR
      iRHSbk  = iRHS
      iVecX   = iVecR
      iRHS    = 7
      iVecR   = 8

      Call PCG_RES(ICONV)
C          CALL PCOLLVEC(IVECX,0)
C          CALL PCOLLVEC(IVECR,0)
      IF (ICONV .NE. 0) THEN
        WRITE (6,'(" Lambda equation did not converge...")')
        WRITE (6,'(" Continue anyway?")')
      END IF
C
      !! Restore contravariant and covariant representations of the
      !! non-variational T-amplitude
      iVecX   = iVecXbk
      iVecR   = iVecRbk
      iRHS    = iRHSbk
      CALL PTRTOC(0,IVECX,IVECC)
      CALL PTRTOC(1,IVECX,IVECC2)
      End If
C     Do iCase = 1, 13
C       write(6,*) "icase=",icase
C       Do iSym = 1, nSym
C         nIN = nINDEP(iSym,iCase)
C         IF(NIN.EQ.0) Cycle
C         nAS = nASUP(iSym,iCase)
C         nIS = nISUP(iSym,iCase)

C         Call RHS_ALLO(nIN,nIS,lg_V1)
C         !! Read the solution vector
C         Call RHS_Read(nIN,nIS,lg_V1,iCase,iSym,iVecR)
C         do i = 1, nin*nis
C           write(6,'(i3,f20.10)') i,work(lg_v1+i-1)
C         end do
C         Call RHS_Free(nIN,nIS,lg_V1)
C       End Do
C     End Do
C
C
C
      RETURN
      END
C
C-----------------------------------------------------------------------
C
      !! RHS_SGMDIA
      SUBROUTINE CASPT2_ResD(Mode,NIN,NIS,lg_W,DIN,DIS)
      IMPLICIT REAL*8 (A-H,O-Z)

#include "rasdim.fh"
#include "caspt2.fh"
#include "output.fh"
#include "WrkSpc.fh"
#include "eqsolv.fh"
      DIMENSION DIN(*),DIS(*)

C Apply the resolvent of the diagonal part of H0 to an RHS array

#include "para_info.fh"
#ifdef _MOLCAS_MPP_
#include "global.fh"
#include "mafdecls.fh"
#endif

#ifdef _MOLCAS_MPP_
      IF (Is_Real_Par()) THEN
        CALL GA_Sync()
        myRank = GA_NodeID()
C-SVC: get the local vertical stripes of the lg_W vector
        CALL GA_Distribution (lg_W,myRank,iLo,iHi,jLo,jHi)
        IF (iLo.NE.0.AND.jLo.NE.0) THEN
          NROW=iHi-iLo+1
          NCOL=jHi-jLo+1
          CALL GA_Access (lg_W,iLo,iHi,jLo,jHi,mW,LDW)
          CALL CASPT2_ResD2(MODE,NROW,NCOL,DBL_MB(mW),LDW,DIN(iLo),
     &                DIS(jLo),SHIFT,SHIFTI)
          CALL GA_Release_Update (lg_W,iLo,iHi,jLo,jHi)
        END IF
        CALL GA_Sync()
C       CALL GAdSUM_SCAL(DOVL)
      ELSE
        CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                   SHIFT,SHIFTI)
      END IF
#else
      CALL CASPT2_ResD2(MODE,NIN,NIS,WORK(lg_W),NIN,DIN,DIS,
     &                 SHIFT,SHIFTI)
#endif

      END
C
C-----------------------------------------------------------------------
C
      !! RESDIA
      SUBROUTINE CASPT2_ResD2(Mode,NROW,NCOL,W,LDW,DIN,DIS,
     &                  SHIFT,SHIFTI)
      IMPLICIT REAL*8 (A-H,O-Z)

      DIMENSION W(LDW,*),DIN(*),DIS(*)

      DO J=1,NCOL
        DO I=1,NROW
          If (Mode.eq.1) Then
            DELTA  = SHIFT+DIN(I)+DIS(J)
            DELINV = DELTA/(DELTA**2+SHIFTI**2)
            !! The following SCAL is the actual residual
            SCAL   = 1.0D+00 - (DIN(I)+DIS(J))*DELINV
C           write(6,*) "residue = ", scal
C           if (abs(residue).ge.1.0d-08) write(6,*) "residue = ", scal
            !! Another scaling is required for lambda
            SCAL   =-SCAL*DELINV
          ELse If (Mode.eq.2) Then
            SCAL   =-SHIFTI/(DIN(I)+DIS(J))
          End If
          W(I,J) = SCAL*W(I,J)
        END DO
      END DO
      END
C
C-----------------------------------------------------------------------
C
      SUBROUTINE PCG_RES(ICONV)
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
C     CALL GETMEM('LISTS','ALLO','INTE',LLISTS,NLSTOT)
C     CALL MKLIST(iWORK(LLISTS))


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
       WRITE(6,*) "RASPT2 with level shift is non-variational,"
       WRITE(6,*) "so the Lambda equation has to be solved"//
     *            " for analytic gradients"
       WRITE(6,*) "Following values are nonsense (or I just don't"//
     *            " know the meaning)"
       WRITE(6,*)
C      WRITE(6,*)'The contributions to the second order'//
C    &     ' correlation energy in atomic units.'
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
      CALL SIGMA_CASPT2(1.0D00,0.0D0,IVECP,IVECT)
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
       WRITE(6,'(25A5)')('-----',I=1,25)
       WRITE(6,*)
C      WRITE(6,*)' FINAL CASPT2 RESULT:'
      END IF
C     CALL POVLVEC(IRHS,IVECX,ECORR)
C     EVJTU=ECORR(0,1)
C     EVJTI=ECORR(0,2)+ECORR(0,3)
C     EATVX=ECORR(0,4)
C     EAIVX=ECORR(0,5)
C     EVJAI=ECORR(0,6)+ECORR(0,7)
C     EBVAT=ECORR(0,8)+ECORR(0,9)
C     EBJAT=ECORR(0,10)+ECORR(0,11)
C     EBJAI=ECORR(0,12)+ECORR(0,13)
C     E2NONV=ECORR(0,0)
C     CALL POVLVEC(IVECX,IVECX,OVLAPS)
C     DENORM=1.0D0+OVLAPS(0,0)
C     REFWGT=1.0D00/DENORM
CPAM Insert: Compute the variational second-order energy.
CPAM Use unshifted H0. Save any shifts, then restore them.
C     SAV=SHIFT
C     SAVI=SHIFTI
C     SHIFT=0.0d0
C     SHIFTI=0.0d0
C     CALL SIGMA_CASPT2(1.0d0,0.0d0,IVECX,IVECT)
C     SHIFT=SAV
C     SHIFTI=SAVI
C     CALL POVLVEC(IVECX,IVECT,OVLAPS)
C     E2CORR=2.0D0*E2NONV+OVLAPS(0,0)
CPAM End of insert.
C     ESHIFT=E2CORR-E2NONV
C     E2TOT=EREF+E2CORR

C     IF(IPRGLB.GT.USUAL.or.iprglb.ne.silent) THEN
C       WRITE(6,*)
C       WRITE(6,*)' Correlation energy /Case, /Symm, and sums:'
C       DO IC=1,13
C        WRITE(6,'(1X,A8,9F12.8)')
C    &      CASES(IC),(ECORR(IS,IC),IS=1,NSYM),ECORR(0,IC)
C       END DO
C       WRITE(6,'(1X,A8,9F12.8)')
C    &    'Summed: ', (ECORR(IS,0),IS=1,NSYM),ECORR(0,0)
C     ENDIF

      IF (IPRGLB.GE.TERSE) THEN
      !  WRITE(6,*)

      !  If (.not.Input % LovCASPT2) Then
C     !     WRITE(6,'(6x,a,f18.10)')'Reference energy:     ',EREF
      !     WRITE(6,'(6x,a,f30.20)')'Reference energy:     ',EREF
C     !     WRITE(6,'(6x,a,f18.10)')'E2 (Non-variational): ',E2NONV
      !     WRITE(6,'(6x,a,f30.20)')'E2 (Non-variational): ',E2NONV
      !     IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
C     !       WRITE(6,'(6x,a,f18.10)')'Shift correction:     ',ESHIFT
      !       WRITE(6,'(6x,a,f30.20)')'Shift correction:     ',ESHIFT
      !     END IF
C     !     WRITE(6,'(6x,a,f18.10)')'E2 (Variational):     ',E2CORR
      !     WRITE(6,'(6x,a,f30.20)')'E2 (Variational):     ',E2CORR
      !     If (.not.Input % FnoCASPT2) Then
C     !        WRITE(6,'(6x,a,f18.10)')'Total energy:         ',E2TOT
      !        write(6,'(6x,a,f30.20)')'Total energy:         ',E2TOT
      !     Else
      !        WRITE(6,'(6x,a,f18.10,a)')'FNO correction:       ',EMP2,
     &!             '   (estimate)   '
      !        WRITE(6,'(6x,a,f13.5)')
      !        E2TOT=E2TOT+EMP2
      !        WRITE(6,'(6x,a,f18.10,a)')'Total energy:         ',E2TOT,
     &!             '   (FNO-CASPT2) '
      !     EndIf
      !     WRITE(6,'(6x,a,f18.10)')'Residual norm:        ',RNORM
C     !     WRITE(6,'(6x,a,f13.5)') 'Reference weight:     ',REFWGT
      !     WRITE(6,'(6x,a,f30.20)') 'Reference weight:     ',REFWGT
      !  Else
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Reference energy:                 ',EREF
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Active-Site E2 (Non-variational): ',E2NONV
      !     IF(SHIFT.NE.0.0d0.or.SHIFTI.ne.0.0d0) THEN
      !       WRITE(6,'(6x,a,f18.10)')
     &!             'Shift correction:                 ',ESHIFT
      !     END IF
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Active-Site E2 (Variational):     ',E2CORR
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Frozen region E2 :                ',EMP2
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Residual norm:                    ',RNORM
      !     WRITE(6,'(6x,a,f13.5)')
     &!             'Reference weight:                 ',REFWGT
      !     WRITE(6,'(6x,a,f13.5)')
      !     E2TOT=E2TOT+EMP2
      !     WRITE(6,'(6x,a,f18.10)')
     &!             'Total energy (LovCASPT2):         ',E2TOT
      !  EndIf
      END IF

* In automatic verification calculations, the precision is lower
* in case of Cholesky calculation.
C     LAXITY=8
C     IF(IfChol) LAXITY=Cho_X_GetTol(LAXITY)
C     Call Add_Info('E_CASPT2',[E2TOT],1,LAXITY)

C     IF(IPRGLB.GE.USUAL) THEN
C      WRITE(6,*)
C      WRITE(6,'(6x,a)')
C    &  'Contributions to the CASPT2 correlation energy'
C      WRITE(6,'(6x,a,F18.10)')
C    &  'Active & Virtual Only:    ',EATVX+EBVAT
C      WRITE(6,'(6x,a,F18.10)')
C    &  'One Inactive Excited:     ',EVJTU+EAIVX+EBJAT
C      WRITE(6,'(6x,a,F18.10)')
C    &  'Two Inactive Excited:     ',EVJTI+EVJAI+EBJAI
C      WRITE(6,*)
C     END IF
C     CALL GETMEM('LISTS','FREE','INTE',LLISTS,NLSTOT)
      CALL QEXIT('PCG')
C
      RETURN
C
      END SUBROUTINE PCG_RES
