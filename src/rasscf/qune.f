************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE QUNE(NCALLS,ENOW,BK,XSX,VL,VM,XQN,XOLD,V1,V2,NDIM,
     &                LUQUNE,TMIN,QNSTEP,QNUPDT,KSDFT)
C
C     RASSCF PROGRAM VERSION IBM-3090: SX SECTION
C
C     PURPOSE: THIS SUBROUTINE MAKES A QUASI NEWTON UPDATE
C              OF THE ROTATION MATRIX X, OR A LINE SEARCH.
C
C

      IMPLICIT REAL*8 (A-H,O-Z)
#include "output_ras.fh"
      Parameter (ROUTINE='XXXXXXXX')

      CHARACTER*2 QNSTEP
      CHARACTER*3 QNUPDT
      Character*16 KSDFT
      DIMENSION BK(NDIM),XSX(NDIM),VL(NDIM),VM(NDIM)
      DIMENSION XQN(NDIM),XOLD(NDIM),V1(NDIM),V2(NDIM)
      SAVE ALPHA,BETA,ELAST,FPLAST,NVEC,NLS
CBOR 000906      PARAMETER (MXVEC=52)
CBOR 000906     DIMENSION ALPHA(MXVEC),BETA(MXVEC)

#include "rasdim.fh"
      DIMENSION ALPHA(mxiter+2),BETA(mxiter+2)
      NCALLS=NCALLS+1

      IF(NCALLS.EQ.1) THEN
       NVEC = 0
* NLS: Nr of consecutive line searches.
       NLS=0
       CALL DCOPY_(mxiter+2,[0.0d0],0,ALPHA,1)
       CALL DCOPY_(mxiter+2,[0.0d0],0,BETA,1)
       IAD=0
       CALL DDAFILE(LUQUNE,1,BK,NDIM,IAD)
       CALL DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
       CALL DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
       ELAST=ENOW
       FPLAST=2.0D0*DDOT_(NDIM,BK,1,XSX,1)
       TMIN=0.D0
       QNSTEP='SX'
       QNUPDT=' NO'
       RETURN
      ENDIF
C -- 1. CALC L AND M ARRAYS
* Read from the beginning of LUQUNE. There is the BLB gradient, the proposed optimal
* (QN) step at the previous iteration, and the actual step taken in that iteration.
* VM should be computed as the difference in QN steps, where the new one is computed
* without the new update vectors (which are still unknown).
      IAD=0
      CALL DDAFILE(LUQUNE,2,V1,NDIM,IAD)
      CALL DCOPY_(NDIM,BK,1,VL,1)
      CALL DAXPY_(NDIM,-1.0D00,V1,1,VL,1)
* VL is the difference between the old and the new BLB gradient array:
* Read old QN step into VM:
      CALL DDAFILE(LUQUNE,2,VM,NDIM,IAD)
      CALL DDAFILE(LUQUNE,2,XOLD,NDIM,IAD)
C -- FP WILL BE USED IN LINE SEARCH ANALYSIS LATER.
      FP=2.0D0*DDOT_(NDIM,BK,1,XOLD,1)
C -- QN VECTOR UPDATED WITHOUT NEWEST UPDATE VECTORS: XQN=-HINV*BK
      CALL DCOPY_(NDIM,XSX,1,XQN,1)
      DO IVEC=1,NVEC
        CALL DDAFILE(LUQUNE,2,V1,NDIM,IAD)
        CALL DDAFILE(LUQUNE,2,V2,NDIM,IAD)
        X=-DDOT_(NDIM,V1,1,BK,1)
        Y=-DDOT_(NDIM,V2,1,BK,1)
        P1=ALPHA(IVEC)*X+BETA(IVEC)*Y
        P2=BETA(IVEC)*X
        CALL DAXPY_(NDIM,P1,V1,1,XQN,1)
        CALL DAXPY_(NDIM,P2,V2,1,XQN,1)
      END DO
* Subtract. VM will now contain the difference in QN steps, where the new one
* is computed without the update (which we still do not know).
      CALL DAXPY_(NDIM,-1.0D00,XQN,1,VM,1)
C -- NOTE: M ARRAY = HK(INV)*LK
C -- 2. DECIDE ON WHICH ROTATION TO USE.
C -- A LINE SEARCH ANALYSIS.
      C0=ELAST
      C1=FPLAST
      C2=3.0D0*(ENOW-ELAST)-2.0D0*FPLAST-FP
      C3=-2.0D0*(ENOW-ELAST)+FPLAST+FP
C     Write(LF,*)('*',I=1,60)
C     Write(LF,*)' SEARCH FOR MINIMUM:'
C     Write(LF,*)' POLYNOMIAL COEFFICIENTS:'
C     Write(LF,'(A,4F16.8)') ' ENOW,ELAST,FP,FPLAST',ENOW,ELAST,FP,FPLAST
C     Write(LF,'(A,4F16.8)') ' C0,C1,C2,C3',C0,C1,C2,C3
      P=3.0D0*C1*C3
      Q=C2**2
C -- NLM=NR OF LOCAL MINIMA
      NLM=0
      TLM =  0.0D0
      IF(ABS(P).GT.0.001D00*Q) THEN
        IF(Q.GT.P) THEN
C         Write(LF,*)' THIS IS A 3RD DEGREE POLY WITH 2 STAT. POINTS.'
          NLM=1
          TLM=(SQRT(Q-P)-C2)/(3.0D0*C3)
C         Write(LF,*)' THERE IS A LOCAL MINIMUM AT'
C         Write(LF,*)' TLM=',TLM
        ELSE
C         Write(LF,*)' THIS IS A MONOTONOUS 3RD DEGREE POLY.'
        END IF
      ELSE IF(ABS(C2).GT.0.001D00*C1) THEN
        IF(C2.GT.0.0D00) THEN
C         Write(LF,*)' THIS IS A 2ND DEGREE POLY WITH A MINIMUM.'
          NLM=1
          TLM=-C1/(2.0D0*C2)
C         Write(LF,*)' THERE IS A LOCAL MINIMUM AT'
C         Write(LF,*)' TLM=',TLM
        ELSE
C         Write(LF,*)' THIS IS A 2ND DEGREE POLY WITH A MAXIMUM.'
        END IF
      END IF
C     IF(NLM.EQ.0)Write(LF,*)' NO LOCAL MINIMUM.'
      T0=-0.5D00
      T1=+2.5D00
      IF(NLM.EQ.1) THEN
        ELM=C0+TLM*(C1+TLM*(C2+TLM*C3))
C       Write(LF,*)' LOCAL MINIMUM AT TLM=',TLM
C       Write(LF,*)'     ENERGY AT TLM IS=',ELM
C --  DISREGARD LOCAL MINIMUM IF TOO FAR AWAY:
        IF(TLM.GT.T1) NLM=0
        IF(TLM.LT.T0) NLM=0
C       IF(NLM.EQ.0) Write(LF,*)' TOO FAR AWAY -- REJECTED.'
      END IF
C -- RULES: 1. IF THERE IS A LOCAL MINIMUM IN THE TRUST REGION
C              -0.4<TLM<1.4, THIS IS USED.
C --        2. ELSE, SELECT THE LOWEST MINIMUM IN -0.5..+2.5.
      IF((NLM.EQ.1).AND.(ABS(TLM-0.5D00).LT.0.9D00)) THEN
        TMIN=TLM
        EMIN=ELM
C       Write(LF,*)' THE LOCAL MINIMUM IS USED.'
      ELSE
        E0=C0+T0*(C1+T0*(C2+T0*C3))
        E1=C0+T1*(C1+T1*(C2+T1*C3))
C       Write(LF,*)' EXTRAP. ENERGIES IN T=-0.5 AND 2.5 ARE'
C       Write(LF,*) E0,E1
        IF(E0.LT.E1) THEN
          TMIN=T0
          EMIN=E0
        ELSE
          TMIN=T1
          EMIN=E1
        END IF
        IF(NLM.EQ.1) THEN
          ELM=C0+TLM*(C1+TLM*(C2+TLM*C3))
          IF(ELM.LT.EMIN) THEN
            TMIN=TLM
            EMIN=ELM
          ELSE
            NLM=0
          END IF
        END IF
      END IF
*      Write(LF,*)'  SELECTED TMIN:',TMIN
*      Write(LF,*)' PREDICTED EMIN:',EMIN
*      Write(LF,*)('*',I=1,60)

* Predicted energy lowering, depending on step taken.
* EPRED_LS should be quite exact, if step length is not too long.
* But EPRED_SX and EPRED_QN are much more uncertain.
      EPRED_LS=EMIN-ENOW
      EPRED_SX=0.5D0*DDOT_(NDIM,BK,1,XSX,1)
      EPRED_QN=0.5D0*DDOT_(NDIM,BK,1,XQN,1)
*      WRITE(LF,*)' Predicted QN energy change:', EPRED_QN
*      WRITE(LF,*)' Predicted SX energy change:', EPRED_SX
*      WRITE(LF,*)' Predicted LS energy change:', EPRED_LS

* Here follows decision whether to update inverse Hessian, or not:
      IF(NLM.EQ.1) THEN
        X=C3*(1.0D00-TMIN)/SQRT(Q-P)
        IF((ABS(X).LT.0.2D00).AND.(TMIN.GT.0.5D00)) THEN
CPAM THEN THE ERROR ARISING FROM NONLINEARITY OF GRADIENT ALONG THE
CPAM SEARCH DIRECTION IS SMALLER THAN ABOUT 25 PERCENT, SO IT IS
CPAM MEANINGFUL INFORMATION ABOUT GRADIENT HESSIAN. ALSO, THE
CPAM INFORMATION GIVEN BY THE LAST GRADIENT IS BETTER THAN THE ONE
CPAM FROM THE NEXT-TO-LAST. THEREFORE, UPDATE INVERSE HESSIAN.
C -- WE ARE PROBABLY IN A QUADRATIC MINIMUM REGION, SO UPDATE
C         Write(LF,*)' QUAD MINIMUM REGION -- UPDATE INVERSE HESSIAN.'
          QNUPDT='YES'
        ELSE
C -- WE ARE OUTSIDE QUADRATIC MINIMUM.
C         Write(LF,*)' NOT IN QUADRATIC MIN. INV HESSIAN NOT UPDATED.'
          QNUPDT=' NO'
        END IF
      ELSE
C -- DEFINITELY NON-QUADRATIC REGION.
*        Write(LF,*)' VERY NON-QUADRATIC REGION.'
*        Write(LF,*)' SHOULD WE SCRAP ALL THE UPDATE VECTORS??'
*        Write(LF,*)' NO, BUT DO NOT UPDATE.'
        QNUPDT=' NO'
*        If(KSDFT(1:3).ne.'SCF') Then
*           NVEC=0
*        End If
      END IF

* I strongly mistrust the above analysis. For the moment, decide to update
* always. The analysis should instead be based on if whether a large rotation
* remains (do not update) or has been done (scrap the vectors).
      IF(QNUPDT.NE.'YES') THEN
*        Write(LF,*)' DECIDE TO UPDATE ANYWAY! (PAM Dec 2006)'
        QNUPDT='YES'
      END IF

      IF(QNUPDT.EQ.'YES') THEN
* Determine new pair of update vectors, and also use it to correct XQN:
* Recall: VL is the difference between the old and the new BLB gradient array:
* VM is the difference between the old, and the new provisional, XQN arrays.
        NVEC=NVEC+1
C METHOD: BFGS
        CALL DCOPY_(NDIM,XOLD,1,V1,1)
        CALL DCOPY_(NDIM,  VM,1,V2,1)
        X=DDOT_(NDIM,V1,1,VL,1)
        Y=DDOT_(NDIM,V2,1,VL,1)
        ALPHA(NVEC)=(1.0D00+Y/X)/X
        BETA(NVEC)=-1.0D00/X
C METHOD: SYMMETRIZED POWELL 1-RANK:
*        CALL DCOPY_(NDIM,VL,1,V1,1)
*        CALL DCOPY_(NDIM,XOLD,1,V2,1)
*        CALL DAXPY_(NDIM,-1.0D00,VM,1,V2,1)
*        X=DDOT_(NDIM,V1,1,V1,1)
*        Y=DDOT_(NDIM,V2,1,V1,1)
*        ALPHA(NVEC)=-Y/(X**2)
*        BETA(NVEC)=1.0D00/X
C --  ADD THE NEWEST UPDATE CONTRIBUTION INTO XQN:
        X=-DDOT_(NDIM,V1,1,BK,1)
        Y=-DDOT_(NDIM,V2,1,BK,1)
        P1=ALPHA(NVEC)*X+BETA(NVEC)*Y
        P2=BETA(NVEC)*X
        CALL DAXPY_(NDIM,P1,V1,1,XQN,1)
        CALL DAXPY_(NDIM,P2,V2,1,XQN,1)
        CALL DDAFILE(LUQUNE,1,V1,NDIM,IAD)
        CALL DDAFILE(LUQUNE,1,V2,NDIM,IAD)
* We have added to XQN the two-rank update,
*  (alpha*|V1><V1| + beta*|v1><V2| + beta |V2><V1|)|BK>
      END IF

      EPRED_QN=0.5D0*DDOT_(NDIM,BK,1,XQN,1)
*      WRITE(LF,*)' Predicted QN energy change (revised):', EPRED_QN

C -- IF LINE SEARCH IS NEARLY CONVERGED, USE THE QN OR SX STEP
C -- TO GET NEW DIRECTION, ELSE CONTINUE LINE SEARCH:
      X=TMIN-1.0D00
      IF((ABS(X).LT.0.4D00).OR.(ABS(EPRED_LS).LT.1.0D-08).or.
     &                                    KSDFT(1:3).ne.'SCF') THEN
*        Write(LF,*)' THE LINE SEARCH MINIMUM IS PREDICTED TO BE'
*        Write(LF,*)' RATHER CLOSE TO THE CURRENT POINT.'
*        Write(LF,*)' THEREFORE, DO NOT CONTINUE BY LS.'
        IF(NVEC.EQ.0) QNSTEP='SX'
        IF(NVEC.GT.0) QNSTEP='QN'
      ELSE IF (NLS.GE.2) THEN
*        Write(LF,*)' WE HAVE ALREADY USED LINE SEARCHES CONSECUTIVELY'
*        Write(LF,*)'    FOR SEVERAL ITERATIONS. WE MUST TRY TO BREAK'
*        Write(LF,*)'    OUT FROM THIS LINE SEARCH: DISABLE IT.'
        IF(NVEC.EQ.0) QNSTEP='SX'
        IF(NVEC.GT.0) QNSTEP='QN'
      ELSE
C Else we should use a line search. But maybe check predicted effect:
C If predicted LS-energy minimum is not an improvement
C then do not use line search.
*        Write(LF,*)' THE LINE-SEARCH ANALYSIS PREDICTS A MINIMUM'
*        Write(LF,*)'    PRETTY FAR FROM WHERE WE GOT WITH THE LAST'
*        Write(LF,*)'    STEP. THEREFORE; WE MAY CONSIDER CONTINUING'
*        Write(LF,*)'    BY LS, ALONG THE SAME DIRECTION AS LAST STEP.'
*        Write(LF,*)' BUT WILL WE GAIN FROM THAT? LETS SEE!'
        IF(EPRED_LS.GT.EPRED_QN) THEN
*          Write(LF,*)' LS REJECTED -- GIVES TOO LITTLE.'
          IF(NVEC.EQ.0) QNSTEP='SX'
          IF(NVEC.GT.0) QNSTEP='QN'
        ELSE IF(EPRED_LS.GT.EPRED_SX) THEN
*          Write(LF,*)' USE SX.'
          QNSTEP='SX'
        ELSE
          QNSTEP='LS'
        END IF
      END IF
      IF (QNSTEP.EQ.'LS') THEN
*        Write(LF,*)' USE LINE SEARCH.'
        X=TMIN-1.0D00
        CALL DYAX(NDIM,X,XOLD,1,XSX,1)
      END IF
      IF (QNSTEP.EQ.'QN') THEN
*          Write(LF,*)' USE QN STEP.'
          CALL DCOPY_(NDIM,XQN,1,XSX,1)
      END IF
*      IF (QNSTEP.EQ.'SX') Write(LF,*)' USE SX STEP.'

* Often, use of qune results in getting caught in a cyclic or caotic
* attractor towards the end. As an attempt to break such a pattern,
* apply a scaling whenever energy goes up:
      IF(ENOW.gt.ELAST) THEN
        SCLFCT=0.7D0
        CALL DSCAL_(NDIM,SCLFCT,XSX,1)
      END IF

* Before committing the finally suggested step (which is now in XSX),
* also apply a step size limitation. The squared 2-norm of the vector XSX
* sum of squares of rotation angles, so the 2-norm is a strict limit on
* largest rotation angle, which we (arbitrarily) limit to 0.5 (say):
      XSXNRM=DNRM2_(NDIM,XSX,1)
      SCLFCT=1.0D0/(1.0D0+2.0D0*XSXNRM)
      CALL DSCAL_(NDIM,SCLFCT,XSX,1)

      ELAST=ENOW
      FPLAST=2.0D0*DDOT_(NDIM,BK,1,XSX,1)
      IAD=0
      CALL DDAFILE(LUQUNE,1,BK,NDIM,IAD)
      CALL DDAFILE(LUQUNE,1,XQN,NDIM,IAD)
      CALL DDAFILE(LUQUNE,1,XSX,NDIM,IAD)
      NLS=NLS+1
      IF(QNSTEP.NE.'LS') NLS=0
      RETURN
      END
