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
      SUBROUTINE DAVCRE(C,HC,HH,CC,E,HD,SC,Q,QQ,S,
     *                  SXSEL,NROOT,ITMAX,NDIM,ITERSX,NSXS)
C
C RASSCF program: version IBM-3090: SX section
C
C Using the Davidson method, in the multiple root version suggested
C by B.LIU. This routine finds the NROOT lowest eigenvalues and
C eigenvectors of a secular problem of dimension NDIM.
C This is a vectorized version adapted to run optimally on IBM 3090
C It is completely CPU bound (except maybe for the construction of
C the sigma vector). Multiple calls to this subroutine may therefore
C be necessary for cases where more iterations are needed than
C is allowed by the core space requirements to store all C and
C sigma vectors.
C
C Externals: HMAT (set up the Davidson H-matrix HH)
C            COVLP (calculate the overlap between two CI vectors)
C            Jacob (full diagonalization routine)
C            DGEMM and other ESSL routines for matrix operations
C Parameters: C  SuperCI-vectors
C             HC SuperCI Sigma vectors
C             HH Davidson's H-matrix
C             CC      "     eigenvectors
C             E       "     eigenvalues
C             HD diagonal elements of SuperCI
C             Q  the Davidson update vectors
C             QQ the norm of all Q-vectors
C             SC scratch area
C
C ********** IBM-3090 Release 88 09 08 *****
C
      IMPLICIT REAL*8 (A-H,O-Z)
#include "rasdim.fh"
#include "warnings.fh"
#include "rasrc.fh"
#include "WrkSpc.fh"
#include "wadr.fh"
#include "output_ras.fh"
#include "fciqmc.fh"
      Parameter (ROUTINE='DAVCRE  ')
      CHARACTER*4 IOUTW,IOUTX
      DIMENSION C(*),HC(*),HH(*),CC(*),E(*),SC(*),
     *          Q(*),QQ(*),S(*)
      CHARACTER*(*) SXSEL
      DIMENSION HD(NROOT+NSXS)
cvv   DATA THRA/1.D-13/,THRLD2/1.D-15/,THRQ/1.D-07/,THRZ/1.D-06/,
cvv   Thrld2 changed to 1.D-14 to avoid numerial unstabillity
      DATA THRA/1.D-13/,THRLD2/5.D-14/,THRQ/1.D-07/,THRZ/1.D-06/,
     &     THRLD1/1.D-08/
C
      Call qEnter('DAVCRE')
C Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
C
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*) 'Super-CI diagonalization. Max iterations: ',ITMAX
      END IF
C
C Set values at loop head and starting guess of the vectors C
C
      Rc_SX = 0
      NTRIAL=NROOT
      NCR=NROOT*NDIM
      CALL VCLR(C,1,NCR)
      II=0
      DO I=1,NROOT
        C(II+I)=1.0D0
        II=II+NDIM
      END DO
C
C  memory allocation for COVLP
C
       CALL GETMEM('SXC1','ALLO','REAL',LC1,NSXS)
       CALL GETMEM('SXC2','ALLO','REAL',LC2,NSXS)
       CALL GETMEM('SXX2','ALLO','REAL',LX,NSXS)
C
C Begin Davidson iterations
C Start by setting up the Davidson HH-matrix
C The dimension of this matrix, NDIMH, is updated in HMAT
C
      NDIMH=0
      ITERSX=0

**************************************************************
*               Head of Super-CI iteration loop:
**************************************************************
100   CONTINUE
      ITERSX=ITERSX+1
      CALL HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
      KDIMH=(NDIMH+NDIMH**2)/2
C
C Echo the HH-matrix before diagonalizing it.
C
      CALL DCOPY_(KDIMH,HH,1,HH(KDIMH+1),1)
C
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Davidson H-matrix in iteration ',ITERSX
        Write(LF,*)' Davidson H-matrix triangular of size =  ', NDIMH
        Write(LF,'(1x,8F14.6)') (HH(I),I=1,KDIMH)
      END IF
      E(1)=HH(1)
      CC(1)=1.d0
      IF(NDIMH.GT.1) THEN
C
C set eigenvector array to identity before JACO call
C
       NDIMH2=NDIMH**2
       CALL VCLR(CC,1,NDIMH2)
       II=-NDIMH
       DO I=1,NDIMH
        II=II+NDIMH+1
        CC(II)=1.0d0
       END DO
C
       CALL Jacob(HH(KDIMH+1),CC,NDIMH,NDIMH)
       IF(SXSEL.eq.'LOWEST  ') THEN
C Root selection here assumes picking the lowest root(s).
         CALL JACORD(HH(KDIMH+1),CC,NDIMH,NDIMH)
         II=0
         DO I=1,NROOT
          II=II+I
          E(I)=HH(II+KDIMH)
         END DO
       ELSE
C Root selection by maximum overlap.
         ISEL=1
         XMX=ABS(CC(1))
         DO I=2,NDIMH
          X=ABS(CC(1+NDIMH*(I-1)))
          IF(X.GT.XMX) THEN
            ISEL=I
            XMX=X
          END IF
         END DO
         IF(ISEL.NE.I) THEN
           DO I=1,NDIMH
             SWAP=CC(I+NDIMH*(ISEL-1))
             CC(I+NDIMH*(ISEL-1))=CC(I)
             CC(I)=-SWAP
           END DO
           SWAP=HH(KDIMH+(ISEL*(ISEL+1))/2)
           HH(KDIMH+(ISEL*(ISEL+1))/2)=HH(KDIMH+1)
           HH(KDIMH+1)=SWAP
         END IF
       END IF
C
       IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Eigenvalues:'
        Write(LF,'(1X,8F14.6)') (E(I),I=1,NROOT)
        Write(LF,*)' Eigenvectors:'
        NST=0
        DO I=1,NROOT
         Write(LF,'(1X,10F11.6)') (CC(J+NST),J=1,NDIMH)
         NST=NST+NDIMH
        END DO
       ENDIF
      ENDIF
C
C Now perform the Davidson update
C First form for each root I, a vector Q(I), where
C Q(I,K)=sum(J,M) CC(I,J,M)*(HC(M,J,K)-E(I)*C(M,J,K))
C here J runs over iterations and M over the number of trial
C vectors in each iterations E(I) is the current estimate of the
C energy of the I:th root, C(M,J) is the M:th trial vector used in
C the J:th iteration, HC being the corresponding sigma vector. CC
C is the matrix of eigenvectors of the Davidson Hamiltonian.
C
C STEP 1: form the matrix CCE=-CC*E (CCE in array SC)
      JI=0
      DO I=1,NROOT
       EI=E(I)
       DO J=1,NDIMH
        JI=JI+1
        S(JI)=-CC(JI)*EI
       END DO
      END DO
C
C Step 2: form the matrix Q=HC*CC
      CALL DGEMM_('N','N',
     &            NDIM,NROOT,NDIMH,
     &            1.0d0,HC,NDIM,
     &            CC,NDIMH,
     &            0.0d0,Q(NDIM+1),NDIM)
C Step 3: add the contribution C*CCE
      CALL DGEMM_('N','N',NDIM,NROOT,NDIMH,
     & 1.D0,C,NDIM,S,NDIMH,1.D0,Q(NDIM+1),NDIM)
C
C Note that the NDIM first positions in Q are untouched so far
C the vector will be moved later
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' The residual vectors Q of size:', NDIM
        Write(LF,'(1x,8F14.10)')(Q(NDIM+I),I=1,NDIM)
      END IF
C
C Check the norms of the Q-vectors for convergence, and obtain
C bounds to the eigenvalues. E(I) is an upper bound to the I:th
C eigenvalue, from McDonald's theorem, while E(I)-NORM(Q(I))
C is a lower bound to the eigenvalue. This is the Weinstein lower
C bound formula in a multiple root form.
C
      ICONVQ=0
      IST=1+NDIM
      DO I=1,NROOT
       CALL COVLP(Q(IST),Q(IST),WORK(LDIA),WORK(LPA),WORK(LSXN),
     &            WORK(LC1),WORK(LC2),WORK(LX),QQ(I))
       IST=IST+NDIM
      END DO
C
      DO I=1,NROOT
          QNORM=SQRT(QQ(I))
          IF(QNORM.LT.THRQ)ICONVQ = ICONVQ + 1
          EI = E(I)
          ENO = EI - QNORM
          IF(NROOT.GT.1) THEN
            IF(I.GT.1) THEN
              IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(20X,F16.8,A,I2,A,F16.8)') ENO,
     &                   ' <  SX energy ',I,'  < ',EI
              END IF
            ELSE
              IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)')ITERSX,
     &               ' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
              END IF
            END IF
          ELSE
            IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)')ITERSX,
     &               ' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
            END IF
          END IF
          IF(IPRLEV.GE.DEBUG) THEN
            Write(LF,*)'  Norm of Q:',QNORM
          END IF
          IST = IST + NDIM
      END DO
CPAM00 End of replacement.
C
C    Reset ICONVQ unless all roots are converged
C
      IF(ICONVQ.LT.NROOT) ICONVQ = 0
C
C    Check the expansion vectors for convergence. This is done
C    by computing the sum of the squares of CC(I,J), for each
C    root I, and for the J trial vectors of the current iteration
C
      NST = NDIMH - NTRIAL + 1
      ICONVA = 0
      DO I = 1,NROOT
          ASQ=DDOT_(NTRIAL,CC(NST),1,CC(NST),1)
          IF(ASQ.LT.THRA) ICONVA = ICONVA + 1
          NST = NST + NDIMH
      END DO
C
C    Reset ICONVA unless all roots are converged.
C
      IF(ICONVA.LT.NROOT) ICONVA = 0
C
C    Branch out of iteration loop if convergence has been
C....    achieved on either the energy or the Davidson expansion
C....    vectors.
C
      IF(ICONVQ.NE.0.OR.ICONVA.NE.0) GO TO 35
C
C Calculate the D-vectors as D(I,K)=Q(I,K)/(E(I)-HD(K,K))
C
      IST=1
      DO I=1,NROOT
       EI=E(I)
       DO K=1,NDIM
        SC(K)=EI-HD(K)
        IF(ABS(SC(K)).LT.THRZ) SC(K)=1.0d0
       END DO
       CALL VDIV(SC,1,Q(IST+NDIM),1,Q(IST),1,NDIM)
       IST=IST+NDIM
      END DO
C Remove any unwanted components. These are signalled by
C huge elements of SX hamiltonian diagonal (set in SXHAM).
      DO I=1,NDIM
       IF(HD(I).GT.1.0D20) THEN
         DO J=1,NROOT
           Q(I+NDIM*(J-1))=0.0D0
         END DO
       END IF
      END DO
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' The correction vectors D(stored in Q):'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM)
      END IF
C
C Q now contains the Davidson correction vectors. These will now
C be orthogonalized to all old vectors and then to each other.
CPAM01 The ON is done in two passes. A vector smaller than THRLD
C after the first ON pass is considered to contribute nothing worthwhile
C to the Davidson procedure, and is discarded. Else, it is normalized.
C If it is not still normalized to within an accuracy THRLD after the
C second pass, then again it is discarded: it means that the rescaling
C after pass 1, together with numerical noise amplification in
C COVLP, is beginning to be noticable. The latter will of course happen
C when some orbital rotation(s) are almost redundant.
C
      IPASS=0
      ICONVL=0
      NTOTDC=NROOT

 23   CONTINUE

C First form the overlap matrix
      IJ=0
      ISTQ=1
      DO I=1,NTOTDC
       ISTC=1
       DO J=1,NDIMH
        IJ=IJ+1
        CALL COVLP(Q(ISTQ),C(ISTC),WORK(LDIA),WORK(LPA),WORK(LSXN),
     &  WORK(LC1),WORK(LC2),WORK(LX),OVL)
        S(IJ)=-OVL
        ISTC=ISTC+NDIM
       END DO
       ISTQ=ISTQ+NDIM
      END DO
      CALL DGEMM_('N','N',NDIM,NTOTDC,NDIMH,
     &           1.D0,C,NDIM,S,NDIMH,
     &           1.D0,Q,NDIM)
C
C The Q-vectors are now orthogonal to all C-vectors.
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)'  D vector orthogonal to all C vectors:'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM)
      END IF
C Now orthogonalize them to one another, rejecting those which
C appear with too small a norm
C
      IST=1
      NTRIAL=0

* Long loop over NTOTDC vectors:
      DO I=1,NTOTDC
       IF(I.NE.1.AND.NTRIAL.NE.0) THEN
C
C  Orthogonalize this vector to the preceeding Q's of this iteration
C
        JST=1
        DO J=1,NTRIAL
         CALL COVLP(Q(IST),Q(JST),WORK(LDIA),WORK(LPA),WORK(LSXN),
     *              WORK(LC1),WORK(LC2),WORK(LX),OVL)
         S(J)=-OVL
         JST=JST+NDIM
        END DO
*        CALL DGEMX(NDIM,NTRIAL,1.D0,Q,NDIM,S,1,Q(IST),1)
        CALL DGEMV_('N',NDIM,NTRIAL,1.D0,Q,NDIM,S,1,1.0D0,Q(IST),1)
       IF(IPRLEV.GE.INSANE) THEN
          Write(LF,*)'  Q vector orthogonal to preceeding Q vectors:'
          Write(LF,'(1x,8F14.10)')(Q(k),k=1,NDIM)
       END IF
       ENDIF
C
C  Normalize this vector and move to trial set if norm large enough
C
       CALL COVLP(Q(IST),Q(IST),WORK(LDIA),WORK(LPA),WORK(LSXN),
     *            WORK(LC1),WORK(LC2),WORK(LX),XNORM)
C Due to large noise amplification in COVLP, the squared-norm
C can actually come out as a negative number.
C Acceptable, only if it is very close to zero. Else, quit.
       IF(XNORM.LT.-1.0D-09)  then
         Write(LF,*)
         Write(LF,*)'      *** Error in subroutine DAVCRE ***'
         Write(LF,*)' The squared norm of a trial vector has been'
         Write(LF,*)' computed to be negative:'
         Write(LF,*)'      XNORM=',XNORM
         Write(LF,*)' This is possible only for some severe malfunction'
         Write(LF,*)' of the rasscf program. Please issue a bug report.'
         Write(LF,*)
         if(.not.iDoNECI) then
           Call qTrace
           Call Quit(_RC_GENERAL_ERROR_)
         else
           Write(LF,*)' non positive-semi definite matrix occurred.'
           Write(LF,*)' Calculation will continue. '
           Write(LF,*)' Divergent results might occur'
           Write(LF,*)' Tests for possible solution on the way...'
         end if
       End iF
       XNORM=sqrt(MAX(0.0D0,XNORM))
       IF(IPRLEV.GE.INSANE) THEN
         Write(LF,'(1X,A,I3,A,I3,A,E16.8)') 'Pass ',IPASS,
     &                 ' New orthogonal vector ',I,' has norm ',XNORM
       END IF

CPAM01 Two different treatments, depending on if this is first or
C second orthonormalization pass:
       IF(IPASS.EQ.0) THEN
C First pass:
         IF(ITERSX.EQ.1 .or. XNORM.GT.THRLD1) THEN
          ISTQ=NTRIAL*NDIM+1
          NTRIAL=NTRIAL+1
          XNORM=1.0D00/(XNORM+1.0D-24)
CPAM01 Note that ISTQ can be (and is!) the same as IST:
          IF(ISTQ.EQ.IST) THEN
           CALL DSCAL_(NDIM,XNORM,Q(IST),1)
          ELSE
           CALL DYAX(NDIM,XNORM,Q(IST),1,Q(ISTQ),1)
          ENDIF
         ENDIF
       ELSE
C Second pass: Demand accurate normalization, else we know that
C poor independence, maybe with strong rounding-error amplification
C in COVLP, has began to erode the orthonormalization of the
C basis vectors, hence the integrity of the Davidson procedure.
         IF(ITERSX.EQ.1 .or. ABS(XNORM-1.0D0).LT.THRLD2) THEN
          ISTQ=NTRIAL*NDIM+1
          NTRIAL=NTRIAL+1
          XNORM=1.0D00/(XNORM+1.0D-24)
CPAM01 Note that ISTQ can be (and is!) the same as IST:
          IF(ISTQ.EQ.IST) THEN
           CALL DSCAL_(NDIM,XNORM,Q(IST),1)
          ELSE
           CALL DYAX(NDIM,XNORM,Q(IST),1,Q(ISTQ),1)
          ENDIF
         ENDIF
       END IF
       IST=IST+NDIM

* End over long loop over NTOTDC vectors I=1..NTOTDC
      END DO
C
C NTRIAL new orthogonal vectors have now been formed. if NTRIAL
C equals zero there are no new linearly independent trial vectors
C and we branch out. If NTRIAL does not equal zero make a second
C orthonormalization pass
C
      IF(NTRIAL.EQ.0) ICONVL=1
      NTOTDC=NTRIAL
      IPASS=IPASS+1
      IF(IPASS.NE.2.AND.NTOTDC.NE.0) GO TO 23
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' The correction vectors, after ON:'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM*NTRIAL)
      END IF
C
C Check if the set of vectors has become linearly dependent
C
      IF(ICONVL.EQ.1) GO TO 36
C
C Move the new orthogonal vectors to C
C
      NST=1+NDIMH*NDIM
      LENGTH=NTRIAL*NDIM
      CALL DCOPY_(LENGTH,Q,1,C(NST),1)
C
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,'(1X,A,I2,A)') ' Adding ',NTRIAL,' new vectors.'
      END IF
      IF(IPRLEV.GE.INSANE) THEN
       IST = NDIMH*NDIM
       DO I = 1,NTRIAL
          Write(LF,'(1X,A,I2,A)') ' New vector ',I,' is:'
          Write(LF,'(1X,10F11.6)') (C(IST+K),K = 1,NDIM)
          IST = IST + NDIM
       END DO
      ENDIF
C
C    At this point, the calculation has neither converged, nor
C    produced linearly dependent trial vectors. Branch back to
C    the head of the iteration loop unless more than ITMAX
C    iterations have been performed.
C
      IF(ITERSX.LT.ITMAX) GO TO 100
C
C Here after ITMAX iterations, and no convergence
C
      IF(IPRLEV.GE.DEBUG) THEN
       IF(ITMAX.LT.MXSXIT) THEN
        Write(LF,*)' Super-CI not converged. Max SX iter increased.'
       ELSE
        Write(LF,*)' Super-CI not converged.'
       END IF
      END IF
      ITMAX=min(ITMAX+2,MXSXIT)
      Rc_SX = 16
      GO TO 34
C
C Here as calculation has converged.
C
35    CONTINUE
      IOUTW = 'y   '
      IF(NROOT.GT.1) IOUTW = 'ies '
      IOUTX = '    '
      IF(NROOT.GT.1) IOUTX = 's   '
      IF(ICONVQ.NE.0.AND.IPRLEV.GE.DEBUG) THEN
         Write(LF,*) ' Convergence on CI energ'//IOUTW
      END IF
      IF(ICONVA.NE.0.AND.IPRLEV.GE.DEBUG) THEN
         Write(LF,*) ' Convergence on Davidson expansion vector'//IOUTX
        Write(LF,*) IOUTX
      END IF
      GO TO 34
C
C Here as trial vectors have become linearly dependent
C
36    CONTINUE
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Trial vector set has become linearly dependent'
      END IF
34    CONTINUE
C
C Compute CI vectors for all roots
C
C CI NEW (K) =  sum(J,M)CC(I,J M)*C M J (K)
C
C Here J runs over the iterations, and M over the
C number of trial vectors used in each iteration J
C
      CALL DGEMM_('N','N',
     &            NDIM,NROOT,NDIMH,
     &            1.0d0,C,NDIM,
     &            CC,NDIMH,
     &            0.0d0,Q,NDIM)
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' Unnormalized final CI vectors (in Q):'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM)
      END IF
C
C Normalize the final CI vectors
C
      IST=1
      DO I=1,NROOT
       CALL COVLP(Q(IST),Q(IST),WORK(LDIA),WORK(LPA),WORK(LSXN),
     * WORK(LC1),WORK(LC2),WORK(LX),XNORM)
       XNORM=1.0D00/XNORM
       CALL DYAX(NDIM,XNORM,Q(IST),1,C(IST),1)
       IST=IST+NDIM
      END DO
C
C Print vector if desired
C
      IF(IPRLEV.GE.INSANE) THEN
       IST = 0
       DO I = 1,NROOT
          Write(LF,*) ' SX-CI vector for root ',I
          Write(LF,'(1X,8F8.4)') (C(IST+K),K = 1,NDIM)
          IST = IST + NDIM
       END DO
      ENDIF
C
C End of diagonalization
C Free memory for COVLP
C
      CALL GETMEM('XXXX','FREE','REAL',LC1,NSXS)
      CALL GETMEM('XXXX','FREE','REAL',LC2,NSXS)
      CALL GETMEM('XXXX','FREE','REAL',LX,NSXS)
      CALL QEXIT('DAVCRE')
      RETURN
      END
