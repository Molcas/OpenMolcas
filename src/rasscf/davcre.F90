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
      SUBROUTINE DAVCRE(C,HC,HH,CC,E,HD,SC,Q,QQ,S,                      &
     &                  SXSEL,NROOT,ITMAX,NDIM,ITERSX,NSXS)
!
! RASSCF program: version IBM-3090: SX section
!
! Using the Davidson method, in the multiple root version suggested
! by B.LIU. This routine finds the NROOT lowest eigenvalues and
! eigenvectors of a secular problem of dimension NDIM.
! This is a vectorized version adapted to run optimally on IBM 3090
! It is completely CPU bound (except maybe for the construction of
! the sigma vector). Multiple calls to this subroutine may therefore
! be necessary for cases where more iterations are needed than
! is allowed by the core space requirements to store all C and
! sigma vectors.
!
! Externals: HMAT (set up the Davidson H-matrix HH)
!            COVLP (calculate the overlap between two CI vectors)
!            Jacob (full diagonalization routine)
!            DGEMM and other ESSL routines for matrix operations
! Parameters: C  SuperCI-vectors
!             HC SuperCI Sigma vectors
!             HH Davidson's H-matrix
!             CC      "     eigenvectors
!             E       "     eigenvalues
!             HD diagonal elements of SuperCI
!             Q  the Davidson update vectors
!             QQ the norm of all Q-vectors
!             SC scratch area
!
! ********** IBM-3090 Release 88 09 08 *****
!
      use fciqmc, only : DoNECI
      use wadr, only: PA, DIA, SXN
      use stdalloc, only: mma_allocate, mma_deallocate
      use PrintLevel, only: DEBUG,INSANE
      use output_ras, only: LF,IPRLOC,RC_SX
      use RASDim, only: MxSXIt
      IMPLICIT None
      INTEGER NROOT,ITMAX,NDIM,ITERSX,NSXS
      Real*8 C((NROOT+NSXS)*NROOT*(ITMAX+1))
      Real*8 HC((NROOT+NSXS)*NROOT*ITMAX)
      Real*8 HH((ITMAX*NROOT)*(ITMAX*NROOT+1))
      Real*8 CC((ITMAX*NROOT)**2)
      Real*8 E((ITMAX*NROOT))
      Real*8 HD(NROOT+NSXS)
      Real*8 SC((NROOT+NSXS))
      Real*8 Q((NROOT+NSXS)*(NROOT+1))
      Real*8 QQ(NROOT)
      Real*8 S(ITMAX*NROOT**2)
      CHARACTER(LEN=*) SXSEL
#include "warnings.h"
      Character(LEN=16):: ROUTINE='DAVCRE  '
      CHARACTER(LEN=4) IOUTW,IOUTX
      Real*8 :: THRA=1.D-13,THRLD2=5.D-14,THRQ=1.D-07,THRZ=1.D-06,      &
     &          THRLD1=1.D-08
      Real*8, Allocatable:: C1(:), C2(:), X(:)
      Integer iPrLev,nTrial,nCR,ii,i,nDimH,kDimH,nDimH2,iSel,nST,j,ji,  &
     &        iConvQ,iST,iConvA,k,iPass,iConvL,nTotDC,ij,iSTQ,iSTC,jST, &
     &        Length
      Real*8  XMX,XX,SWAP,Ei,QNorm,ENO,ASQ,Ovl,XNorm
      Real*8, External:: DDot_
!
! Local print level (if any)
      IPRLEV=IPRLOC(1)
      IF(IPRLEV.ge.DEBUG) THEN
        WRITE(LF,*)' Entering ',ROUTINE
      END IF
!
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*) 'Super-CI diagonalization. Max iterations: ',ITMAX
      END IF
!
! Set values at loop head and starting guess of the vectors C
!
      Rc_SX = 0
      NTRIAL=NROOT
      NCR=NROOT*NDIM
      CALL FZERO(C,NCR)
      II=0
      DO I=1,NROOT
        C(II+I)=1.0D0
        II=II+NDIM
      END DO
!
!  memory allocation for COVLP
!
       CALL mma_allocate(C1,NSXS,Label='C1')
       CALL mma_allocate(C2,NSXS,Label='C2')
       CALL mma_allocate(X,NSXS,Label='X')
!
! Begin Davidson iterations
! Start by setting up the Davidson HH-matrix
! The dimension of this matrix, NDIMH, is updated in HMAT
!
      NDIMH=0
      ITERSX=0

!*************************************************************
!               Head of Super-CI iteration loop:
!*************************************************************
100   CONTINUE
      ITERSX=ITERSX+1
      CALL HMAT(C,HC,HH,HD,NDIM,NDIMH,NTRIAL)
      KDIMH=(NDIMH+NDIMH**2)/2
!
! Echo the HH-matrix before diagonalizing it.
!
      CALL DCOPY_(KDIMH,HH,1,HH(KDIMH+1),1)
!
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Davidson H-matrix in iteration ',ITERSX
        Write(LF,*)' Davidson H-matrix triangular of size =  ', NDIMH
        Write(LF,'(1x,8F14.6)') (HH(I),I=1,KDIMH)
      END IF
      E(1)=HH(1)
      CC(1)=1.d0
      IF(NDIMH.GT.1) THEN
!
! set eigenvector array to identity before JACO call
!
       NDIMH2=NDIMH**2
       CALL FZERO(CC,NDIMH2)
       II=-NDIMH
       DO I=1,NDIMH
        II=II+NDIMH+1
        CC(II)=1.0d0
       END DO
!
       CALL Jacob(HH(KDIMH+1),CC,NDIMH,NDIMH)
       IF(SXSEL.eq.'LOWEST  ') THEN
! Root selection here assumes picking the lowest root(s).
         CALL JACORD(HH(KDIMH+1),CC,NDIMH,NDIMH)
         II=0
         DO I=1,NROOT
          II=II+I
          E(I)=HH(II+KDIMH)
         END DO
       ELSE
! Root selection by maximum overlap.
         ISEL=1
         XMX=ABS(CC(1))
         DO I=2,NDIMH
          XX=ABS(CC(1+NDIMH*(I-1)))
          IF(XX.GT.XMX) THEN
            ISEL=I
            XMX=XX
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
!
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
!
! Now perform the Davidson update
! First form for each root I, a vector Q(I), where
! Q(I,K)=sum(J,M) CC(I,J,M)*(HC(M,J,K)-E(I)*C(M,J,K))
! here J runs over iterations and M over the number of trial
! vectors in each iterations E(I) is the current estimate of the
! energy of the I:th root, C(M,J) is the M:th trial vector used in
! the J:th iteration, HC being the corresponding sigma vector. CC
! is the matrix of eigenvectors of the Davidson Hamiltonian.
!
! STEP 1: form the matrix CCE=-CC*E (CCE in array SC)
      JI=0
      DO I=1,NROOT
       EI=E(I)
       DO J=1,NDIMH
        JI=JI+1
        S(JI)=-CC(JI)*EI
       END DO
      END DO
!
! Step 2: form the matrix Q=HC*CC
      CALL DGEMM_('N','N',                                              &
     &            NDIM,NROOT,NDIMH,                                     &
     &            1.0d0,HC,NDIM,                                        &
     &            CC,NDIMH,                                             &
     &            0.0d0,Q(NDIM+1),NDIM)
! Step 3: add the contribution C*CCE
      CALL DGEMM_('N','N',NDIM,NROOT,NDIMH,                             &
     & 1.D0,C,NDIM,S,NDIMH,1.D0,Q(NDIM+1),NDIM)
!
! Note that the NDIM first positions in Q are untouched so far
! the vector will be moved later
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' The residual vectors Q of size:', NDIM
        Write(LF,'(1x,8F14.10)')(Q(NDIM+I),I=1,NDIM)
      END IF
!
! Check the norms of the Q-vectors for convergence, and obtain
! bounds to the eigenvalues. E(I) is an upper bound to the I:th
! eigenvalue, from McDonald's theorem, while E(I)-NORM(Q(I))
! is a lower bound to the eigenvalue. This is the Weinstein lower
! bound formula in a multiple root form.
!
      ICONVQ=0
      IST=1+NDIM
      DO I=1,NROOT
       CALL COVLP(Q(IST),Q(IST),DIA,PA,SXN,                             &
     &            C1,C2,X,QQ(I))
       IST=IST+NDIM
      END DO
!
      DO I=1,NROOT
          QNORM=SQRT(QQ(I))
          IF(QNORM.LT.THRQ)ICONVQ = ICONVQ + 1
          EI = E(I)
          ENO = EI - QNORM
          IF(NROOT.GT.1) THEN
            IF(I.GT.1) THEN
              IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(20X,F16.8,A,I2,A,F16.8)') ENO,               &
     &                   ' <  SX energy ',I,'  < ',EI
              END IF
            ELSE
              IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)')ITERSX,      &
     &               ' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
              END IF
            END IF
          ELSE
            IF(IPRLEV.GE.DEBUG) THEN
                Write(LF,'(1X,I2,A,4X,F16.8,A,I2,A,F16.8)')ITERSX,      &
     &               ' SX iteration ',ENO,' <  SX energy ',I,'  < ',EI
            END IF
          END IF
          IF(IPRLEV.GE.DEBUG) THEN
            Write(LF,*)'  Norm of Q:',QNORM
          END IF
          IST = IST + NDIM
      END DO
!PAM00 End of replacement.
!
!    Reset ICONVQ unless all roots are converged
!
      IF(ICONVQ.LT.NROOT) ICONVQ = 0
!
!    Check the expansion vectors for convergence. This is done
!    by computing the sum of the squares of CC(I,J), for each
!    root I, and for the J trial vectors of the current iteration
!
      NST = NDIMH - NTRIAL + 1
      ICONVA = 0
      DO I = 1,NROOT
          ASQ=DDOT_(NTRIAL,CC(NST),1,CC(NST),1)
          IF(ASQ.LT.THRA) ICONVA = ICONVA + 1
          NST = NST + NDIMH
      END DO
!
!    Reset ICONVA unless all roots are converged.
!
      IF(ICONVA.LT.NROOT) ICONVA = 0
!
!    Branch out of iteration loop if convergence has been
!....    achieved on either the energy or the Davidson expansion
!....    vectors.
!
      IF(ICONVQ.NE.0.OR.ICONVA.NE.0) GO TO 35
!
! Calculate the D-vectors as D(I,K)=Q(I,K)/(E(I)-HD(K,K))
!
      IST=1
      DO I=1,NROOT
       EI=E(I)
       DO K=1,NDIM
        SC(K)=EI-HD(K)
        IF(ABS(SC(K)).LT.THRZ) SC(K)=1.0d0
       END DO
       Q(IST:IST+NDIM-1) = Q(IST+NDIM:IST+2*NDIM-1)/SC(1:NDIM)
       IST=IST+NDIM
      END DO
! Remove any unwanted components. These are signalled by
! huge elements of SX hamiltonian diagonal (set in SXHAM).
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
!
! Q now contains the Davidson correction vectors. These will now
! be orthogonalized to all old vectors and then to each other.
!PAM01 The ON is done in two passes. A vector smaller than THRLD
! after the first ON pass is considered to contribute nothing worthwhile
! to the Davidson procedure, and is discarded. Else, it is normalized.
! If it is not still normalized to within an accuracy THRLD after the
! second pass, then again it is discarded: it means that the rescaling
! after pass 1, together with numerical noise amplification in
! COVLP, is beginning to be noticable. The latter will of course happen
! when some orbital rotation(s) are almost redundant.
!
      IPASS=0
      ICONVL=0
      NTOTDC=NROOT

 23   CONTINUE

! First form the overlap matrix
      IJ=0
      ISTQ=1
      DO I=1,NTOTDC
       ISTC=1
       DO J=1,NDIMH
        IJ=IJ+1
        CALL COVLP(Q(ISTQ),C(ISTC),DIA,PA,SXN,                          &
     &             C1,C2,X,OVL)
        S(IJ)=-OVL
        ISTC=ISTC+NDIM
       END DO
       ISTQ=ISTQ+NDIM
      END DO
      CALL DGEMM_('N','N',NDIM,NTOTDC,NDIMH,                            &
     &           1.D0,C,NDIM,S,NDIMH,                                   &
     &           1.D0,Q,NDIM)
!
! The Q-vectors are now orthogonal to all C-vectors.
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)'  D vector orthogonal to all C vectors:'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM)
      END IF
! Now orthogonalize them to one another, rejecting those which
! appear with too small a norm
!
      IST=1
      NTRIAL=0

! Long loop over NTOTDC vectors:
      DO I=1,NTOTDC
       IF(I.NE.1.AND.NTRIAL.NE.0) THEN
!
!  Orthogonalize this vector to the preceding Q's of this iteration
!
        JST=1
        DO J=1,NTRIAL
         CALL COVLP(Q(IST),Q(JST),DIA,PA,SXN,                           &
     &              C1,C2,X,OVL)
         S(J)=-OVL
         JST=JST+NDIM
        END DO
!        CALL DGEMX(NDIM,NTRIAL,1.D0,Q,NDIM,S,1,Q(IST),1)
        CALL DGEMV_('N',NDIM,NTRIAL,1.D0,Q,NDIM,S,1,1.0D0,Q(IST),1)
       IF(IPRLEV.GE.INSANE) THEN
          Write(LF,*)'  Q vector orthogonal to preceding Q vectors:'
          Write(LF,'(1x,8F14.10)')(Q(k),k=1,NDIM)
       END IF
       ENDIF
!
!  Normalize this vector and move to trial set if norm large enough
!
       CALL COVLP(Q(IST),Q(IST),DIA,PA,SXN,                             &
     &            C1,C2,X,XNORM)
! Due to large noise amplification in COVLP, the squared-norm
! can actually come out as a negative number.
! Acceptable, only if it is very close to zero. Else, quit.
       IF(XNORM.LT.-1.0D-09)  then
         Write(LF,*)
         Write(LF,*)'      *** Error in subroutine DAVCRE ***'
         Write(LF,*)' The squared norm of a trial vector has been'
         Write(LF,*)' computed to be negative:'
         Write(LF,*)'      XNORM=',XNORM
         Write(LF,*)' This is possible only for some severe malfunction'
         Write(LF,*)' of the rasscf program. Please issue a bug report.'
         Write(LF,*)
         if (.not. DoNECI) then
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
         Write(LF,'(1X,A,I3,A,I3,A,ES16.8)') 'Pass ',IPASS,             &
     &                 ' New orthogonal vector ',I,' has norm ',XNORM
       END IF

!PAM01 Two different treatments, depending on if this is first or
! second orthonormalization pass:
       IF(IPASS.EQ.0) THEN
! First pass:
         IF(ITERSX.EQ.1 .or. XNORM.GT.THRLD1) THEN
          ISTQ=NTRIAL*NDIM+1
          NTRIAL=NTRIAL+1
          XNORM=1.0D00/(XNORM+1.0D-24)
!PAM01 Note that ISTQ can be (and is!) the same as IST:
          IF(ISTQ.EQ.IST) THEN
           CALL DSCAL_(NDIM,XNORM,Q(IST),1)
          ELSE
           CALL DYAX(NDIM,XNORM,Q(IST),1,Q(ISTQ),1)
          ENDIF
         ENDIF
       ELSE
! Second pass: Demand accurate normalization, else we know that
! poor independence, maybe with strong rounding-error amplification
! in COVLP, has began to erode the orthonormalization of the
! basis vectors, hence the integrity of the Davidson procedure.
         IF(ITERSX.EQ.1 .or. ABS(XNORM-1.0D0).LT.THRLD2) THEN
          ISTQ=NTRIAL*NDIM+1
          NTRIAL=NTRIAL+1
          XNORM=1.0D00/(XNORM+1.0D-24)
!PAM01 Note that ISTQ can be (and is!) the same as IST:
          IF(ISTQ.EQ.IST) THEN
           CALL DSCAL_(NDIM,XNORM,Q(IST),1)
          ELSE
           CALL DYAX(NDIM,XNORM,Q(IST),1,Q(ISTQ),1)
          ENDIF
         ENDIF
       END IF
       IST=IST+NDIM

! End over long loop over NTOTDC vectors I=1..NTOTDC
      END DO
!
! NTRIAL new orthogonal vectors have now been formed. if NTRIAL
! equals zero there are no new linearly independent trial vectors
! and we branch out. If NTRIAL does not equal zero make a second
! orthonormalization pass
!
      IF(NTRIAL.EQ.0) ICONVL=1
      NTOTDC=NTRIAL
      IPASS=IPASS+1
      IF(IPASS.NE.2.AND.NTOTDC.NE.0) GO TO 23
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' The correction vectors, after ON:'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM*NTRIAL)
      END IF
!
! Check if the set of vectors has become linearly dependent
!
      IF(ICONVL.EQ.1) GO TO 36
!
! Move the new orthogonal vectors to C
!
      NST=1+NDIMH*NDIM
      LENGTH=NTRIAL*NDIM
      CALL DCOPY_(LENGTH,Q,1,C(NST),1)
!
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
!
!    At this point, the calculation has neither converged, nor
!    produced linearly dependent trial vectors. Branch back to
!    the head of the iteration loop unless more than ITMAX
!    iterations have been performed.
!
      IF(ITERSX.LT.ITMAX) GO TO 100
!
! Here after ITMAX iterations, and no convergence
!
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
!
! Here as calculation has converged.
!
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
!
! Here as trial vectors have become linearly dependent
!
36    CONTINUE
      IF(IPRLEV.GE.DEBUG) THEN
        Write(LF,*)' Trial vector set has become linearly dependent'
      END IF
34    CONTINUE
!
! Compute CI vectors for all roots
!
! CI NEW (K) =  sum(J,M)CC(I,J M)*C M J (K)
!
! Here J runs over the iterations, and M over the
! number of trial vectors used in each iteration J
!
      CALL DGEMM_('N','N',                                              &
     &            NDIM,NROOT,NDIMH,                                     &
     &            1.0d0,C,NDIM,                                         &
     &            CC,NDIMH,                                             &
     &            0.0d0,Q,NDIM)
      IF(IPRLEV.GE.INSANE) THEN
        Write(LF,*)' Unnormalized final CI vectors (in Q):'
        Write(LF,'(1x,8F14.10)')(Q(I),I=1,NDIM)
      END IF
!
! Normalize the final CI vectors
!
      IST=1
      DO I=1,NROOT
       CALL COVLP(Q(IST),Q(IST),DIA,PA,SXN,                             &
     &            C1,C2,X,XNORM)
       XNORM=1.0D00/XNORM
       CALL DYAX(NDIM,XNORM,Q(IST),1,C(IST),1)
       IST=IST+NDIM
      END DO
!
! Print vector if desired
!
      IF(IPRLEV.GE.INSANE) THEN
       IST = 0
       DO I = 1,NROOT
          Write(LF,*) ' SX-CI vector for root ',I
          Write(LF,'(1X,8F8.4)') (C(IST+K),K = 1,NDIM)
          IST = IST + NDIM
       END DO
      ENDIF
!
! End of diagonalization
! Free memory for COVLP
!
      CALL mma_deallocate(C1)
      CALL mma_deallocate(C2)
      CALL mma_deallocate(X)
      END SUBROUTINE DAVCRE
