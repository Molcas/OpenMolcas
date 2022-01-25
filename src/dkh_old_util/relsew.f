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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************
! $Id: relsew.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de
!   Modified: 2005 Jesper Wisborg Krogh, Jesper.Krogh@teokem.lu.se
!
!
      SUBROUTINE SCFCLI(idbg,epsilon,S,H,V,PVP,N,ISIZE,VELIT,           &
     &                  TMP1,TMP2,TMP3,TMPA,TMPB,TMPC,EW,E,AA,RR,TT,    &
     &                  TMP4,TMPD,TMPE,TMPF,TWRK4,IDIM)
      USE DKH_Info, ONLY: IRELAE
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(ISIZE),H(ISIZE),V(ISIZE),PVP(ISIZE)
      DIMENSION TMP1(IDIM*(IDIM+1)/2),TMP2(ISIZE),TMP3(ISIZE),          &
     &          TMPA(N,N),TMPB(N,N),TMPC(N,N),                          &
     &          EW(N),E(N),AA(N),RR(N),TT(N)
      DIMENSION TMP4(IDIM*(IDIM+1)/2)
      DIMENSION TMPD(IDIM,IDIM)
      DIMENSION TMPE(N,N)
      DIMENSION TMPF(N,N)
      DIMENSION TWRK4(N*200)
!
!
!      write(6,*) ' in SCFCLI', N, iSize
      PREA=1d0/(VELIT*VELIT)
      CON2=PREA+PREA
      CON=1.D0/PREA
#ifdef _DEBUGPRINT_
!
!     CALCULATE DETERMINANT
!
      Call Square(S,TMPB,n,1,n)
!      do i=1,n
!         write(6,'(5f10.5)') (TMPB(i,j),j=1,n)
!      enddo
      icontr=-1
      dtol=1.D-14
      CALL dcopiv(TMPB,TMPB,n,1,n,dtol,det,iex,icontr,TMP2)
      if(idbg.gt.0)WRITE (idbg,2016) icontr,det,iex
2016  FORMAT('  relsew| DCOPIV rc=',I2,', |S|=',D20.6,'x 10**(',I4,') ')
      IF (icontr.NE.0) THEN
         WRITE (6,2016) icontr,det,iex
         WRITE (6,2012) dtol
2012     FORMAT('  relsew|****** '/,                                    &
     & '        |****** WARNING - OVERLAP MATRIX SINGULAR '/,           &
     & '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,    &
     & '        |******'//)
         CALL errex_rel(' relsew| singular overlap matrix')
      ENDIF
#endif
!
!     SCHMIDT-ORTHOGONALIZE
!
      CALL Sogr(idbg,N,S,TMPE,TMP2,TMPC,EW)
!
! ** TMPE CONTAINS TRANSFORMATION TO ORTHOGONAL AO-BASIS
!--------------------------------------------------------------------
!     NON-RELATIVISTIC KINETIC ENERGY
!--------------------------------------------------------------------
      CALL Diagr(H,N,TMPF,EW,TMPE,TMPB,TMP2)
      if(idbg.gt.0)WRITE (idbg,556)
  556 FORMAT(//,7X,'- NREL. ENERG.  -  DIVIDED BY C ',                  &
     &             '- REL.  ENERG.  -  MOMENTUM    -',                  &
     &             ' TERMS OF POWER SERIES (LOW ENERGY ONLY)'//)
      DO 4 I=1,N
         IF (ew(i).LT.0.D0) THEN
             WRITE (6,*) ' scfcli| ew(',i,') = ',ew(i)
         CALL errex_rel('kinetic energy eigenvalue less than zero')
         ENDIF
!
!     IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO
!
         RATIO=EW(I)/VELIT
!
!     CALCULATE RELATIVISTIC ENERGY AND MOMENTUM
!
         SR=sqrt(2.D0*EW(I))
         TT(I)=EW(I)
         IF (RATIO.GT.0.02D0) GOTO 11
         TV1=EW(I)
         TV2=-TV1*EW(I)*PREA/2.D0
         TV3=-TV2*EW(I)*PREA
         TV4=-TV3*EW(I)*PREA*1.25D0
         EW(I)=TV1+TV2+TV3+TV4
         if(idbg.gt.0)WRITE (idbg,100) I,TV1,RATIO,EW(I),SR,TV2,TV3,TV4
  100    FORMAT(1X,I4,7(1X,D15.8))
         GOTO 12
  11     TV1=EW(I)
         EW(I)=CON*(sqrt(1.D0+CON2*EW(I))-1.D0)
         if(idbg.gt.0)WRITE (idbg,100) I,TV1,RATIO,EW(I),SR
  12     CONTINUE
         E(I)=EW(I)+CON
4     CONTINUE
!---------------------------------------------------------------------
!     CALCULATE REVERSE TRANSFORMATION
!---------------------------------------------------------------------
      Call DGEMM_('N','N',n,n,n,1.0D0,TMPE,n,TMPF,n,0.0D0,TMPB,n)
#ifdef MOLPRO
      Call Square(TMPC,S,N,N)
#else
      Call Square(S,TMPC,N,1,N)
#endif
      Call DGEMM_('N','N',n,n,n,1.0D0,TMPC,n,TMPB,n,0.0D0,TMPA,n)
! ** TMPC  CONTAINS OVERLAP MATRIX IN FULL
!
      IF(IRELAE.NE.21.AND.IRELAE.NE.22.AND.IRELAE.NE.23) THEN
!
         Call dCopy_(iSize,[0.0D0],0,H,1)
         Do K = 1,N
            IJ = 0
            Do I = 1,N
               Do J = 1,I
                  IJ = IJ + 1
                  H(IJ) = H(IJ) + TMPA(I,K)*TMPA(J,K)*EW(K)
               End Do
            End Do
         End Do
!
      ELSE
         Call dCopy_(iSize,[0.0D0],0,H,1)
      ENDIF
!
!     CALCULATE KINEMATICAL FACTORS
!
      IF(IRELAE.NE.11) THEN
!
         DO 362 I=1,N
            AA(I)=sqrt((CON+E(I)) / (2.D0*E(I)))
            RR(I)=sqrt(CON)/(CON+E(I))
362      CONTINUE
!
      ELSE IF(IRELAE.EQ.11) THEN
!
         DO I=1,N
            AA(I)=(sqrt(1.0D0+CON*TT(I)*2.0D0/((CON+E(I))*(CON+E(I))))) &
     &                                /(CON+E(I))
! O OPERATOR
            RR(I)=sqrt(CON)/(CON+E(I))
! Q OPERATOR
         ENDDO
!
      ENDIF
!
!     POTENTIAL
!
!
!     BEYOND THIS POINT, TMPC IS USED AS SCRATCH ARRAY
!
!
!    TRANSFORM V TO T-BASIS
!
      CALL TrSmr2(V,TMPE,TMP3,N,TMPB,TMPF,TMPC)
!ulf
      if(idbg.gt.0)CALL PRMAT(IDBG,V,N,0,'v oper  ')
!
!    MULTIPLY
!
      IF(IRELAE.EQ.0.OR.IRELAE.EQ.1.OR.IRELAE.EQ.2.OR.                  &
     &   IRELAE.EQ.3) THEN
!
         IJ=0
         DO 2005 I=1,N
            DO 2006 J=1,I
               IJ=IJ+1
               TMP2(IJ)=TMP3(IJ)
               TMP3(IJ)=TMP3(IJ)*AA(I)*AA(J)
2006        CONTINUE
2005     CONTINUE
         IF (IRELAE .eq. 3) Then
#ifdef MOLPRO
           Call Square(TMPD,TMP3,N,N)
#else
           Call Square(TMP3,TMPD,N,1,N)
#endif
         ENDIF
!
      ELSE IF(IRELAE.EQ.11) THEN
!
         IJ=0
         DO I=1,N
             DO J=1,I
                  IJ=IJ+1
                  TMP2(IJ)=TMP3(IJ)
                  TMP3(IJ)=VELIT*TMP3(IJ)*(sqrt(RR(I)*RR(J))            &
     &                      *AA(I)/AA(J)+sqrt(RR(J)*RR(I))*AA(J)/AA(I))
             ENDDO
         ENDDO
!
      ELSE IF(IRELAE.EQ.22) THEN   ! ZORA(FP)
!
         IJ=0
         DO I=1,N
            DO J=1,I
               IJ=IJ+1
               TMP3(IJ)=TMP3(IJ)*AA(I)*AA(J)
            ENDDO
         ENDDO
!
      ENDIF
!
      CALL TrSmtr(TMP3,TMPA,H,1.0D0,N,TMPB,TMPC)
!ulf
!
!     PVP INTEGRALS AND TRANSFORM THEM TO T-BASIS
!
      if(idbg.gt.0)CALL PRMAT(IDBG,pvp,N,0,'raw pvp integrals  ')
      CALL TrSmr2(PVP,TMPE,PVP,N,TMPB,TMPF,TMPC)
!
!    MULTIPLY
!
      IF(IRELAE.EQ.0.OR.IRELAE.EQ.1.OR.IRELAE.EQ.2.OR.                  &
     &   IRELAE.EQ.3) THEN
!
         IJ=0
         DO 3005 I=1,N
            DO 3006 J=1,I
               IJ=IJ+1
               TMP3(IJ)=PVP(IJ)
               PVP(IJ)=PVP(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
               IF (IRELAE .eq. 3) Then
                  TMPD(I,J)=TMPD(I,J)+PVP(IJ)
                  TMPD(J,I)=TMPD(I,J)
               End If
3006        CONTINUE
3005     CONTINUE
!
      ELSE IF(IRELAE.EQ.11) THEN
!
         IJ=0
         DO I=1,N
            DO J=1,I
               IJ=IJ+1
               TMP3(IJ)=PVP(IJ)
               PVP(IJ)=PVP(IJ)*(RR(I)*RR(J)*AA(I)/AA(J)                 &
     &                      +RR(J)*RR(I)*AA(J)/AA(I))*0.5D0
            ENDDO
         ENDDO
!
      ELSE IF(IRELAE.EQ.21.OR.IRELAE.EQ.22.OR.IRELAE.EQ.23) THEN
!
         IJ=0
         DO I=1,N
            DO J=1,I
               IJ=IJ+1
               PVP(IJ)=-PVP(IJ)*0.5D0*PREA
               IF(I.EQ.J) PVP(IJ)=PVP(IJ)+TT(I)*2.0D0
            ENDDO
         ENDDO
#ifdef MOLPRO
         Call Square(TMPD,PVP,N,N)
#else
         Call Square(PVP,TMPD,N,1,N)
#endif
!
!     ----- inverse operator -----
!
         EPS=1.d-13
         CALL MINVD(TMPD,N,N,EPS,ILL)
         IJ=0
         DO I=1,N
            DO J=1,I
               IJ=IJ+1
               PVP(IJ)=TMPD(I,J)
               PVP(IJ)=PVP(IJ)*TT(I)*TT(J)*2.0D0
               IF(IRELAE.EQ.22) PVP(IJ)=PVP(IJ)*AA(I)*AA(J)
            ENDDO
         ENDDO
!
      ENDIF
!
      CALL TrSmtr(PVP,TMPA,H,1.0D0,N,TMPB,TMPC)
!ulf
!
      IF(IRELAE.EQ.23) THEN
         IJ=0
#ifdef MOLPRO
         Call Square(TMPE,PVP,N,N)
#else
         Call Square(PVP,TMPE,N,1,N)
#endif
         DO I=1,N
            DO J=1,I
               IJ=IJ+1
               TMPD(I,J)=PVP(IJ)/TT(J)*0.5D0*PREA
               TMPD(J,I)=PVP(IJ)/TT(I)*0.5D0*PREA
            ENDDO
         ENDDO
         M=N
         Call dCopy_(N*N,[0.0D0],0,TMPC,1)
         Call dCopy_(N,[1.0D0],0,TMPC,N+1)
         Call DGEMM_('N','N',N,N,N,1.0D0,TMPD,M,TMPE,M,1.0D0,TMPC,M)
! ----- modified overlap is incorporated into PVP
         Call DGEMM_('N','N',N,N,N,1.0D0,TMPA,N,TMPC,N,0.0D0,TMPB,N)
         Call dGemm_tri('N','T',N,N,N,1.0D0,TMPB,N,TMPA,N,0.0D0,PVP,N)
      ENDIF
      IF(IRELAE.EQ. 1.OR.IRELAE.EQ.11.OR.                               &
     &   IRELAE.EQ.21.OR.IRELAE.EQ.22.OR.IRELAE.EQ.23)                  &
     &   GOTO 1000
!
      If (IRELAE .eq. 3) Then
! --- KEEP T-BASIS VEXT INTO TMP1 FOR HIGHER-ORDER DK
         Call dCopy_(N*(N+1)/2,TMP2,1,TMP1,1)
! --- KEEP T-BASIS PVP INTO TMP4 FOR HIGHER-ORDER DK
         Call dCopy_(N*(N+1)/2,TMP3,1,TMP4,1)
      End If
!
!     CALCULATE Even2r OPERATOR
!
      CALL Even2r(idbg,N,TMP2,TMP3,E,AA,RR,TT,TMPE,TMPB,TMPC,           &
     &            TMPF)
!
!    TRANSFORM BACK
!
!ulf
      if(idbg.gt.0)CALL PRMAT(IDBG,TMP3,n,0,'ev2 orig')
      CALL TrSmtr(TMP3,TMPA,H,1.0D0,N,TMPB,TMPC)
!ulf
      IF(IRELAE.EQ.0.OR.IRELAE.EQ.2) GOTO 1000   ! DK2
!
!     ----- CALCULATE Even3r OPERATOR -----
!
      CALL Even3r(idbg,N,TMP2,TMP3,E,AA,RR,TT,TMPB,                     &
     &            TMPD,TMP1,TMP4,TMPE,TMPF)
!
!    TRANSFORM BACK
!
!ulf
      if(idbg.gt.0)CALL PRMAT(IDBG,TMP3,n,0,'ev2 orig')
      CALL TrSmtr(TMP3,TMPA,H,1.0D0,N,TMPB,TMPC)
!ulf
!
      IF(IRELAE.EQ.3) GOTO 1000   ! DK3
!
!     More to come!
!
 1000 CONTINUE
!
!ulf
      if(idbg.gt.0)CALL PRMAT(IDBG,h,n,0,'h   oper')
      CALL Sogr(idbg,N,S,TMPE,TMP2,TMPC,EW)
      CALL Diagr(H,N,TMPF,EW,TMPE,TMPB,TMP2)
      if(idbg.gt.0)CALL PRMAT(IDBG,h,n,0,'h   oper(final)')
      if(idbg.gt.0)WRITE (idbg,*) '--- EIGENVALUES OF H MATRIX ---'
      if(idbg.gt.0)WRITE (idbg,'(4D20.12)') EW
      RETURN
! Avoid unused argument warnings
      IF (.FALSE.) THEN
        CALL Unused_real(epsilon)
        CALL Unused_real_array(TWRK4)
      END IF
      END
      SUBROUTINE TrSmr(A,B,C,N,H,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N*(N+1)/2),B(N,N),C(N*(N+1)/2),H(N,N),W(N,N)
!
!     TRANSFORM SYMMETRIC MATRIX A BY UNITARY TRANSFORMATION
!     IN B. RESULT IS IN C
!
#ifdef MOLPRO
      Call Square(W,A,n,n)
#else
      Call Square(A,W,n,1,n)
#endif
      Call DGEMM_('T','N',n,n,n,1.0D0,B,n,W,n,0.0D0,H,n)
      Call dGemm_tri('N','N',n,n,n,1.0D0,H,n,B,n,0.0D0,C,n)
      RETURN
      END
!
      SUBROUTINE TrSmr2(A,B,C,N,H,G,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N*(N+2)/2),B(N,N),C(N*(N+1)/2),H(N,N),W(N,N),G(N,N)
!
!     Performs the equivalent operations as two calls to TrSmr
!     RESULT IS IN C = G^T * (B^T * A * B) * G
!     and where C and A are triangular packed.
!
#ifdef MOLPRO
      Call Square(W,A,n,n)
#else
      Call Square(A,W,n,1,n)
#endif
      Call DGEMM_('T','N',n,n,n,1.0D0,B,n,W,n,0.0D0,H,n)
      Call DGEMM_('N','N',n,n,n,1.0D0,H,n,B,n,0.0D0,W,n)
      Call DGEMM_('T','N',n,n,n,1.0D0,G,n,W,n,0.0D0,H,n)
      Call dGemm_tri('N','N',n,n,n,1.0D0,H,n,G,n,0.0D0,C,n)
      RETURN
      END
!
      SUBROUTINE TrSmtr(A,B,C,FACTOR,N,H,W)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(N*(N+1)/2),B(N,N),C(N*(N+1)/2),H(N,N),W(N,N)
!
!     B*A*BT
!
#ifdef MOLPRO
      Call Square(W,A,n,n)
#else
      Call Square(A,W,n,1,n)
#endif
      Call DGEMM_('N','N',n,n,n,1.0D0,B,n,W,n,0.0D0,H,n)
      Call dGemm_tri('N','T',n,n,n,1.0D0,H,n,B,n,Factor,C,n)
      RETURN
      END
      SUBROUTINE Sogr(idbg,N,SS,SINV,P,G,A1)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION SS(*),P(*),G(*),A1(*),SINV(N,N)
!
!     SUBROUTINE TO CALCULATE TRANSFORMATION TO SCHMIDT-
!     ORTHOGONALIZED BASIS.
!     N              DIMENSION OF MATRICES. ISIZE=N*(N+1)/2
!     SS(ISIZE)      ORIGINAL OVERLAP MATRIX (LOWER TRIANGULAR)
!                    WILL NOT BE DESTROYED
!     P (ISIZE)      OUTPUT TRANSFORMATION MATRIX
!     G (ISIZE)      SCRATCH
!     A1(N)          SCRATCH
!
      INTEGER ierr
      If (iDbg.gt.0) Call PrMat(idbg,SS,n,0,'SS')
      ierr=0
      JL=0
      IQ=0
      DO 349 J=1,N
         IL=JL
         JQ=IQ
         S1KK=SS(IQ+J)
         G(IL+J)=1.D0
         IF(J.EQ.1)GO TO 341
         J1=J-1
         JL=0
         DO 342 K=1,J1
            LG=JQ
            ETOT=0.D0
            DO 343 L=1,K
               LG=LG+1
               JL=JL+1
               ETOT=ETOT+SS(LG)*G(JL)
  343       CONTINUE
            S1KK=S1KK-ETOT*ETOT
            A1(K)=ETOT
  342    CONTINUE
         JF=1
         JL=IL
         DO 344 K=1,J1
            SUM=0.D0
            JL=JL+1
            JF=JF+K-1
            IH=JF
            DO 345 L=K,J1
               IH=IH+L-1
               SUM=SUM+A1(L)*G(IH)
  345       CONTINUE
            G(JL)=-SUM
  344    CONTINUE
  341    CONTINUE
         IF (s1kk .LE. 1.D-16) THEN
            WRITE (6,*) '    Sogr| j=',j,' s1kk=',s1kk
            ierr=ierr+1
         ENDIF
         S1KK=1.D0/sqrt(S1KK)
         JL=IL
         DO 340 K=1,J
            JL=JL+1
            IQ=IQ+1
            G(JL)=G(JL)*S1KK
            P(IQ)=G(JL)
  340    CONTINUE
  349 CONTINUE
      IJ=0
      DO 1 I=1,N
         DO 2 J=1,I
            IJ=IJ+1
            SINV(I,J)=0.D0
            SINV(J,I)=P(IJ)
2        CONTINUE
1     CONTINUE
      IF (ierr.GT.0)  CALL errex_rel('function has negative norm')
      If (iDbg.gt.0) Call PrMat(idbg,P,n,0,'P')
      RETURN
      END
!***** NAME AddMar
      SUBROUTINE AddMar(N,S,OVE)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION S(N),OVE(N)
      DO 1 I=1,N
         OVE(I)=OVE(I)+S(I)
1     CONTINUE
      RETURN
      END
!
      SUBROUTINE Diagr(A,N,EIG,EW,SINV,AUX,AUXI)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(*),AUX(N,N),SINV(N,N),EIG(N,N),EW(*),AUXI(*)
#include "WrkSpc.fh"
      If (n.eq.0) Return
#ifdef MOLPRO
      Call Square(Aux,A,n,n)
#else
      Call Square(A,Aux,n,1,n)
#endif
      Call DGEMM_('N','N',n,n,n,1.0D0,Aux,n,Sinv,n,0.0D0,Eig,n)
      Call dGemm_tri('T','N',n,n,n,1.0D0,SINV,n,EIG,n,0.0D0,AUXI,n)
!
      Call dCopy_(n*n,[0.0D0],0,Eig,1)
      Call dCopy_(n,[1.0D0],0,Eig,n+1)
      Call dCopy_(N*(N+1)/2,AUXI,1,AUX,1)
!      Call NIDiag(AUXI,EIG,N,N)
      Call NIDiag_New(AUXI,EIG,N,N)
      Call vEig(N,AUXI,EW)
      Call JacOrd2(EW, Eig, n, n)
!
      RETURN
      END
      SUBROUTINE Even2r(idbg,N,V,G,E,A,R,TT,AUXF,AUXG,AUXH,W1W1)

!  DELC$
!     EVEN2 - BERND HESS - V 1.0 - 5.2.86
!     CALCULATE EVEN2 OPERATORS
!
!
!     N       DIMENSION OF MATRICES
!     V       POTENTIAL MATRIX
!     G       MATRIX OF PVP OPERATOR. WILL CONTAIN EVEN2 OPERATORS
!             ON OUTPUT
!     E       RELATIVISTIC ENERGY (DIAGONAL)
!     A       A-FACTORS (DIAGONAL)
!     R       R-FACTORS (DIAGONAL)
!     TT      NONREL. KINETIC ENERGY (DIAGONAL)
!     AUXF,AUXG,AUXH  SCRATCH ARAYS
!
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(N*(N+1)/2),G(N*(N+1)/2),E(N),R(N),A(N),TT(N),         &
     &          AUXF(N,N),AUXG(N,N),AUXH(N,N)
      DIMENSION W1W1(N,N)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL PRMAT(IDBG,V,N,0,'V       ')
      if(idbg.gt.0)CALL PRMAT(IDBG,G,N,0,'G       ')
      if(idbg.gt.0)CALL PRMAT(IDBG,E,N,1,'E       ')
      if(idbg.gt.0)CALL PRMAT(IDBG,A,N,1,'A       ')
      if(idbg.gt.0)CALL PRMAT(IDBG,R,N,1,'R       ')
      if(idbg.gt.0)CALL PRMAT(IDBG,TT,N,1,'TT      ')
#endif
      M=N
      Call dCopy_(N*N,[0.0D0],0,AUXH,1)
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            V(IJ)=V(IJ)/(E(I)+E(J))
            G(IJ)=G(IJ)/(E(I)+E(J))
         End Do
      End Do
      DO J=1,N
         IJ = J*(J-1)/2 + 1
         DO I=J,N
            IJ = IJ + I-1
            AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*A(J)
            AUXG(I,J)=R(I)*V(IJ)*A(J)
         End Do
      End Do
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*A(I)
            AUXG(J,I)=R(J)*V(IJ)*A(I)
         End Do
      End Do
!
!     ARQA ARQA
!
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
#ifdef _DEBUGPRINT_
      IF (IE.NE.0) Call SysHalt('relsew')
!ulf
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   1',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
         End Do
      End Do
!
!     ARQA AQRA
!
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   2',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(J,I)=A(J)*V(IJ)*A(I)*A(I)*R(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXF(I,J)=A(I)*V(IJ)*A(J)*A(J)*R(J)
         End Do
      End Do
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
         End Do
      End Do
!
!     AQRA ARQA
!
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   3',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=G(IJ)*A(I)*R(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=G(IJ)*A(J)*R(J)
         End Do
      End Do
!
!     AQRA AQRA
!
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!
!     KEEP W1*W1 FOR HIGHER-ORDER DK
!
      Call dCopy_(N*N,AuxH,1,W1W1,1)
!
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(IDBG,'W*W     ',AUXH,N)
#endif
!
!     1/2 EW*W + 1/2 W*WE
!
      DO 610 I=1,N
         DO 611 J=1,N
            AUXH(I,J)=0.5D0*( AUXH(I,J)*E(I) + AUXH(I,J)*E(J) )
611      CONTINUE
610   CONTINUE
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH SYM',auxh,n)
#endif
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)*E(I)*A(I)
            AUXG(J,I)=R(J)*V(IJ)*A(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)*E(J)*A(J)
            AUXG(I,J)=R(I)*V(IJ)*A(J)
         End Do
      End Do
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   5',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=-0.5D0/TT(J)*G(IJ)*A(I)*R(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=-0.5D0/TT(I)*G(IJ)*A(J)*R(J)
         End Do
      End Do
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   6',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(J,I)=A(J)*V(IJ)*R(I)*A(I)*E(I)*A(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXF(I,J)=A(I)*V(IJ)*R(J)*A(J)*E(J)*A(J)
         End Do
      End Do
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=-2.D0*TT(J)*R(J)*V(IJ)*A(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=-2.D0*TT(I)*R(I)*V(IJ)*A(J)
         End Do
      End Do
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   7',auxh,n)
#endif
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXG(J,I)=G(IJ)*A(I)*R(I)
         End Do
      End Do
      DO J=1,N
         IJ = INT(DBLE(J*(J-1))*0.5D0 + 1D0)
         DO I=J,N
            IJ = IJ + I-1
            AUXG(I,J)=G(IJ)*A(J)*R(J)
         End Do
      End Do
      CALL CpLabr(AUXF,AUXG,N,N,N,M,M,AUXH,M,IE)
!ulf
#ifdef _DEBUGPRINT_
      if(idbg.gt.0)CALL prsq(idbg,'AUXH   8',auxh,n)
#endif
!
!     SYMMETRISIEREN
!
      IJ=0
      DO 430 I=1,N
         DO 431 J=1,I
            IJ=IJ+1
            G(IJ)=-0.5D0*(AUXH(I,J)+AUXH(J,I))
 431     CONTINUE
 430  CONTINUE
!
!CC   CALL PRM('OUTPUT  ',G,N)
      RETURN
#ifndef _DEBUGPRINT_
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(idbg)
#endif
      END
      SUBROUTINE CpLabr (A,B,L,M,N,IA,IB,C,IC,IER)
      implicit real*8(a-h,o-z)
      REAL*8   A(IA,M),B(IB,N),C(IC,N)
#ifdef _DEBUGPRINT_
      IF (IA .GE. L .AND. IB .GE. M .AND. IC .GE. L) GO TO 5
      IER=129
      GO TO 9000
    5 CONTINUE
#endif
      IER = 0
      Call DGEMM_('N','N',L,M,N,1.0D0,A,IA,B,IB,1.0D0,C,IC)
#ifdef _DEBUGPRINT_
 9000 CONTINUE
#endif
      RETURN
      END
      subroutine errex_rel(char)
      character*(*) char
      write(6,'(a50)') char
      Call Abend
      end
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE CZERO2(XX,N,M,NDIM)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION XX(NDIM,*)
      DATA ZERO /0.0D+00/
!
      DO J=1,N
         DO I=1,M
            XX(I,J)=ZERO
         ENDDO
      ENDDO
!
      RETURN
      END
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE Even3r(idbg,N,V,G,E,A,R,TT,AUXF,                       &
     &                  AUXH,VEXTT,PVPT,RE1R,W1W1)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION V(N*(N+1)/2),G(N*(N+1)/2),E(N),R(N),A(N),TT(N)
      DIMENSION AUXF(N,N)                          ! Scratch
      DIMENSION AUXH(N,N)                          ! Input/Scratch
      DIMENSION RE1R(N,N)                          ! Scratch
      DIMENSION VEXTT(N*(N+1)/2),PVPT(N*(N+1)/2)   ! Input
      DIMENSION W1W1(N,N)                          ! Input/Scratch
!
!     ----- CONSTRUCT RE1R -----
!
!     write(6,*) 'Hello from Even3R'
      Call xflush(6)
      CALL DKRE1R(A,R,E,TT,V,G,RE1R,VEXTT,PVPT,N)
!
      M=N
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            V(IJ)=VEXTT(IJ)/(E(I)+E(J))
            G(IJ)=PVPT(IJ)/(E(I)+E(J))
         ENDDO
      ENDDO
!
!     ----- W1*W1*E1 -----
!
!     ------- 1/2 E1*W1*W1 + 1/2 W1*W1*E1
!
      DO J=1,N
         DO I=1,N
            AUXF(I,J)=0.5D0*AUXH(I,J)
         ENDDO
      ENDDO
      Call dCopy_(N*N,[0.0D0],0,AUXH,1)
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
      CALL CpLabr(AUXF,W1W1,N,N,N,M,M,AUXH,M,IE)
!
!     ----- W1*E1*W1 TERM -----
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5D0
            AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5D0
         ENDDO
      ENDDO
      CALL CZERO2(W1W1,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,W1W1,M,IE)
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=-A(I)*V(IJ)*A(J)
            AUXF(J,I)=-A(J)*V(IJ)*A(I)
         ENDDO
      ENDDO
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)*V(IJ)*A(J)
            AUXF(J,I)=A(J)*V(IJ)*A(I)
         ENDDO
      ENDDO
      CALL CZERO2(W1W1,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,W1W1,M,IE)
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=-A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5D0
            AUXF(J,I)=-A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5D0
         ENDDO
      ENDDO
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)*V(IJ)*A(J)
            AUXF(J,I)=A(J)*V(IJ)*A(I)
         ENDDO
      ENDDO
      CALL CZERO2(W1W1,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,W1W1,M,IE)
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)*V(IJ)*A(J)
            AUXF(J,I)=A(J)*V(IJ)*A(I)
         ENDDO
      ENDDO
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)*R(I)*G(IJ)*A(J)/R(J)/TT(J)*0.5D0
            AUXF(J,I)=A(J)*R(J)*G(IJ)*A(I)/R(I)/TT(I)*0.5D0
         ENDDO
      ENDDO
      CALL CZERO2(W1W1,N,N,N)
      CALL CpLabr(AUXF,RE1R,N,N,N,M,M,W1W1,M,IE)
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            AUXF(I,J)=A(I)/R(I)*G(IJ)*A(J)*R(J)/TT(I)*0.5D0
            AUXF(J,I)=A(J)/R(J)*G(IJ)*A(I)*R(I)/TT(J)*0.5D0
         ENDDO
      ENDDO
      CALL CpLabr(W1W1,AUXF,N,N,N,M,M,AUXH,M,IE)
!
!     ----- SYMMETRISIEREN -----
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            G(IJ)=0.5D0*(AUXH(I,J)+AUXH(J,I))
         ENDDO
      ENDDO
!
      RETURN
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_integer(idbg)
      END
!
!
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE DKRE1R(A,R,E,TT,V,G,RE1R,VEXTT,PVPT,N)
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!
!     ----- CONSTRUCT RE1R FOR DK3 -----
!
      IMPLICIT REAL*8(A-H,O-Z)
!
      DIMENSION E(N),A(N),R(N),TT(N)
      DIMENSION V(N*(N+1)/2),G(N*(N+1)/2)
      DIMENSION VEXTT(N*(N+1)/2),PVPT(N*(N+1)/2)
      DIMENSION RE1R(N,N)
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            V(IJ)=VEXTT(IJ)
            G(IJ)=PVPT(IJ)
         ENDDO
      ENDDO
!
!     ------- MULTIPLY
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            V(IJ)=V(IJ)*A(I)*A(J)*R(I)*R(I)*R(J)*R(J)*TT(I)*TT(J)*4.0D0
            RE1R(I,J)=V(IJ)
            RE1R(J,I)=RE1R(I,J)
         ENDDO
      ENDDO
!
      IJ=0
      DO I=1,N
         DO J=1,I
            IJ=IJ+1
            G(IJ)=G(IJ)*A(I)*A(J)*R(I)*R(J)
            RE1R(I,J)=RE1R(I,J)+G(IJ)
            RE1R(J,I)=RE1R(I,J)
         ENDDO
      ENDDO
!
      RETURN
! Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real_array(E)
      END
!
!#NUMPAC#MINVD               REVISED ON 1984-11-30
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      SUBROUTINE MINVD(A,KA,N,EPS,ILL)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(KA,N)
      INTEGER MX(1000)
      IF(N.LT.1.OR.N.GT.1000.OR.N.GT.KA.OR.EPS.LE.0.) GO TO 250
!-----LU DECOMPOSITION--------------------------------------------------
      JM1=0  ! dummy initialize
      M  =0  ! dummy initialize
      NM1=N-1
      DO 90 J=1,N
         IF(J.EQ.1) GO TO 30
         JM1=J-1
         DO 20 I=1,JM1
            M=MX(I)
            S=A(M,J)
            A(M,J)=A(I,J)
            IF(I.EQ.1) GO TO 21
            IM1=I-1
            DO 10 K=1,IM1
               S=A(I,K)*A(K,J)+S
   10       CONTINUE
   21       A(I,J)=S
   20    CONTINUE
   30    AM=0.
         DO 60 I=J,N
            S=A(I,J)
            IF(J.EQ.1) GO TO 50
            DO 40 K=1,JM1
               S=A(I,K)*A(K,J)+S
   40       CONTINUE
            A(I,J)=S
   50       AA=abs(S)
            IF(AA.LE.AM) GO TO 60
            AM=AA
            M=I
   60    CONTINUE
         IF(AM.LT.EPS) GO TO 240
         MX(J)=M
         IF(M.EQ.J) GO TO 80
         DO 70 K=1,J
            W=A(M,K)
            A(M,K)=A(J,K)
            A(J,K)=W
   70    CONTINUE
   80    IF(J.EQ.N) GO TO 100
         JP1=J+1
         W=-A(J,J)
         DO 91 I=JP1,N
            A(I,J)=A(I,J)/W
   91    CONTINUE
   90 CONTINUE
  100 IF(N.LE.2) GO TO 130
!-----INPLACE INVERSION OF L-COMPONENT----------------------------------
      DO 120 I=3,N
         IM1=I-1
         IM2=I-2
         DO 121 J=1,IM2
            S=A(I,J)
            JP1=J+1
            DO 110 K=JP1,IM1
                S=A(I,K)*A(K,J)+S
  110       CONTINUE
            A(I,J)=S
  121    CONTINUE
  120 CONTINUE
!-----INPLACE INVERSION OF U-COMPONENT----------------------------------
  130 A(1,1)=1./A(1,1)
      IF(N.EQ.1) GO TO 230
      DO 150 J=2,N
         A(J,J)=1./A(J,J)
         P=-A(J,J)
         JM1=J-1
         DO 151 I=1,JM1
            S=0.
            DO 140 K=I,JM1
               S=A(I,K)*A(K,J)+S
  140       CONTINUE
            A(I,J)=S*P
  151    CONTINUE
  150 CONTINUE
!-----INPLACE MULTIPLICATION OF L AND U COMPONENT-----------------------
      DO 190 J=1,NM1
         JP1=J+1
         DO 170 I=1,J
            S=A(I,J)
            DO 160 K=JP1,N
               S=A(I,K)*A(K,J)+S
  160       CONTINUE
            A(I,J)=S
  170    CONTINUE
         DO 191 I=JP1,N
            S=0.
            DO 180 K=I,N
                S=A(I,K)*A(K,J)+S
  180       CONTINUE
            A(I,J)=S
  191    CONTINUE
  190 CONTINUE
!------INTERCHANGE OF COLUMNS-------------------------------------------
      J=NM1
  200 M=MX(J)
      IF(M.EQ.J) GO TO 220
      DO 210 I=1,N
         W=A(I,M)
         A(I,M)=A(I,J)
         A(I,J)=W
  210 CONTINUE
  220 J=J-1
      IF(J.GE.1) GO TO 200
  230 ILL=0
      RETURN
  240 ILL=J
      RETURN
  250 ILL=30000
      RETURN
      END
