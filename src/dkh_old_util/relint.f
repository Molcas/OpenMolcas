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
      SUBROUTINE BESSKA(AVAL,XVAL,KA,KA1)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 AVAL,XVAL,EPS,KA,KA1
C
C     SUBROUTINE TO CALCULATE MODIFIED BESSEL FUNCTION OF THE
C     THIRD KIND. FORTRAN SUBPROGRAM ADAPTED FROM ALGOL PROGRAM
C     IN: N.M.TEMME, J.COMP.PHYS 19 (1975) 324
C
C     MODIFICATIONS: RESULT IS EXP(X) TIMES BESSEL FUNCTION
C                    INSTEAD OF BESSEL FUNCTION ITSELF.
C     RELATIVE ACCURACY FIXED TO 5*1.D-14
C
      REAL*8 A1,B,C,D,E,F,G,H,P,PI,Q,S,A,X,RCPG,SINHX
culf      INTEGER*4 N,NA
      INTEGER N,NA
      LOGICAL REC,REV
*
      NA=0
      A=AVAL
      X=XVAL
      EPS=5.D-14
      PI=4.0D0*ATAN(1.0D0)
      REV=A.LT.-5D0
      IF (REV) A=-A-1.D0
      REC=A.GE.0.5D0
      IF (.NOT. REC) GOTO 1
      NA=INT(A+0.5D0)
      A=A-DBLE(NA)
1     CONTINUE
      IF (A.NE.-.5D0) GOTO 2
      F=sqrt(PI/X*.5D0)
      G=F
      GOTO 3
2     IF (X.GE.1.0D0) GOTO 4
      B=X*.5D0
      D=-LOG(B)
      E=A*D
      C=A*PI
      IF (abs(C).GE.1.D-15) GOTO 5
      C=1.D0
      GOTO 6
5     C=C/SIN(C)
6     CONTINUE
      IF (abs(E).GE.1.D-15) GOTO 7
      S=1.D0
      GOTO 8
7     S=SINHX(E)/E
8     CONTINUE
      E=EXP(E)
      A1=(E+1.0D0/E)*.5D0
      G=RCPG(A,P,Q)*E
      KA=C*(P*A1+Q*S*D)
      F=KA
      E=A*A
      P=0.5D0*G*C
      Q=0.5D0/G
      C=1.D0
      D=B*B
      KA1=P
      N=1
      GOTO 10
9     N=N+1
10    F=(F*DBLE(N)+P+Q)/(DBLE(N*N)-E)
      C=C*D/DBLE(N)
      P=P/(DBLE(N)-A)
      Q=Q/(DBLE(N)+A)
      G=C*(P-DBLE(N)*F)
      H=C*F
      KA=KA+H
      KA1=KA1+G
      IF (H/KA+abs(G)/KA.GT.EPS) GOTO 9
      F=KA*EXP(X)
      G=KA1*EXP(X)/B
      GOTO 3
4     C=0.25D0-A*A
      G=1.D0
      F=0.D0
      E=X*COS(A*PI)/PI/EPS
      N=1
      GOTO 13
12    N=N+1
13    H=(2.D0*(DBLE(N)+X)*G
     &   -(DBLE(N)-1.D0+C/DBLE(N))*F)/(DBLE(N)+1)
      F=G
      G=H
      IF (H*DBLE(N).LT.E) GOTO 12
      P=F/G
      Q=P
      B=X+X
      E=B-2.D0
15    P=(DBLE(N)-1.D0+C/DBLE(N)) / (E+(DBLE(N)+1.D0)*(2.D0-P))
      Q=P*(Q+1.D0)
      N=N-1
      IF (N.GT.0) GOTO 15
      F=sqrt(PI/B)/(1.D0+Q)
      G=F*(A+X+0.5D0-P)/X
3     CONTINUE
      IF (.NOT.REC) GOTO 16
      X=2.D0/X
      N=1
17    H=F+(A+DBLE(N))*X*G
      F=G
      G=H
      N=N+1
      IF (N.LE.NA) GOTO 17
16    IF (.NOT. REV) GOTO 18
      KA1=F
      KA=G
      GOTO 19
18    KA=F
      KA1=G
19    RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
      SUBROUTINE CHECK_DKH(ICORE,K,TEXT)
      CHARACTER*8 TEXT
      DATA MSK /0/, KSK /5000000/
C
C     GENUEGEND KERNSPEICHER VORHANDEN?
C
      IF (ICORE.EQ.0) GOTO 4
      L=IABS(ICORE)-K
      IF (KSK.GT.L) KSK=L
      IF (L .GE. 0) GOTO 1
      L=-L
      L=(L+1023)/1024
      WRITE (6,100) TEXT,K,ICORE,L
100   FORMAT(' ***'/,' *** NEED MORE CORE AT ',A8/,
     *       ' *** NEED ',I10,' IBM+ WORDS'/,
     *       ' *** HAVE ',I10,' IBM+ WORDS'/,
     *       ' *** NEED ',I10,' K WORDS MORE'/,
     *       ' ***')
      Call Abend
C
C     PRINTOUT,ENDE ODER RETURN?
C
1     IC=IABS(ICORE)
      IF (MSK.LT.IC) MSK=IC
      IF (ICORE.GT.0) GOTO 3
      L=MSK-KSK
      IC=(L+1023)/1024
      WRITE (6,101) L,IC
101   FORMAT(' DYNAMIC CORE USED SO FAR ',I10,' WORDS (',I10,
     *'K)')
3     RETURN
C
C     ENDE DER RECHNUNG
C
4     IC=(MSK+1023)/1024
      KSK=MSK-KSK
      L=(KSK+1023)/1024
      WRITE (6,102) MSK,IC,KSK,L
102   FORMAT(' DYNAMIC STORAGE: ',I8,' WORDS (',I5,'K ) - USED ',
     *I8,' WORDS (',I5,'K)'//,
     *' E N D   O F   C A L C U L A T I O N')
      Call Finish(0)
      END
*                                                                      *
************************************************************************
*                                                                      *
      SUBROUTINE CPLAB (A,B,L,M,N,IA,IB,C,IC,IER)
      implicit real*8(a-h,o-z)
Culf      INTEGER            L,M,N,IA,IB,IC,IER
Culf      real*8   A(IA,M),B(IB,N),C(IC,N)
Culf      real*8   TEMP
      DIMENSION   A(IA,M),B(IB,N),C(IC,N)
      IF (IA .GE. L .AND. IB .GE. M .AND. IC .GE. L) GO TO 5
      IER=129
      GO TO 9000
    5 IER = 0
      DO 15 I=1,L
         DO 15 J=1,N
            TEMP=0.0D0
            DO 10 K=1,M
               TEMP=A(I,K)*B(K,J)+TEMP
   10       CONTINUE
            C(I,J)=C(I,J)+TEMP
   15 CONTINUE
      GO TO 9005
 9000 CONTINUE
 9005 RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
      SUBROUTINE KBR(BETA,X,RESULT)
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     PROGRAM TO CALCULATE CONTINUED FRACTION
C     1/1+ A1*X/1+ A2*X/1+ ...
C
C     V 1.0 - 12.5.86 - BERND HESS
C
C     A(2N+1) = N+BETA
C     A(2N)   = N-1/2
C
      THR=1.D-12
      THR1=1.D-8
      SUM=1.D0
      I=100
1     OLD=SUM
      SUM=1.D0
      J=I
3     SUM=X*(DBLE(J)-0.5D0)/SUM+1.D0
      SUM=X*(DBLE(J)+BETA )/SUM+1.D0
      J=J-1
      IF (J.GT.0) GOTO 3
      SUM=1.D0/(1.D0-0.5D0*X/SUM)
      IF (I.GT.2000) WRITE (6,991) I,BETA,X,SUM
  991 FORMAT(' KBR - I,BETA,X,SUM ',I5,F10.3,D20.10,D30.20)
      I=I+100
      IF (I.LT.300 .OR. abs(OLD-SUM).GT.THR.AND. I.LT.3000) GOTO 1
      RESULT=SUM
      IF (I.LT.3000) RETURN
      DEL=OLD-SUM
      WRITE (6,999) DEL,THR
999   FORMAT(' CONTINUED FRACTION DEL=',D20.10,' LARGER THAN THR=',
     *D20.10)
      IF (DEL.GT.THR1) Call Abend
      RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
      real*8 FUNCTION DCOF(E,LX,IX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     FUNCTION DCOF CALCULATES COEFFICIENT D NEEDED FOR FT
C     DCOF(E,LX,IX)= ((-E)**IX * FAK(LX))/(FAK(LX-2*IX)*FAK(IX))
C     V 1.0 - 12.3.86 - BERND HESS
C
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GA(20),IMAX
      SAVE /CRELOP/
      I=LX-2*IX
      D=1.D0
      IF (IX.EQ.0) GOTO 2
      DO 3 K=1,IX
3     D=-D*E
2     CONTINUE
      D=D*FAK(LX+1)/(FAK(I+1)*FAK(IX+1))
      DCOF=D
      RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
      real*8 FUNCTION RCPG(XVAL,ODD,EVEN)
      REAL*8 XVAL,ODD,EVEN
C
C     PROCEDURE TO CALCULATE RECIPROCAL OF GAMMA FUNCTION
C
      REAL*8 X,ALFA,BETA,X2,B(12)
      INTEGER I
      DATA B /-.28387 65422 76024D0, -.07685 28408 44786D0,
     *         .00170 63050 71096D0,  .00127 19271 36655D0,
     *         .00007 63095 97586D0, -.00000 49717 36704D0,
     *        -.00000 08659 20800D0, -.00000 00331 26120D0,
     *         .00000 00017 45136D0,  .00000 00002 42310D0,
     *         .00000 00000 09161D0, -.00000 00000 00170D0/
      X=XVAL
      X2=X*X*8.D0
      ALFA=-.00000 00000 00001D0
      BETA=0.D0
      I=12
1     BETA=-(ALFA*2.D0+BETA)
      ALFA=-BETA*X2-ALFA+B(I)
      I=I-2
      IF (I.GE.2) GOTO 1
      EVEN=(0.5D0*BETA+ALFA)*X2-ALFA+0.92187 02936 50453D0
      ALFA=-0.00000 00000 00034D0
      BETA=0.D0
      I=11
2     BETA=-(ALFA*2.D0+BETA)
      ALFA=-BETA*X2-ALFA+B(I)
      I=I-2
      IF (I.GE.1) GOTO 2
      ODD=(ALFA+BETA)*2.D0
      RCPG=ODD*X+EVEN
      RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
      real*8 FUNCTION SINHX(XVAL)
      REAL*8 XVAL
      REAL*8 AX,Y,X
      X=XVAL
      AX=abs(X)
      IF (AX.GE.0.3D0) GOTO 1
      IF (AX.GE.0.1D0) GOTO 2
      Y=X*X
      GOTO 3
2     Y=X*X/9.D0
3     CONTINUE
      X=(((1.D0/5040D0*Y+1.D0/120.D0)*Y+1.D0/6.D0)*Y+1.D0)*X
      IF (AX.GE.0.1D0) GOTO 5
      SINHX=X
      RETURN
5     SINHX=X*(1.D0+4.D0*X*X/27.D0)
      RETURN
1     AX=EXP(AX)
      SINHX=sign( (0.5D0*(AX-1.D0/AX)),X )
      RETURN
      END
*                                                                      *
************************************************************************
*                                                                      *
C      SQROPY - PROGRAM TO CALCULATE RELATIVISTIC ONE-ELECTRON
C      OPERATORS IN A BASIS SET OF GAUSSIAN FUNCTIONS
C      V 1.0 - 12.3.86 - BERND HESS
C      V 1.1 - 21.5.87 - BERND HESS - NORMALIZED FUNCTIONS
C      SUBROUTINE "RELOP" MUST BE CALLED TO INITITALIZE COMMON
C      BLOCK BEFORE FIRST EXECUTION OF THIS ROUTINE
C      INTEGRALS ARE CALCULATED BETWEEN   NORMALIZED FUNCTIONS
      real*8 FUNCTION SQROPY(AA,AB,LA,MA,NA,LB,MB,NB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     SUBROUTINE SQROPY CALCULATES MATRIXELEMENT BETWEEN FUNCTION
C     XA**LA * YA**MA * ZA**NA * EXP(-AA*RA**2) AND
C     XB**LB * YB**MB * ZB**NB * EXP(-AB*RB**2)
C     WITH XA=X ,YA=Y ZA=Z
C     V 1.1 - 1.6.86 - BERND HESS
C
      COMMON /CRELOP/ PI,ZWP,SQPI,VELIT,PREA,CSQ,ZWPH32,FAK(26),
     *ZWPH12,BCO(210),GAM(20),IMAX
C@    WRITE (6,922) LA,MA,NA,LB,MB,NB,AA,AB
C@922 FORMAT(' SQROPY ',3I2,2X,3I2,2X,3D14.6)
C
C     WE ALSO NEED OVERLAP. CALCULATE IT IN POSITION SPACE
C
      SQROPY=0.D0
      S=0.D0
      IMAX=LA+MA+NA+LB+MB+NB+3
      IF (IMAX.GT.20) GOTO 1100
      II=LA+LB
      JJ=MA+MB
      KK=NA+NB
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
C
C     AND USE ANGULAR PART TO DETERMINE IF IT VANISHES BY SYMMETRY
C
      IF (ANG.EQ.0.D0) THEN
       SQROPY=0.D0
       RETURN
      ENDIF
      EX=-0.5D0*DBLE(IMAX)
      OVL=0.5D0*ANG*GAM(IMAX)*((AA+AB)**EX)
C
C     CALCULATE SQR INTEGRAL IN MOMENTUM SPACE
C
      EA=1.D0/(4.D0*AA)
      EB=1.D0/(4.D0*AB)
      DELTA=EA+EB
      LLA=LA+MA+NA
      LLB=LB+MB+NB
      C=((2.D0*AA)**LLA * (2.D0*AB)**LLB)
      C=C*(4.D0*AA*AB)**1.5D0
      C=1.D0/C
      LAX=LA/2+1
      MAX=MA/2+1
      NAX=NA/2+1
      LBX=LB/2+1
      MBX=MB/2+1
      NBX=NB/2+1
C
C     COMPUTE VALUE OF BESSEL FUNCTION
C
      ARG=0.5D0*DELTA*CSQ
      CALL BESSKA(0.D0,ARG,X0,X1)
C#    WRITE (6,901) ARG,X0,X1
C#901 FORMAT(' ARG,X0,X1 ',3D20.10)
C
C     COMPUTE ARGUMENT OF CONTINUED FRACTION
C
      ARG=1.D0/(CSQ*DELTA)
C
C     AND RUN THE TEDIOUS SUM OVER SIMPLE FOURIER INTEGRALS ...
C
      DO 1 IA=1,LAX
      D1=DCOF(AA,LA,IA-1)
      DO 2 JA=1,MAX
      D2=D1*DCOF(AA,MA,JA-1)
      DO 3 KA=1,NAX
      D3=D2*DCOF(AA,NA,KA-1)
      DO 4 IB=1,LBX
      D4=D3*DCOF(AB,LB,IB-1)
      DO 5 JB=1,MBX
      D5=D4*DCOF(AB,MB,JB-1)
      DO 6 KB=1,NBX
      D6=D5*DCOF(AB,NB,KB-1)
      II=LA-2*IA+LB-2*IB+4
      JJ=MA-2*JA+MB-2*JB+4
      KK=NA-2*KA+NB-2*KB+4
C
C     ANGULAR PART OF FOURIER INTEGRAL
C
      IMAX=II+JJ+KK+3
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      IF (ANG.EQ.0.D0) GOTO 6
C
C     RADIAL PART
C
      U=X1*D6
      IF (IMAX.EQ.3) GOTO 60
      JMAX=IMAX-1
      DO 70 I=3,JMAX,2
      BET=0.5D0*DBLE(I)
      CALL KBR(BET,ARG,R)
70    U=U*R*BET/DELTA
60    CONTINUE
      U=0.25D0*VELIT*U/DELTA
      S=S+U*ANG
6     CONTINUE
5     CONTINUE
4     CONTINUE
3     CONTINUE
2     CONTINUE
1     CONTINUE
      IF (MOD(LLA-LLB,4).EQ.2) S=-S
      S=S*C
      S=CSQ*(S-OVL)
C
C     NORMALIZATION

C
      II=LA+LA
      JJ=MA+MA
      KK=NA+NA
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DBLE(II+JJ+KK+3)
      OV1=0.5D0*ANG*GAM(II+JJ+KK+3)*((AA+AA)**EX)
C#    WRITE (6,*) ' SRQ  OV1',LA,MA,NA,AA,ANG,sqrt(1/OV1)
      II=LB+LB
      JJ=MB+MB
      KK=NB+NB
      ANG=THETA(II+JJ,KK)*PHI(JJ,II)
      EX=-0.5D0*DBLE(II+JJ+KK+3)
      OV2=0.5D0*ANG*GAM(II+JJ+KK+3)*((AB+AB)**EX)
C#    WRITE (6,*) ' SRQ  OV2',LB,MB,NB,AB,ANG,sqrt(1/OV2)
      SQROPY=S/sqrt(OV1*OV2)
      RETURN
C
C     ANGULAR MOMENTUM TOO LARGE
C
1100  WRITE (6,1101) LA,MA,NA,LB,MB,NB,AA,AB
1101  FORMAT(' *** ANGULAR MOMENTUM TOO LARGE ***'/,
     *6I3,3X,2D17.10)
      Call Abend
      END
