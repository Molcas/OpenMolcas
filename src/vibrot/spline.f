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
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  THIS ROUTINE PERFORMS A SPLINE FIT TO A SET OF POINTS (XFIT,YFIT).  C
C  IT INTERPOLATES A SET OF POINTS (XINT,YINT), AND LOCATES ALL        C
C  EXTREME VALUES (XEXT,YEXT) AND DETERMINES THE TYPE (IEXT).          C
C  VARIOUS TYPES OF BOUNDARY CONDITIONS CAN BE CHOOSEN (IBOUND),       C
C  AT PRESENT IBOUND=1 AND 2 ARE IMPLEMENTED.                          C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      SUBROUTINE SPLINE(XFIT,YFIT,NFIT,XINT,YINT,NINT,XEXT,YEXT,IEXT,
     *   NEXT,IBOUND)
      IMPLICIT REAL*8 (A-H,O-Z)
#include "WrkSpc.fh"
      Parameter (ndim=5000)
C NFIT, NINT ARE INPUT ONLY
      DIMENSION XFIT(NFIT),YFIT(NFIT),XINT(NINT),YINT(NINT)
C NEXT IS INPUT (DIMENSION BOUND) AND OUTPUT (NR OF EXTREMA).
      DIMENSION XEXT(NEXT),YEXT(NEXT),IEXT(NEXT),T0(2)
C      CHARACTER*12 TEXT(3)
      DATA ZERO,ONE,TWO,THREE/0.0D0,1.0D0,2.0D0,3.0D0/
C      DATA TEXT/'MAX. POINT','SADDLE POINT','MIN. POINT'/
      IF(NFIT.GT.NDIM) Then
        Write(6,*)'SPLINE Error: NFIT.gt.NDIM.'
        Write(6,'(1x,a,2I8)')'NFIT,NDIM:',NFIT,NDIM
        Call Abend
      End If
      NOS=NFIT-1
      Call GetMem('MatA','Allo','Real',indexA,nfit*nfit)
      Call GetMem('VecH','Allo','Real',indexH,ndim)
      Call GetMem('VecD','Allo','Real',indexD,ndim)
      Call GetMem('VecC','Allo','Real',indexC,ndim)
      Call GetMem('VecB','Allo','Real',indexB,ndim)
      DO 10 I=0,nfit*nfit-1
         Work(indexA+I)=ZERO
10    CONTINUE
      DO 20 I=1,NOS
         Work(indexH+I-1)=XFIT(I+1)-XFIT(I)
         Work(indexD+I-1)=(YFIT(I+1)-YFIT(I))/Work(indexH+I-1)
20    CONTINUE
      if(IBOUND.eq.1) goto 100
      if(IBOUND.eq.2) goto 200
      Write(6,*)'SPLINE Error: IBOUND should be 1 or 2.'
      Write(6,'(1x,a,2I8)')'But IBOUND=',IBOUND
      Call Abend
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  IBOUND=1:                                                           C
C  BOTH END POINTS ARE EXTRAPOLATED WITH STRAIGHT LINES.               C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
100   CONTINUE
      Work(indexA)=TWO*Work(indexH)
      Work(indexA+nfit)=Work(indexH)
      Work(indexB)=THREE*Work(indexH)*Work(indexD)
      Work(indexA+nfit*nfit-nfit-1)=Work(indexH+NOS-1)
      Work(indexA+nfit*nfit-1)=TWO*Work(indexH+NOS-1)
      Work(indexB+nfit-1)=THREE*Work(indexH+NOS-1)*Work(indexD+NOS-1)
      DO 110 I=2,NOS
         Work(indexA+(I-2)*nfit+I-1)=Work(indexH+I-1)
         Work(indexA+(I-1)*nfit+I-1)=
     *       TWO*(Work(indexH+I-1)+Work(indexH+I-2))
         Work(indexA+I*nfit+I-1)=Work(indexH+I-2)
         Work(indexB+I-1)=THREE*(Work(indexH+I-2)*Work(indexD+I-1)+
     *        Work(indexH+I-1)*Work(indexD+I-2))
110   CONTINUE
      Call Dminv(nfit,nfit,Work(indexA))
      DO 120 I=1,NFIT
         SUM=ZERO
         DO 121 J=1,NFIT
            SUM=SUM+Work(indexA+(J-1)*nfit+I-1)*Work(indexB+J-1)
121      CONTINUE
         Work(indexC+I-1)=SUM
120   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PERFORM INTERPOLATION (EXTRAPOLATION).                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 130 IP=1,NINT
         IF(XINT(IP).LT.XFIT(1)) THEN
            YINT(IP)=YFIT(1)+Work(indexC)*(XINT(IP)-XFIT(1))
         ELSE IF(XINT(IP).GT.XFIT(NFIT)) THEN
            YINT(IP)=YFIT(NFIT)+
     *               Work(indexC+nfit-1)*(XINT(IP)-XFIT(NFIT))
         ELSE
            NS=1
131         IF(XINT(IP).GT.XFIT(NS+1)) THEN
               NS=NS+1
               GOTO 131
            END IF
            T=(XINT(IP)-XFIT(NS))/Work(indexH+NS-1)
            T1=ONE-T
            YINT(IP)=T*YFIT(NS+1)+T1*YFIT(NS)+
     *         Work(indexH+NS-1)*T*T1*
     *         (T1*(Work(indexC+NS-1)-Work(indexD+NS-1))-
     *         T*(Work(indexC+NS)-Work(indexD+NS-1)))
         END IF
130   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  LOCATE ALL EXTREMUM POINTS.                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NEXT=0
      DO 150 NS=1,NOS
         A2=THREE*(Work(indexC+NS-1)+Work(indexC+NS)-
     *      TWO*Work(indexD+NS-1))
         A1=TWO*(THREE*Work(indexD+NS-1)-
     *      TWO*Work(indexC+NS-1)-Work(indexC+NS))
         A0=Work(indexC+NS-1)
         B1=TWO*A2
         B0=A1
         IF(abs(A2).GT.(abs(A1)+abs(A0))*1.0D-3) THEN
***** A2 LARGE *****
            P=-A1/A2/TWO
            Q=A0/A2
            S=P*P-Q
            IF(abs(S).LT.1.0D-10) THEN
***** DOUBLE ROOT *****
               IT=2
               T=P
               IF(T.GE.ZERO.AND.T.LE.ONE) THEN
                  X=XFIT(NS)+T*Work(indexH+NS-1)
                  Y=YFIT(NS)+Work(indexH+NS-1)*
     *              T*(A0+T*(A1/TWO+T*A2/THREE))
C                 WRITE(6,2001) TEXT(IT),X,Y
C2001  FORMAT(1X,A,'  X=',E12.5,' ,  Y=',E12.5)
                  NEXT=NEXT+1
                  XEXT(NEXT)=X
                  YEXT(NEXT)=Y
                  IEXT(NEXT)=IT
               END IF
            ELSE IF(S.GT.ZERO) THEN
***** TWO ROOTS *****
               T0(1)=P+sqrt(S)
               T0(2)=P-sqrt(S)
               DO 151 I=1,2
                  T=T0(I)
                  IF(T.GE.ZERO.AND.T.LE.ONE) THEN
                     SE=B0+T*B1
                     X=XFIT(NS)+T*Work(indexH+NS-1)
                     Y=YFIT(NS)+
     *                 Work(indexH+NS-1)*T*(A0+T*(A1/TWO+T*A2/THREE))
                     IT=1
                     IF(SE.GT.ZERO) IT=3
C                    WRITE(6,2001) TEXT(IT),X,Y
                     NEXT=NEXT+1
                     XEXT(NEXT)=X
                     YEXT(NEXT)=Y
                     IEXT(NEXT)=IT
                  END IF
151            CONTINUE
            END IF
         ELSE IF(abs(A1).GE.0.9D0*abs(A0)) THEN
***** A2 SMALL *****
            T=-A0/A1
            T=-(A0+A2*T*T)/A1
            T=-(A0+A2*T*T)/A1
            T=-(A0+A2*T*T)/A1
            IF(T.GE.ZERO.AND.T.LE.ONE) THEN
               IT=1
               IF(B0.GT.ZERO) IT=3
               X=XFIT(NS)+T*Work(indexH+NS-1)
               Y=YFIT(NS)+work(indexH+NS-1)*T*(A0+T*(A1/TWO+T*A2/THREE))
C              WRITE(6,2001) TEXT(IT),X,Y
               NEXT=NEXT+1
               XEXT(NEXT)=X
               YEXT(NEXT)=Y
               IEXT(NEXT)=IT
            END IF
         END IF
150   CONTINUE
      Call GetMem('MatA','Free','Real',indexA,nfit*nfit)
      Call GetMem('VecH','Free','Real',indexH,ndim)
      Call GetMem('VecD','Free','Real',indexD,ndim)
      Call GetMem('VecC','Free','Real',indexC,ndim)
      Call GetMem('VecB','Free','Real',indexB,ndim)
      RETURN
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  IBOUND=2:                                                           C
C  EXTRAPOLATION TO PLUS INFINITY IS DONE WITH A STRAIGHT              C
C  HORIZONTAL LINE.                                                    C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
200   CONTINUE
#if 0
!IFG IBOUND=2 implemented, but never used
      Work(indexA+nfit*nfit-1)=ONE
      Work(indexB+nfit-1)=ZERO
      Work(indexA+(NOS-1)*nfit+NOS-1)=TWO
      Work(indexA+(nfit-1)*nfit+NOS-1)=ONE
      Work(indexB+NOS-1)=THREE*Work(indexD+NOS-1)
      DO 210 I=2,NOS
         Work(indexA+(I-2)*nfit+I-2)=Work(indexH+I-1)
         Work(indexA+(I-1)*nfit+I-2)=
     *       TWO*(Work(indexH+I-1)+Work(indexH+I-2))
         Work(indexA+I*nfit+I-2)=Work(indexH+I-2)
         Work(indexB+I-2)=THREE*(Work(indexH+I-2)*Work(indexD+I-1)+
     *          Work(indexH+I-1)*Work(indexD+I-2))
210   CONTINUE
      Call Dminv(nfit,nfit,Work(indexA))
      DO 220 I=1,NFIT
         SUM=ZERO
         DO 221 J=1,NFIT
            SUM=SUM+Work(indexA+(J-1)*nfit+I-1)*Work(indexB+J-1)
221      CONTINUE
         Work(indexC+I-1)=SUM
220   CONTINUE
      S1=Work(indexC)
      S2=(Work(indexC)+TWO*Work(indexC+1)-
     *   THREE*Work(indexD))/Work(indexH)
      C2=-S2/S1
      C1=S2/C2/C2/EXP(-C2*XFIT(1))
      C0=YFIT(1)-C1*EXP(-C2*XFIT(1))
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  PERFORM INTERPOLATION (EXTRAPOLATION).                              C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      DO 230 IP=1,NINT
         IF(XINT(IP).LT.XFIT(1)) THEN
            YINT(IP)=C0+C1*EXP(-C2*XINT(IP))
         ELSE IF(XINT(IP).GT.XFIT(NFIT)) THEN
            YINT(IP)=YFIT(NFIT)+
     *               Work(indexC+nfit-1)*(XINT(IP)-XFIT(NFIT))
         ELSE
            NS=1
231         IF(XINT(IP).GT.XFIT(NS+1)) THEN
               NS=NS+1
               GOTO 231
            END IF
            T=(XINT(IP)-XFIT(NS))/Work(indexH+NS-1)
            T1=ONE-T
            YINT(IP)=T*YFIT(NS+1)+T1*YFIT(NS)+
     *         Work(indexH+NS-1)*T*T1*
     *         (T1*(Work(indexC+NS-1)-Work(indexD+NS-1))-
     *         T*(Work(indexC+NS)-Work(indexD+NS-1)))
         END IF
230   CONTINUE
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C                                                                      C
C  LOCATE ALL EXTREMUM POINTS.                                         C
C                                                                      C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      NEXT=0
      DO 250 NS=1,NOS
         A2=THREE*(Work(indexC+NS-1)+Work(indexC+NS)-
     *      TWO*Work(indexD+NS-1))
         A1=TWO*(THREE*Work(indexD+NS-1)-
     *      TWO*Work(indexC+NS-1)-Work(indexC+NS))
         A0=Work(indexC+NS-1)
         B1=TWO*A2
         B0=A1
         IF(abs(A2).GT.(abs(A1)+abs(A0))*1.0D-3) THEN
***** A2 LARGE *****
            P=-A1/A2/TWO
            Q=A0/A2
            S=P*P-Q
            IF(abs(S).LT.1.0D-10) THEN
***** DOUBLE ROOT *****
               IT=2
               T=P
               IF(T.GE.ZERO.AND.T.LE.ONE) THEN
                  X=XFIT(NS)+T*Work(indexH+NS-1)
                  Y=YFIT(NS)+
     *              Work(indexH+NS-1)*T*(A0+T*(A1/TWO+T*A2/THREE))
C                 WRITE(6,2001) TEXT(IT),X,Y
                  NEXT=NEXT+1
                  XEXT(NEXT)=X
                  YEXT(NEXT)=Y
                  IEXT(NEXT)=IT
               END IF
            ELSE IF(S.GT.ZERO) THEN
***** TWO ROOTS *****
               T0(1)=P+sqrt(S)
               T0(2)=P-sqrt(S)
               DO 251 I=1,2
                  T=T0(I)
                  IF(T.GE.ZERO.AND.T.LE.ONE) THEN
                     SE=B0+T*B1
                     X=XFIT(NS)+T*Work(indexH+NS-1)
                     Y=YFIT(NS)+
     *                 Work(indexH+NS-1)*T*(A0+T*(A1/TWO+T*A2/THREE))
                     IT=1
                     IF(SE.GT.ZERO) IT=3
C                    WRITE(6,2001) TEXT(IT),X,Y
                     NEXT=NEXT+1
                     XEXT(NEXT)=X
                     YEXT(NEXT)=Y
                     IEXT(NEXT)=IT
                  END IF
251            CONTINUE
            END IF
         ELSE IF(abs(A1).GE.0.9D0*abs(A0)) THEN
***** A2 SMALL *****
            T=-A0/A1
            T=-(A0+A2*T*T)/A1
            T=-(A0+A2*T*T)/A1
            T=-(A0+A2*T*T)/A1
            IF(T.GE.ZERO.AND.T.LE.ONE) THEN
               IT=1
               IF(B0.GT.ZERO) IT=3
               X=XFIT(NS)+T*Work(indexH+NS-1)
               Y=YFIT(NS)+Work(indexH+NS-1)*T*(A0+T*(A1/TWO+T*A2/THREE))
C              WRITE(6,2001) TEXT(IT),X,Y
               NEXT=NEXT+1
               XEXT(NEXT)=X
               YEXT(NEXT)=Y
               IEXT(NEXT)=IT
            END IF
         END IF
250   CONTINUE
      Call GetMem('MatA','Free','Real',indexA,nfit*nfit)
      Call GetMem('VecH','Free','Real',indexH,ndim)
      Call GetMem('VecD','Free','Real',indexD,ndim)
      Call GetMem('VecC','Free','Real',indexC,ndim)
      Call GetMem('VecB','Free','Real',indexB,ndim)
      RETURN
#endif
      End
