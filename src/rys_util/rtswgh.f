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
      SUBROUTINE RTSWGH(TARR,NT,U2,WGH,NRYS)
      use vRys_RW
      IMPLICIT REAL*8 (A-H,O-Z)
#include "itmax.fh"
#include "real.fh"
#include "print.fh"
#include "abtab.fh"
CMAW
#include "FMM.fh"
      DIMENSION TARR(NT),U2(NRYS,NT),WGH(NRYS,NT)
      DIMENSION ROOT(MXSIZ1,MXSIZ1),RYS(0:MXSIZ1),RYSD(0:MXSIZ1)
      DIMENSION ALPHA(0:MXSIZ1),BETA(0:MXSIZ1)
      DIMENSION BINV(0:MXSIZ1)
      parameter (coef1=-1.0d0/120d0, coef2=-5d0*coef1, coef3=-2d0*coef2)
      parameter (coef4=-coef3, coef5=-coef2, coef6=-coef1)
c
      iRout = 78
      iPrint = nPrint(iRout)
      RYSD(0)=ZERO
c
      IF (NRYS.gt.maxdeg) THEN
        CALL WarningMessage(2,' Too many requested Rys roots.')
        CALL AbEnd()
      ENDIF
c
      DO 5000 IT=1,NT
        T=TARR(IT)
c Use asymptotic formulae if outside table.
CMAW    if(t.gt.TVALUE(NTAB2-NTAB1-1)) then
*
*  For the FMM we use the asymptotic limit to compute the
*  multipole-component of the integrals
*
      If( (t.gt.TVALUE(NTAB2-NTAB1-1)) .OR. asymptotic_Rys ) Then
        do 25 iroot=1,nRYS
          tmp = 1.0D0/T
          U2(iroot,IT)=HerR2(iHerR2(nRys)+iroot-1)*tmp
          WGH(iroot,IT)=HerW2(iHerW2(nRys)+iroot-1)*SQRT(tmp)
 25     Continue
        Go To 5000
      end if
c translate to tabulation function for equidist. interp.
c xn=interpol. variable.
c Ex: T=0.0--0.05 gives xn=0.0--1.0 (approx.)
c itab= Start tab index for interp. Ex above: itab=1.
      xn=5d0*T + 200d0*T/(14d0+T)
      nx=Int(xn)
      itab=nx-1-NTAB1
      p=xn-dble(nx)
      a2=(p+2d0)
      a3=(p+1d0)*a2
      a4=(p    )*a3
      a5=(p-1d0)*a4
      a6=(p-2d0)*a5
      b5=(p-3d0)
      b4=(p-2d0)*b5
      b3=(p-1d0)*b4
      b2=(p    )*b3
      b1=(p+1d0)*b2
      c1=coef1*   b1
      c2=coef2*a2*b2
      c3=coef3*a3*b3
      c4=coef4*a4*b4
      c5=coef5*a5*b5
      c6=coef6*a6
      ALPHA(0)=c1*ATAB(0,itab)+c2*ATAB(0,itab+1)+
     &         c3*ATAB(0,itab+2)+c4*ATAB(0,itab+3)+
     &         c5*ATAB(0,itab+4)+c6*ATAB(0,itab+5)
      do 20 k=1,NRYS
        ALPHA(k)=c1*ATAB(k,itab)+c2*ATAB(k,itab+1)+
     &         c3*ATAB(k,itab+2)+c4*ATAB(k,itab+3)+
     &         c5*ATAB(k,itab+4)+c6*ATAB(k,itab+5)
        BETA(k)=c1*BTAB(k,itab)+c2*BTAB(k,itab+1)+
     &         c3*BTAB(k,itab+2)+c4*BTAB(k,itab+3)+
     &         c5*BTAB(k,itab+4)+c6*BTAB(k,itab+5)
        BINV(K)=ONE/BETA(K)
  20  continue
      rys(0)=c1*p0(itab)+c2*p0(itab+1)+
     &       c3*p0(itab+2)+c4*p0(itab+3)+
     &       c5*p0(itab+4)+c6*p0(itab+5)
        ROOT(1,1)=ALPHA(0)
        x1=(ALPHA(0)+ALPHA(1))/2d0
        x2=(ALPHA(0)-ALPHA(1))/2d0
        x3=SQRT(x2**2+BETA(1)**2)
        ROOT(1,2)=x1-x3
        ROOT(2,2)=x1+x3
C LOOP OVER DEGREE OF RYS POLY
        DO 2000 IDEG=3,NRYS
C ESTIMATE POSITION OF U2 OF THIS DEGREE:
          ROOT(1,IDEG)=(DBLE(IDEG)-0.5D00)*ROOT(1,IDEG-1)/DBLE(IDEG)
          ROOT(IDEG,IDEG)=
     &      ONE-(DBLE(IDEG)-0.5D00)*(ONE-ROOT(IDEG-1,IDEG-1))/DBLE(IDEG)
          DO 80 IROOT=2,IDEG-1
            R1=ROOT(IROOT,IDEG-1)
            R2=ROOT(IROOT-1,IDEG-1)
            x2=(DBLE(IROOT)-0.5D00)/DBLE(IDEG)
            x1=ONE-x2
            R=x1*R1+x2*R2
            ROOT(IROOT,IDEG)=R
  80      CONTINUE
C         IF(IDEG.EQ.NRYS) ITER=0
          RYSD(1)=RYS(0)*BINV(1)
          DO 1000 IROOT=1,IDEG
            Z=ROOT(IROOT,IDEG)
C Compute the correction coefficient:
            corr=zero
            do 91 J=1,iroot-1
              corr=corr+one/(root(iroot,ideg)-root(j,ideg))
  91        continue
            do 92 J=iroot+1,ideg
              corr=corr+one/(root(iroot,ideg)-root(j,ideg))
  92        continue
C COMPUTE RYS AND FIRST DERIVATIVE, DO NEWTON-RAPHSON:
  99        CONTINUE
            RYS(1)=(Z-ALPHA(0))*RYSD(1)
            ZZ=(Z-ALPHA(1))
            RYSD(2)=(ZZ*RYSD(1)+RYS(1))*BINV(2)
            RYS(2) =(ZZ*RYS(1)-BETA(1)*RYS(0))*BINV(2)
            DO 110 K=2,IDEG-1
              ZZ=Z-ALPHA(K)
              BK= BETA(K)
              RYSD(K+1)=(RYS(K)+ZZ*RYSD(K)-BK*RYSD(K-1))*BINV(K+1)
              RYS(K+1) =(       ZZ*RYS(K) -BK*RYS (K-1))*BINV(K+1)
 110        CONTINUE
            DELTA=-RYS(IDEG)/(RYSD(IDEG)-CORR*RYS(IDEG))
            Z=Z+DELTA
C           IF(IDEG.EQ.NRYS) ITER=ITER+1
            IF(ABS(DELTA).GT.1.0d-08) goto 99
            ROOT(IROOT,IDEG)=Z
 1000     CONTINUE
 2000   CONTINUE
C     IF(NRYS.GT.2)
C    &      write(*,'(1x,a,f8.2)')' Avg. iter/root:',(iter*one)/NRYS
      DO 3010 IROOT=1,NRYS
        Z=ROOT(IROOT,NRYS)
C COMPUTE RYS VALUES AND ADD SQUARES TO GET WGH:
        SUM=RYS(0)**2
        IF(NRYS.EQ.1) goto 3021
        RYS(1)=(Z-ALPHA(0))*RYS(0)*BINV(1)
        SUM=SUM+RYS(1)**2
        IF(NRYS.EQ.2) goto 3021
        ZZ=(Z-ALPHA(1))
        RYS(2) =(ZZ*RYS(1)-BETA(1)*RYS(0))*BINV(2)
        SUM=SUM+RYS(2)**2
        IF(NRYS.EQ.3) goto 3021
        DO 3020 K=2,NRYS-2
          ZZ=Z-ALPHA(K)
          BK= BETA(K)
          RYS(K+1) =(ZZ*RYS(K) -BK*RYS (K-1))*BINV(K+1)
          SUM=SUM+RYS(K+1)**2
 3020   CONTINUE
 3021   CONTINUE
        WGH(IROOT,IT)=ONE/SUM
        U2(IROOT,IT)=ROOT(IROOT,NRYS)
 3010 CONTINUE
C
 5000 CONTINUE
C
#ifdef _DEBUG_
      If (iPrint.ge.99) Then
         Call RecPrt(' Roots',' ',U2,nRys,nT)
         Call RecPrt(' Weights',' ',Wgh,nRys,nT)
      End If
#endif
      RETURN
      END
