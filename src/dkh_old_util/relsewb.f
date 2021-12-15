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
* Copyright (C) 1995, Bernd Artur Hess                                 *
************************************************************************
C $Id: relsewb.r,v 1.4 1995/05/08 14:08:53 hess Exp $
C calculate relativistic operators
C   Bernd Artur Hess, hess@uni-bonn.de
C
C
      SUBROUTINE SCFCLI2(idbg,epsilon,S,H,V,PVP,N,ISIZE,VELIT,
     *BU,P,G,EV2,EIG,SINV,REVT,AUX,OVE,EW,E,AA,RR,TT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION S(ISIZE),H(ISIZE),V(ISIZE),PVP(ISIZE)
      DIMENSION BU(ISIZE),P(ISIZE),G(ISIZE),EV2(ISIZE),
     *EIG(N,N),SINV(N,N),REVT(N,N),AUX(N,N),OVE(N,N),
     *EW(N),E(N),AA(N),RR(N),TT(N)
C
*
*
      TOL=1.D-14
      PREA=1/(VELIT*VELIT)
      CON2=PREA+PREA
      CON=1.D0/PREA
      DO 90 I=1,ISIZE
      BU(I)=S(I)
90    CONTINUE
      ii=0
      DO i=1,n
      DO j=1,i
      ii=ii+1
      aux(I,J)=s(ii)
      aux(J,I)=s(ii)
      ENDDO
      ENDDO
*      write(*,*) 'OVERLAP MATRIX'
*      do i=1,n
*         write(6,'(5f10.5)') (aux(i,j),j=1,n)
*      enddo
C
C     CALCULATE DETERMINANT
C
      icontr=-1
      dtol=tol
      CALL dcopiv(aux,aux,n,1,n,dtol,det,iex,icontr,p)
*     if(idbg.gt.0)WRITE (idbg,2016) icontr,det,iex
*2016  FORMAT(' relsewb| DCOPIV rc=',I2,', |S|=',D20.6,'x 10**(',I4,') ')
      IF (icontr.NE.0) THEN
*     WRITE (6,2016) icontr,det,iex
*      WRITE (6,2012) dtol
*2012   FORMAT('  relsewb|****** '/,
*     * '        |****** WARNING - OVERLAP MATRIX SINGULAR '/,
*     * '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,
*     * '        |******'//)
      CALL errex_rel(' relsewb| singular overlap matrix')
      ENDIF
C
C     SCHMIDT-ORTHOGONALIZE
C
*
      CALL Sogr(iDbg,N,BU,SINV,P,OVE,EW)
*
      CALL Square(BU,OVE,N,1,N)
C      CALL FilMar(N,BU,OVE)

C
C ** SINV CONTAINS TRANSFORMATION TO ORTHOGONAL AO-BASIS
C ** OVE  CONTAINS OVERLAP MATRIX IN FULL
C--------------------------------------------------------------------
C     NON-RELATIVISTIC KINETIC ENERGY
C--------------------------------------------------------------------
*
      CALL Diagr(H,N,EIG,EW,SINV,AUX,BU)
*
*     if(idbg.gt.0)WRITE (idbg,556)
* 556 FORMAT(//,7X,'- NREL. ENERG.  -  DIVIDED BY C ',
*    *             '- REL.  ENERG.  -  MOMENTUM    -',
*    *             ' TERMS OF POWER SERIES (LOW ENERGY ONLY)'//)
      DO 4 I=1,N
      IF (ew(i).LT.0.D0) THEN
*      WRITE (*,*) ' scfcli2| ew(',i,') = ',ew(i)
       CALL errex_rel('kinetic energy eigenvalue less than zero')
      ENDIF
C
C     IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO
C
      RATIO=EW(I)/VELIT
C
C     CALCULATE RELATIVISTIC ENERGY AND MOMENTUM
C
      SR=sqrt(2.D0*EW(I))
      TT(I)=EW(I)
      IF (RATIO.GT.0.02D0) GOTO 11
      TV1=EW(I)
      TV2=-TV1*EW(I)*PREA/2.D0
      TV3=-TV2*EW(I)*PREA
      TV4=-TV3*EW(I)*PREA*1.25D0
      EW(I)=TV1+TV2+TV3+TV4
*     if(idbg.gt.0)WRITE (idbg,100) I,TV1,RATIO,EW(I),SR,TV2,TV3,TV4
* 100 FORMAT(1X,I4,7(2X,D14.8))
      GOTO 12
  11  TV1=EW(I)
      EW(I)=CON*(sqrt(1.D0+CON2*EW(I))-1.D0)
*     if(idbg.gt.0)WRITE (idbg,100) I,TV1,RATIO,EW(I),SR
  12  CONTINUE
      E(I)=EW(I)+CON
4     CONTINUE
C---------------------------------------------------------------------
C     CALCULATE REVERSE TRANSFORMATION
C---------------------------------------------------------------------
      DO 3 I=1,N
      DO 30 J=1,N
      AUX(I,J)=0.D0
      DO 31 K=I,N
      AUX(I,J)=AUX(I,J)+SINV(I,K)*EIG(K,J)
   31 CONTINUE
   30 CONTINUE
    3 CONTINUE
      DO 6 I=1,N
      DO 60 J=1,N
      REVT(I,J)=0.D0
      DO 5 K=1,N
      REVT(I,J)=REVT(I,J)+OVE(I,K)*AUX(K,J)
    5 CONTINUE
   60 CONTINUE
    6 CONTINUE
      IJ=0
      DO 8 I=1,N
      DO 80 J=1,I
      IJ=IJ+1
      H(IJ)=0.D0
      DO 7 K=1,N
      H(IJ)=H(IJ)+REVT(I,K)*REVT(J,K)*EW(K)
    7 CONTINUE
   80 CONTINUE
    8 CONTINUE
C
C     CALCULATE KINEMATICAL FACTORS
C
      DO 362 I=1,N
      AA(I)=sqrt((CON+E(I)) / (2.D0*E(I)))
      RR(I)=sqrt(CON)/(CON+E(I))
362   CONTINUE
C
C     POTENTIAL
C
C
C     BEYOND THIS POINT, OVE IS USED AS SCRATCH ARRAY
C
C
C    TRANSFORM V TO T-BASIS
C
      CALL TrSmr(V,SINV,G,N,AUX,OVE)
*
      CALL TrSmr(G,EIG,BU,N,AUX,OVE)

culf
      if(idbg.gt.0)CALL PRMAT(IDBG,V,N,0,'v oper  ')
C
C    MULTIPLY
C
      IJ=0
      DO 2005 I=1,N
      DO 2006 J=1,I
      IJ=IJ+1
      P(IJ)=BU(IJ)
      BU(IJ)=BU(IJ)*AA(I)*AA(J)
2006  CONTINUE
2005  CONTINUE
*
      CALL TrSmtr(BU,REVT,G,0.0D0,N,AUX,OVE)
*
culf
      if(idbg.gt.0)CALL PRMAT(IDBG,g,N,0,'g oper  ')
*
      CALL AddMar(ISIZE,G,H)
C
C     PVP INTEGRALS AND TRANSFORM THEM TO T-BASIS
C
      if(idbg.gt.0)CALL PRMAT(IDBG,pvp,N,0,'raw pvp integrals  ')
      CALL TrSmr(PVP,SINV,G,N,AUX,OVE)
*
      CALL TrSmr(G,EIG,BU,N,AUX,OVE)
C
C    MULTIPLY
C
      IJ=0
      DO 3005 I=1,N
      DO 3006 J=1,I
      IJ=IJ+1
      G(IJ)=BU(IJ)
      BU(IJ)=BU(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
3006  CONTINUE
3005  CONTINUE
*
      CALL TrSmtr(BU,REVT,EV2,0.0D0,N,AUX,OVE)
culf
      if(idbg.gt.0)CALL PRMAT(IDBG,ev2,n,0,'pvp oper')
*
      CALL AddMar(ISIZE,EV2,H)
C
C     CALCULATE Even2r OPERATOR
C
*     if(idbg.gt.0)CALL Even2r(idbg,N,P,G,E,AA,RR,TT,EIG,AUX,OVE)
*
C
C    TRANSFORM BACK
C
culf
      if(idbg.gt.0)CALL PRMAT(IDBG,g,n,0,'ev2 orig')
*     CALL TrSmtr(G,REVT,EV2,0.0D0,N,AUX,OVE)
culf
*     if(idbg.gt.0)CALL PRMAT(IDBG,ev2,n,0,'ev2 oper')
*     CALL AddMar(ISIZE,EV2,H)
culf
      if(idbg.gt.0)CALL PRMAT(IDBG,h,n,0,'h   oper')
*     CALL Sogr(iDbg,N,S,SINV,P,OVE,EW)
*     CALL Diagr(H,N,EIG,EW,SINV,AUX,0)
*     write(*,*) 'END OF SCFCLI2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
C
*     if(idbg.gt.0)WRITE (idbg,*) '--- EIGENVALUES OF H MATRIX ---'
*     if(idbg.gt.0)WRITE (idbg,'(4D20.12)') EW
      RETURN
c Avoid unused argument warnings
      IF (.FALSE.) CALL Unused_real(epsilon)
      END
