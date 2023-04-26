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
!
c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
      SUBROUTINE SCECORas(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,
     1                                    NCN,V,SDRDY,BMAX,VLIM,DGDV2)
c** Subroutine calculates (approximate!) semiclassical estimate of
c  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential
c  {V(i),i=1,NDP} (in 'internal BCFT units' {V[cm-1]*BFCT}), and uses
c  those results to estimate energy of level  KVLEV
c** If the 'clever' semiclassical procedure fails - try a brute force
c  step-by-step search, using alternately INNER & OUTER well starting
c** BMAX is internal barrier maximum energy for double-well case,
c   and very large negative number for single-well potential
c** On return, negative DGDV2 signals error!  No phase integrals found
c
      INTEGER I,II,I1,I2,I3,I4,IV1,IV2,INNER,ICOR,JROT,KV,KVB,KVLEV,
     1  KVDIF,NDP,NCN,IDIF,BRUTE,IB,IWR
      REAL*8 EO,DE0,RH,BFCT,ARG2,ARG3,EINT,VPH1,VPH2,DGDV1,DGDV2,DGDVM,
     1  DGDV2P,DGDVB,DGDVBP,EBRUTE,DEBRUTE,DE1,DE2,Y1,Y2,Y3,RT,ANS1,
     2  ANS2,XDIF,VLIM,BMAX,Pi,Pi2,PNCN,PP1,V(NDP),SDRDY(NDP)
      SAVE BRUTE,EBRUTE,DEBRUTE,DGDVB,Pi,Pi2
      DATA DGDVB/-1.d0/,KVB/-1/,Pi/3.1415926454d0/,Pi2/6.283185308d0/
c
      ARG3 = 0
      I1 = 0
      II = 0
      DGDVBP=-1.d0
      DGDV2= -1.d0
      EINT= EO*BFCT
      IF(KVLEV.EQ.0) DGDVB= -1.d0
      KVDIF= KVLEV- KV
      IF(ICOR.EQ.1) BRUTE= 0
      I3= NDP
!     PNCN= DFLOAT(NCN-2)/DFLOAT(NCN+2)
      PNCN= DBLE(NCN-2)/DBLE(NCN+2)
      PP1= 1.d0/pNCN + 1.d0
c*** For Quasibound levels, first search inward to classically forbidden
      IF(EO.GT.VLIM) THEN
          PNCN= 1.d0
          PP1= 2.d0
          DO  I= NDP,1,-1
              I3= I
              IF(V(I).GT.EINT) GOTO 8
              ENDDO
          ENDIF
c*** First, search inward for outermost turning point
    8 DO  I= I3,1,-1
          I4= I
          IF(V(I).LT.EINT) GOTO 10
          ENDDO
c*** If never found an 'outer' turning point (e.g., above qbdd. barier)
c  then simply return with negative  DGDV2  as error flag
      RETURN
c... Now collect vibrational phase and its energy deriv. over outer well
   10 Y1= EINT- V(I4+1)
      Y2= EINT- V(I4)
      Y3= EINT- V(I4-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH2= 0.5d0*ARG2 + ANS2*SDRDY(I4)**2/RH
      DGDV2= 0.5d0/ARG2 + ANS1*SDRDY(I4)**2/RH
      DO  I= I4-2,1,-1
c... now, collect (v+1/2) and dv/dG integrals to next turning point ...
          II= I
          IF(V(I).GT.EINT) GO TO 12
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH2= VPH2+ ARG2*SDRDY(I)**2
          DGDV2= DGDV2+ SDRDY(I)**2/ARG2
          ENDDO
   12 I3= II+1
      Y1= EINT- V(I3-1)
      Y2= EINT- V(I3)
      Y3= EINT- V(I3+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH2= (VPH2 - ARG2 - 0.5d0*ARG3 + ANS2*SDRDY(I3)**2/RH)/Pi
      DGDV2= DGDV2 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1*SDRDY(I3)**2/RH
      DGDV2= Pi2/(BFCT*DGDV2)
c*** Next, search for innermost turning point
      DO  I= 1,NDP
          I1= I
          IF(V(I).LT.EINT) GOTO 20
c... then collect vibrational phase and its energy deriv. over outer well
          ENDDO
c
   20 IF(I1.EQ.1) THEN
          WRITE(6,602) JROT,EO
!         STOP
          CALL ABEND()
          ENDIF
      IF(I1.GE.I3) THEN
c*** For single-well potential or above barrier of double-well potential
          IF(IWR.GE.2) WRITE(6,600) ICOR,KV,JROT,EO,VPH2-0.5d0,DGDV2
          IF((KV.NE.(KVLEV-1)).AND.(DGDVB.GT.0.d0)) THEN
c... If got wrong level (KV not one below KVLEV) and NOT first call ...
              IF((EO-BMAX).GT.(2.d0*DGDV2)) THEN
c... 'Normal' case: use B-S plot area to estimate correct energy
!                 DE0= KVDIF*(DGDV2- 0.5d0*(DGDV2-DGDVB)/DFLOAT(KV-KVB))
                  DE0= KVDIF*(DGDV2- 0.5d0*(DGDV2-DGDVB)/DBLE(KV-KVB))
                  EO= EO+ DE0
                  KV= KVB
                  KVLEV= KV+1
                  RETURN
                ELSE
c... but close to barrier in double-well potential, switch to 'BRUTE'
                  BRUTE=BRUTE+ 1
                  DGDV1= DGDV2
                  XDIF= SIGN(1,KVDIF)
                  GOTO 54
                ENDIF
              ENDIF
          IF(KVLEV.EQ.0) THEN
c*** If looking for v=0, just use local DVDG2 to estimate energy correction
              EO= EO + KVDIF*DGDV2
              RETURN
              ENDIF
          IF(KV.EQ.0) THEN
c** Normally:  use B-S plot considerations to estimate next level energy
c... use harmonic estimate for v=1 energy
              EO= EO+ DGDV2
            ELSE
c... estimate Delta(G) based on linear Birge-Sponer
              DE0= 0.5d0*(3.d0*DGDV2 - DGDVB)
              IF((2.d0*DGDV2).GT.DGDVB) THEN
c... if linear Birge-Sponer predicts another level, then use it
                  EO= EO+ DE0
                ELSE
c... otherwise, use N-D theory extrapolation for next level...
                  DGDV2P= DGDV2**PNCN
                  DE0= (DGDV2P+DGDV2P-DGDVBP)
                  IF(DE0.GT.0.d0) THEN
                      DE0= (DE0**PP1- DGDV2P**PP1)/(PP1*(DGDV2P-DGDVBP))
                      EO= EO+ DE0
                    ELSE
c... but if NDT predicts no more levels, quit, and (optionally) print
                      IF(IWR.GT.0) WRITE(6,604) KV,EO
  604 FORMAT(10x,'Find highest bound level is   E(v=',i3,')=',1PD18.10)
                      RETURN
                    ENDIF
                ENDIF
            ENDIF
          DGDVB= DGDV2
          DGDVBP= DGDVB**PNCN
          KVB= KV
          INNER= 0
          RETURN
          ENDIF
c
c*** For a double-well potential, collect vibrational phase and its
c   energy derivative over the inner well
      Y1= EINT- V(I1-1)
      Y2= EINT- V(I1)
      Y3= EINT- V(I1+1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      ARG2= DSQRT(Y3)
      VPH1= 0.5d0*ARG2 + ANS2*SDRDY(I1)**2/RH
      DGDV1= 0.5d0/ARG2 + ANS1*SDRDY(I1)**2/RH
      DO  I= I1+2,NDP
c... now, collect integral and count nodes outward to next turning point ...
          IF(V(I).GT.EINT) GO TO 22
          ARG3= ARG2
          ARG2= DSQRT(EINT - V(I))
          VPH1= VPH1+ ARG2*SDRDY(I)**2
          DGDV1= DGDV1+ SDRDY(I)**2/ARG2
          ENDDO
   22 I2= I-1
      Y1= EINT- V(I2+1)
      Y2= EINT- V(I2)
      Y3= EINT- V(I2-1)
      CALL LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
      VPH1= (VPH1 - ARG2 - 0.5d0*ARG3 + ANS2*SDRDY(I2)**2/RH)/Pi
      DGDV1= DGDV1 -1.d0/ARG2 - 0.5d0/ARG3 + ANS1*SDRDY(I2)**2/RH
      DGDV1= Pi2/(BFCT*DGDV1)
      DGDVM= DGDV1*DGDV2/(DGDV1+DGDV2)
      WRITE(6,*) DGDVM ! make it "referenced"
      IF(KVDIF.EQ.0) THEN
c** If already at level sought, return
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
          RETURN
          ENDIF
c
c** Check whether looking for higher or lower level ...
      IDIF= SIGN(1,KVDIF)
      XDIF= IDIF
      IF((ICOR.GE.6).AND.((IABS(KVDIF).EQ.1).OR.(BRUTE.GT.0))) GOTO 50
c*** 'Conventional' semiclassical search for neared INNER or OUTER well level
      IF(INNER.LE.0) THEN
c... and current energy EO is for an outer-well level ...
          DE2= DGDV2*XDIF
          IV1= INT(VPH1+ 0.5d0)
!         DE1= (DFLOAT(IV1) + 0.5d0 - VPH1)*DGDV1*XDIF
          DE1= (DBLE(IV1) + 0.5d0 - VPH1)*DGDV1*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   30     IF(DABS(DE1).LT.DABS(DE2)) THEN
              INNER= 1
              EO= EO+ DE1
              DE1= DGDV1*XDIF
            ELSE
              INNER= 0
              EO= EO+ DE2
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 30
          ENDIF
      IF(INNER.GT.0) THEN
c... and current energy EO is for an inner-well level ...
          DE1= DGDV1*XDIF
          IV2= INT(VPH2+ 0.5d0)
!         DE2= (DFLOAT(IV2) + 0.5d0 - VPH2)*DGDV2*XDIF
          DE2= (DBLE(IV2) + 0.5d0 - VPH2)*DGDV2*XDIF
          IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   40     IF(DABS(DE2).LT.DABS(DE1)) THEN
              INNER= 0
              EO= EO+ DE2
              DE2= DGDV2*XDIF
            ELSE
              INNER= 1
              EO= EO+ DE1
            ENDIF
          KVDIF= KVDIF-IDIF
          IF(KVDIF.EQ.0) THEN
              RETURN
              ENDIF
          GOTO 40
          ENDIF
   50 BRUTE= BRUTE+ 1
c*** Now .. Brute force search for desired level !
      IF(IWR.GE.2) WRITE(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,
     1                                           ICOR,VPH2-0.5d0,DGDV2
   54 IF(BRUTE.EQ.1) THEN
c... in first brute-force step, use previous energy with opposite INNER
          EBRUTE= EO
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          DEBRUTE= DMIN1(DGDV1,DGDV2)*XDIF*0.3d0
          RETURN
          ENDIF
      IB= BRUTE/2
c... in subsequent even steps, lower EO by DEBRUTE/10 for same INNER
      IF((IB+IB).EQ.BRUTE) THEN
          EBRUTE= EBRUTE+ DEBRUTE
          EO= EBRUTE
          RETURN
        ELSE
c... in subsequent odd steps, lower repeat previous EO with INNER changed
          IF(INNER.EQ.0) THEN
              INNER= 1
            ELSE
              INNER= 0
            ENDIF
          EO= EBRUTE
          RETURN
        ENDIF
c     RETURN
  600 FORMAT(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',f10.2,
     1 '  v(SC)=',F8.3,'  dGdv=',f8.3)
  602 FORMAT(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2 )
  610 FORMAT(' Double well   E(v=',i3,', J=',I3,')=',f9.3,
     1 ':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,
     2 ' (ICOR=',I2,')',8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2 )
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
