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
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE MATXEL(KV1,JROT1,IOMEG1,EO1,KV2,JROT2,IOMEG2,IRFN,EO2, &
     &  NBEG,NEND,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
!     USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
!** Subroutine to calculate matrix elements of powers of the distance
!  coordinate between vib. eigenfunction WF1(i) for v=KV1, J=JROT1 of
!  potential-1 & WF2(I), corresponding to KV2 & JROT2 of potentl.-2
      INTEGER I,J,IOMEG1,IOMEG2,IOMUP,IOMLW,IRFN,JROT1,JROT2,JUP,JLW,   &
     & KV1,KV2,KVUP,KVLW,LXPCT,NBEG,NEND,MORDR,NDIMR
      REAL*8  ZMAT(0:20),WF1(NEND),WF2(NEND),RFN(NEND),DM(0:MORDR),     &
     &  DRDY2(NDIMR)
      REAL*8  AEINST,DEG,DME,DSM,EO1,EO2,ELW,FCF,FREQ,OMUP,RH,RI,SJ,    &
     &    ZJUP
      CHARACTER*1  DJ(-3:3)
      DATA DJ/'N','O','P','Q','R','S','T'/
      ZMAT(0)= 0.D0
      IF(MORDR.GE.1) THEN
          DO  J= 1,MORDR
              ZMAT(J)= 0.D0
              ENDDO
          ENDIF
      IF(IRFN.NE.-4) THEN
!** For regular power series or function matrix elements ...
          DO  I=NBEG,NEND
              DSM= WF2(I)*WF1(I) * DRDY2(I)
              ZMAT(0)= ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J= 1,MORDR
                      DSM= DSM*RI
                      ZMAT(J)= ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
        ELSE
!** For partial derivative matrix elements ...
          DO  I= NBEG+1,NEND-1
              DSM= WF1(I)*(WF2(I+1)- WF2(I-1)) * DRDY2(I)
              ZMAT(0)= ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J= 1,MORDR
                      DSM= DSM*RI
                      ZMAT(J)= ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
          DO  J= 0,MORDR
              ZMAT(J)= ZMAT(J)/(2.d0*RH)
              ENDDO
        ENDIF
      DME= 0.D0
      FCF= (ZMAT(0)*RH)**2
      IF(MORDR.GE.0) THEN
          DO  J= 0,MORDR
              ZMAT(J)= ZMAT(J)*RH
              DME= DME+DM(J)*ZMAT(J)
              ENDDO
          ENDIF
      FREQ= EO2-EO1
      ELW= DMIN1(EO1,EO2)
!** Now calculate the Honl-London Factor for the particular transition
!   Factors updated as per Hansson & Watson JMS (2005).
      SJ= 0.D0
      KVUP= KV1
      KVLW= KV2
      JUP= JROT1
      JLW= JROT2
      IOMUP= MAX(IOMEG1,0)
      IOMLW= MAX(IOMEG2,0)
      IF(EO2.GT.EO1) THEN
          KVUP= KV2
          KVLW= KV1
          JUP= JROT2
          JLW= JROT1
          IOMUP= MAX(IOMEG2,0)
          IOMLW= MAX(IOMEG1,0)
          ENDIF
      ZJUP= JUP
      OMUP= IOMUP
      DEG= 2*JUP+ 1
      IF((JLW.LT.IOMLW).OR.(JUP.LT.IOMUP)) GO TO 50
      IF(IOMUP.EQ.IOMLW) THEN
!** Factors for  DELTA(LAMBDA) = 0  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+ OMUP)*(JUP- IOMUP)/ZJUP
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))                               &
     &                       SJ= DEG*OMUP**2/(ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1)) SJ= (ZJUP+1.D0+OMUP)*(JUP+1-IOMUP)/        &
     &                                                     (ZJUP+1.D0)
          ENDIF
      IF(IOMUP.EQ.(IOMLW+1)) THEN
!** Factors for  DELTA(LAMBDA) = +1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+OMUP)*(JUP-1+IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))                               &
     &       SJ= (ZJUP+OMUP)*(JUP+1-IOMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1))                                            &
     &       SJ= (JUP+1-IOMUP)*(ZJUP+2.D0-OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF
      IF(IOMUP.LT.IOMLW) THEN
!** Factors for  DELTA(LAMBDA) = -1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (JUP-IOMUP)*(JUP-1-IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))                               &
     &      SJ= (JUP-IOMUP)*(ZJUP+1.D0+OMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1))                                            &
     &           SJ= (ZJUP+1.D0+OMUP)*(ZJUP+2.D0+OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF
!... finally, include Hansson-Watson  w0/w1  term in Honl-London factor
      IF((MIN(IOMUP,IOMLW).EQ.0).and.(IOMUP.NE.IOMLW)) SJ= SJ+SJ
!
!** For FREQ in  cm-1  and dipole moment in  debye , AEINST is the
!  absolute Einstein radiative emission rate (s-1) , using the
!  rotational intensity factors for sigma-sigma transitions.
   50 CONTINUE
      AEINST = DABS(3.1361891D-7 *DABS(FREQ)**3*DME**2 * SJ/DEG)
      IF(LXPCT.GT.0) THEN
          WRITE(6,600) KV1,JROT1,EO1,KV2,JROT2,EO2
          IF(IABS(IRFN).LE.9) WRITE(6,602) (J,ZMAT(J),J= 0,MORDR)
          WRITE(6,604) FCF,DME,FREQ,AEINST
          WRITE(6,606)
          ENDIF
      IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5).AND.(SJ.GT.0.D0)) THEN
          IF(IABS(JUP-JLW).LE.3) WRITE(8,801) DJ(JUP-JLW),JLW,KVUP,     &
     &                                    KVLW,ELW,FREQ,AEINST,FCF,DME
!... Special printout for Hui/LeRoy N2 Quadrupole paper [JCP 1XX (2007)]
!cc       E00= 1175.7693d0
!c        WRITE(11,811) -FREQ,KVUP,JUP,KVLW,JLW,-FREQ,ELW-FREQ-E00,
!c   1                                      ELW-E00,DME**2
!c811 FORMAT(F12.4,2I4,I6,I4,3f12.4,1PD15.6)
          IF(IABS(JUP-JLW).GT.3) WRITE(8,802) JUP-JLW,JLW,KVUP,         &
     &                                    KVLW,ELW,FREQ,AEINST,FCF,DME
          ENDIF
      IF(IABS(LXPCT).GE.5)                                              &
!    1         WRITE(7,701) KVUP,JUP,KVLW,JLW,(ZMAT(J),J=0,MORDR)
     &         WRITE(7,701) KVUP,JUP,KVLW,JLW,FREQ,(ZMAT(J),J=0,MORDR)
      RETURN
  600 FORMAT(' Coupling   E(v=',I3,', J=',I3,')=',F12.4,'   to   E(v=', &
     & I3,', J=',I3,')=',F12.4)
  602 FORMAT(5x,'Moment matrix elements:',2('   <X**',I2,'>=',F14.10:), &
     &  1x/(3x,3('   <X**',I2,'>=',F14.10:),1x))
  604 FORMAT(' FCF=',1PD11.4,'   <M>=',D12.5,'   d(E)=',0PF10.2,        &
     &  '   A(Einst)=',1PD11.4,' s-1')
  606 FORMAT(1X,79('+'))
  701 FORMAT(4I4,F12.4,4F12.8:/(4X,6F12.8))
! 701 FORMAT(4I4,6F12.8:/(16X,6F12.8))
  801 FORMAT(1x,A1,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
  802 FORMAT(i2,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
      END
