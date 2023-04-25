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

!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
subroutine MATXEL(KV1,JROT1,IOMEG1,EO1,KV2,JROT2,IOMEG2,IRFN,EO2,NBEG,NEND,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
!** Subroutine to calculate matrix elements of powers of the distance
!  coordinate between vib. eigenfunction WF1(i) for v=KV1, J=JROT1 of
!  potential-1 & WF2(I), corresponding to KV2 & JROT2 of potentl.-2

!use stdalloc, only: mma_allocate, mma_deallocate

integer I, J, IOMEG1, IOMEG2, IOMUP, IOMLW, IRFN, JROT1, JROT2, JUP, JLW, KV1, KV2, KVUP, KVLW, LXPCT, NBEG, NEND, MORDR, NDIMR
real*8 ZMAT(0:20), WF1(NEND), WF2(NEND), RFN(NEND), DM(0:MORDR), DRDY2(NDIMR)
real*8 AEINST, DEG, DME, DSM, EO1, EO2, ELW, FCF, FREQ, OMUP, RH, RI, SJ, ZJUP
character*1 DJ(-3:3)
data DJ/'N','O','P','Q','R','S','T'/

ZMAT(0) = 0.d0
if (MORDR >= 1) then
  do J=1,MORDR
    ZMAT(J) = 0.d0
  end do
end if
if (IRFN /= -4) then
  ! For regular power series or function matrix elements ...
  do I=NBEG,NEND
    DSM = WF2(I)*WF1(I)*DRDY2(I)
    ZMAT(0) = ZMAT(0)+DSM
    RI = RFN(I)
    if (MORDR >= 1) then
      do J=1,MORDR
        DSM = DSM*RI
        ZMAT(J) = ZMAT(J)+DSM
      end do
    end if
  end do
else
  ! For partial derivative matrix elements ...
  do I=NBEG+1,NEND-1
    DSM = WF1(I)*(WF2(I+1)-WF2(I-1))*DRDY2(I)
    ZMAT(0) = ZMAT(0)+DSM
    RI = RFN(I)
    if (MORDR >= 1) then
      do J=1,MORDR
        DSM = DSM*RI
        ZMAT(J) = ZMAT(J)+DSM
      end do
    end if
  end do
  do J=0,MORDR
    ZMAT(J) = ZMAT(J)/(2.d0*RH)
  end do
end if
DME = 0.d0
FCF = (ZMAT(0)*RH)**2
if (MORDR >= 0) then
  do J=0,MORDR
    ZMAT(J) = ZMAT(J)*RH
    DME = DME+DM(J)*ZMAT(J)
  end do
end if
FREQ = EO2-EO1
ELW = dmin1(EO1,EO2)
! Now calculate the Honl-London Factor for the particular transition
! Factors updated as per Hansson & Watson JMS (2005).
SJ = 0.d0
KVUP = KV1
KVLW = KV2
JUP = JROT1
JLW = JROT2
IOMUP = max(IOMEG1,0)
IOMLW = max(IOMEG2,0)
if (EO2 > EO1) then
  KVUP = KV2
  KVLW = KV1
  JUP = JROT2
  JLW = JROT1
  IOMUP = max(IOMEG2,0)
  IOMLW = max(IOMEG1,0)
end if
ZJUP = JUP
OMUP = IOMUP
DEG = 2*JUP+1
if ((JLW < IOMLW) .or. (JUP < IOMUP)) go to 50
if (IOMUP == IOMLW) then
  ! Factors for  DELTA(LAMBDA) = 0  transitions of spin singlets
  if (JUP == (JLW+1)) SJ = (ZJUP+OMUP)*(JUP-IOMUP)/ZJUP
  if ((JUP == JLW) .and. (JUP > 0)) SJ = DEG*OMUP**2/(ZJUP*(ZJUP+1.d0))
  if (JUP == (JLW-1)) SJ = (ZJUP+1.d0+OMUP)*(JUP+1-IOMUP)/(ZJUP+1.d0)
end if
if (IOMUP == (IOMLW+1)) then
  ! Factors for  DELTA(LAMBDA) = +1  transitions of spin singlets
  if (JUP == (JLW+1)) SJ = (ZJUP+OMUP)*(JUP-1+IOMUP)/(2.d0*ZJUP)
  if ((JUP == JLW) .and. (JUP > 0)) SJ = (ZJUP+OMUP)*(JUP+1-IOMUP)*DEG/(2.d0*ZJUP*(ZJUP+1.d0))
  if (JUP == (JLW-1)) SJ = (JUP+1-IOMUP)*(ZJUP+2.d0-OMUP)/(2.d0*ZJUP+2.d0)
end if
if (IOMUP < IOMLW) then
  ! Factors for  DELTA(LAMBDA) = -1  transitions of spin singlets
  if (JUP == (JLW+1)) SJ = (JUP-IOMUP)*(JUP-1-IOMUP)/(2.d0*ZJUP)
  if ((JUP == JLW) .and. (JUP > 0)) SJ = (JUP-IOMUP)*(ZJUP+1.d0+OMUP)*DEG/(2.d0*ZJUP*(ZJUP+1.d0))
  if (JUP == (JLW-1)) SJ = (ZJUP+1.d0+OMUP)*(ZJUP+2.d0+OMUP)/(2.d0*ZJUP+2.d0)
end if
!... finally, include Hansson-Watson  w0/w1  term in Honl-London factor
if ((min(IOMUP,IOMLW) == 0) .and. (IOMUP /= IOMLW)) SJ = SJ+SJ

! For FREQ in  cm-1  and dipole moment in  debye , AEINST is the
! absolute Einstein radiative emission rate (s-1) , using the
! rotational intensity factors for sigma-sigma transitions.
50 continue
AEINST = dabs(3.1361891D-7*dabs(FREQ)**3*DME**2*SJ/DEG)
if (LXPCT > 0) then
  write(6,600) KV1,JROT1,EO1,KV2,JROT2,EO2
  if (iabs(IRFN) <= 9) write(6,602) (J,ZMAT(J),J=0,MORDR)
  write(6,604) FCF,DME,FREQ,AEINST
  write(6,606)
end if
if ((iabs(LXPCT) == 4) .or. (iabs(LXPCT) == 5) .and. (SJ > 0.d0)) then
  if (iabs(JUP-JLW) <= 3) write(8,801) DJ(JUP-JLW),JLW,KVUP,KVLW,ELW,FREQ,AEINST,FCF,DME
  !... Special printout for Hui/LeRoy N2 Quadrupole paper [JCP 1XX (2007)]
  !E00 = 1175.7693d0
  !write(11,811) -FREQ,KVUP,JUP,KVLW,JLW,-FREQ,ELW-FREQ-E00,ELW-E00,DME**2
  !811 format(F12.4,2I4,I6,I4,3f12.4,1PD15.6)
  if (iabs(JUP-JLW) > 3) write(8,802) JUP-JLW,JLW,KVUP,KVLW,ELW,FREQ,AEINST,FCF,DME
end if
if (iabs(LXPCT) >= 5) write(7,701) KVUP,JUP,KVLW,JLW,FREQ,(ZMAT(J),J=0,MORDR)
!if (iabs(LXPCT) >= 5) write(7,701) KVUP,JUP,KVLW,JLW,(ZMAT(J),J=0,MORDR)

return

600 format(' Coupling   E(v=',I3,', J=',I3,')=',F12.4,'   to   E(v=',I3,', J=',I3,')=',F12.4)
602 format(5x,'Moment matrix elements:',2('   <X**',I2,'>=',F14.10:),1x/(3x,3('   <X**',I2,'>=',F14.10:),1x))
604 format(' FCF=',1PD11.4,'   <M>=',D12.5,'   d(E)=',0PF10.2,'   A(Einst)=',1PD11.4,' s-1')
606 format(1X,79('+'))
701 format(4I4,F12.4,4F12.8:/(4X,6F12.8))
!701 format(4I4,6F12.8:/(16X,6F12.8))
801 format(1x,A1,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
802 format(i2,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))

end subroutine MATXEL
