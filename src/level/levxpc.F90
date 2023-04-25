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
subroutine LEVXPC(KV,JR,EPR,GAMA,NPP,WF,RFN,V,VLIM,YH,DREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Calculates expectation values of the kinetic energy and of X**IP
!  (IP=1,MORDR), denoted XPTKE and XPCTR(IP), respectively, for level
!  v=KV, J=JR, E=EPR(cm-1), using wave function WF(i), (i=NBEG,NEND).
!** Assumes units of length are (Angstroms) .
!** Division by BFCT converts potential V(I) to units (cm-1).
!** If (|LXPCT| = 2  or  4), "punch" (WRITE(7,XXX)) results to channel-7
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: KV, JR, NPP, NBEG, NEND, LXPCT, MORDR, IRFN
real(kind=wp) :: EPR, GAMA, WF(NPP), RFN(NPP), V(NPP), VLIM, YH, DREF, DM(0:20), BFCT
integer(kind=iwp) :: I, IPNCH, ITRY, K
real(kind=wp) :: DER, DMR, DRT, DS, EINN, RR, RXPCT, SF2, SS2, XPCTR(0:11), XPTKE
real(kind=wp), allocatable :: DRDY2(:), FAS(:), RVB(:), SDRDY(:), VBZ(:), YVB(:)
integer(kind=iwp), parameter :: NDIMR = 200001

call mma_allocate(RVB,NDIMR,LABEL='RVB')
call mma_allocate(YVB,NDIMR,LABEL='YVB')
call mma_allocate(DRDY2,NDIMR,LABEL='DRDY2')
call mma_allocate(FAS,NDIMR,LABEL='FAS')
call mma_allocate(SDRDY,NDIMR,LABEL='SDRDY')
call mma_allocate(VBZ,NDIMR,LABEL='VBZ')
EINN = BFCT*EPR
IPNCH = 0
if ((abs(LXPCT) == 2) .or. (abs(LXPCT) >= 4)) IPNCH = 1
! MORDR is the highest-power expectation value considered.
if (MORDR > 11) MORDR = 11
ITRY = 20
if (((IRFN == -1) .or. ((IRFN >= 1) .and. (IRFN <= 9))) .and. (DREF <= Zero)) ITRY = 0
! Start by calculating contributions at end points
2 continue
SS2 = WF(NBEG)**2*DRDY2(NBEG)
SF2 = WF(NEND)**2*DRDY2(NEND)
XPTKE = Half*(SS2*(EINN-V(NBEG))+SF2*(EINN-V(NEND)))
if (MORDR > 0) then
  XPCTR(0) = One/YH
  do K=1,MORDR
    SS2 = SS2*RFN(NBEG)
    SF2 = SF2*RFN(NEND)
    XPCTR(K) = Half*(SS2+SF2)
  end do
end if
if (IRFN > -4) then
  ! For regular expectation values of a radial function ...
  do I=NBEG+1,NEND-1
    DS = WF(I)**2*DRDY2(I)
    XPTKE = XPTKE+DS*(EINN-V(I))
    if (MORDR > 0) then
      RR = RFN(I)
      do K=1,MORDR
        DS = DS*RR
        XPCTR(K) = XPCTR(K)+DS
      end do
    end if
  end do
else
  ! For expectation values involving partial derivative operator ...
  do K=0,MORDR
    XPCTR(K) = Zero
  end do
  do I=NBEG+1,NEND-1
    DS = WF(I)**2*DRDY2(I)
    XPTKE = XPTKE+DS*(EINN-V(I))
    DS = WF(I)*(WF(I+1)-WF(I-1))*DRDY2(I)
    if (MORDR > 0) then
      RR = RFN(I)
      do K=1,MORDR
        DS = DS*RR
        XPCTR(K) = XPCTR(K)+DS
      end do
    end if
  end do
  do K=0,MORDR
    XPCTR(K) = XPCTR(K)/(Two*YH)
  end do
end if
XPTKE = XPTKE*YH/BFCT
if (MORDR < 0) go to 99
DMR = Zero
do K=0,MORDR
  XPCTR(K) = XPCTR(K)*YH
  DMR = DMR+DM(K)*XPCTR(K)
end do
if ((LXPCT == 1) .or. (abs(LXPCT) == 2)) then
  if (EPR <= VLIM) write(u6,600) KV,JR,EPR,DMR,XPTKE
  if (EPR > VLIM) write(u6,602) KV,JR,EPR,DMR,XPTKE,GAMA
  if (abs(IRFN) <= 9) write(u6,604) (K,XPCTR(K),K=1,MORDR)
  if (IPNCH >= 1) write(7,701) KV,JR,EPR,GAMA,XPTKE,DMR,(XPCTR(K),K=1,MORDR)
end if
if (ITRY > 19) go to 99
! If appropriate, iteratively correct DREF value till distance
! coordinate expectation value is identically zero.
if (IRFN == -1) then
  ! For Dunham expansion parameter, define revised function here
  DREF = XPCTR(1)
  DRT = DREF
  write(u6,603) ITRY,DRT,DREF
  do I=1,NPP
    RVB(I) = RVB(I)/DREF-One
  end do
  ITRY = 99
  go to 2
end if
! For Surkus-type expansion parameter, define revised function
ITRY = ITRY+1
if (ITRY == 1) then
  RXPCT = XPCTR(1)
  DREF = Zero
  DRT = RXPCT
else
  DER = -IRFN/(Two*DREF)
  DRT = -XPCTR(1)/DER
end if
DREF = DREF+DRT
write(u6,603) ITRY,DRT,DREF
! Redefine Surkus-type distance variable RFN using new DREF
do I=1,NPP
  RFN(I) = (RVB(I)**IRFN-DREF**IRFN)/(RVB(I)**IRFN+DREF**IRFN)
end do
if (abs(DRT/DREF) >= 1.0e-12_wp) go to 2
call mma_deallocate(RVB)
call mma_deallocate(YVB)
call mma_deallocate(DRDY2)
call mma_deallocate(FAS)
call mma_deallocate(SDRDY)
call mma_deallocate(VBZ)

99 continue
return

600 format(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,'   <KE>=',F11.3)
602 format(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,'   <KE>=',F11.3/'   Tunneling predissociation  Width(FWHM)=', &
           G13.6,'    <X**',I2,'>=',F13.8)
604 format((8x,3('   <X**',I2,'>=',F13.8:)))
603 format(' On iteration #',I2,'  change DREF by',1PD10.2,'  to   DREF=',0PF13.10,' [Angstroms]')
701 format(2I4,F11.3,G11.4,F11.3,3(F12.7)/(5X,6F12.7))

end subroutine LEVXPC
