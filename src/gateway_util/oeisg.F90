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

subroutine OEISG(REL,SREL,TREL,UREL,ZETA,ZN,NSYM,NBAS1,unrel,tnrel,hcorr,iprint,VEXTT,PVPT,EVN1,RE1R,AUXI,W1W1,iScratch,ScratchZ, &
                 ScratchSQ,Scratch)

use Constants, only: One, Two, Eight, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NSYM, NBAS1, iprint
real(kind=wp), intent(out) :: REL(nBAS1*(nBAS1+1)/2), SREL(nBAS1*(nBAS1+1)/2), TREL(nBAS1*(nBAS1+1)/2), UREL(nBAS1*(nBAS1+1)/2), &
                              unrel(nBAS1*(nBAS1+1)/2), tnrel(nBAS1*(nBAS1+1)/2), hcorr(nBAS1*(nBAS1+1)/2), VEXTT(*), PVPT(*), &
                              EVN1(NBAS1,NBAS1), RE1R(NBAS1,NBAS1), AUXI(NBAS1,NBAS1), W1W1(NBAS1,NBAS1), &
                              ScratchZ(nBAS1*(nBAS1+1)/2,3), ScratchSQ(NBAS1*NBAS1,5), Scratch(NBAS1,5)
integer(kind=iwp), intent(out) :: iScratch(NBAS1*(NBAS1+1)/2)
real(kind=wp), intent(in) :: ZETA(nBAS1), ZN
integer(kind=iwp) :: I, I1, I2, I3, IJ, IM1, iparm, J, J1, J2, J3, L, NBSZ, np, NPQ, NPQ1, NPQ2, NQ
real(kind=wp) :: facto(16), TERM1, TERM2, VP, VPQ, VPQ1, VPQM2, VPQP2, VQ, W_P, WQ, ZP, ZPQ, ZQ
real(kind=wp), parameter :: SQRT8PI = sqrt(Eight/Pi)
real(kind=wp), external :: EXTC

! ONE CONFIGURATION
! ONE ELECTRON INTEGRAL PROGRAM
! COMPUTES THE MATRICES S,U,T WITH SLATER (NFLAG1=0) OR GAUSSIAN
!  (NFLAG1=1) ORBITALS AS  BASIS

! initialize
call RELOP()
! coefficients
FACTO(1) = One
FACTO(2) = One
do I=3,16
  IM1 = I-1
  FACTO(I) = IM1*FACTO(I-2)
end do
iparm = 2
!write(u6,*) ' MX100=',MX100
!write(u6,*) ' nsym',nsym
if (iprint >= 10) then
  write(u6,*) ' symmetry',nsym
  write(u6,*) ' number of basis functions',nbas1
  write(u6,*) ' charge',zn
  do i=1,nbas1
    write(u6,*) zeta(i)
  end do
end if
L = NSYM

IJ = 0
np = L
NQ = NP
do I=1,NBAS1
  ZP = ZETA(I)
  W_P = Two*(NP-L)/ZP
  do J=1,I
    IJ = IJ+1
    ZQ = ZETA(J)
    ZPQ = Half*(ZP+ZQ)
    NPQ = NP+NQ+1
    WQ = Two*(NQ-L)/ZQ
    NPQ1 = NPQ-1
    NPQ2 = NPQ-2
    VPQ = FACTO(NPQ1)/ZPQ**(Half*NPQ)
    VP = FACTO(2*NP)/ZP**(NP+Half)
    VQ = FACTO(2*NQ)/ZQ**(NQ+Half)
    if (NPQ1 > 2) then
      VPQM2 = FACTO(NPQ2-1)/ZPQ**(Half*NPQ2)
    else
      VPQM2 = One
    end if
    VPQ1 = FACTO(NPQ2)*SQRT8PI/ZPQ**(Half*NPQ1)
    VPQP2 = FACTO(NPQ+1)/ZPQ**(Half*(NPQ+2))
    TERM2 = VPQP2-VPQ*(W_P+WQ)+W_P*WQ*VPQM2
    TERM1 = One/sqrt(VP*VQ)
    SREL(IJ) = TERM1*VPQ
    UREL(IJ) = TERM1*VPQ1
    UNREL(IJ) = TERM1*VPQ1
    TNREL(IJ) = Half*ZP*ZQ*TERM1*TERM2
    !TNRE = Half*ZP*ZQ*TERM1*TERM2
    if (IPARM > 0) then
      I1 = np-1
      I2 = 0
      I3 = 0
      J1 = nq-1
      J2 = 0
      J3 = 0

      if (np == 3) then
        ! PVP FOR XY-FUNCTIONS

        REL(IJ) = EXTC(L,ZP,ZQ,1,1,0,1,1,0)

      else if (np == 4) then
        ! PVP FOR XYZ-FUNCTIONS

        REL(IJ) = EXTC(L,ZP,ZQ,1,1,1,1,1,1)

      else
        ! PVP FOR X**(NP-1)-FUNCTIONS

        REL(IJ) = EXTC(L,ZP,ZQ,I1,I2,I3,J1,J2,J3)

      end if

      !TREL(IJ) = SQROPY(ZP,ZQ,I1,I2,I3,J1,J2,J3)
      TREL(IJ) = Half*ZP*ZQ*TERM1*TERM2
    else
      TREL(IJ) = Half*ZP*ZQ*TERM1*TERM2
      !TNRE = Half*ZP*ZQ*TERM1*TERM2
    end if
  end do
end do

! CONVERT TO RELATIVISTIC INTEGRALS

NBSZ = NBAS1*(NBAS1+1)/2
!  OVERLAP  POTENTIAL  KIN.ENERGY  PVP  MULT  BU  P  G  EIG  SINV  REVT  AUX  OVE  EW  E  AA  RR  TT
call AT34R(NBAS1,NBSZ,ZN,SREL,UREL,TREL,REL,iScratch,ScratchZ(:,1),ScratchZ(:,1),ScratchZ(:,3),ScratchSQ(:,1),ScratchSQ(:,2), &
           ScratchSQ(:,3),ScratchSQ(:,4),ScratchSQ(:,5),Scratch(:,1),Scratch(:,2),Scratch(:,3),Scratch(:,4),Scratch(:,5),iprint, &
           VEXTT,PVPT,EVN1,RE1R,AUXI,W1W1)

! TRANSFER RELATIVISTIC KINETIC ENERGY AND POTENTIAL INTEGRALS

if (iprint >= 10) then
  write(u6,*) ' matrices'
  write(u6,*) l,nbas1
  write(u6,*) ' srel'
  write(u6,12) (srel(j),j=1,(nbas1*(nbas1+1))/2)
  write(u6,*) ' trel'
  write(u6,12) (trel(j),j=1,(nbas1*(nbas1+1))/2)
  write(u6,*) ' urel'
  write(u6,12) (urel(j),j=1,(nbas1*(nbas1+1))/2)
  write(u6,*) ' tnrel'
  write(u6,12) (tnrel(j),j=1,(nbas1*(nbas1+1))/2)
  write(u6,*) ' unrel'
  write(u6,12) (unrel(j),j=1,(nbas1*(nbas1+1))/2)
  write(u6,*) ' rel'
  write(u6,12) (rel(j),j=1,(nbas1*(nbas1+1))/2)
end if
ij = 0
do i=1,nbas1
  do j=1,i
    ij = ij+1
    hcorr(ij) = (trel(ij)-zn*urel(ij))-(tnrel(ij)-zn*unrel(ij))
    trel(ij) = trel(ij)-tnrel(ij)
    urel(ij) = -(zn*urel(ij)-zn*unrel(ij))
  end do
end do
if (iprint >= 20) then
  write(u6,*) ' full correction metrix'
  write(u6,*) l,nbas1
  write(u6,12) (hcorr(j),j=1,(nbas1*(nbas1+1))/2)
end if

return

12 format(4f18.14)

end subroutine OEISG
