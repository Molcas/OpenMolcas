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

subroutine AT34R(N,ISIZE,CHARGE,SMAT,V,H,EV2,MULT,BU,P,G,EIG,SINV,REVT,AUX,OVE,EW,E,AA,RR,TT,iprint,VEXTT,PVPT,EVN1,RE1R,AUXI,W1W1)
! INPUT: SMAT  OVERLAP MATRIX
!        V     POTENTIAL
!        H     RELATIVISTIC KINETIC ENERGY
!        EV2   PVP INTEGRALS

use DKH_Info, only: cLightAU, IRELMP
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, iprint
integer(kind=iwp), intent(inout) :: ISIZE
real(kind=wp), intent(in) :: CHARGE, SMAT(ISIZE)
real(kind=wp), intent(inout) :: V(ISIZE), EV2(ISIZE)
real(kind=wp), intent(inout) :: H(ISIZE), BU(ISIZE), P(ISIZE), G(ISIZE), EIG(N,N), SINV(N,N), REVT(N,N), AUX(N,N), OVE(N,N), &
                                EW(N), E(N), AA(N), RR(N), TT(N), VEXTT(ISIZE), PVPT(ISIZE), EVN1(N,N), RE1R(N,N), AUXI(N,N), &
                                W1W1(N,N)
integer(kind=iwp), intent(out) :: MULT(N)
integer(kind=iwp) :: I, IJ, J, K
real(kind=wp) :: CON, CON2, CR, PREA, RATIO, TV1, TV2, TV3, TV4, VELIT

!call PRMAT(u6,SMAT,N,0,'SMAT    ')
VELIT = cLightAU
ISIZE = N*(N+1)/2
PREA = 1/(VELIT*VELIT)
CON2 = PREA+PREA
CON = One/PREA
MULT(1) = 0
do I=1,N
  MULT(I+1) = MULT(I)+I
end do

! SCHMIDT-ORTHOGONALIZE OVERLAP MATRIX

call SOG(N,SMAT,SINV,P,OVE,EW)
call SQUARE(SMAT,OVE,1,N,N)

!-----------------------------------------------------------------------
! MATRIX REPRESENTATION CALCULATED FROM NONRELATIVISTIC T MATRIX
!-----------------------------------------------------------------------
call DIAG_DKH(H,N,EIG,EW,SINV,AUX,0)
if (iprint >= 10) then
  write(u6,*) ' eigenvalues in at34r'
  write(u6,*) (ew(i),i=1,n)
end if
do I=1,N

  ! IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATION

  TT(I) = EW(I)
  RATIO = EW(I)/VELIT
  if (RATIO > 0.02_wp) then
    EW(I) = CON*(sqrt(One+CON2*EW(I))-One)
  else
    TV1 = EW(I)
    TV2 = -TV1*EW(I)*PREA*Half
    TV3 = -TV2*EW(I)*PREA
    TV4 = -TV3*EW(I)*PREA*1.25_wp
    EW(I) = TV1+TV2+TV3+TV4
  end if
  E(I) = EW(I)+CON
end do
!-----------------------------------------------------------------------
! CALCULATE REVERSE TRANSFORMATION
!-----------------------------------------------------------------------

! CALCULATE TRANSFORMATION MATRICES

do I=1,N
  do J=1,N
    AUX(I,J) = Zero
    do K=I,N
      AUX(I,J) = AUX(I,J)+SINV(I,K)*EIG(K,J)
    end do
  end do
end do
do I=1,N
  do J=1,N
    REVT(I,J) = Zero
    do K=1,N
      REVT(I,J) = REVT(I,J)+OVE(I,K)*AUX(K,J)
    end do
  end do
end do
IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    H(IJ) = Zero
    do K=1,N
      H(IJ) = H(IJ)+REVT(I,K)*REVT(J,K)*EW(K)
    end do
  end do
end do
if (IRELMP /= 11) then
  do I=1,N
    AA(I) = sqrt((CON+E(I))/(Two*E(I)))
    RR(I) = sqrt(CON)/(CON+E(I))
  end do
else if (IRELMP == 11) then  ! RESC
  do I=1,N
    AA(I) = (sqrt(One+CON*TT(I)*Two/((CON+E(I))*(CON+E(I)))))/(CON+E(I))  ! O OPERATOR
    RR(I) = sqrt(CON)/(CON+E(I))  ! Q OPERATOR
  end do
end if

! BEYOND THIS POINT, OVE IS USED AS SCRATCH ARRAY

! TRANSFORM V TO T-BASIS

call TRSM_DKH(V,SINV,G,N,AUX,OVE)
call TRSM_DKH(G,EIG,BU,N,AUX,OVE)

!MULTIPLY

if (IRELMP /= 11) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      P(IJ) = -BU(IJ)*CHARGE
      VEXTT(IJ) = P(IJ) ! KEEP T-BASIS VEXT INTO VEXTT FOR HIGHER-ORDER DK
      BU(IJ) = P(IJ)*AA(I)*AA(J)
      EVN1(I,J) = BU(IJ)
      EVN1(J,I) = EVN1(I,J)
    end do
  end do

else if (IRELMP == 11) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      P(IJ) = -BU(IJ)*CHARGE
      BU(IJ) = VELIT*P(IJ)*(sqrt(RR(I)*RR(J))*AA(I)/AA(J)+sqrt(RR(J)*RR(I))*AA(J)/AA(I))
    end do
  end do

end if

call TRSMT(BU,REVT,V,N,AUX,OVE)

! PVP INTEGRALS

call TRSM_DKH(EV2,SINV,G,N,AUX,OVE)
call TRSM_DKH(G,EIG,BU,N,AUX,OVE)

! MULTIPLY

if (IRELMP /= 11) then
  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      G(IJ) = -BU(IJ)*CHARGE
      PVPT(IJ) = G(IJ) ! KEEP T-BASIS PVP INTO PVPT FOR HIGHER-ORDER DK
      BU(IJ) = G(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
      EVN1(I,J) = EVN1(I,J)+BU(IJ)
      EVN1(J,I) = EVN1(I,J)
    end do
  end do
else if (IRELMP == 11) then
  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      G(IJ) = -BU(IJ)*CHARGE
      BU(IJ) = G(IJ)*(RR(I)*RR(J)*AA(I)/AA(J)+RR(J)*RR(I)*AA(J)/AA(I))*Half
    end do
  end do
end if
call TRSMT(BU,REVT,EV2,N,AUX,OVE)
!call PRMAT(6,EV2,N,0,'PVPFULL ')
V(1:ISIZE) = V(1:ISIZE)+EV2(1:ISIZE)

if ((IRELMP /= 1) .and. (IRELMP /= 11)) then

  ! CALCULATE EVEN2 OPERATOR

  !call PRMAT(u6,TT,N,1,'TT      ')
  !call PRMAT(u6,E,N,1,'E       ')
  call EVEN2(N,P,G,E,AA,RR,TT,EIG,AUX,OVE,W1W1)

  ! TRANSFORM BACK

  call TRSMT(G,REVT,EV2,N,AUX,OVE)
  V(1:ISIZE) = V(1:ISIZE)+EV2(1:ISIZE)

  if ((IRELMP /= 0) .and. (IRELMP /= 2)) then  ! DK2

    ! CALCULATE EVEN3 OPERATOR

    call EVEN3(N,P,G,E,AA,RR,TT,EIG,AUX,OVE,EVN1,VEXTT,PVPT,RE1R,W1W1,AUXI)

    ! TRANSFORM BACK FOR DK3

    call TRSMT(G,REVT,EV2,N,AUX,OVE)
    V(1:ISIZE) = V(1:ISIZE)+EV2(1:ISIZE)

    if (IRELMP /= 3) then  ! DK3

      ! More to come here

    end if
  end if
end if

CR = 1/CHARGE
do I=1,ISIZE
  V(I) = -V(I)*CR
end do
!call PRMAT(u6,BU,N,0,'TOTAL H ')

return

end subroutine AT34R
