!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1995, Bernd Artur Hess                                 *
!***********************************************************************

subroutine SCFCLI2(idbg,S,H,V,PVP,N,ISIZE,VELIT,BU,P,G,EV2,EIG,SINV,REVT,AUX,OVE,EW,E,AA,RR,TT)
! $Id: relsewb.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de

use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: idbg, N, ISIZE
real(kind=wp), intent(in) :: S(ISIZE), V(ISIZE), PVP(ISIZE), VELIT
real(kind=wp), intent(inout) :: H(ISIZE)
real(kind=wp), intent(out) :: BU(ISIZE), P(ISIZE), G(ISIZE), EIG(N,N), EV2(ISIZE), SINV(N,N), REVT(N,N), AUX(N,N), OVE(N,N), &
                              EW(N), E(N), AA(N), RR(N), TT(N)
integer(kind=iwp) :: I, icontr, iex, ii, IJ, j, K
real(kind=wp) :: CON, CON2, det, dtol, PREA, RATIO, TOL, TV1, TV2, TV3, TV4

TOL = 1.0e-14_wp
PREA = 1/(VELIT*VELIT)
CON2 = PREA+PREA
CON = One/PREA
do I=1,ISIZE
  BU(I) = S(I)
end do
ii = 0
do i=1,n
  do j=1,i
    ii = ii+1
    aux(I,J) = s(ii)
    aux(J,I) = s(ii)
  end do
end do
!write(u6,*) 'OVERLAP MATRIX'
!do i=1,n
!  write(u6,'(5f10.5)') (aux(i,j),j=1,n)
!end do

! CALCULATE DETERMINANT

icontr = -1
dtol = tol
call dcopiv(aux,aux,n,1,n,dtol,det,iex,icontr,p)
!if (idbg > 0) write(idbg,2016) icontr,det,iex
if (icontr /= 0) then
  !write(u6,2016) icontr,det,iex
  !write(u6,2012) dtol
  write(u6,'(a)') ' relsewb| singular overlap matrix'
  call Abend()
end if

! SCHMIDT-ORTHOGONALIZE

call Sogr(iDbg,N,BU,SINV,P,OVE,EW)

call Square(BU,OVE,N,1,N)
!call FilMar(N,BU,OVE)

! ** SINV CONTAINS TRANSFORMATION TO ORTHOGONAL AO-BASIS
! ** OVE  CONTAINS OVERLAP MATRIX IN FULL
!-----------------------------------------------------------------------
! NON-RELATIVISTIC KINETIC ENERGY
!-----------------------------------------------------------------------

call Diagr(H,N,EIG,EW,SINV,AUX,BU)

!if (idbg > 0) write(idbg,556)
do I=1,N
  if (ew(i) < Zero) then
    !write(u6,*) ' scfcli2| ew(',i,') = ',ew(i)
    write(u6,'(a)') 'kinetic energy eigenvalue less than zero'
    call Abend()
  end if

  ! IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO

  RATIO = EW(I)/VELIT

  ! CALCULATE RELATIVISTIC ENERGY AND MOMENTUM

  !SR = sqrt(Two*EW(I))
  TT(I) = EW(I)
  if (RATIO > 0.02_wp) then
    TV1 = EW(I)
    EW(I) = CON*(sqrt(One+CON2*EW(I))-One)
    !if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR
  else
    TV1 = EW(I)
    TV2 = -TV1*EW(I)*PREA*Half
    TV3 = -TV2*EW(I)*PREA
    TV4 = -TV3*EW(I)*PREA*1.25_wp
    EW(I) = TV1+TV2+TV3+TV4
    !if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR,TV2,TV3,TV4
  end if
  E(I) = EW(I)+CON
end do
!-----------------------------------------------------------------------
! CALCULATE REVERSE TRANSFORMATION
!-----------------------------------------------------------------------
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

! CALCULATE KINEMATICAL FACTORS

do I=1,N
  AA(I) = sqrt((CON+E(I))/(Two*E(I)))
  RR(I) = sqrt(CON)/(CON+E(I))
end do

! POTENTIAL

! BEYOND THIS POINT, OVE IS USED AS SCRATCH ARRAY

! TRANSFORM V TO T-BASIS

call TrSmr(V,SINV,G,N,AUX,OVE)

call TrSmr(G,EIG,BU,N,AUX,OVE)

!ulf
if (idbg > 0) call PRMAT(IDBG,V,N,0,'v oper  ')

! MULTIPLY

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    P(IJ) = BU(IJ)
    BU(IJ) = BU(IJ)*AA(I)*AA(J)
  end do
end do

call TrSmtr(BU,REVT,G,Zero,N,AUX,OVE)

!ulf
if (idbg > 0) call PRMAT(IDBG,g,N,0,'g oper  ')

H(:) = H+G

! PVP INTEGRALS AND TRANSFORM THEM TO T-BASIS

if (idbg > 0) call PRMAT(IDBG,pvp,N,0,'raw pvp integrals  ')
call TrSmr(PVP,SINV,G,N,AUX,OVE)

call TrSmr(G,EIG,BU,N,AUX,OVE)

! MULTIPLY

IJ = 0
do I=1,N
  do J=1,I
    IJ = IJ+1
    G(IJ) = BU(IJ)
    BU(IJ) = BU(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
  end do
end do

call TrSmtr(BU,REVT,EV2,Zero,N,AUX,OVE)
!ulf
if (idbg > 0) call PRMAT(IDBG,ev2,n,0,'pvp oper')

H(:) = H+EV2

! CALCULATE Even2r OPERATOR

!if (idbg > 0) call Even2r(idbg,N,P,G,E,AA,RR,TT,EIG,AUX,OVE)

! TRANSFORM BACK

!ulf
if (idbg > 0) call PRMAT(IDBG,g,n,0,'ev2 orig')
!call TrSmtr(G,REVT,EV2,Zero,N,AUX,OVE)
!ulf
!if (idbg > 0) call PRMAT(IDBG,ev2,n,0,'ev2 oper')
!H(:) = H+EV2
!ulf
if (idbg > 0) call PRMAT(IDBG,h,n,0,'h   oper')
!call Sogr(iDbg,N,S,SINV,P,OVE,EW)
!call Diagr(H,N,EIG,EW,SINV,AUX,0)
!write(u6,*) 'END OF SCFCLI2 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

!if (idbg > 0) write(idbg,*) '--- EIGENVALUES OF H MATRIX ---'
!if (idbg > 0) write(idbg,'(4D20.12)') EW

return

!100 format(1X,I4,7(2X,D14.8))
!556 format(//,7X,'- NREL. ENERG.  -  DIVIDED BY C - REL.  ENERG.  -  MOMENTUM    - TERMS OF POWER SERIES (LOW ENERGY ONLY)'//)
!2012 format('  relsewb|****** '/,'        |****** WARNING - OVERLAP MATRIX SINGULAR '/, &
!            '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,'        |******'//)
!2016 format(' relsewb| DCOPIV rc=',I2,', |S|=',D20.6,'x 10**(',I4,') ')

end subroutine SCFCLI2
