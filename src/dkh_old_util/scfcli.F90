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
! Copyright (C) 1986,1995, Bernd Artur Hess                            *
!               2005, Jesper Wisborg Krogh                             *
!***********************************************************************

subroutine SCFCLI(idbg,eps,S,H,V,PVP,N,ISIZE,VELIT,TMP1,TMP2,TMP3,TMPA,TMPB,TMPC,EW,E,AA,RR,TT,TMP4,TMPD,TMPE,TMPF,TMP5,I_DIM)
! $Id: relsew.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de
!   Modified: 2005 Jesper Wisborg Krogh, Jesper.Krogh@teokem.lu.se

use DKH_Info, only: IRELAE
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: idbg, N, ISIZE, I_DIM
real(kind=wp), intent(out) :: eps, TMP1(I_DIM*(I_DIM+1)/2), TMP2(ISIZE), TMP3(ISIZE), TMPA(N,N), TMPB(N,N), TMPC(N,N), EW(N), &
                              E(N), AA(N), RR(N), TT(N), TMP4(I_DIM*(I_DIM+1)/2), TMPD(I_DIM,I_DIM), TMPE(N,N), TMPF(N,N), &
                              TMP5(ISIZE)
real(kind=wp), intent(in) :: S(ISIZE), V(ISIZE), VELIT
real(kind=wp), intent(inout) :: H(ISIZE), PVP(ISIZE)
integer(kind=iwp) :: I, IJ, ILL, J, K, M
real(kind=wp) :: CON, CON2, PREA, RATIO, SR, TV1, TV2, TV3, TV4
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: icontr, iex
real(kind=wp) :: det, dtol
#endif

!write(u6,*) ' in SCFCLI', N, iSize
PREA = One/(VELIT*VELIT)
CON2 = PREA+PREA
CON = One/PREA

#ifdef _DEBUGPRINT_
! CALCULATE DETERMINANT

call Square(S,TMPB,n,1,n)
!do i=1,n
!  write(u6,'(5f10.5)') (TMPB(i,j),j=1,n)
!end do
icontr = -1
dtol = 1.0e-14_wp
call dcopiv(TMPB,TMPB,n,1,n,dtol,det,iex,icontr,TMP2)
if (idbg > 0) write(idbg,2016) icontr,det,iex
if (icontr /= 0) then
  write(u6,2016) icontr,det,iex
  write(u6,2012) dtol
  write(u6,'(a)') ' relsew| singular overlap matrix'
  call Abend()
end if
#endif

! SCHMIDT-ORTHOGONALIZE

call Sogr(idbg,N,S,TMPE,TMP2,TMPC,EW)

! ** TMPE CONTAINS TRANSFORMATION TO ORTHOGONAL AO-BASIS
!-----------------------------------------------------------------------
! NON-RELATIVISTIC KINETIC ENERGY
!-----------------------------------------------------------------------
call Diagr(H,N,TMPF,EW,TMPE,TMPB,TMP2)
if (idbg > 0) write(idbg,556)
do I=1,N
  if (ew(i) < Zero) then
    write(u6,*) ' scfcli| ew(',i,') = ',ew(i)
    write(u6,'(a)') 'kinetic energy eigenvalue less than zero'
    call Abend()
  end if

  ! IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO

  RATIO = EW(I)/VELIT

  ! CALCULATE RELATIVISTIC ENERGY AND MOMENTUM

  SR = sqrt(Two*EW(I))
  TT(I) = EW(I)
  if (RATIO > 0.02_wp) then
    TV1 = EW(I)
    EW(I) = CON*(sqrt(One+CON2*EW(I))-One)
    if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR
  else
    TV1 = EW(I)
    TV2 = -TV1*EW(I)*PREA*Half
    TV3 = -TV2*EW(I)*PREA
    TV4 = -TV3*EW(I)*PREA*1.25_wp
    EW(I) = TV1+TV2+TV3+TV4
    if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR,TV2,TV3,TV4
  end if
  E(I) = EW(I)+CON
end do
!-----------------------------------------------------------------------
! CALCULATE REVERSE TRANSFORMATION
!-----------------------------------------------------------------------
call DGEMM_('N','N',n,n,n,One,TMPE,n,TMPF,n,Zero,TMPB,n)
call Square(S,TMPC,N,1,N)
call DGEMM_('N','N',n,n,n,One,TMPC,n,TMPB,n,Zero,TMPA,n)
! ** TMPC  CONTAINS OVERLAP MATRIX IN FULL

if ((IRELAE /= 21) .and. (IRELAE /= 22) .and. (IRELAE /= 23)) then

  H(:) = Zero
  do K=1,N
    IJ = 0
    do I=1,N
      do J=1,I
        IJ = IJ+1
        H(IJ) = H(IJ)+TMPA(I,K)*TMPA(J,K)*EW(K)
      end do
    end do
  end do

else
  H(:) = Zero
end if

! CALCULATE KINEMATICAL FACTORS

if (IRELAE /= 11) then

  do I=1,N
    AA(I) = sqrt((CON+E(I))/(Two*E(I)))
    RR(I) = sqrt(CON)/(CON+E(I))
  end do

else if (IRELAE == 11) then

  do I=1,N
    AA(I) = (sqrt(One+CON*TT(I)*Two/((CON+E(I))*(CON+E(I)))))/(CON+E(I)) ! O OPERATOR
    RR(I) = sqrt(CON)/(CON+E(I)) ! Q OPERATOR
  end do

end if

! POTENTIAL

! BEYOND THIS POINT, TMPC IS USED AS SCRATCH ARRAY

! TRANSFORM V TO T-BASIS

call TrSmr2(V,TMPE,TMP3,N,TMPB,TMPF,TMPC)
!ulf
if (idbg > 0) call PRMAT(IDBG,V,N,0,'v oper  ')

! MULTIPLY

if ((IRELAE == 0) .or. (IRELAE == 1) .or. (IRELAE == 2) .or. (IRELAE == 3)) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMP2(IJ) = TMP3(IJ)
      TMP3(IJ) = TMP3(IJ)*AA(I)*AA(J)
    end do
  end do
  if (IRELAE == 3) then
    call Square(TMP3,TMPD,N,1,N)
  end if

else if (IRELAE == 11) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMP2(IJ) = TMP3(IJ)
      TMP3(IJ) = VELIT*TMP3(IJ)*(sqrt(RR(I)*RR(J))*AA(I)/AA(J)+sqrt(RR(J)*RR(I))*AA(J)/AA(I))
    end do
  end do

else if (IRELAE == 22) then  ! ZORA(FP)

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMP3(IJ) = TMP3(IJ)*AA(I)*AA(J)
    end do
  end do

end if

call TrSmtr(TMP3,TMPA,H,One,N,TMPB,TMPC)
!ulf

! PVP INTEGRALS AND TRANSFORM THEM TO T-BASIS

if (idbg > 0) call PRMAT(IDBG,pvp,N,0,'raw pvp integrals  ')
TMP5(:) = PVP(:)
call TrSmr2(TMP5,TMPE,PVP,N,TMPB,TMPF,TMPC)

! MULTIPLY

if ((IRELAE == 0) .or. (IRELAE == 1) .or. (IRELAE == 2) .or. (IRELAE == 3)) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMP3(IJ) = PVP(IJ)
      PVP(IJ) = PVP(IJ)*AA(I)*RR(I)*AA(J)*RR(J)
      if (IRELAE == 3) then
        TMPD(I,J) = TMPD(I,J)+PVP(IJ)
        TMPD(J,I) = TMPD(I,J)
      end if
    end do
  end do

else if (IRELAE == 11) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMP3(IJ) = PVP(IJ)
      PVP(IJ) = PVP(IJ)*(RR(I)*RR(J)*AA(I)/AA(J)+RR(J)*RR(I)*AA(J)/AA(I))*Half
    end do
  end do

else if ((IRELAE == 21) .or. (IRELAE == 22) .or. (IRELAE == 23)) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      PVP(IJ) = -PVP(IJ)*Half*PREA
      if (I == J) PVP(IJ) = PVP(IJ)+TT(I)*Two
    end do
  end do
  call Square(PVP,TMPD,N,1,N)

  ! inverse operator

  EPS = 1.0e-13_wp
  call MINVD(TMPD,N,N,EPS,ILL)
  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      PVP(IJ) = TMPD(I,J)
      PVP(IJ) = PVP(IJ)*TT(I)*TT(J)*Two
      if (IRELAE == 22) PVP(IJ) = PVP(IJ)*AA(I)*AA(J)
    end do
  end do

end if

call TrSmtr(PVP,TMPA,H,One,N,TMPB,TMPC)
!ulf

if (IRELAE == 23) then
  IJ = 0
  call Square(PVP,TMPE,N,1,N)
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMPD(I,J) = PVP(IJ)/TT(J)*Half*PREA
      TMPD(J,I) = PVP(IJ)/TT(I)*Half*PREA
    end do
  end do
  M = N
  call unitmat(TMPC,N)
  call DGEMM_('N','N',N,N,N,One,TMPD,M,TMPE,M,One,TMPC,M)
  ! modified overlap is incorporated into PVP
  call DGEMM_('N','N',N,N,N,One,TMPA,N,TMPC,N,Zero,TMPB,N)
  call dGemm_tri('N','T',N,N,N,One,TMPB,N,TMPA,N,Zero,PVP,N)
end if

if ((IRELAE /= 1) .and. (IRELAE /= 11) .and. (IRELAE /= 21) .and. (IRELAE /= 22) .and. (IRELAE /= 23)) then

  if (IRELAE == 3) then
    ! KEEP T-BASIS VEXT INTO TMP1 FOR HIGHER-ORDER DK
    TMP1(1:N*(N+1)/2) = TMP2(1:N*(N+1)/2)
    ! KEEP T-BASIS PVP INTO TMP4 FOR HIGHER-ORDER DK
    TMP4(1:N*(N+1)/2) = TMP3(1:N*(N+1)/2)
  end if

  ! CALCULATE Even2r OPERATOR

  call Even2r(idbg,N,TMP2,TMP3,E,AA,RR,TT,TMPE,TMPB,TMPC,TMPF)

  ! TRANSFORM BACK

  !ulf
  if (idbg > 0) call PRMAT(IDBG,TMP3,n,0,'ev2 orig')
  call TrSmtr(TMP3,TMPA,H,One,N,TMPB,TMPC)
  !ulf

  if ((IRELAE /= 0) .and. (IRELAE /= 2)) then  ! DK2

    ! CALCULATE Even3r OPERATOR

    call Even3r(N,TMP2,TMP3,E,AA,RR,TT,TMPB,TMPD,TMP1,TMP4,TMPE,TMPF)

    ! TRANSFORM BACK

    !ulf
    if (idbg > 0) call PRMAT(IDBG,TMP3,n,0,'ev2 orig')
    call TrSmtr(TMP3,TMPA,H,One,N,TMPB,TMPC)
    !ulf

    if (IRELAE /= 3) then  ! DK3

      ! More to come!

    end if
  end if
end if

!ulf
if (idbg > 0) call PRMAT(IDBG,h,n,0,'h   oper')
call Sogr(idbg,N,S,TMPE,TMP2,TMPC,EW)
call Diagr(H,N,TMPF,EW,TMPE,TMPB,TMP2)
if (idbg > 0) then
  call PRMAT(IDBG,h,n,0,'h   oper(final)')
  write(idbg,*) '--- EIGENVALUES OF H MATRIX ---'
  write(idbg,'(4ES20.12)') EW
end if

return

100 format(1X,I4,7(1X,ES15.8))
556 format(//,7X,'- NREL. ENERG.  -  DIVIDED BY C - REL.  ENERG.  -  MOMENTUM    - TERMS OF POWER SERIES (LOW ENERGY ONLY)'//)
#ifdef _DEBUGPRINT_
2012 format('  relsew|****** '/,'        |****** WARNING - OVERLAP MATRIX SINGULAR '/, &
            '        |****** PIVOTAL ELEMENT LESS THAN ',ES20.4,' FOUND'/,'        |******'//)
2016 format('  relsew| DCOPIV rc=',I2,', |S|=',ES20.6,'x 10**(',I4,') ')
#endif

end subroutine SCFCLI
