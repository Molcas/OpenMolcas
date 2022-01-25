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

subroutine SCFCLI(idbg,epsilon,S,H,V,PVP,N,ISIZE,VELIT,TMP1,TMP2,TMP3,TMPA,TMPB,TMPC,EW,E,AA,RR,TT,TMP4,TMPD,TMPE,TMPF,TWRK4,IDIM)
! $Id: relsew.r,v 1.4 1995/05/08 14:08:53 hess Exp $
! calculate relativistic operators
!   Bernd Artur Hess, hess@uni-bonn.de
!   Modified: 2005 Jesper Wisborg Krogh, Jesper.Krogh@teokem.lu.se

use DKH_Info, only: IRELAE

implicit real*8(A-H,O-Z)
dimension S(ISIZE), H(ISIZE), V(ISIZE), PVP(ISIZE)
dimension TMP1(IDIM*(IDIM+1)/2), TMP2(ISIZE), TMP3(ISIZE), TMPA(N,N), TMPB(N,N), TMPC(N,N), EW(N), E(N), AA(N), RR(N), TT(N)
dimension TMP4(IDIM*(IDIM+1)/2)
dimension TMPD(IDIM,IDIM)
dimension TMPE(N,N)
dimension TMPF(N,N)
dimension TWRK4(N*200)

!write(6,*) ' in SCFCLI', N, iSize
PREA = 1d0/(VELIT*VELIT)
CON2 = PREA+PREA
CON = 1.d0/PREA
#ifdef _DEBUGPRINT_

! CALCULATE DETERMINANT

call Square(S,TMPB,n,1,n)
!do i=1,n
!  write(6,'(5f10.5)') (TMPB(i,j),j=1,n)
!end do
icontr = -1
dtol = 1.D-14
call dcopiv(TMPB,TMPB,n,1,n,dtol,det,iex,icontr,TMP2)
if (idbg > 0) write(idbg,2016) icontr,det,iex
2016 format('  relsew| DCOPIV rc=',I2,', |S|=',D20.6,'x 10**(',I4,') ')
if (icontr /= 0) then
  write(6,2016) icontr,det,iex
  write(6,2012) dtol
2012 format('  relsew|****** '/,'        |****** WARNING - OVERLAP MATRIX SINGULAR '/, &
            '        |****** PIVOTAL ELEMENT LESS THAN ',D20.4,' FOUND'/,'        |******'//)
  call errex_rel(' relsew| singular overlap matrix')
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
556 format(//,7X,'- NREL. ENERG.  -  DIVIDED BY C - REL.  ENERG.  -  MOMENTUM    - TERMS OF POWER SERIES (LOW ENERGY ONLY)'//)
do I=1,N
  if (ew(i) < 0.d0) then
    write(6,*) ' scfcli| ew(',i,') = ',ew(i)
    call errex_rel('kinetic energy eigenvalue less than zero')
  end if

  ! IF T SUFFICIENTLY SMALL, USE SERIES EXPANSION TO AVOID CANCELLATIO

  RATIO = EW(I)/VELIT

  ! CALCULATE RELATIVISTIC ENERGY AND MOMENTUM

  SR = sqrt(2.d0*EW(I))
  TT(I) = EW(I)
  if (RATIO > 0.02d0) goto 11
  TV1 = EW(I)
  TV2 = -TV1*EW(I)*PREA/2.d0
  TV3 = -TV2*EW(I)*PREA
  TV4 = -TV3*EW(I)*PREA*1.25d0
  EW(I) = TV1+TV2+TV3+TV4
  if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR,TV2,TV3,TV4
100 format(1X,I4,7(1X,D15.8))
  goto 12
11 TV1 = EW(I)
  EW(I) = CON*(sqrt(1.d0+CON2*EW(I))-1.d0)
  if (idbg > 0) write(idbg,100) I,TV1,RATIO,EW(I),SR
12 continue
  E(I) = EW(I)+CON
end do
!-----------------------------------------------------------------------
! CALCULATE REVERSE TRANSFORMATION
!-----------------------------------------------------------------------
call DGEMM_('N','N',n,n,n,1.0d0,TMPE,n,TMPF,n,0.0d0,TMPB,n)
#ifdef MOLPRO
call Square(TMPC,S,N,N)
#else
call Square(S,TMPC,N,1,N)
#endif
call DGEMM_('N','N',n,n,n,1.0d0,TMPC,n,TMPB,n,0.0d0,TMPA,n)
! ** TMPC  CONTAINS OVERLAP MATRIX IN FULL

if ((IRELAE /= 21) .and. (IRELAE /= 22) .and. (IRELAE /= 23)) then

  call dCopy_(iSize,[0.0d0],0,H,1)
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
  call dCopy_(iSize,[0.0d0],0,H,1)
end if

! CALCULATE KINEMATICAL FACTORS

if (IRELAE /= 11) then

  do I=1,N
    AA(I) = sqrt((CON+E(I))/(2.d0*E(I)))
    RR(I) = sqrt(CON)/(CON+E(I))
  end do

else if (IRELAE == 11) then

  do I=1,N
    AA(I) = (sqrt(1.0d0+CON*TT(I)*2.0d0/((CON+E(I))*(CON+E(I)))))/(CON+E(I)) ! O OPERATOR
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
#   ifdef MOLPRO
    call Square(TMPD,TMP3,N,N)
#   else
    call Square(TMP3,TMPD,N,1,N)
#   endif
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

call TrSmtr(TMP3,TMPA,H,1.0d0,N,TMPB,TMPC)
!ulf

! PVP INTEGRALS AND TRANSFORM THEM TO T-BASIS

if (idbg > 0) call PRMAT(IDBG,pvp,N,0,'raw pvp integrals  ')
call TrSmr2(PVP,TMPE,PVP,N,TMPB,TMPF,TMPC)

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
      PVP(IJ) = PVP(IJ)*(RR(I)*RR(J)*AA(I)/AA(J)+RR(J)*RR(I)*AA(J)/AA(I))*0.5d0
    end do
  end do

else if ((IRELAE == 21) .or. (IRELAE == 22) .or. (IRELAE == 23)) then

  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      PVP(IJ) = -PVP(IJ)*0.5d0*PREA
      if (I == J) PVP(IJ) = PVP(IJ)+TT(I)*2.0d0
    end do
  end do
# ifdef MOLPRO
  call Square(TMPD,PVP,N,N)
# else
  call Square(PVP,TMPD,N,1,N)
# endif

  ! inverse operator

  EPS = 1.d-13
  call MINVD(TMPD,N,N,EPS,ILL)
  IJ = 0
  do I=1,N
    do J=1,I
      IJ = IJ+1
      PVP(IJ) = TMPD(I,J)
      PVP(IJ) = PVP(IJ)*TT(I)*TT(J)*2.0d0
      if (IRELAE == 22) PVP(IJ) = PVP(IJ)*AA(I)*AA(J)
    end do
  end do

end if

call TrSmtr(PVP,TMPA,H,1.0d0,N,TMPB,TMPC)
!ulf

if (IRELAE == 23) then
  IJ = 0
# ifdef MOLPRO
  call Square(TMPE,PVP,N,N)
# else
  call Square(PVP,TMPE,N,1,N)
# endif
  do I=1,N
    do J=1,I
      IJ = IJ+1
      TMPD(I,J) = PVP(IJ)/TT(J)*0.5d0*PREA
      TMPD(J,I) = PVP(IJ)/TT(I)*0.5d0*PREA
    end do
  end do
  M = N
  call dCopy_(N*N,[0.0d0],0,TMPC,1)
  call dCopy_(N,[1.0d0],0,TMPC,N+1)
  call DGEMM_('N','N',N,N,N,1.0d0,TMPD,M,TMPE,M,1.0d0,TMPC,M)
  ! modified overlap is incorporated into PVP
  call DGEMM_('N','N',N,N,N,1.0d0,TMPA,N,TMPC,N,0.0d0,TMPB,N)
  call dGemm_tri('N','T',N,N,N,1.0d0,TMPB,N,TMPA,N,0.0d0,PVP,N)
end if
if ((IRELAE == 1) .or. (IRELAE == 11) .or. (IRELAE == 21) .or. (IRELAE == 22) .or. (IRELAE == 23)) goto 1000

if (IRELAE == 3) then
  ! KEEP T-BASIS VEXT INTO TMP1 FOR HIGHER-ORDER DK
  call dCopy_(N*(N+1)/2,TMP2,1,TMP1,1)
  ! KEEP T-BASIS PVP INTO TMP4 FOR HIGHER-ORDER DK
  call dCopy_(N*(N+1)/2,TMP3,1,TMP4,1)
end if

! CALCULATE Even2r OPERATOR

call Even2r(idbg,N,TMP2,TMP3,E,AA,RR,TT,TMPE,TMPB,TMPC,TMPF)

! TRANSFORM BACK

!ulf
if (idbg > 0) call PRMAT(IDBG,TMP3,n,0,'ev2 orig')
call TrSmtr(TMP3,TMPA,H,1.0d0,N,TMPB,TMPC)
!ulf
if ((IRELAE == 0) .or. (IRELAE == 2)) goto 1000   ! DK2

! CALCULATE Even3r OPERATOR

call Even3r(idbg,N,TMP2,TMP3,E,AA,RR,TT,TMPB,TMPD,TMP1,TMP4,TMPE,TMPF)

! TRANSFORM BACK

!ulf
if (idbg > 0) call PRMAT(IDBG,TMP3,n,0,'ev2 orig')
call TrSmtr(TMP3,TMPA,H,1.0d0,N,TMPB,TMPC)
!ulf

if (IRELAE == 3) goto 1000   ! DK3

! More to come!

1000 continue

!ulf
if (idbg > 0) call PRMAT(IDBG,h,n,0,'h   oper')
call Sogr(idbg,N,S,TMPE,TMP2,TMPC,EW)
call Diagr(H,N,TMPF,EW,TMPE,TMPB,TMP2)
if (idbg > 0) call PRMAT(IDBG,h,n,0,'h   oper(final)')
if (idbg > 0) write(idbg,*) '--- EIGENVALUES OF H MATRIX ---'
if (idbg > 0) write(idbg,'(4D20.12)') EW

return
! Avoid unused argument warnings
if (.false.) then
  call Unused_real(epsilon)
  call Unused_real_array(TWRK4)
end if

end subroutine SCFCLI
