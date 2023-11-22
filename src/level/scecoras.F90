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
subroutine SCECORas(KV,KVLEV,JROT,INNER,ICOR,IWR,EO,RH,BFCT,NDP,NCN,V,SDRDY,BMAX,VLIM,DGDV2)
!** Subroutine calculates (approximate!) semiclassical estimate of
!  dG/dv for level  v= KV  with energy  EO [cm-1]  on potential
!  {V(i),i=1,NDP} (in 'internal BCFT units' {V[cm-1]*BFCT}), and uses
!  those results to estimate energy of level  KVLEV
!** If the 'clever' semiclassical procedure fails - try a brute force
!  step-by-step search, using alternately INNER & OUTER well starting
!** BMAX is internal barrier maximum energy for double-well case,
!   and very large negative number for single-well potential
!** On return, negative DGDV2 signals error!  No phase integrals found

use Constants, only: Zero, One, Two, Three, Half, Pi
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: KV, KVLEV, INNER
integer(kind=iwp), intent(in) :: JROT, ICOR, IWR, NDP, NCN
real(kind=wp), intent(inout) :: EO
real(kind=wp), intent(in) :: RH, BFCT, V(NDP), SDRDY(NDP), BMAX, VLIM
real(kind=wp), intent(out) :: DGDV2
integer(kind=iwp) :: BRUTE, I, I1, I2, I3, I4, IB, IDIF, II, IV1, IV2, KVB = -1, KVDIF
real(kind=wp) :: ANS1, ANS2, ARG2, ARG3, DE0, DE1, DE2, DGDV1, DGDV2P, DGDVB = -One, DGDVBP, EINT, PNCN, PP1, RT, VPH1, VPH2, &
                 XDIF, Y1, Y2, Y3
logical(kind=iwp) :: Found, Skip
real(kind=wp), save :: DEBRUTE, EBRUTE

ARG3 = 0
I1 = 0
II = 0
DGDVBP = -One
DGDV2 = -One
EINT = EO*BFCT
if (KVLEV == 0) DGDVB = -One
KVDIF = KVLEV-KV
if (ICOR == 1) BRUTE = 0
I3 = NDP
PNCN = real(NCN-2,kind=wp)/real(NCN+2,kind=wp)
PP1 = One/pNCN+One
! For Quasibound levels, first search inward to classically forbidden
if (EO > VLIM) then
  PNCN = One
  PP1 = Two
  do I=NDP,1,-1
    I3 = I
    if (V(I) > EINT) exit
  end do
end if
! First, search inward for outermost turning point
Found = .false.
do I=I3,1,-1
  I4 = I
  if (V(I) < EINT) then
    Found = .true.
    exit
  end if
end do
! If never found an 'outer' turning point (e.g., above qbdd. barier)
! then simply return with negative  DGDV2  as error flag
if (.not. Found) return
!... Now collect vibrational phase and its energy deriv. over outer well
Y1 = EINT-V(I4+1)
Y2 = EINT-V(I4)
Y3 = EINT-V(I4-1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
ARG2 = sqrt(Y3)
VPH2 = Half*ARG2+ANS2*SDRDY(I4)**2/RH
DGDV2 = Half/ARG2+ANS1*SDRDY(I4)**2/RH
do I=I4-2,1,-1
  !... now, collect (v+1/2) and dv/dG integrals to next turning point ...
  II = I
  if (V(I) > EINT) exit
  ARG3 = ARG2
  ARG2 = sqrt(EINT-V(I))
  VPH2 = VPH2+ARG2*SDRDY(I)**2
  DGDV2 = DGDV2+SDRDY(I)**2/ARG2
end do
I3 = II+1
Y1 = EINT-V(I3-1)
Y2 = EINT-V(I3)
Y3 = EINT-V(I3+1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
VPH2 = (VPH2-ARG2-Half*ARG3+ANS2*SDRDY(I3)**2/RH)/Pi
DGDV2 = DGDV2-One/ARG2-Half/ARG3+ANS1*SDRDY(I3)**2/RH
DGDV2 = Two*Pi/(BFCT*DGDV2)
! Next, search for innermost turning point
do I=1,NDP
  I1 = I
  if (V(I) < EINT) exit
end do
!... then collect vibrational phase and its energy deriv. over outer well
if (I1 == 1) then
  write(u6,602) JROT,EO
  !stop
  call ABEND()
end if
Skip = .false.
if (I1 >= I3) then
  ! For single-well potential or above barrier of double-well potential
  if (IWR >= 2) write(u6,600) ICOR,KV,JROT,EO,VPH2-Half,DGDV2
  if ((KV /= KVLEV-1) .and. (DGDVB > Zero)) then
    !... If got wrong level (KV not one below KVLEV) and NOT first call ...
    if ((EO-BMAX) > (Two*DGDV2)) then
      !... 'Normal' case: use B-S plot area to estimate correct energy
      DE0 = KVDIF*(DGDV2-Half*(DGDV2-DGDVB)/real(KV-KVB,kind=wp))
      EO = EO+DE0
      KV = KVB
      KVLEV = KV+1
      return
    else
      !... but close to barrier in double-well potential, switch to 'BRUTE'
      BRUTE = BRUTE+1
      DGDV1 = DGDV2
      XDIF = sign(1,KVDIF)
      Skip = .true.
    end if
  end if
  if (.not. Skip) then
    if (KVLEV == 0) then
      ! If looking for v=0, just use local DVDG2 to estimate energy correction
      EO = EO+KVDIF*DGDV2
      return
    end if
    if (KV == 0) then
      ! Normally:  use B-S plot considerations to estimate next level energy
      !... use harmonic estimate for v=1 energy
      EO = EO+DGDV2
    else
      !... estimate Delta(G) based on linear Birge-Sponer
      DE0 = Half*(Three*DGDV2-DGDVB)
      if ((Two*DGDV2) > DGDVB) then
        !... if linear Birge-Sponer predicts another level, then use it
        EO = EO+DE0
      else
        !... otherwise, use N-D theory extrapolation for next level...
        DGDV2P = DGDV2**PNCN
        DE0 = (DGDV2P+DGDV2P-DGDVBP)
        if (DE0 > Zero) then
          DE0 = (DE0**PP1-DGDV2P**PP1)/(PP1*(DGDV2P-DGDVBP))
          EO = EO+DE0
        else
          !... but if NDT predicts no more levels, quit, and (optionally) print
          if (IWR > 0) write(u6,604) KV,EO
          return
        end if
      end if
    end if
    DGDVB = DGDV2
    DGDVBP = DGDVB**PNCN
    KVB = KV
    INNER = 0
    return
  end if
end if

! For a double-well potential, collect vibrational phase and its
! energy derivative over the inner well
Y1 = EINT-V(I1-1)
Y2 = EINT-V(I1)
Y3 = EINT-V(I1+1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
ARG2 = sqrt(Y3)
VPH1 = Half*ARG2+ANS2*SDRDY(I1)**2/RH
DGDV1 = Half/ARG2+ANS1*SDRDY(I1)**2/RH
do I=I1+2,NDP
  !... now, collect integral and count nodes outward to next turning point ...
  if (V(I) > EINT) exit
  ARG3 = ARG2
  ARG2 = sqrt(EINT-V(I))
  VPH1 = VPH1+ARG2*SDRDY(I)**2
  DGDV1 = DGDV1+SDRDY(I)**2/ARG2
end do
I2 = I-1
Y1 = EINT-V(I2+1)
Y2 = EINT-V(I2)
Y3 = EINT-V(I2-1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
VPH1 = (VPH1-ARG2-Half*ARG3+ANS2*SDRDY(I2)**2/RH)/Pi
DGDV1 = DGDV1-One/ARG2-Half/ARG3+ANS1*SDRDY(I2)**2/RH
DGDV1 = Two*Pi/(BFCT*DGDV1)
if (KVDIF == 0) then
  ! If already at level sought, return
  if (IWR >= 2) write(u6,610) KV,JROT,EO,VPH1-Half,DGDV1,KVLEV,ICOR,VPH2-Half,DGDV2
  return
end if

! Check whether looking for higher or lower level ...
IDIF = sign(1,KVDIF)
XDIF = IDIF
if ((ICOR < 6) .or. ((abs(KVDIF) /= 1) .and. (BRUTE <= 0))) then
  ! 'Conventional' semiclassical search for neared INNER or OUTER well level
  if (INNER <= 0) then
    !... and current energy EO is for an outer-well level ...
    DE2 = DGDV2*XDIF
    IV1 = int(VPH1+Half)
    DE1 = (real(IV1,kind=wp)+Half-VPH1)*DGDV1*XDIF
    if (IWR >= 2) write(u6,610) KV,JROT,EO,VPH1-Half,DGDV1,KVLEV,ICOR,VPH2-Half,DGDV2
    do
      if (abs(DE1) < abs(DE2)) then
        INNER = 1
        EO = EO+DE1
        DE1 = DGDV1*XDIF
      else
        INNER = 0
        EO = EO+DE2
      end if
      KVDIF = KVDIF-IDIF
      if (KVDIF == 0) return
    end do
  end if
  if (INNER > 0) then
    !... and current energy EO is for an inner-well level ...
    DE1 = DGDV1*XDIF
    IV2 = int(VPH2+Half)
    DE2 = (real(IV2,kind=wp)+Half-VPH2)*DGDV2*XDIF
    if (IWR >= 2) write(u6,610) KV,JROT,EO,VPH1-Half,DGDV1,KVLEV,ICOR,VPH2-Half,DGDV2
    do
      if (abs(DE2) < abs(DE1)) then
        INNER = 0
        EO = EO+DE2
        DE2 = DGDV2*XDIF
      else
        INNER = 1
        EO = EO+DE1
      end if
      KVDIF = KVDIF-IDIF
      if (KVDIF == 0) return
    end do
  end if
end if
if (.not. Skip) then
  BRUTE = BRUTE+1
  ! Now .. Brute force search for desired level !
  if (IWR >= 2) write(u6,610) KV,JROT,EO,VPH1-Half,DGDV1,KVLEV,ICOR,VPH2-Half,DGDV2
end if
if (BRUTE == 1) then
  !... in first brute-force step, use previous energy with opposite INNER
  EBRUTE = EO
  if (INNER == 0) then
    INNER = 1
  else
    INNER = 0
  end if
  DEBRUTE = min(DGDV1,DGDV2)*XDIF*0.3_wp
  return
end if
IB = BRUTE/2
!... in subsequent even steps, lower EO by DEBRUTE/10 for same INNER
if ((IB+IB) == BRUTE) then
  EBRUTE = EBRUTE+DEBRUTE
  EO = EBRUTE
  return
else
  !... in subsequent odd steps, lower repeat previous EO with INNER changed
  if (INNER == 0) then
    INNER = 1
  else
    INNER = 0
  end if
  EO = EBRUTE
  return
end if

!return

600 format(' Single well  ICOR=',I2,':  E(v=',i3,',J=',I3,')=',f10.2,'  v(SC)=',F8.3,'  dGdv=',f8.3)
602 format(/' *** ERROR ***  V(1) < E(J=',i3,')=',f10.2)
604 format(10x,'Find highest bound level is   E(v=',i3,')=',ES18.10)
610 format(' Double well   E(v=',i3,', J=',I3,')=',f9.3,':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,' (ICOR=',I2,')', &
           8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2)

end subroutine SCECORas
