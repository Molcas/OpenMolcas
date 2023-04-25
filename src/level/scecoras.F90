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

integer I, II, I1, I2, I3, I4, IV1, IV2, INNER, ICOR, JROT, KV, KVB, KVLEV, KVDIF, NDP, NCN, IDIF, BRUTE, IB, IWR
real*8 EO, DE0, RH, BFCT, ARG2, ARG3, EINT, VPH1, VPH2, DGDV1, DGDV2, DGDVM, DGDV2P, DGDVB, DGDVBP, EBRUTE, DEBRUTE, DE1, DE2, Y1, &
       Y2, Y3, RT, ANS1, ANS2, XDIF, VLIM, BMAX, Pi, Pi2, PNCN, PP1, V(NDP), SDRDY(NDP)
save BRUTE, EBRUTE, DEBRUTE, DGDVB, Pi, Pi2
data DGDVB/-1.d0/,KVB/-1/,Pi/3.1415926454d0/,Pi2/6.283185308d0/

ARG3 = 0
I1 = 0
II = 0
DGDVBP = -1.d0
DGDV2 = -1.d0
EINT = EO*BFCT
if (KVLEV == 0) DGDVB = -1.d0
KVDIF = KVLEV-KV
if (ICOR == 1) BRUTE = 0
I3 = NDP
!PNCN = dfloat(NCN-2)/dfloat(NCN+2)
PNCN = dble(NCN-2)/dble(NCN+2)
PP1 = 1.d0/pNCN+1.d0
! For Quasibound levels, first search inward to classically forbidden
if (EO > VLIM) then
  PNCN = 1.d0
  PP1 = 2.d0
  do I=NDP,1,-1
    I3 = I
    if (V(I) > EINT) goto 8
  end do
end if
! First, search inward for outermost turning point
8 continue
do I=I3,1,-1
  I4 = I
  if (V(I) < EINT) goto 10
end do
! If never found an 'outer' turning point (e.g., above qbdd. barier)
! then simply return with negative  DGDV2  as error flag
return
!... Now collect vibrational phase and its energy deriv. over outer well
10 continue
Y1 = EINT-V(I4+1)
Y2 = EINT-V(I4)
Y3 = EINT-V(I4-1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
ARG2 = dsqrt(Y3)
VPH2 = 0.5d0*ARG2+ANS2*SDRDY(I4)**2/RH
DGDV2 = 0.5d0/ARG2+ANS1*SDRDY(I4)**2/RH
do I=I4-2,1,-1
  !... now, collect (v+1/2) and dv/dG integrals to next turning point ...
  II = I
  if (V(I) > EINT) go to 12
  ARG3 = ARG2
  ARG2 = dsqrt(EINT-V(I))
  VPH2 = VPH2+ARG2*SDRDY(I)**2
  DGDV2 = DGDV2+SDRDY(I)**2/ARG2
end do
12 continue
I3 = II+1
Y1 = EINT-V(I3-1)
Y2 = EINT-V(I3)
Y3 = EINT-V(I3+1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
VPH2 = (VPH2-ARG2-0.5d0*ARG3+ANS2*SDRDY(I3)**2/RH)/Pi
DGDV2 = DGDV2-1.d0/ARG2-0.5d0/ARG3+ANS1*SDRDY(I3)**2/RH
DGDV2 = Pi2/(BFCT*DGDV2)
! Next, search for innermost turning point
do I=1,NDP
  I1 = I
  if (V(I) < EINT) goto 20
end do
!... then collect vibrational phase and its energy deriv. over outer well
20 continue
if (I1 == 1) then
  write(6,602) JROT,EO
  !stop
  call ABEND()
end if
if (I1 >= I3) then
  ! For single-well potential or above barrier of double-well potential
  if (IWR >= 2) write(6,600) ICOR,KV,JROT,EO,VPH2-0.5d0,DGDV2
  if ((KV /= KVLEV-1) .and. (DGDVB > 0.d0)) then
    !... If got wrong level (KV not one below KVLEV) and NOT first call ...
    if ((EO-BMAX) > (2.d0*DGDV2)) then
      !... 'Normal' case: use B-S plot area to estimate correct energy
      !DE0 = KVDIF*(DGDV2-0.5d0*(DGDV2-DGDVB)/dfloat(KV-KVB))
      DE0 = KVDIF*(DGDV2-0.5d0*(DGDV2-DGDVB)/dble(KV-KVB))
      EO = EO+DE0
      KV = KVB
      KVLEV = KV+1
      return
    else
      !... but close to barrier in double-well potential, switch to 'BRUTE'
      BRUTE = BRUTE+1
      DGDV1 = DGDV2
      XDIF = sign(1,KVDIF)
      goto 54
    end if
  end if
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
    DE0 = 0.5d0*(3.d0*DGDV2-DGDVB)
    if ((2.d0*DGDV2) > DGDVB) then
      !... if linear Birge-Sponer predicts another level, then use it
      EO = EO+DE0
    else
      !... otherwise, use N-D theory extrapolation for next level...
      DGDV2P = DGDV2**PNCN
      DE0 = (DGDV2P+DGDV2P-DGDVBP)
      if (DE0 > 0.d0) then
        DE0 = (DE0**PP1-DGDV2P**PP1)/(PP1*(DGDV2P-DGDVBP))
        EO = EO+DE0
      else
        !... but if NDT predicts no more levels, quit, and (optionally) print
        if (IWR > 0) write(6,604) KV,EO
        604 format(10x,'Find highest bound level is   E(v=',i3,')=',1PD18.10)
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

! For a double-well potential, collect vibrational phase and its
! energy derivative over the inner well
Y1 = EINT-V(I1-1)
Y2 = EINT-V(I1)
Y3 = EINT-V(I1+1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
ARG2 = dsqrt(Y3)
VPH1 = 0.5d0*ARG2+ANS2*SDRDY(I1)**2/RH
DGDV1 = 0.5d0/ARG2+ANS1*SDRDY(I1)**2/RH
do I=I1+2,NDP
  !... now, collect integral and count nodes outward to next turning point ...
  if (V(I) > EINT) go to 22
  ARG3 = ARG2
  ARG2 = dsqrt(EINT-V(I))
  VPH1 = VPH1+ARG2*SDRDY(I)**2
  DGDV1 = DGDV1+SDRDY(I)**2/ARG2
end do
22 continue
I2 = I-1
Y1 = EINT-V(I2+1)
Y2 = EINT-V(I2)
Y3 = EINT-V(I2-1)
call LEVQAD(Y1,Y2,Y3,RH,RT,ANS1,ANS2)
VPH1 = (VPH1-ARG2-0.5d0*ARG3+ANS2*SDRDY(I2)**2/RH)/Pi
DGDV1 = DGDV1-1.d0/ARG2-0.5d0/ARG3+ANS1*SDRDY(I2)**2/RH
DGDV1 = Pi2/(BFCT*DGDV1)
DGDVM = DGDV1*DGDV2/(DGDV1+DGDV2)
write(6,*) DGDVM ! make it "referenced"
if (KVDIF == 0) then
  ! If already at level sought, return
  if (IWR >= 2) write(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,ICOR,VPH2-0.5d0,DGDV2
  return
end if

! Check whether looking for higher or lower level ...
IDIF = sign(1,KVDIF)
XDIF = IDIF
if ((ICOR >= 6) .and. ((iabs(KVDIF) == 1) .or. (BRUTE > 0))) goto 50
! 'Conventional' semiclassical search for neared INNER or OUTER well level
if (INNER <= 0) then
  !... and current energy EO is for an outer-well level ...
  DE2 = DGDV2*XDIF
  IV1 = int(VPH1+0.5d0)
  !DE1 = (dfloat(IV1)+0.5d0-VPH1)*DGDV1*XDIF
  DE1 = (dble(IV1)+0.5d0-VPH1)*DGDV1*XDIF
  if (IWR >= 2) write(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,ICOR,VPH2-0.5d0,DGDV2
  30 continue
  if (dabs(DE1) < dabs(DE2)) then
    INNER = 1
    EO = EO+DE1
    DE1 = DGDV1*XDIF
  else
    INNER = 0
    EO = EO+DE2
  end if
  KVDIF = KVDIF-IDIF
  if (KVDIF == 0) then
    return
  end if
  goto 30
end if
if (INNER > 0) then
  !... and current energy EO is for an inner-well level ...
  DE1 = DGDV1*XDIF
  IV2 = int(VPH2+0.5d0)
  !DE2 = (dfloat(IV2)+0.5d0-VPH2)*DGDV2*XDIF
  DE2 = (dble(IV2)+0.5d0-VPH2)*DGDV2*XDIF
  if (IWR >= 2) write(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,ICOR,VPH2-0.5d0,DGDV2
  40 continue
  if (dabs(DE2) < dabs(DE1)) then
    INNER = 0
    EO = EO+DE2
    DE2 = DGDV2*XDIF
  else
    INNER = 1
    EO = EO+DE1
  end if
  KVDIF = KVDIF-IDIF
  if (KVDIF == 0) then
    return
  end if
  goto 40
end if
50 continue
BRUTE = BRUTE+1
! Now .. Brute force search for desired level !
if (IWR >= 2) write(6,610) KV,JROT,EO,VPH1-0.5d0,DGDV1,KVLEV,ICOR,VPH2-0.5d0,DGDV2
54 continue
if (BRUTE == 1) then
  !... in first brute-force step, use previous energy with opposite INNER
  EBRUTE = EO
  if (INNER == 0) then
    INNER = 1
  else
    INNER = 0
  end if
  DEBRUTE = dmin1(DGDV1,DGDV2)*XDIF*0.3d0
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
610 format(' Double well   E(v=',i3,', J=',I3,')=',f9.3,':   v1(SC)=',F7.3,'   dGdv1=',f8.2/8x,'seeking  v=',I3,' (ICOR=',I2,')', &
           8x,':   v2(SC)=',F7.3,'   dGdv2=',f8.2)

end subroutine SCECORas
