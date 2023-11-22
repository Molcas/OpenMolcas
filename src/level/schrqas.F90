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
subroutine SCHRQas(KV,JROT,EO,GAMA,VMAX,VLIM,V,WF,BFCT,EEPS,YMIN,YH,NPP,NBEG,NEND,INNODE,INNER,IWR,LPRWF)
!*****  subroutine SCHRQ, last updated  21 August 2007 *****
!   !!!!   Form modified for handling Stolyarov radial variable  !!!!!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** SCHRQ solves radial Schrodinger equation in dimensionless form
!  d^2WF/dy^2 = - [(E-V(R))*BFCT*(r')^2 - F(y)]*WF(R) ,  where WF(I) is
!  the wave function and  y  the reduced radial vble.  y_p(r).
!** Integrate by Numerov method over NPP mesh points with increment
!  H=YH across range beginning at YMIN .
!** Input trial energy EO, eigenvalue convergence criterion EEPS
!  potential asymptote VLIM, and all returned energies (EO, GAMA & VMAX)
!  have units (cm-1).
!** On entry, the input potential V(I) must include the centrifugal
!  term, the factor:  'BFCT'=2*mu*(YH/hbar)**2 (1/cm-1) as well as the
!  Stolyarov conversion factors (r')^2 and F(y).
!  BFCT is also internally incorporated into EO, VLIM & EEPS.
!* Note that these reduced quantities (& the internal eigenvalue E)
!  contain a factor of the squared integration increment  YH**2 .
!  This saves arithmetic work in the innermost loop of the algorithm.
!** For energy in (cm-1), BFCT=ZMU(u)*H(Angst)**2/16.85762920 (1/cm-1)
!** INNODE > 0  specifies that wavefx. initiates at YMIN with a node
!     (normal default case);  INNODE <= 0  specifies  zero slope  at
!     YMIN (for finding symmetric eigenfunctions of symmetric potential
!     with potential mid-point @ YMIN).
!** INNER specifies wave function matching condition: INNER = 0  makes
!     matching of inward & outward solutions occur at outermost turning
!     point;  INNER > 0 makes matching occur at innermost turning point.
! * Normally use  INNER=0 ,  but to find inner-well levels of double
!     minimum potential, set  INNER > 0 .
!-----------------------------------------------------------------------
!** Output vibrational quantum number KV, eigenvalue EO, normalized
!  wave function WF(I), and range, NBEG <= I <= NEND  over
!  which WF(I) is defined. *** Have set  WF(I)=0  outside this range.
!* (NBEG,NEND), defined by requiring  abs(WF(I)) < RATST=1.0e-9  outside.
!** If(LPRWF /= 0) write every LPRWF-th value of wavefunction WF(I) to
!   a file on channel-10 (i.e., WRITE(10,XXX)), starting at YVB(NBEG)
!   with step size  |LPRWF|*YH.
!** For energies above the potential asymptote VLIM, locate quasibound
!  levels using Airy function boundary condition and return the level
!  width GAMA and barrier height VMAX, as well as EO.
!** ERROR condition on return is  KV < 0 ; usually KV=-1, but return
!  KV=-2 if error appears to arise from too low trial energy.
!** If(IWR /= 0) print error & warning descriptions
!  If (IWR > 0) also print final eigenvalues & node count.
!  If (IWR >= 2) also show end-of-range wave function amplitudes
!  If (IWR >= 3) print also intermediate trial eigenvalues, etc.
!** If input KV >= 998 , tries to find highest bound level, and
!  trial energy should be only slightly less than VLIM.
!-----------------------------------------------------------------------
!++ "SCHRQ" calls subroutines "QBOUND" and "WIDTH", and the latter
!++ calls "LEVQAD" .
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use LEVEL_COMMON, only: ARV, DRDY2, PRV, RVB, SDRDY, VBZ, YVB
use Constants, only: Zero, One, Two, Ten, Twelve, Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: KV
integer(kind=iwp), intent(in) :: JROT, NPP, INNODE, INNER, IWR, LPRWF
real(kind=wp), intent(inout) :: EO
real(kind=wp), intent(out) :: GAMA, VMAX, WF(NPP)
real(kind=wp), intent(in) :: VLIM, V(NPP), BFCT, EEPS, YMIN, YH
integer(kind=iwp), intent(out) :: NBEG, NEND
integer(kind=iwp) :: I, IBEGIN, ICOR, IT, ITER, ITP1, ITP1P, ITP2, ITP3, J, J1, J2, JPSIQ, KKV, KVIN, M, MS, MSAVE, NPR
real(kind=wp) :: DE, DEP, DEPRN, DF, DOLD, DSOC, DXPW, E, EPS, F, GB, GI, H, HT, PPROD, PROD, RATIN, RATOUT, REND, RR, SB, SI, SM, &
                 SN, SRTGB, SRTGI, VMX, VPR, WKBTST, XPR, Y1, Y2, Y3, YIN, YM, YMINN, YOUT
logical(kind=iwp) :: Conv, Found, Skip
integer(kind=iwp), parameter :: NDN = 10
real(kind=wp), parameter :: RATST = 1.0e-9_wp, XPW = log(1.0e10_wp)

!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING
!write(u6,*) 'After entering schrqas we have:'
!write(u6,*) 'KV=',KV
!write(u6,*) 'JROT=',JROT
!write(u6,*) 'EO=',EO
!write(u6,*) 'GAMA=',GAMA
!write(u6,*) 'VMAX=',VMAX
!write(u6,*) 'VLIM=',VLIM
!do I=1,3
!  write(u6,*) 'V=',V(I)
!  write(u6,*) 'WF=',WF(I)
!end do
!write(u6,*) 'BFCT=',BFCT
!write(u6,*) 'EEPS=',EEPS
!write(u6,*) 'YMIN=',YMIN
!write(u6,*) 'YH=',YH
!write(u6,*) 'NPP=',NPP
!write(u6,*) 'NBEG=',NBEG
!write(u6,*) 'NEND=',NEND
!write(u6,*) 'INNODE=',INNODE
!write(u6,*) 'INNER=',INNER
!write(u6,*) 'IWR=',IWR
!write(u6,*) 'LPRWF=',LPRWF
YOUT = 0
YM = 0
YIN = 0
ITP1P = 0
DXPW = XPW/real(NDN,kind=wp)
ICOR = 0
KVIN = KV
KV = -1
YMINN = YMIN-YH
GAMA = Zero
VMAX = VLIM
VMX = VMAX*BFCT
H = YH
HT = One/Twelve
E = EO*BFCT
EPS = EEPS*BFCT
DSOC = VLIM*BFCT
DE = Zero
RATIN = Zero
RATOUT = Zero
if (IWR > 2) then
  if (KVIN >= 998) then
    write(u6,610) EO
  else
    write(u6,601) KVIN,JROT,EO,INNER
  end if
  write(u6,602)
end if
NEND = NPP
GI = Zero
GB = Zero
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!write(u6,*) 'NEND=',NEND
!write(u6,*) 'V(1)=',V(1)
!write(u6,*) 'V(NEND)=',V(NEND)
!write(u6,*) 'E=',E
!write(u6,*) 'DSOC=',DSOC
! Start iterative loop; try to converge for up to 15 iterations.
! Actually allow only 10 because garble "finds" v=10 with 12 iterations.
Conv = .false.
do IT=1,10
  ITER = IT
  do
    ! OPTIONALLY write when debugging:
    !write(u6,*) 'INNER=',INNER,'If >0, SKIP'
    if (INNER <= 0) then
      !if (E > DSOC) then
      ! For quasibound levels, initialize wave function in "QBOUND"
      !end if
      if (ITER <= 2) then
        ! For  E < DSOC  begin inward integration by using JWKB to estimate
        ! optimum (minimum) inward starting point which will still give
        ! RATOUT < RATST = exp(-XPW) (ca. 1.0e-10) [not needed after 1st 2 ITER]
        NEND = NPP-1
        GB = VBZ(NEND)-E
        ! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
        !write(u6,*) 'VBZ(NEND)=',VBZ(NEND)
        !write(u6,*) 'GB=',GB
        ! ... first do rough inward search for outermost turning point
        Found = .False.
        do M=NEND-NDN,1,-NDN
          ITP2 = M
          GI = VBZ(M)-E
          !write(u6,*) 'VBZ(M)=',VBZ(M)
          !write(u6,*) 'E=',E
          !write(u6,*) 'GI=',GI,'If <= 0, EXIT'
          if ((GI <= Zero) .and. (E <= Zero)) then
            Found = .true.
            exit
          end if
          !if (GI <= Zero) then
          !  Found = .true.
          !  exit
          !end if
          GB = GI
        end do
        if (.not. Found) then
          if (IWR /= 0) write(u6,611) JROT,EO
          KV = -1
          return
        end if
        SM = GB/(GI-GB)
        SM = Half*(One+SM)*sqrt(GB)
        ITP2 = ITP2+2*NDN
        ! OPTIONALLY write when debugging:
        !write(u6,*) 'ITP2=',ITP2,'If < NEND, EXIT'
        if (ITP2 < NEND) then
          ! ... now integrate exponent till JWKB wave fx. would be negligible
          do M=ITP2,NPP-1,NDN
            NEND = M
            SM = SM+sqrt(VBZ(M)-E)*SDRDY(M)**2
            ! OPTIONALLY write when debugging:
            !write(u6,*) 'SM=',SM,'If SM > DXPW, EXIT'
            if (SM > DXPW) exit
          end do
        end if
      end if
      ! Now, checking that {[V-E](r')**2 + FAS} small enuf that Numerov,
      ! stable, and if necessary, step inward till  {[V-E](r')**2 - F} < 10
      do
        GB = V(NEND)-E*DRDY2(NEND)
        ! OPTIONALLY PRINT THESE VARIABLES WHEN DEBUGGING:
        !write(u6,*) 'NEND=',NEND
        !write(u6,*) 'V(NEND)=',V(NEND)
        !write(u6,*) 'DRDY2(NEND)=',DRDY2(NEND)
        !write(u6,*) 'GB=',GB
        if (GB <= Ten) exit
        ! If potential has [V-E] so high that H is (locally) much too large,
        ! then shift outer starting point inward & use WKB starting condition.
        ! [extremely unlikely condition w. WKB initialization]
        NEND = NEND-1
        ! OPTIONALLY write when debugging:
        !write(u6,*) 'NEND=',NEND,'If <=1, EXIT'
        if (NEND <= 1) then
          if (IWR /= 0) write(u6,613)
          KV = -1
          return
        end if
      end do
      if ((ITER <= 1) .and. (IWR >= 2) .and. (NEND < NPP-1)) write(u6,6609) JROT,EO,NEND,YVB(NEND)
      if (NEND == NPP-1) then
        !!! Initialize with node if at end of range  (YMAX= 1)
        NEND = NPP
        SB = Zero
        Y1 = Zero
        SI = One
        GI = GB
      else
        ! For truly bound state initialize wave function as 1-st order WKB
        ! solution increasing inward
        GB = V(NEND)-E*DRDY2(NEND)
        GI = V(NEND-1)-E*DRDY2(NEND-1)
        MS = NEND-1
        if (GI < Zero) then
          ! ERROR condition if  E > V(R)  at outer end of integration range.
          XPR = YMINN+MS*H
          VPR = V(MS)/BFCT
          if (IWR /= 0) write(u6,608) EO,MS,VPR,XPR,IT
          return
        end if
        ! Below is an even stronger condition to return. Basically print an error if the level above 0cm=-1
        ! Comment the three lines below (IF statement) if you want to allow levels above dissociation:
        !write(u6,*) 'EO=',EO
        !write(u6,*) 'GI=',GI
        !if (EO > Zero) then
        !  write(u6,*) 'Level is not bound!'
        !  XPR = YMINN+MS*H
        !  VPR = V(MS)/BFCT
        !  if (IWR /= 0) write(u6,608) EO,MS,VPR,XPR,IT
        !  return
        !end if
        SRTGB = sqrt(VBZ(NEND)-E)
        SRTGI = sqrt(VBZ(NEND-1)-E)
        SB = One
        SI = SB*sqrt(SRTGB/SRTGI)*exp((SRTGB+SRTGI)*Half*(RVB(NEND)-RVB(NEND-1))/YH)
        if (SB > SI) then
          ! WOOPS - JWKB gives inward DEcreasing solution, so initialize with node
          if (IWR /= 0) write(u6,618) JROT,EO,SB/SI
          SI = One
          SB = Zero
          GI = V(NEND-1)-E*DRDY2(NEND-1)
        end if
      end if
      M = NEND-1
      Y1 = (One-HT*GB)*SB
      Y2 = (One-HT*GI)*SI
      WF(NEND) = SB
      WF(NEND-1) = SI
      MS = NEND
      IBEGIN = 3
      if (INNER > 0) IBEGIN = ITP1+2
      ! Actual inward integration loop starts here
      Found = .false.
      do I=IBEGIN,NEND
        M = M-1
        Y3 = Y2+Y2-Y1+GI*SI
        GI = V(M)-E*DRDY2(M)
        SB = SI
        SI = Y3/(One-HT*GI)
        WF(M) = SI
        if (abs(SI) >= 1.0e17_wp) then
          ! Renormalize to prevent overflow of  WF(I)  in classically
          !forbidden region where  (V(I) > E)
          SI = One/SI
          WF(M:MS) = WF(M:MS)*SI
          !MS = M
          Y2 = Y2*SI
          Y3 = Y3*SI
          SB = SB*SI
          SI = One
        end if
        Y1 = Y2
        Y2 = Y3
        ! Test for outermost maximum of wave function.
        ! ... old matching condition - turning point works OK & is simpler.
        !if ((INNER == 0) .and. (abs(SI) <= abs(SB))) then
        !  Found = .true.
        !  exit
        !end if
        !** Test for outer well turning point
        !write(u6,*) 'GI=',GI
        if ((INNER == 0) .and. (GI < Zero)) then
          Found = .true.
          exit
        end if
      end do
      if ((.not. Found) .and. (INNER == 0)) then
        ! Error mode ... inward propagation finds no turning point
        KV = -2
        if (IWR /= 0) write(u6,616) KV,JROT,EO
        KV = -1
        return
      end if
      ! Scale outer part of wave function before proceding
      SI = One/SI
      YIN = Y1*SI
      MSAVE = M
      RR = YMINN+MSAVE*H
      RATOUT = WF(NEND-1)*SI
      WF(MSAVE:NEND) = WF(MSAVE:NEND)*SI
    end if
    !-------------------------------------------------------------------
    !** Set up to prepare for outward integration **********************
    NBEG = 2
    if (INNODE <= 0) then
      ! Option to initialize with zero slope at beginning of the range
      SB = One
      GB = V(1)-E*DRDY2(1)
      Y1 = SB*(One-HT*GB)
      Y2 = Y1+GB*SB*Half
      GI = V(2)-E*DRDY2(2)
      SI = Y2/(One-HT*GI)
    else
      ! Initialize outward integration with a node at beginning of range
      do
        GB = V(NBEG)-E*DRDY2(NBEG)
        if (GB <= Ten) exit
        ! If potential has [V(i)-E] so high that H is (locally) much too
        ! large, then shift inner starting point outward.
        NBEG = NBEG+1
        if (NBEG >= NPP) then
          if (IWR /= 0) write(u6,613)
          KV = -1
          return
        end if
      end do
      if (NBEG == 2) NBEG = 1
      if ((ITER <= 1) .and. (IWR /= 0)) then
        if (NBEG > 1) write(u6,609) JROT,EO,NBEG,YVB(NBEG)
        if (GB <= Zero) write(u6,604) JROT,EO,NBEG,V(NBEG)/BFCT
      end if
      ! Initialize outward wave function with a node:  WF(NBEG) = 0.
      SB = Zero
      SI = One
      GI = V(NBEG+1)-E*DRDY2(NBEG+1)
      Y1 = SB*(One-HT*GB)
      Y2 = SI*(One-HT*GI)
    end if

    WF(NBEG) = SB
    WF(NBEG+1) = SI
    if (INNER > 0) MSAVE = NPP
    ! Actual outward integration loops start here
    Found = .false.
    do I=NBEG+2,MSAVE
      Y3 = Y2+Y2-Y1+GI*SI
      GI = V(I)-E*DRDY2(I)
      SI = Y3/(One-HT*GI)
      WF(I) = SI
      if (abs(SI) >= 1.0e17_wp) then
        ! Renormalize to prevent overflow of  WF(I)  in classically forbidden
        ! region where  V(I) > E
        SI = One/SI
        WF(NBEG:I) = WF(NBEG:I)*SI
        Y2 = Y2*SI
        Y3 = Y3*SI
        SI = One
      end if
      Y1 = Y2
      Y2 = Y3
      ITP1P = I
      ! Exit from this loop at onset of classically allowed region
      if (GI <= Zero) then
        Found = .true.
        exit
      end if
    end do
    Skip = .false.
    if (.not. Found) then
      MS = MSAVE
      if ((INNER == 0) .and. (GB <= Zero)) then
        Skip = .true.
      else
        if (IWR /= 0) write(u6,612) KVIN,JROT,EO,MSAVE
        KV = -1
        return
      end if
    end if
    if (.not. Skip) then
      ! ITP1 is last point of AS-forbidden region & ITP1P 1'st point in allowed
      ITP1 = ITP1P-1
      MS = ITP1
      if (INNER <= 0) then
        do I=ITP1P+1,MSAVE
          Y3 = Y2+Y2-Y1+GI*SI
          GI = V(I)-E*DRDY2(I)
          SI = Y3/(One-HT*GI)
          WF(I) = SI
          if (abs(SI) > 1.0e17_wp) then
            ! Renormalize to prevent overflow of  WF(I) , as needed.
            SI = One/SI
            WF(NBEG:I) = WF(NBEG:I)*SI
            Y2 = Y2*SI
            Y3 = Y3*SI
            SI = One
          end if
          Y1 = Y2
          Y2 = Y3
        end do
        MS = MSAVE
      end if
    end if
    ! Finished outward integration.  Normalize w.r.t. WF(MSAVE)
    SI = One/SI
    YOUT = Y1*SI
    YM = Y2*SI
    RATIN = WF(NBEG+1)*SI
    WF(NBEG:MS) = WF(NBEG:MS)*SI
    if (INNER <= 0) exit
  end do
  !----- Finished numerical integration ... now correct trial energy
  ! DF*H  is the integral of  (WF(I))**2 dR
  DF = sum(DRDY2(NBEG:NEND)*WF(NBEG:NEND)**2)
  !** Add edge correction to DF assuming wave function dies off as simple
  !  exponential past R(NEND);  matters only if WF(NEND) unusually large.
  !if ((E <= DSOC) .and. (WF(NEND) /= 0)) then
  !
  !!! huh ... how do I fix this for AS ??? - or is it no longer necessary ??
  !
  !if  ((KVIN >= -10) .and. (WF(NEND-1)/WF(NEND) > ONe)) DF = DF+WF(NEND)**2/(Two*log(WF(NEND-1)/WF(NEND)))
  !
  !. note that by construction, at this point  WF(MSAVE)= 1.0
  F = -YOUT-YIN+Two*YM+GI
  DOLD = DE
  if (abs(F) <= 1.0e30_wp) then
    DE = F/DF
  else
    F = 9.9e30_wp
    DF = F
    DE = abs(0.01_wp*(DSOC-E))
  end if
  if (IWR > 2) then
    DEPRN = DE/BFCT
    ! RATIN & RATOUT  are wave fx. amplitude at inner/outer ends of range
    ! relative to its value at outermost extremum.
    write(u6,603) IT,EO,F,DF,DEPRN,MSAVE,RR,RATIN,RATOUT,NBEG,ITP1
  end if
  ! Test trial eigenvalue for convergence
  if (abs(DE) <= abs(EPS)) then
    Conv = .true.
    exit
  end if
  E = E+DE
  ! KV >= 998  Option ... Search for highest bound level.  Adjust new
  ! trial energy downward if it would have been above dissociation.
  if ((KVIN >= 998) .and. (E > VMX)) E = VMX-Two*(VMX-E+DE)
  EO = E/BFCT
  if ((IT > 4) .and. (abs(DE) >= abs(DOLD)) .and. (DOLD*DE <= Zero)) then
    ! Adjust energy increment if having convergence difficulties.  Not
    ! usually needed except for some quasibounds extremely near  VMAX .
    ICOR = ICOR+1
    DEP = DE/BFCT
    if (IWR /= 0) write(u6,617) IT,DEP
    DE = Half*DE
    E = E-DE
    EO = E/BFCT
  end if
end do
!** End of iterative loop which searches for eigenvalue ************
!-------------------------------------------------------------------*
if (.not. Conv) then
  ! Convergence fails, so return in error condition
  E = E-DE
  EO = E/BFCT
  DEPRN = DE/BFCT
  if (IWR /= 0) write(u6,620) KVIN,JROT,ITER,DEPRN
  KV = -1
  return
end if
if (IWR /= 0) then
  if (IWR >= 3) write(u6,619)
  if ((abs(RATIN) > RATST) .and. (INNODE > 0) .and. (YMIN > Zero)) write(u6,614) JROT,EO,RATIN
  if ((E < DSOC) .and. (abs(RATOUT) > RATST)) then
    WKBTST = Half*abs(V(NEND)-V(NEND-1))/sqrt((V(NEND)-E)**3)
    if (WKBTST > 1.0e-3_wp) write(u6,615) JROT,EO,RATOUT,RATST,WKBTST
  end if
end if
KKV = 0
! Perform node count on converged solution
PROD = WF(ITP1)*WF(ITP1-1)
J1 = ITP1+1
J2 = NEND-1
do J=J1,J2
  PPROD = PROD
  PROD = WF(J)*WF(J-1)
  if ((PPROD <= Zero) .and. (PROD > Zero)) KKV = KKV+1
end do
KV = KKV

!write(12,699) kv,jrot,EO,nend
!699 format('   v=',i3,'    J='i3,'   E=',f10.3,'   NEND=',i6)

! Normalize & find interval (NBEG,NEND) where WF(I) is non-negligible
SN = One/sqrt(H*DF)
WF(NBEG:NEND) = WF(NBEG:NEND)*SN
if (ITP1 > 1) then
  J = ITP1P
  do I=1,ITP1
    J = J-1
    if (abs(WF(J)) < RATST) exit
  end do
  NBEG = J
  if (NBEG > 1) then
    J = J-1
    WF(1:J) = Zero
  end if
end if
! Move NEND inward to where wavefunction "non-negligible"
J = NEND-1
do I=NBEG,NEND
  if (abs(WF(J)) > RATST) exit
  J = J-1
end do
NEND = J+1
if (NEND < NPP) then
  ! Zero out wavefunction array at distances past NEND
  WF(NEND+1:NPP) = Zero
end if
if (LPRWF < 0) then
  ! If desired, write every |LPRWF|-th point of wave function to a file
  ! on channel-10, starting at mesh point # NBEG for radial distance
  ! YVB(NBEG), with the NPR values separated by mesh step  JPSIQ*YH
  JPSIQ = -LPRWF
  NPR = 1+(NEND-NBEG)/JPSIQ
  ! Write every JPSIQ-th point of the wave function for level  v=KV
  ! J=JROT , beginning at mesh point NBEG & distance RSTT where
  write(10,701) KV,JROT,EO,NPR,YVB(NBEG),YH*JPSIQ,NBEG,JPSIQ
  write(10,702) (YVB(I),WF(I),I=NBEG,NEND,JPSIQ)
end if
if (IWR == 1) write(u6,607) KV,JROT,EO
if (IWR >= 2) then
  REND = ARV*((One+YVB(NEND-1))/(One-YVB(NEND-1)))**(One/PRV)
  RATIN = RATIN*SDRDY(NBEG+1)/SDRDY(MSAVE)
  RATOUT = RATOUT*SDRDY(NEND-1)/SDRDY(MSAVE)
  write(u6,607) KV,JROT,EO,ITER,RR,RATIN,NBEG,REND,RATOUT,NEND-1
end if
! For quasibound levels, calculate width in subroutine "WIDTH"
if ((E > DSOC) .and. (KVIN > -10)) call WIDTHas(KV,JROT,E,EO,DSOC,VBZ,WF,SDRDY,VMX,YMIN,H,BFCT,IWR,ITP1,ITP2,ITP3,INNER,NPP,GAMA)

return

601 format(/' Solve for  v=',I3,'   J=',I3,'   ETRIAL=',ES15.7,'  INNER=',i2,'   WF(1st)  WF(NEND)')
602 format('ITER    ETRIAL',8X,'F(E)      DF(E)     D(E)',5X,'M    yp(M)   /WF(M)    /WF(M)  NBEG  ITP1'/1X,96('-'))
603 format(I3,ES15.7,3es10.2,I6,F7.3,2ES9.1,I5,I6)
604 format('   NOTE:  for  J=',I3,'   EO=',F12.4,' >= V(',i3,')=',F12.4)
607 format('E(v=',I3,',J=',I3,')=',F15.8,I3,' iterations  yp(M)=',F6.3,'  WF(NBEG)/WF(M)=',ES8.1,'   NBEG=',i5/40x,'R(NEND)=', &
           f9.2,'   WF(NEND)/WF(M)=',ES8.1,'   NEND=',i5)
608 format(' *** SCHRQ Error:  E=',F9.2,' > V(',I5,')=',F9.2,'  at  Rmax=',F6.2,'  for  IT=',I2)
6609 format(' *** For  J=',I3,'   E=',ES15.7,"  integration can't start till inside"/21x,'mesh point',I6,' (yp=',f8.4, &
            '),  so YMAX larger than needed')
609 format(' *** For  J=',I3,'   E=',ES15.7,"  integration can't start till past"/23x,'mesh point',I5,' (yp=',f6.2, &
           '),  so YMIN smaller than needed')
610 format(/' Seek highest bound level:   ETRIAL =',ES9.2,17x,'WF(1st)  WF(NEND)')
611 format(' *** SCHRQ inward search at   J=',i3,'   E=',f11.2,' finds no classical region')
612 format(/' *** ERROR *** for   v =',I3,'   J =',I3,'   E =',F12.4,'  Innermost turning point not found by   M = MSAVE =',I5)
613 format(/' *** ERROR in potential array ... V(I) everywhere too big to integrate with given  increment')
614 format(' *** CAUTION *** For  J=',I3,'  E=',G15.8/16x,'WF(first)/WF(Max)=',ES9.2,'  suggests  YMIN  may be too large')
615 format(' ** CAUTION ** For  J=',I3,'  E=',ES13.6,'  WF(NEND)/WF(Max)=',ES8.1,' >',ES8.1/4X,' initialization quality test ', &
           ES8.1,' > 1.D-3   so RMAX may be too small')
616 format(' ** WARNING *** For  v=',I2,', J=',I3,' at  E=',G14.7, &
           ':  inward propagation finds no turning point ... v=-2 means: Trial energy is too low (!), or potential is too weak')
617 format(' *** SCHRQ has a convergence problem, so for  IT=',I2,'  cut  DE=',ES10.2,'  in HALF')
618 format(' *** For  J=',I3,'  E=',F9.2,'  JWKB start gives  SB/SI=',ES10.3,'  so use a node.')
619 format(1X,96('-'))
620 format(' *** CAUTION for  v=',I3,'  J=',I3,"  SCHRQ doesn't converge by  ITER=",I2,'  DE=',ES9.2)
701 format(/2x,'Level  v=',I3,'   J=',I3,'   E=',F12.4,' ,  wave function at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8, &
           '   NBEG=',I4,'   |LPRWF|=',I3)
702 format((4(f10.6,f11.7)))

end subroutine SCHRQas
