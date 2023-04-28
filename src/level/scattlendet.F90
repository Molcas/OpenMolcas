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
subroutine SCATTLEN(JROT,SL,VLIM,V,WF,BFCT,YMIN,YH,NPP,CNN,NCN,IWR,LPRWF)
!***** Subroutine SCATTLEN, last updated 30 April 2011 ****
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** SCATTLEN solves the radial Schrodinger equation in dimensionless
!  form  d^2WF/dy^2 = - [(VLIM-V(R))*(r')^2 - F(y)]*WF(y) ,  where WF(I)
!  is the wave function,  y  the reduced radial vble. y_p(r), and  VLIM
!  the energy at the potential asymptote, to determine the scattering
!  length  SL.
!** Integrate by Numerov method over NPP mesh points with increment
!  H=YH across range from YMIN  to  YMAX= 1. Then uses interpoaltion
!  over last 5 \phi(y) values to the \phi'(y=1) in order to generate
!  SL from log-derivative equation.  After completing the calculation
!  using mesh YH, repeat it with mesh 2*YH and thn use Richardson
!  extrapolation (RE) to improve the final result.  Scheme used requires
!  PRV= 1.0.
!** Input potential asymptote VLIM has have units (cm-1).
!** On entry, the input potential V(I) must include the centrifugal
!  term and the factor:  'BFCT'=2*mu*(2*pi*YH/hbar)**2  (1/cm-1)
!   = ZMU[u]*YH[Angst]**2/16.85762920 (1/cm-1)  which is also
!  (internally) incorporated into VLIM.  Note that inclusion of the
!  squared integration increment YH**2 saves arithmetic in the
!  innermost loop of the algorithm.
!-----------------------------------------------------------------------
!** Output scattering length SL [Angst] normalized wave function WF(I)
!  and range, NBEG <=  I <=  NEND  over which WF(I) is defined. Define
!  WF(I)=0  outside the range (NBEG,NEND), which is defined by requiring
!  abs(WF(I)) < RATST=1.0e-9  outside.
!** If(LPRWF > 0) print [WRITE(u6,xxx)] wavefx WF(I) every LPRWF-th point.
!* If(LPRWF < 0) every |LPRWF|-th point of the wave function to Channel
!      10 starting at R(NBEG)
!** If(IWR /= 0) print error & warning descriptions
!  If (IWR >= 2) also show end-of-range wave function amplitudes
!  If (IWR >= 3) print also intermediate trial eigenvalues, etc.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use LEVEL_COMMON, only: ARV, DRDY2, PRV, RVB, SDRDY, YVB
use Constants, only: Zero, One, Two, Three, Four, Six, Eight, Ten, Twelve
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: JROT
integer(kind=iwp), intent(in) :: NPP, NCN, IWR, LPRWF
real(kind=wp), intent(out) :: SL, WF(NPP)
real(kind=wp), intent(in) :: VLIM, V(NPP), BFCT, YMIN, YH, CNN
integer(kind=iwp) :: I, ITP1, ITP1P, JPSIQ, NBEG, NNH, NODE, NPR
real(kind=wp) :: C4BAR, DSOC, GI, GN, HT, PHIp1, PHIp2, PHIp3, PHIp4, RATIN, RINC, RSTT, SB, SI, SL2, SLcor, sumSL, sumVV, WF0, &
                 WF1, WF2, WF3, WF4, Y1, Y2, Y3, YMINN, Z4
logical(kind=iwp) :: Found
real(kind=wp), parameter :: RATST = 1.0e-9_wp

WF4 = 0
if (abs(PRV-One) > Zero) then
  ! Scattering length calculation assumes  PRV=1  s.th.  FAS= 0.0
  write(u6,620) PRV
  SL = Zero
  return
end if
Z4 = Four
YMINN = YMIN-YH
HT = One/Twelve
DSOC = VLIM*BFCT
RATIN = Zero
NBEG = 1
C4BAR = Zero
if (NCN == 4) C4BAR = BFCT*CNN/(Two*ARV)**2
! Begin by checking that Numerov is stable at innermost end of range ...
do
  GN = V(NBEG)-DSOC*DRDY2(NBEG)
  if (GN <= Ten) exit
  ! If potential has [V(i)-E] so high that H is (locally) too ;arge,
  ! then shift inner starting point outward.
  NBEG = NBEG+1
  if (NBEG >= NPP) then
    if (IWR /= 0) write(u6,600)
    ! ERROR condition if  E > V(R)  at outer end of integration range.
    ! Return in error mode
    JROT = -1
    return
  end if
end do
if (IWR /= 0) then
  if (NBEG > 1) write(u6,602) JROT,NBEG,YVB(NBEG)
  if (GN <= Zero) write(u6,604) JROT,NBEG,V(NBEG)/BFCT
end if
NNH = (NPP-NBEG)/2
if ((NPP-NBEG) > (2*NNH)) then
  ! If necessary, adjust NBEG by 1 to ensure interval has an even
  ! no. mesh points in order to simplify RE correction step
  NBEG = NBEG+1
  GN = V(NBEG)-DSOC*DRDY2(NBEG)
end if
! Initialize outward wave function with a node:  WF(NBEG) = 0.
SB = Zero
SI = One
GI = V(NBEG+1)-DSOC*DRDY2(NBEG+1)
Y1 = SB*(One-HT*GN)
Y2 = SI*(One-HT*GI)
WF(NBEG) = SB
WF(NBEG+1) = SI
NODE = 0
!sumSL = SI*(GI/SDRDY(NBEG+1))*(One+YVB(NBEG+1))/(One-YVB(NBEG+1))
sumSL = SI*GI*(One+YVB(NBEG+1))
! Actual outward integration loops start here
Found = .false.
do I=NBEG+2,NPP
  Y3 = Y2+Y2-Y1+GI*SI
  GI = V(I)-DSOC*DRDY2(I)
  SI = Y3/(One-HT*GI)
  WF(I) = SI
  !sumSL = sumSL+(GI*SI/SDRDY(I))*(One+YVB(I))/(One-YVB(I))
  sumSL = sumSL+GI*SI*(One+YVB(I))
  if (abs(SI) >= 1.0e36_wp) then
    ! Renormalize to prevent overflow of  WF(I)  in classically forbidden
    ! region where  V(I) > E
    SI = One/SI
    sumSL = sumSL*SI
    WF(NBEG:I) = WF(NBEG:I)*SI
    Y2 = Y2*SI
    Y3 = Y3*SI
    SI = One
  end if
  ITP1 = I
  ! Exit from this loop at onset of classically allowed region
  if (GI <= Zero) then
    Found = .true.
    exit
  end if
  Y1 = Y2
  Y2 = Y3
end do
if (.not. Found) then
  if (IWR /= 0) write(u6,606) JROT,NPP
  ! ERROR condition if  E > V(R)  at outer end of integration range.
  ! Return in error mode
  JROT = -1
  return
end if
ITP1P = ITP1+1
do I=ITP1P,NPP-1 ! Now - integrate automatically to second-last mesh point ...
  Y1 = Y2
  Y2 = Y3
  Y3 = Y2+Y2-Y1+GI*SI
  GI = V(I)-DSOC*DRDY2(I)
  SB = SI
  SI = Y3/(One-HT*GI)
  sumSL = sumSL+GI*SI*(One+YVB(I))
  ! perform node count ...
  if ((SI*SB <= Zero) .and. (abs(SI) > Zero)) NODE = NODE+1
  WF(I) = SI
end do
! Finally ... complete integration to very last mesh point at  y= 1,
Y1 = Y2
Y2 = Y3
Y3 = Y2+Y2-Y1+GI*SI
if (NCN > 4) GI = Zero
if (NCN == 4) GI = -C4BAR
SB = SI
SI = Y3/(One-HT*GI)
WF(NPP) = SI
! Now generate a value for \phi'(y=1) from the WF values using
! "Newton's formula for forward interpolation" as described in
! Sect. 1.4 of K. Smith "Calculation of Atomic Collision Processes"`
PHIp1 = (WF(NPP)-WF(NPP-1))/YH
PHIp2 = PHIp1+(WF(NPP)-Two*WF(NPP-1)+WF(NPP-2))/(Two*YH)
PHIp3 = PHIp2+(WF(NPP)-Three*WF(NPP-1)+Three*WF(NPP-2)-WF(NPP-3))/(Three*YH)
PHIp4 = PHIp3+(WF(NPP)-Four*WF(NPP-1)+Six*WF(NPP-2)-Four*WF(NPP-3)+WF(NPP-4))/(Four*YH)
SL = ARV*(Two*PHIp4/WF(NPP)-One)
write(u6,608) SL,PHIp4/SI,PHIp1,PHIp2,PHIp3,PHIp4
!=======================================================================

! If desired, calculate partial derivatives of scattering length
! w.r.t. parameters.
! DF*H  is the integral of  (WF(I))**2 dR
!if (NPARM > 0) then
!DADPARM(1:NPARM) = Zero
!do I= NBEG,NPP
!  DF = DRDY2(I)*WF(I)**2
!  DADPARM(1:NPARM) = DADPARM(1:NPARM)+DF*DVDP(I,1:NPARM)
!end do
!DADPARM(1:NPARM) = DADPARM(1:NPARM)*YH

if ((abs(RATIN) > RATST) .and. (YMIN > Zero)) write(u6,614) JROT,RATIN
if (LPRWF < 0) then
  ! If desired, write every |LPRWF|-th point of the wave function
  ! to a file on channel-10, starting at the NBEG-th mesh point.
  JPSIQ = -LPRWF
  NPR = 1+(NPP-NBEG)/JPSIQ
  RINC = YH*JPSIQ
  RSTT = YMINN+NBEG*YH
  ! Write every JPSIQ-th point of the wave function, beginning at mesh
  ! point NBEG & distance RSTT where
  ! the NPR values written separated by mesh step RINC=JPSIQ*YH
  write(10,701) JROT,NPR,RSTT,RINC,NBEG,JPSIQ
  write(10,702) (YVB(I),WF(I),I=NBEG,NPP,JPSIQ)
end if

! Now ... re-do SL calculation with twice the step size to allow
! Richardson Extraoplation correction extimation ...
! Initialize outward wave function with a node:  WF(NBEG) = 0.
SB = Zero
SI = One
GN = Z4*(V(NBEG)-DSOC*DRDY2(NBEG))
GI = Z4*(V(NBEG+2)-DSOC*DRDY2(NBEG+2))
Y1 = SB*(One-HT*GN)
Y2 = SI*(One-HT*GI)
WF1 = SB
WF0 = SI
!sumSL = SI*(GI/SDRDY(NBEG+1))*(One+YVB(NBEG+1))/(One-YVB(NBEG+1))
sumSL = SI*GI*(One+YVB(NBEG+2))
! Actual outward integration loops start here
Found = .false.
do I=NBEG+4,NPP,2
  WF1 = WF0
  Y3 = Y2+Y2-Y1+GI*SI
  GI = (V(I)-DSOC*DRDY2(I))*Z4
  SI = Y3/(One-HT*GI)
  WF0 = SI
  !sumSL = sumSL+(GI*SI/SDRDY(I))*(One+YVB(I))/(One-YVB(I))
  sumSL = sumSL+GI*SI*(One+YVB(I))
  if (abs(SI) >= 1.0e36_wp) then
    ! Renormalize to prevent overflow of  WF(I)  in classically forbidden
    ! region where  V(I) > E
    SI = One/SI
    WF1 = WF1*SI
    sumSL = sumSL*SI
    Y2 = Y2*SI
    Y3 = Y3*SI
    SI = One
    WF0 = SI
  end if
  ITP1 = I
  ! Exit from this loop at onset of classically allowed region
  if (GI <= Zero) then
    Found = .true.
    exit
  end if
  Y1 = Y2
  Y2 = Y3
end do
if (.not. Found) then
  if (IWR /= 0) write(u6,606) JROT,NPP
  ! ERROR condition if  E > V(R)  at outer end of integration range.
  ! Return in error mode
  JROT = -1
  return
end if
ITP1P = ITP1+2
WF2 = WF1
WF3 = WF2
do I=ITP1P,NPP,2
  ! Now - integrate automatically to second-last mesh point ...
  Y1 = Y2
  Y2 = Y3
  Y3 = Y2+Y2-Y1+GI*SI
  GI = (V(I)-DSOC*DRDY2(I))*Z4
  if (I == NPP) then
    if (NCN > 4) GI = Zero
    if (NCN == 4) GI = -C4BAR*Z4
  end if
  !if (NCN > 4) GI = Zero
  !... HEY ... RJ should go & figure out how to treat the C4 case!
  !C4bar = BFCT*C4/(4*ARV**2)   ???
  WF4 = WF3
  WF3 = WF2
  WF2 = WF1
  WF1 = WF0
  SI = Y3/(One-HT*GI)
  WF0 = SI
  sumSL = sumSL+GI*SI*(One+YVB(I))
end do
! Now generate a value for  \phi'(y=1) from the WF values using
! "Newton's formula for forward interpolation" as described in
! Sect. 1.4 of K. Smith "Calculation of Atomic Collision Processes"`
PHIp1 = (WF0-WF1)/(Two*YH)
PHIp2 = PHIp1+(WF0-Two*WF1+WF2)/(Four*YH)
PHIp3 = PHIp2+(WF0-Three*WF1+Three*WF2-WF3)/(Six*YH)
PHIp4 = PHIp3+(WF0-Four*WF1+Six*WF2-Four*WF3+WF4)/(Eight*YH)
!... SL2  is scattering length associated with mesh of  2*YH
SL2 = ARV*(Two*PHIp4/WF0-One)
write(u6,608) SL2,PHIp4/SI,PHIp1,PHIp2,PHIp3,PHIp4
! Finally - user Ricardson expraolation of results for mesh  YH  and
! 2*YH  to obtain final optimum  SL estimate!
SLcor = SL+(SL-SL2)/15.0_wp
write(u6,610) YH,SLCOR,SL2,SL
write(8,610) YH,SLCOR,SL2,SL
!write(u6,612) NODE-1
SL = SLcor
! Now ... use second-last mesh point to normalize wavefunction to
! correspond to asymptotic normalization  \psi(r) \sim r .
SI = (RVB(NPP-1)-SL)/(WF(NPP-1)*SDRDY(NPP-1))
SUMVV = Zero
do I=NBEG,NPP-1
  WF(I) = WF(I)*SI
  ! ... and calculate expectation values of  V(r)  in cm-1
  SUMVV = SUMVV+DRDY2(I)*V(I)*WF(I)**2
end do
SUMVV = SUMVV/YH
write(u6,616) SUMVV,BFCT

write(u6,612) NODE-1

return

600 format(/' *** ERROR in potential array ... V(I) everywhere too big to integrate with given  increment')
602 format(' *** For  J=',I3,"  integration can't start till past"/23x,'mesh point',I5,' (yp=',0pf6.2, &
           '),  so YMIN smaller than needed')
604 format('   NOTE:  for  J=',I3,'   V(',i3,')=',F12.4,' <= 0.0')
606 format(/' *** ERROR *** for   J =',I3,'  Innermost turning point not found by   M = MSAVE =',I5)
608 format(/' Calculate  SL=',1PD21.13,'   log-derivative(y=1)=',D20.12/'     with slope convergence:',D21.13/(28x,D21.13))
610 format(/' YH=',f10.7,'  gives  SL(RE)=',1PD21.13,':  SL2=',D21.13/55x,'SL=',D21.13)
612 format(/' Last bound level of this potential is   v=',i3////)
614 format(' *** CAUTION *** For  J=',I3,'   WF(first)/WF(Max)=',D9.2,'  suggests  YMIN  may be too large')
616 format(' Expectation value of  V(r) is:',1PD17.8,'   BFCT=',D17.8)
620 format(/' *** ERROR in scattlen ***  Input  PRV=',F7.3,'  /= 1')
701 format(/2x,'For   J=',I3,',  wave function at',I6,' points.'/7x,'R(1-st)=',F12.8,'   mesh=',F12.8,'   NBEG=',I4,'   |LPRWF|=', &
           I3)
702 format((1X,4(0Pf9.5,1PD13.5)))

end subroutine SCATTLEN
