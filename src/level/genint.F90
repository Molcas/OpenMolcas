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
subroutine GENINT(LNPT,NPP,XX,YY,NUSE,IR2,NTP,XI,YI,VLIM,ILR,NCN,CNN)
!** GENINT produces a smooth function YY(i) at the NPP input distances
!  XX(i) by performing numerical interpolation over the range of the
!  NTP input function values YI(j) at the distances XI(j), and using
!  analytic functions to extrapolate beyond their range to with an
!  exponential at short range and a form specified by ILR, NCN & CNN
!** ILR specifies how to extrapolate beyond largest given turning pts
!   If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
!   If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
!   If ILR = 1 : fit last two points to:  VLIM - A/R**B .
!* If(ILR >= 2) fit last turning points to:  VLIM - sum(of ILR
!  inverse-power terms beginning with  1/R**NCN). *** If CNN /= 0 ,
!  leading coefficient fixed at  CNN ; otherwise get it from points too.
!* Assume read-in CNN value has units:  ((cm-1)(Angstroms)**'NCN').
!  If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
!  If ILR > 3 : this factor is  1/R .
!=== Calls subroutines PLYINTRP, SPLINT & SPLINE ==================
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use Constants, only: Zero, One, Two
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LNPT, NPP, NUSE, IR2, NTP
real(kind=wp), intent(in) :: XX(NPP), XI(NTP), YI(NTP), VLIM
real(kind=wp), intent(out) :: YY(NPP)
integer(kind=iwp), intent(inout) :: ILR, NCN
real(kind=wp), intent(inout) :: CNN
integer(kind=iwp) :: I, IDER, IFXCN, J, MBEG, MF(10), MFIN, MI(10), MINNER, NCN2, NCN4, NN, NORD, NUST
real(kind=wp) :: ADCSR, DCSR, DUMM(20), DX1, DX2, DX3, EX1, EX2, EX3, PDCSR, VRAT, X1, X2, X3, XJ(20), Y1, Y2, Y3, YJ(20)
logical(kind=iwp) :: DoPrint
integer(kind=iwp), save :: IMX1, ISR, JMAX, JR2, NMX
real(kind=wp), save :: ALR, ASR, BLR, BSR, CLR, CSR

Y3 = 0
Y2 = 0
Y1 = 0
X3 = 0
X1 = 0
NUST = NUSE/2
if (NUSE <= 0) NUST = 2
IDER = 0
! Determine if/where need to begin extrapolation beyond input data
! XX(MI(J))  is the 1-st mesh point past turning point  XI(J) .
! XX(MF(J))  is the last mesh point before turning point  XI(NTP+1-J)
outer: do J=1,NUST
  MI(J) = 1
  MF(J) = 0
  do I=1,NPP
    if (XX(I) <= XI(J)) MI(J) = I+1
    if (XX(I) >= XI(NTP+1-J)) exit outer
    MF(J) = I
  end do
end do outer
if (NUST < 2) then
  MFIN = MI(1)-1
else
  MFIN = MI(2)-1
end if
if (LNPT > 0) then
  !-----------------------------------------------------------------------
  ! For a new case determine analytic functions for extrapolating beyond
  ! the range of the input points (if necessary) on this or later calls.
  ! Try to fit three innermost turning points to  V(R)=A+B*EXP(-C*R).
  ! If unsatisfactory, extrapolate inward with inverse power function
  if (IR2 <= 0) then
    YJ(1:4) = YI(1:4)
  else
    YJ(1:4) = YI(1:4)/XI(1:4)**2
  end if
  X1 = XI(1)
  X2 = XI(2)
  X3 = XI(3)
  Y1 = YJ(1)
  Y2 = YJ(2)
  Y3 = YJ(3)
  if ((Y1-Y2)*(Y2-Y3) <= Zero) then
    ! If 3 innermost points not monotonic, use A+B/X inward extrapoln.
    ISR = 0
    write(u6,600)
  else
    ! Use cubic through innermost points to get initial trial exponent
    ! from ratio of derivatives,  Y''/Y'
    IDER = 2
    ISR = 4
    call PLYINTRP(XI,YJ,ISR,X2,XJ,ISR,IDER)
    CSR = XJ(3)/XJ(2)
    DCSR = abs(CSR*X2)
    if (DCSR > 1.5e2_wp) then
      ! If exponential causes overflows, use inverse power inward extrapoln.
      ISR = 0
      write(u6,602) CSR
    else
      ! Prepare parameters for inward exponential extrapolation
      VRAT = (Y3-Y2)/(Y1-Y2)
      DX1 = X1-X2
      DX3 = X3-X2
      EX2 = One
      ADCSR = 1.0e99_wp
      ! Now iterate (with actual point) to get exact exponent coefficient
      DoPrint = .true.
      do J=1,15
        PDCSR = ADCSR
        EX1 = exp(CSR*DX1)
        EX3 = exp(CSR*DX3)
        DCSR = (VRAT-(EX3-EX2)/(EX1-EX2))/((X3*EX3-X2-(X1*EX1-X2)*(EX3-EX2)/(EX1-EX2))/(EX1-EX2))
        ADCSR = abs(DCSR)
        if (((ADCSR > PDCSR) .and. (ADCSR < 1.0e-8_wp)) .or. (ADCSR < 1.0e-12_wp)) then
          DoPrint = .false.
          exit
        end if
        CSR = CSR+DCSR
      end do
      if (DoPrint) write(u6,604) DCSR
      BSR = (Y1-Y2)/(EX1-EX2)
      ASR = Y2-BSR*EX2
      BSR = BSR*exp(-CSR*X2)
      write(u6,606) X2,ASR,BSR,CSR
    end if
  end if
  if (ISR <= 0) then
    if ((X1*X2) <= Zero) then
      ! If 1'st two mesh points of opposite sign, extrapolate linearly
      ISR = -1
      ASR = Y2
      BSR = (Y2-Y1)/(X2-X1)
      CSR = X2
      write(u6,608) X2,ASR,BSR,CSR
    else
      ! For inward extrapolation as inverse power through 1'st two points ..
      BSR = (Y1-Y2)*X1*X2/(X2-X1)
      ASR = Y1-BSR/X1
      CSR = X2
      write(u6,610) X2,ASR,BSR
    end if
  end if
end if
600 format('  ** CAUTION ** Exponential inward extrapolation fails'/16x,'since first 3 points not monotonic, ... so ...')
602 format(' *** CAUTION ** inward extrapolation exponent coefficient C=',ES12.4/10x,'could cause overflows, ... so ...')
604 format(' *** CAUTION ** after 15 tries inward extrap. exponent coefft change is',ES9.1)
606 format(' Extrapolate to   X <=',F7.4,'  with'/'   Y=',F13.3,SP,ES15.6,' * exp(',SS,ES13.6,'*X)')
608 format(' Extrapolate to   X <=',F8.4,'   with'/'   Y=',F13.3,SP,ES16.7,' * [X - (',SS,F8.4,')]')
610 format(' Extrapolate to  X <=',F8.4,'   with   Y=',F12.3,SP,ES15.6,')/X**1')

if (MFIN > 0) then
  ! If needed, calculate function in inner extrapolation region
  if (ISR > 0) then
    ! ... either as an exponential
    do I=1,MFIN
      EX1 = CSR*XX(I)
      if (abs(EX1) > 1.0e2_wp) EX1 = 1.0e2_wp*sign(One,EX1)
      YY(I) = ASR+BSR*exp(EX1)
    end do
  else if (ISR == 0) then
    ! ... or if that fails, as an inverse power
    YY(1:MFIN) = ASR+BSR/XX(1:MFIN)
  else if (ISR < 0) then
    ! ... or if X changes sign, extrapolate inward linearly
    YY(1:MFIN) = ASR+BSR*(XX(1:MFIN)-CSR)
  end if
end if
! End of inward extrapolation procedure
!-----------------------------------------------------------------------
MINNER = MFIN
if (NUST > 2) then
  ! If(NUSE > 5) minimize spurious behaviour by interpolating with
  ! order less than NUSE on intervals near inner end of range
  do J=3,NUST
    NORD = 2*(J-1)
    MBEG = MI(J-1)
    MFIN = MI(J)-1
    do I=MBEG,MFIN
      call PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
      YY(I) = DUMM(1)
    end do
  end do
end if
! Main interpolation step begins here
!=======================================================================
MBEG = MI(NUST)
MFIN = MF(NUST)
if (MFIN >= MBEG) then
  if (NUSE <= 0) then
    ! Either ... use cubic spline for main interpolation step
    call SPLINT(LNPT,NTP,XI,YI,MBEG,MFIN,XX,YY)
  else
    ! ... or use piecewise polynomials for main interpolation step
    do I=MBEG,MFIN
      call PLYINTRP(XI,YI,NTP,XX(I),DUMM,NUSE,IDER)
      YY(I) = DUMM(1)
    end do
  end if
end if
if (MFIN < NPP) then
  if (NUST <= 2) then
    ! If(NUSE > 5) minimize spurious behaviour by interpolating with
    ! order less than NUSE on intervals near outer end of range
    MBEG = MF(NUST)+1
  else
    NN = NUST-2
    do J=1,NN
      NORD = 2*(NUST-J)
      MBEG = MF(NUST-J+1)+1
      MFIN = MF(NUST-J)
      if (MFIN >= MBEG) then
        do I=MBEG,MFIN
          call PLYINTRP(XI,YI,NTP,XX(I),DUMM,NORD,IDER)
          YY(I) = DUMM(1)
        end do
      end if
    end do
  end if
end if
MBEG = MFIN+1
! In (IR2 > 0) option, now remove X**2 from the interpolated function
if (IR2 > 0) YY(MINNER+1:MFIN) = YY(MINNER+1:MFIN)/XX(MINNER+1:MFIN)**2
! Print test of smoothness at join with analytic inward extrapolation
!if (LNPT > 0) then
!  MST = max(MINNER-4,1)
!  MFN = MST+8
!  if (MFN > NPP) MFN = NPP
!  if (MFN > MFIN) MFN = MFIN
!  if (MINNER > 0) write(u6,611) X2,((XX(I),YY(I),I=J,MFN,3),J=MST,MST+2)
!  611 format('     Verify smoothness of inner join at   X=',F9.5/(3X,3(F10.5,G15.7)))
!end if
!-----------------------------------------------------------------------
! To extrapolate potential beyond range of given turning points ...
if (LNPT > 0) then
  ! On first entry, calculate things needed for extrapolation constants
  Y1 = YI(NTP)
  Y2 = YI(NTP-1)
  Y3 = YI(NTP-2)
  X1 = XI(NTP)
  X2 = XI(NTP-1)
  X3 = XI(NTP-2)
  if (IR2 > 0) then
    Y1 = Y1/X1**2
    Y2 = Y2/X2**2
    Y3 = Y3/X3**2
  end if
end if
! Check inverse-power tail power ...
if (NCN <= 0) NCN = 6
if (ILR < 0) then
  if (LNPT > 0) then
    ! For  ILR < 0  use  Y = VLIM - ALR * exp[-CLR*(X - BLR)**2]
    EX1 = log((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
    EX2 = log((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
    BLR = (X1+X2-(X2+X3)*EX1/EX2)/(Two-Two*EX1/EX2)
    CLR = -EX1/(X1+X2-Two*BLR)
    ALR = (VLIM-Y1)*exp(CLR*(X1-BLR)**2)
    write(u6,614) X2,VLIM,ALR,CLR,BLR
    if (CLR < Zero) then
      ! ... but replace it by an inverse power of exponent constant negative
      write(u6,612)
      ILR = 1
    end if
  end if
  if (ILR /= 1) YY(MBEG:NPP) = VLIM-ALR*exp(-CLR*(XX(MBEG:NPP)-BLR)**2)
else if (ILR == 0) then
  ! For ILR <= 0  use  Y = VLIM - ALR * X**p * exp(-CLR*X)
  if (LNPT > 0) then
    EX1 = log((VLIM-Y1)/(VLIM-Y2))/(X1-X2)
    EX2 = log((VLIM-Y2)/(VLIM-Y3))/(X2-X3)
    DX1 = log(X1/X2)/(X1-X2)
    DX2 = log(X2/X3)/(X2-X3)
    BLR = (EX1-EX2)/(DX1-DX2)
    CLR = BLR*DX1-EX1
    ALR = (VLIM-Y1)*exp(CLR*X1)/X1**BLR
    write(u6,616) X2,VLIM,ALR,BLR,CLR
    if (CLR < Zero) then
      ! ... but replace it by an inverse power of exponent constant negative
      write(u6,612)
      ILR = 1
    end if
  end if
  if (ILR /= 1) YY(MBEG:NPP) = VLIM-ALR*XX(MBEG:NPP)**BLR*exp(-CLR*XX(MBEG:NPP))
end if
if (ILR == 1) then
  ! For  ILR=1 ,  use     Y = VLIM + ALR/X**BLR
  if (LNPT > 0) then
    BLR = log((VLIM-Y2)/(VLIM-Y1))/log(X1/X2)
    ALR = (Y1-VLIM)*X1**BLR
    NCN = nint(BLR)
    write(u6,618) X2,VLIM,ALR,BLR,NCN
  end if
  YY(MBEG:NPP) = VLIM+ALR/XX(MBEG:NPP)**BLR
else
  ! Set constants for long-range extrapolation
  IFXCN = 0
  if ((CNN > Zero) .or. (CNN < Zero)) IFXCN = 1
  NCN2 = NCN+2
  select case (ILR)
    case (2)
      ! For ILR=2 ,  use   Y = VLIM - CNN/X**NCN - BSR/X**(NCN+2)
      ! If CNN held fixed need ILR > 2  to prevent discontinuity
      if (LNPT > 0) then
        if (IFXCN <= 0) CNN = ((VLIM-Y1)*X1**NCN2-(VLIM-Y2)*X2**NCN2)/(X1**2-X2**2)
        ALR = CNN
        BLR = (VLIM-Y1)*X1**NCN2-CNN*X1**2
        write(u6,620) X2,VLIM,CNN,NCN,BLR,NCN2
      end if
      YY(MBEG:NPP) = VLIM-(ALR+BLR/XX(MBEG:NPP)**2)/XX(MBEG:NPP)**NCN
    case (3)
      ! For ILR=3 , use   Y = VLIM - (CN + CN2/X**2 + CN4/X**4)/X**NCN
      if (LNPT > 0) then
        NCN4 = NCN+4
        if (IFXCN > 0) then
          ALR = CNN
          BLR = (((VLIM-Y1)*X1**NCN-ALR)*X1**4-((VLIM-Y2)*X2**NCN-ALR)*X2**4)/(X1**2-X2**2)
          CLR = ((VLIM-Y1)*X1**NCN-ALR-BLR/X1**2)*X1**4
        else
          EX1 = X1**2
          EX2 = X2**2
          EX3 = X3**2
          DX1 = (VLIM-Y1)*X1**NCN4
          DX2 = (VLIM-Y2)*X2**NCN4
          DX3 = (VLIM-Y3)*X3**NCN4
          BLR = (DX1-DX2)/(EX1-EX2)
          ALR = (BLR-(DX2-DX3)/(EX2-EX3))/(EX1-EX3)
          BLR = BLR-ALR*(EX1+EX2)
          CLR = DX1-(ALR*EX1+BLR)*EX1
        end if
        write(u6,622) X2,VLIM,ALR,NCN,BLR,NCN2,CLR,NCN4
      end if
      YY(MBEG:NPP) = VLIM-(ALR+(BLR/XX(MBEG:NPP)**2+CLR/XX(MBEG:NPP)**4))/XX(MBEG:NPP)**NCN
    case (4)
      ! For ILR >= 4,   Y = VLIM-SUM(BB(K)/X**K) , (K=NCN,NMX=NCN+ILR-1)
      if (LNPT > 0) then
        if (NCN <= 0) NCN = 1
        IMX1 = ILR-1
        NMX = NCN+IMX1
        JR2 = 0
        if (IR2 > 0) JR2 = 2
        IDER = 0
        JMAX = ILR
        if (IFXCN > 0) JMAX = IMX1
        write(u6,624) X2,ILR,NCN,VLIM
        if (IFXCN > 0) write(u6,626) NCN,CNN
      end if
      ! Actually extrapolate with polynomial fitted to the last JMAX
      ! values of  (VLIM - YI(I))*XI(I)**NMX  , & then convert back to  YY(I).
      if (MBEG <= NPP) then
        J = NTP-JMAX
        do I=1,JMAX
          J = J+1
          XJ(I) = XI(J)
          YJ(I) = (VLIM-YI(J)/XI(J)**JR2)*XI(J)**NMX
          if (IFXCN > 0) YJ(I) = YJ(I)-CNN*XI(J)**IMX1
        end do
        do I=MBEG,NPP
          call PLYINTRP(XJ,YJ,JMAX,XX(I),DUMM,JMAX,IDER)
          YY(I) = DUMM(1)
          if (IFXCN > 0) YY(I) = YY(I)+CNN*XX(I)**IMX1
          YY(I) = VLIM-YY(I)/XX(I)**NMX
        end do
      end if
  end select
end if
! Finished extrapolation section.

! Test smoothness at outer join to analytic extrapolation function
!if ((LNPT > 0) .and. (MBEG <= NPP)) then
!  MST = MBEG-5
!  IF (MST < 1) MST = 1
!  MFN = MST+8
!  if (MFN > NPP) MFN = NPP
!  write(u6,627) X2,((XX(I),YY(I),I=J,MFN,3),J=MST,MST+2)
!  627 format('     Verify smoothness of outer join at   X=',F9.5/(3X,3(F10.5,G15.7)))
!end if

return

612 format('  *** BUT *** since exponent has positive coefficient, switch form ...')
614 format(' Function for  X >=',F8.4,'   generated as'/'   Y=',F12.4,' - (',ES13.6,') * exp{-',F10.6,' * (r -',F9.6,')**2}')
616 format(' Function for  X >=',F8.4,'   generated as'/'   Y=',F12.4,' - (',ES13.6,') * r**',F10.6,'  * exp{-(',F11.6,'*r)}')
618 format(' Extrapolate to  X >=',F8.4,'  using'/'   Y=',F12.4,SP,ES15.6,'/X**(',SS,ES13.6,')] ,  yielding   NCN=',I3)
620 format(' Extrapolate to  X >=',F8.4,'  using'/'   Y=',F12.4,' - [',ES13.6,'/X**',I1,SP,ES14.6,'/X**',SS,I1,']')
622 format(' Extrapolate to  X >=',F8.4,'  using'/'   Y=',F12.4,' - [',ES13.6,'/X**',I1,SP,ES14.6,'/X**',SS,I1,SP,ES14.6,'/X**', &
           SS,I2,']')
624 format(' Function for  X >=',F7.3,'  generated by',I3,'-point inverse-power interpolation'/'   with leading term  1/r**',I1, &
           '  relative to dissociation limit   YLIM=',F11.3)
626 format('   and (dimensionless) leading coefficient fixed as   C',I1,'=',G15.8)

end subroutine GENINT
