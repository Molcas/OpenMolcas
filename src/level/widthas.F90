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
!** Subroutine to calculates quasibound level tunneling lifetime/width
!** For relevant theory see Le Roy & Liu [J.Chem.Phys.69,3622-31(1978)]
!  and Connor & Smith [Mol.Phys. 43, 397 (1981)] and Huang & Le Roy
!  [J.Chem.Phys. 119, 7398 (2003); Erratum, ibid, 127, xxxx (2007)]
!** Final level width calculation from Eq.(4.5) of Connor & Smith.
!------------------ Corrected: 12 March 2007 --------------------------
subroutine WIDTHas(KV,JROT,E,EO,DSOC,V,S,SDRDY,VMX,RMIN,H,BFCT,IWR,ITP1,ITP2,ITP3,INNER,NPP,GAMA)
!++ "WIDTH" calls subroutine "LEVQAD" ++++++++++++++++++++++++++++++++++
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

use Constants, only: Zero, One, Two, Four, Ten, Half, Pi, cLight
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KV, JROT, IWR, ITP2, ITP3, INNER, NPP
real(kind=wp), intent(in) :: E, EO, DSOC, V(NPP), S(NPP), SDRDY(NPP), RMIN, H, BFCT
real(kind=wp), intent(inout) :: VMX
integer(kind=iwp), intent(out) :: ITP1
real(kind=wp), intent(out) :: GAMA
integer(kind=iwp) :: I, IMM, IRM, ITP1P, ITP1P1, ITP2M, ITP2M2, ITP2P1, ITP2P2, KVI, KVO, M, M2, NN, NST
real(kind=wp) :: ANS1, ANS2, ARG, COR, D1, D2, D3, DFI, DSGB, DSGN, DWEB, EMSC, EMV, G1, G2, G3, GAMALG, H2, HBW, HBWB, OMEGJC, &
                 PMX, RMINN, RMX, RT, SM, TAU, TAULG, TI, TUN0, U1, U2, VMAX, XJ, XX
logical(kind=iwp) :: Found
! PSI0 = polygamma(1/2) = -gamma-ln(4) [gamma = Euler-Mascheroni constant]
real(kind=wp), parameter :: FACT = 0.5e-2_wp/(Pi*cLight), PSI0 = -1.96351002602142347944_wp
character(len=*), parameter :: LWELL(2) = ['INNER','OUTER']

IMM = 0
PMX = 0
RMINN = RMIN-H
H2 = H*H
! First - locate innermost turning point ...
Found = .false.
do I=1,ITP2
  ITP1 = I
  if (V(I) < E) then
    Found = .true.
    exit
  end if
end do
if (.not. Found) then
  GAMA = Zero
  return
end if
! ITP1 is first mesh point to right of innermost turning point.
loop_1: do
  ITP1P = ITP1+1
  ITP1P1 = ITP1P+1
  IRM = ITP1-1
  ! Calculate JWKB tunneling probability from quadrature over barrier
  ! (ITP2 is first point inside barrier - as determined in QBOUND)
  ITP2P1 = ITP2+1
  ITP2P2 = ITP2+2
  ! ITP2M is the last mesh point before the 2-nd turning point.
  ITP2M = ITP2-1
  ITP2M2 = ITP2-2
  G1 = V(ITP2M)-E
  G2 = V(ITP2)-E
  G3 = V(ITP2P1)-E
  ! Quadrature over barrier starts here.
  call LEVQAD(G1,G2,G3,H,RT,ANS1,ANS2)
  SM = ANS2*SDRDY(ITP2)**2/H
  !SM = ANS2/H
  if (G3 >= Zero) then
    SM = SM+Half*sqrt(G3)*SDRDY(ITP2)**2
    PMX = VMX
    M2 = ITP2P2
    loop_2: do
      Found = .false.
      do I=M2,ITP3
        M = I
        G3 = V(I)-E
        if (V(I) > PMX) PMX = V(I)
        if (G3 < Zero) then
          Found = .true.
          exit
        end if
        SM = SM+sqrt(G3)*SDRDY(I)**2
      end do
      if (.not. Found) then
        if (V(M) > V(M-1)) then
          if (IWR /= 0) write(u6,602) KV,JROT
          return
        end if
        RMX = RMINN+M*H
        U1 = sqrt(G3/(V(M)-DSOC))
        U2 = sqrt((E-DSOC)/(V(M)-DSOC))
        SM = SM-Half*sqrt(G3)+(log((One+U1)/U2)-U1)*RMX*sqrt(V(M)-DSOC)/H
        XJ = (sqrt(One+Four*(V(M)-DSOC)*(RMX/H)**2)-One)*Half
        if (IWR /= 0) write(u6,603) JROT,EO,XJ,RMX
        exit loop_1
      else if (M < ITP3) then
        ! If encounter a double-humped barrier, take care here.
        if (IWR /= 0) write(u6,609) KV,JROT,EO,M
        KVO = 0
        DSGN = sign(One,S(M-1))
        ! Find the effective quantum number for the outer well
        do I=M,ITP3
          DSGB = DSGN
          DSGN = sign(One,S(I))
          if ((DSGN*DSGB) < Zero) KVO = KVO+1
        end do
        KVI = KV-KVO
        if (INNER == 0) then
          ! For levels of outer well, get correct width by changing ITP1
          ITP1 = M
          if (IWR > 0) write(u6,610) KVO,LWELL(2)
        else
          if (IWR > 0) write(u6,610) KVI,LWELL(1)
          ! For "inner-well" levels, locate outer barrier
          do I=M,ITP3
            M2 = I
            G3 = V(I)-E
            if (G3 >= Zero) cycle loop_2
          end do
          exit loop_1
        end if
      else
        G1 = V(M)-E
        G2 = V(M-1)-E
        G3 = V(M-2)-E
        call LEVQAD(G1,G2,G3,H,RT,ANS1,ANS2)
        SM = SM-Half*sqrt(G3)*SDRDY(M-2)**2-sqrt(G2)*SDRDY(M-1)+ANS2*SDRDY(M)**2/H+ANS2/H
        exit loop_1
      end if
    end do loop_2
  end if
end do loop_1
EMSC = -SM/Pi
if (INNER > 0) VMX = PMX
VMAX = VMX/BFCT
! Tunneling factors calculated here ** TUN0 is simple WKB result
! as in Child's eqs.(57c) & (59).
! .....  EPSRJ= -2.* Pi* EMSC
TUN0 = Half*exp(Two*Pi*EMSC)
! ... for permeability calculate Connor-Smith's Eq.(3.7) \omega=OMEGJC
OMEGJC = sqrt(One+Two*TUN0)-One
! ... alternate calculation to give better precision for small TUN0
if (TUN0 < 1.0e-5_wp) OMEGJC = TUN0*(ONe-Half*TUN0*(One-TUN0))
OMEGJC = Four*OMEGJC/(OMEGJC+Two)
! Quadrature for JWKB calculation of vibrational spacing in well HBW
D1 = E-V(IRM)
D2 = E-V(ITP1)
D3 = E-V(ITP1P)
call LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
SM = ANS1*SDRDY(ITP1)**2/H
!SM = ANS1/H
if (D3 >= Zero) then
  SM = SM+Half/sqrt(D3)
  Found = .false.
  do I=ITP1P1,ITP2M2
    IMM = I
    EMV = E-V(I)
    if (EMV < Zero) then
      Found = .true.
      exit
    end if
    SM = SM+SDRDY(I)**2/sqrt(EMV)
  end do
  if (Found) then
    ! If encounter a double-minimum well, take care here.
    D1 = EMV
    D2 = E-V(IMM-1)
    D3 = E-V(IMM-2)
    if (IWR /= 0) write(u6,605) KV,JROT,EO
  else
    D1 = E-V(ITP2)
    D2 = E-V(ITP2M)
    D3 = E-V(ITP2M2)
  end if
  !call LEVQAD(D1,D2,D3,R1,R2,R3,H,RT,ANS1,ANS2)
  call LEVQAD(D1,D2,D3,H,RT,ANS1,ANS2)
  SM = SM-Half*SDRDY(IMM-2)/sqrt(D3)+ANS1*SDRDY(IMM-1)**2/H
end if
!SM = SM-Half*SDRDY(IMM-2)/sqrt(D3)+ANS1/H
! Get HBW in same energy units (1/cm) associated with BFCT
HBW = Two*Pi/(BFCT*SM)
! HBW fix up suggested by Child uses his eqs.(48)&(62) for HBW
! Derivative of complex gamma function argument calculated as
! per eq.(6.1.27) in Abramowitz and Stegun.
NST = int(abs(EMSC)*1.0e2_wp)
NST = max(NST,4)
ARG = PSI0
do I=0,NST
  NN = I
  XX = I+Half
  TI = One/(XX*((XX/EMSC)**2+One))
  ARG = ARG+TI
  if (abs(TI) < 1.0e-10_wp) exit
end do
! ... and use continuum approximation for tail of summation (???)
COR = Half*(EMSC/(NN+One))**2
ARG = ARG+COR-COR**2
! Now use WKL's Weber fx. approx for (?) derivative of barrier integral ..
DWEB = (EO-VMAX)*BFCT/(H2*EMSC)
DFI = (log(abs(EMSC))-ARG)*BFCT/(H2*DWEB)
HBWB = One/(One/HBW+DFI/(Two*Pi))
! Width from formula (4.5) of  Connor & Smith, Mol.Phys.43,397(1981)
! [neglect time delay integral past barrier in their Eq.(4.16)].
if (EMSC > -25.0_wp) then
  GAMA = (HBWB/(Two*Pi))*OMEGJC
  TAU = Zero
  if (GAMA > 1.0e-60_wp) TAU = FACT/GAMA
  ! GAM0 = TUN0*HBW/Pi  is the simple WKB width GAMMA(0) discussed by
  ! Le Roy & Liu in J.C.P.69,3622(1978).
  if (IWR > 0) write(u6,601) TAU,GAMA,HBWB,VMAX
else
  GAMALG = log10(HBWB/(Two*Pi))+Two*Pi*EMSC/log(Ten)
  TAULG = log10(FACT)-GAMALG
  if (IWR > 0) write(u6,611) TAULG,GAMALG,HBWB,VMAX
end if

return

601 format('    Lifetime=',1PD10.3,'(s)   Width=',D10.3,'   dG/dv=',0PF7.2,'   V(max)=',F9.2)
602 format(' *** WARNING ***  For   v =',I3,'   J =',I3,'   cannot calculate width since barrier maximum beyond range')
603 format(' *** For  J=',I3,'  E=',F9.2,'  R(3-rd) beyond range so approx. tunneling calc. uses'/8X, &
           'pure centrifugal potential with  J(app)=',F7.2,'  for  R > R(max)=',F7.2)
605 format(' **** CAUTION *** Width estimate only qualitative, as have a double-minimum well for   E(v=',I3,', J=',I3,')=',F15.7/ &
           15X,'a more stable result may be obtained by searching for the quasibound levels using option: INNER > 0 .')
609 format(' *** CAUTION - Permeability estimate not exact as have a double-humped barrier:  E(v=',I3,', J=',I3,') =',G15.8,I6)
610 format(16X,'(NOTE: this has the node count of a   v=',I3,2X,A5,'-well level')
611 format(12X,'Log10(lifetime/sec)=',F10.5,' ;   Log10(width/cm-1)=',F10.5,'   Spacing=',G12.5,'   V(max)=',G14.7,'(cm-1)')

end subroutine WIDTHas
