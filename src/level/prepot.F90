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
subroutine PREPOT(LNPT,NPP,OMEGA,RR,RM2,VLIM,VV,CNN,NCN,PPAR,QPAR,NSR,NLR,DSCM,REQ,RREF,PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
!** Driver subroutine of package to read parameters and/or generate
!  values of a potential VV(I) at the NPP input distances RR(I).
!====================== Version of  21 Apr 2009 ========================
!**** Subroutine Input:
!----------------------
!  LNPT  is an integer specifying the operational mode:
!      *  LNPT > 0  : for a new case for which all potential-defining
!                     parameters are read in and a description printed
!      *  LNPT <= 0 : if potential points are to be generated in exactly
!                     the same manner as on preceding call, but at
!                     different distances RR(I) (no reads or writes)
!  IAN1 & IAN2 are the atomic numbers and IMN1 & IMN2 the mass numbers
!        of atoms #1 & 2, used (if needed) to specify isotope masses for
!        calculating adiabatic and/or non-adiabatic BOB correction fx.
!  NPP (integer) is the number of input distances  RR(i) (in Angstroms)
!        at which potential values  VV(i) (in cm-1) are to be generated
!  RR  (real array) is set of NPP distances where potential calculated
!  RM2 (real array) on input is the (centrifugal) array of  1/RR(i)**2
!----------------------
!**** Subroutine Output:
!----------------------
!  OMEGA   is the (integer) electronic angular momentum projection q.no.
!  VLIM (cm-1)  is the absolute energy at the potential asymptote
!  VV (real 1D array)  is the set of function values generated (in cm-1)
!  RM2 values returned are (if appropriate) be modified to include BOB
!      corrections to the (centrifugal) potential  1/RR(i)**2
!  NCN is an integer power defining the asymptotically-dominant
!      inverse-power long-range potential tail:  CNN/R**NCN
!  CNN is limiting long-range coefficient in units  cm-1(Angst)^{NCN}
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+ Calls GENINT (which calls PLYINTRP, SPLINT & SPLINE) ,  or POTGEN ++
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Set maximum array dimension for the input function values to be
!  interpolated over & extrapolated beyond

use Constants, only: Zero, Ten
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(inout) :: LNPT, NCN
integer(kind=iwp), intent(in) :: NPP, OMEGA, PPAR, QPAR, NSR, NLR, MMLR(3), NCMM, IVSR, IDSTT
real(kind=wp), intent(in) :: RR(NPP), VLIM, DSCM, REQ, RREF, PARM(4), CMM(3), RHOAB
real(kind=wp), intent(inout) :: RM2(NPP), CNN
real(kind=wp), intent(out) :: VV(NPP)
integer(kind=iwp) :: I, INPTS, J, JWR, LPPOT, LWR, NLIN, NROW, NTP
real(kind=wp) :: D1V(3), D1VB(3), D2V(3), EFACT, RFACT, RH, RWR(3), RWRB(3), VSHIFT, VWR(3), VWRB(3)
! Save variables needed for 'subsequent' LNPT <= 0 calls
integer(kind=iwp), parameter :: NTPMX = 1600
integer(kind=iwp), save :: ILR, IR2, NUSE
real(kind=wp), save :: XI(NTPMX), YI(NTPMX)

D2V(:) = Zero
D1V(:) = Zero
RWR(:) = Zero
RWRB(:) = Zero
VWR(:) = Zero
VWRB(:) = Zero
D1VB(:) = Zero
EFACT = Zero
RFACT = Zero
VSHIFT = Zero
!2 continue
!call READ_INPUT(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS,NTP,LPPOT,IOMEG,VLIM,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ, &
!                RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)

! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!write(u6,*) 'prepot has the following at the beginning:'
!write(u6,*) 'IAN1 = ',IAN1
!write(u6,*) 'IMN1 = ',IMN1
!write(u6,*) 'IAN2 = ',IAN2
!write(u6,*) 'IMN2 = ',IMN2
!write(u6,*) 'CHARGE = ',CHARGE
!write(u6,*) 'NUMPOT = ',NUMPOT
!write(u6,*) 'RH = ',RH
!write(u6,*) 'RMIN = ',RMIN
!write(u6,*) 'PRV = ',PRV
!write(u6,*) 'ARV = ',ARV
!write(u6,*) 'EPS = ',EPS
!write(u6,*) 'NTP = ',NTP
!write(u6,*) 'LPPOT = ',LPPOT
!write(u6,*) 'IOMEG1(now OMEGA) = ',OMEGA
!write(u6,*) 'VLIM = ',VLIM
!write(u6,*) 'IPOTL = ',IPOTL
!write(u6,*) 'PPAR = ',PPAR
!write(u6,*) 'QPAR = ',QPAR
!write(u6,*) 'NSR = ',NSR
!write(u6,*) 'NLR = ',NLR
!write(u6,*) 'IBOB = ',IBOB
!write(u6,*) 'DSCM = ',DSCM
!write(u6,*) 'REQ = ',REQ
!write(u6,*) 'RREF = ',RREF
!write(u6,*) 'NCMM = ',NCMM
!write(u6,*) 'IVSR = ',IVSR
!write(u6,*) 'IDSTT = ',IDSTT
!write(u6,*) 'RHOAB = ',RHOAB
!write(u6,*) 'MMLR = ',MMLR
!write(u6,*) 'CMM = ',CMM
!write(u6,*) 'PARM = ',PARM
!write(u6,*) 'NLEV1 = ',NLEV1
!write(u6,*) 'AUTO1 = ',AUTO1
!write(u6,*) 'LCDC = ',LCDC
!write(u6,*) 'LXPCT = ',LXPCT
!write(u6,*) 'NJM = ',NJM
!write(u6,*) 'JDJR = ',JDJR
!write(u6,*) 'IWF = ',IWF
!write(u6,*) 'LPRWF = ',LPRWF
!write(u6,*) ''
LPPOT = 0
! WHEN COMPILING WITH CMAKE_BUILD_TYPE=GARBLE, NTP GETS CORRUPTED AT
! SOME POINT. FOR NOW WE SUPPORT ONLY ANALYTIC POTENTIALS, SO WE CAN
! RESET NTP TO -1:
NTP = -1
write(u6,*) 'NTP = ',NTP
if (LNPT > 0) then
  !** If NTP > 0    define potential by interpolation over & extrapolation
  !          beyond the NTP read-in turning points using subroutine GENINT
  !   If NTP <= 0   generate a (fully analytic) potential in POTGEN.
  !** If LPPOT > 0  at every |LPPOT|-th point, print potential and
  !        derivatives-by-differences. ***  If  LPPOT < 0  write potential
  !        at every |LPPOT|-th point to channel-8 in a compact format **
  !  OMEGA  is the (integer) total elextronic angular momentum projection
  !         quantum number (required for proper rotational intensities)
  !** VLIM [cm-1]   is the energy associated with the potential asymptote.
  !---------------------------------------------------------------------
  !read(u5,*) NTP,LPPOT,OMEGA,VLIM
  !---------------------------------------------------------------------
  write(u6,600) OMEGA,VLIM
  if (NTP > 0) then
    !** For a pointwise potential (NTP > 0), now read points & parameters
    !  controlling how the interpolation/extrapolation is to be done.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    !** NTP (read above) is number of turning points (XI,YI) to be read in.
    !** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
    !    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE <= 0)
    !    interpolate with cubic spline instead of local polynomials.
    !** If IR2 > 0   interpolate over  YI*XI**2 ; otherwise on  YI  itself
    !     [IR2 > 0 usually improves interpolation for steep repulsive wall]
    !** ILR specifies how to extrapolate beyond largest input distance XI(i)
    !  If ILR < 0   fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
    !  If ILR = 0   fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
    !  If ILR = 1   fit last two points to:  VLIM - A/R**B .
    !** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
    !  inverse-power terms beginning with  1/R**NCN}. *** If CNN /= 0 ,
    !  leading coefficient fixed at  CNN ; otherwise get it from points too.
    !* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
    !* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
    !* If ILR > 3 : successive higher power terms differ by factor  1/R
    !-------------------------------------------------------------------
    !read(u5,*) NUSE,IR2,ILR,NCN,CNN
    !-------------------------------------------------------------------
    ! NTPMX = 1600, I'm setting NTP to 1599 because for some reason, this:
    ! gitlab.com/hpqc-labs/OpenMolcas/-/jobs/3544803772/artifacts/browse/
    NTP = 1599
    if (NTP > NTPMX) then
      write(u6,602) NTP,NTPMX
      !stop
      call ABEND()
    end if
    if (NUSE > 0) write(u6,604) NUSE,NTP
    if (NUSE <= 0) write(u6,606) NTP
    if (IR2 > 0) write(u6,608)
    if ((ILR > 1) .and. (abs(CNN) > Zero)) write(u6,610) CNN,NCN
    !** Read in turning points to be interpolated over
    !** RFACT & EFACT are factors required to convert units of input turning
    !       points (XI,YI) to Angstroms & cm-1, respectively (may be = 1.0)
    !** Turning points (XI,YI) must be ordered with increasing XI(I)
    !** Energy VSHIFT [cm-1] is added to the input potential points to
    !   make their absolute energy consistent with VLIM (often VSHIFT=Te).
    !-----------------------------------------------------------------------
    !read(u5,*) RFACT,EFACT,VSHIFT
    !read(u5,*) (XI(I),YI(I),I=1,NTP)
    !-----------------------------------------------------------------------
    write(u6,612) VSHIFT,RFACT,EFACT
    NROW = (NTP+2)/3
    do J=1,NROW
      if (EFACT <= Ten) then
        write(u6,614) (XI(I),YI(I),I=J,NTP,NROW)
      else
        write(u6,616) (XI(I),YI(I),I=J,NTP,NROW)
      end if
    end do
    write(u6,624)
    YI(1:NTP) = YI(1:NTP)*EFACT+VSHIFT
    XI(1:NTP) = XI(1:NTP)*RFACT
    if (IR2 > 0) YI(1:NTP) = YI(1:NTP)*XI(1:NTP)**2
    if ((abs(YI(NTP)-YI(NTP-1)) <= 0) .and. (XI(NTP) < RR(NPP))) write(u6,618)
  end if
end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
if (NTP > 0) then
  call GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,NCN,CNN)
else
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** If (NTP <= 0) PREPOT uses subroutine POTGEN to generate a fully
  !  analytic potential defined by the following read-in parameters.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
  !  and for all cases except IPOTL=2, the potential asymptote energy is
  !  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
  !  potential minimum and  DSCM  the leading (quadratic) potential coeft.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** IPOTL  specifies the type of potential function to be generated.
  !** PPAR, QPAR, NSR, NLR & NCMM  integers characterize chosen potential
  !** IBOB   specifies whether (if > 0) or not (if <= 0) atomic mass
  !      dependent Born-Oppenheimer breakdown corrections will be included
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !** If IPOTL=1  generate an L.J.(NSR,NLR) potential.
  !** If IPOTL=2  use Seto's modification of Surkus' GPEF expansion in
  !       z = [R**PPAR - Re**PPAR]/[a*R**PPAR + b*Re**PPAR] where
  !       a=PARM(NLR+1) & b=PARM(NLR+2), which incorporates Dunham, SPF,
  !       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
  !       where  c_0[cm-1] is read in as DSCM and the first NLR parameters
  !       PARM(i)'s are the  c_i  (i > 0).  [PPAR is dummy parameter here]
  !  * For Dunham case:  PPAR=1, PARM(NLR+1)= 0.0, PARM(NLR+2)= 1.0
  !  * For SPF case:  PPAR=1, PARM(NLR+1)= 1.0, PARM(NLR+2)= 0.0
  !  * For Ogilvie-Tipping:  PPAR=1, PARM(NLR+1)= 0.5 = PARM(NLR+2)
  !  * NOTE that for Surkus PPAR < 0 case:  z(PPAR,a,b)= z(|PPAR|,-b,-a)
  !      Generate & return the  D_e  value implied by these coefficients.
  !** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
  !      with exponent factor 'beta' defined as a power series of order
  !      max{NLR,NSR} with (max{NLR,NSR}+1) coefficients PARM(i) in vble
  !      y_{PPAR}= (R**PPAR - RREF**PPAR)/(R**PPAR + RREF**PPAR)
  !      where  PPAR >= 1  and inputing  RREF <= 0  sets  RREF= REQ
  !    * For conventional "simple" Morse potential,  NSR=NLR=0 & PPAR dummy
  !*  Special option #1: set  PPAR=0   to produce Wei Hua's 4-parameter
  !      modified Morse function with  b= PARM(1)  and C= PARM(2).
  !** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
  !         with exponent coefficient function
  !              beta(r)= yp*beta_{infty} + [1-yp]*Sum(beta_i{yq}^i
  !         &  yp= y_{PPAR}= (R**PPAR - RREF**PPAR)/(R**PPAR + RREF**PPAR)
  !      where exponent Sum is polynomial of order max{NSR,NRL} in
  !            yq= y_{QPAR}= (R**QPAR - RREF**QPAR)/(R**QPAR + RREF**QPAR)
  !         with NVARB= [max{NSR,NRL}+1] coeffts PARM(j)    and
  !      long-range defined by NCMM inverse-power terms CMM(i)/r^{MMLR(i)}
  !** If IPOTL=5  generate a Double-Exponential Long-Range (DELR)
  !       potential [JCP 119, 7398 (2003)] with additive long-range part
  !       defined by a sum of NCMM damped inverse-power terms, & exponent
  !       polynomial radial variable defined as for the EMO case (IPOTL=3)
  !** If IPOTL=6  generate generalized HFD({m_i},i=1,NCMM) potential.
  !       PARM(1-3) are the parameters defining the HFD damping function
  !       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
  !       PARM(4) the quadratic coefficient in the exponent, and
  !       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
  !              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2]
  !** If IPOTL=7  use Tiemann polynomial potential of order NLR with NLR+1
  !     expansion coefficients a(i) attached to an inverse-power long-range
  !     tail defined by NCMM read-in coefficients plus one additional term,
  !     and an 1/R^{12} (or exponential) inner wall.  NVARB= NLR+4.
  !---------------------------------------------------------------------
  !read(u5,*) IPOTL,PPAR,QPAR,NSR,NLR,NCMM,IBOB
  !read(u5,*) DSCM,REQ,RREF
  !if (IPOTL >= 4) read(u5,*) (MMLR(I),CMM(I),I=1,NCMM)
  !if (NVARB > 0) read(u5,*) (PARM(I),I=1,NVARB)
  !if (IBOB > 0) then
  !  read(u5,*) MN1R,MN2R,PAD,QAD,NU1,NU2,PNA,NT1,NT2
  !  if (NU1 >= 0) read(u5,*) U1INF,(U1(I),I=0,NU1)
  !  if (NU2 >= 0) read(u5,*) U2INF,(U2(I),I=0,NU2)
  !  if (NT1 >= 0) read(u5,*) T1INF,(T1(I),I=0,NT1)
  !  if (NT2 >= 0) read(u5,*) T2INF,(T2(I),I=0,NT2)
  !end if
  !---------------------------------------------------------------------
  NCN = 99
  ! When debugging, you can print De and the first 3 values of V(R):
  !write(u6,*) 'DSCM=',DSCM
  !do I=1,3
  !  write(u6,*) RR(I)
  !end do
  write(u6,*) ''
  write(u6,*) 'Exiting prepot'
  write(u6,*) 'Entering potgen'
  write(u6,*) ''
  ! VV is not yet defined.
  call POTGEN(LNPT,NPP,VLIM,RR,RM2,VV,NCN,CNN,PPAR,QPAR,NSR,NLR,DSCM,REQ,RREF,PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
  !call POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,NCN,CNN)
  write(u6,*) 'Returned from potgen!'
end if
if (LPPOT /= 0) then
  ! If desired, on the first pass (i.e. if LNPT > 0) print the potential
  RH = RR(2)-RR(1)
  INPTS = abs(LPPOT)
  if (LPPOT < 0) then
    ! Option to write resulting function compactly to channel-8.
    !RMIN = RR(1)
    NLIN = NPP/INPTS+1
    write(8,800) NLIN,VLIM
    write(8,802) (RR(I),VV(I),I=1,NPP,INPTS)
  else
    ! Option to print potential & its 1-st three derivatives, the latter
    ! calculated by differences, assuming equally spaced RR(I) values.
    RWRB(:) = Zero
    VWRB(:) = Zero
    D1V(:) = Zero
    write(u6,620)
    NLIN = NPP/(2*INPTS)+1
    RH = INPTS*RH
    do I=1,NLIN
      LWR = 1+INPTS*(I-1)
      do J=1,2
        JWR = LWR+(J-1)*NLIN*INPTS
        if (JWR <= NPP) then
          RWR(J) = RR(JWR)
          VWR(J) = VV(JWR)
          D1V(J) = (VWR(J)-VWRB(J))/(RWR(J)-RWRB(J))
          VWRB(J) = VWR(J)
          D2V(J) = (D1V(J)-D1VB(J))/(RWR(J)-RWRB(J))
          RWRB(J) = RWR(J)
          D1VB(J) = D1V(J)
        else
          RWR(J) = Zero
          VWR(J) = Zero
        end if
        if (I <= 2) then
          D2V(J) = Zero
          if (I == 1) D1V(J) = Zero
        end if
      end do
      write(u6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J=1,2)
    end do
  end if
end if
if (LNPT > 0) write(u6,624)

return

600 format(' State has  OMEGA=',i2,'   and energy asymptote:   Y(lim)=',F12.4,'(cm-1)')
602 format(/' **** ERROR in dimensioning of arrays required',' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
604 format(' Perform',I3,'-point piecewise polynomial interpolation over',I5,' input points')
606 format(' Perform cubic spline interpolation over the',I5,' input points')
608 format(' Interpolation actually performed over modified input array:   Y(I) * r(I)**2')
610 format(' Beyond read-in points extrapolate to limiting asymptotic behaviour:'/20x,'Y(r)  =  Y(lim) - (',ES16.7,')/r**',I2)
612 format(' To make input points Y(i) consistent with  Y(lim),  add  Y(shift)=',F12.4/' Scale input points:  (distance)*', &
           ES16.9,'    (energy)*',ES16.9/13x,'to get required internal units  [Angstroms & cm-1 for potentials]'/ &
           3('      r(i)         Y(i)  ')/3(3X,11('--')))
614 format((3(F13.8,F12.4)))
616 format((3(F12.6,F13.8)))
618 format(/' !!! CAUTION !!! Last two mesh point  YI  values are equal'/17x, &
           'so extrapolation to large  r  will be unreliable !!!'/)
620 format(/'  Function and first 2 derivatives by differences'/2('     r       Y(r)     d1Y/dr1    d2Y/dr2')/2(2X,19('--')))
622 format(2(F8.3,F11.3,ES11.3,ES10.2))
!622 format(2(F7.2,F12.5,ES11.3,ES10.2))
624 format(1x,38('--'))
800 format(/I7,' function values with asymptotic value:',F14.6)
802 format((1X,1(F12.8,F14.6)))

end subroutine PREPOT
