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
!
!***********************************************************************
! Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
!***********************************************************************
      SUBROUTINE PREPOT(LNPT,IAN1,IAN2,IMN1,IMN2,NPP,OMEGA,RR,RM2,VLIM, &
     &  VV,CNN,NCN,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,     &
     &  MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
!** Driver subroutine of package to read parameters and/or generate
!  values of a potential VV(I) at the NPP input distances RR(I).
!====================== Version of  21 Apr 2009 ========================
!**** Subroutine Input:
!----------------------
!  LNPT  is an integer specifying the operational mode:
!      *  LNPT > 0  : for a new case for which all potential-defining
!                     parameters are read in and a description printed
!      *  LNPT.le.0 : if potential points are to be generated in exactly
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
!  interpolated over & extrapolated beyong
      INTEGER NTPMX
      PARAMETER (NTPMX= 1600)
      INTEGER I,J,IAN1,IAN2,IMN1,IMN2,INPTS,ILR,IR2,JWR,LNPT,LPPOT,LWR, &
     &  NCN,NLIN,NPP,NROW,NTP,NUSE,OMEGA,NSR,NLR,IBOB,NCMM,IVSR,IDSTT,  &
     &  MMLR(3),PPAR,QPAR
      REAL*8  RFACT,EFACT,RH,RMIN,VLIM,VSHIFT,VV(NPP),RR(NPP),RM2(NPP), &
     &  XI(NTPMX),YI(NTPMX),RWR(20),RWRB(3),VWR(20),VWRB(3),D1V(3),     &
     &  D1VB(3),D2V(3),CNN,DSCM,REQ,RREF,RHOAB,CMM(3),PARM(4)
!
!** Save variables needed for 'subsequent' LNPT.le.0 calls
      SAVE ILR,IR2,LPPOT,NTP,NUSE
      SAVE VSHIFT,XI,YI
      DO I=1,3
       D2V(I) = 0.d0
       D1V(I) = 0.d0
       RWR(I) = 0.d0
       RWRB(I) = 0.d0
       VWR(I) = 0.d0
       VWRB(I) = 0.d0
       D1VB(I) = 0.d0
      ENDDO
      EFACT = 0.d0
      RFACT = 0.d0
      VSHIFT = 0.d0
!   2 CALL READ_INPUT(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,
!    1 ARV,EPS,NTP,LPPOT,IOMEG,VLIM,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,
!    2 REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV,AUTO1,LCDC,
!    3 LXPCT,NJM,JDJR,IWR,LPRWF)
!
! OPTIONALLY WRITE THESE VARIABLES WHEN DEBUGGING:
!      WRITE(6,*) 'prepot.f has the following at the beginning:'
!      WRITE(6,*) 'IAN1 = ',IAN1
!      WRITE(6,*) 'IMN1 = ',IMN1
!      WRITE(6,*) 'IAN2 = ',IAN2
!      WRITE(6,*) 'IMN2 = ',IMN2
!!     WRITE(6,*) 'CHARGE = ',CHARGE
!!     WRITE(6,*) 'NUMPOT = ',NUMPOT
!!     WRITE(6,*) 'RH = ',RH
!!     WRITE(6,*) 'RMIN = ',RMIN
!!     WRITE(6,*) 'PRV = ',PRV
!!     WRITE(6,*) 'ARV = ',ARV
!!     WRITE(6,*) 'EPS = ',EPS
!!     WRITE(6,*) 'NTP = ',NTP
!!     WRITE(6,*) 'LPPOT = ',LPPOT
!      WRITE(6,*) 'IOMEG1(now OMEGA) = ',OMEGA
!!     WRITE(6,*) 'VLIM = ',VLIM
!      WRITE(6,*) 'IPOTL = ',IPOTL
!      WRITE(6,*) 'PPAR = ',PPAR
!      WRITE(6,*) 'QPAR = ',QPAR
!      WRITE(6,*) 'NSR = ',NSR
!      WRITE(6,*) 'NLR = ',NLR
!      WRITE(6,*) 'IBOB = ',IBOB
!      WRITE(6,*) 'DSCM = ',DSCM
!      WRITE(6,*) 'REQ = ',REQ
!      WRITE(6,*) 'RREF = ',RREF
!      WRITE(6,*) 'NCMM = ',NCMM
!      WRITE(6,*) 'IVSR = ',IVSR
!      WRITE(6,*) 'IDSTT = ',IDSTT
!      WRITE(6,*) 'RHOAB = ',RHOAB
!      WRITE(6,*) 'MMLR = ',MMLR
!      WRITE(6,*) 'CMM = ',CMM
!      WRITE(6,*) 'PARM = ',PARM
!!     WRITE(6,*) 'NLEV1 = ',NLEV1
!!     WRITE(6,*) 'AUTO1 = ',AUTO1
!!     WRITE(6,*) 'LCDC = ',LCDC
!!     WRITE(6,*) 'LXPCT = ',LXPCT
!!     WRITE(6,*) 'NJM = ',NJM
!!     WRITE(6,*) 'JDJR = ',JDJR
!!     WRITE(6,*) 'IWF = ',IWF
!!     WRITE(6,*) 'LPRWF = ',LPRWF
!      WRITE(6,*) ''
      LPPOT= 0
! WHEN COMPILING WITH CMAKE_BUILD_TYPE=GARBLE, NTP GETS CORRUPTED AT
! SOME POINT. FOR NOW WE SUPPORT ONLY ANALYTIC POTENTIALS, SO WE CAN
! RESET NTP TO -1:
      WRITE(6,*) 'NTP = ',NTP
      NTP=-1
      IF(LNPT.GT.0) THEN
!** If NTP > 0    define potential by interpolation over & extrapolation
!          beyond the NTP read-in turning points using subroutine GENINT
!   If NTP.le.0   generate a (fully analytic) potential in POTGEN.
!** If LPPOT > 0  at every |LPPOT|-th point, print potential and
!        derivatives-by-differences. ***  If  LPPOT < 0  write potential
!        at every |LPPOT|-th point to channel-8 in a compact format **
!  OMEGA  is the (integer) total elextronic angular momentum projection
!         quantum number (required for proper rotational intensities)
!** VLIM [cm-1]   is the energy associated with the potential asymptote.
!-----------------------------------------------------------------------
!         READ(5,*) NTP, LPPOT, OMEGA, VLIM
!-----------------------------------------------------------------------
          WRITE(6,600) OMEGA,VLIM
          IF(NTP.GT.0) THEN
!** For a pointwise potential (NTP > 0), now read points & parameters
!  controlling how the interpolation/extrapolation is to be done.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** NTP (read above) is number of turning points (XI,YI) to be read in.
!** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
!    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
!    interpolate with cubic spline instead of local polynomials.
!** If IR2 > 0   interpolate over  YI*XI**2 ; otherwise on  YI  itself
!     [IR2 > 0 usually improves interpolation for steep repulsive wall]
!** ILR specifies how to extrapolate beyond largest input distance XI(i)
!  If ILR < 0   fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
!  If ILR = 0   fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
!  If ILR = 1   fit last two points to:  VLIM - A/R**B .
!** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
!  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
!  leading coefficient fixed at  CNN ; otherwise get it from points too.
!* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
!* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
!* If ILR > 3 : successive higher power terms differ by factor  1/R
!-----------------------------------------------------------------------
!             READ(5,*) NUSE, IR2, ILR, NCN, CNN
!-----------------------------------------------------------------------
! NTPMX = 1600, I'm setting NTP to 1599 because for some reason, this:
! gitlab.com/hpqc-labs/OpenMolcas/-/jobs/3544803772/artifacts/browse/
              NTP = 1599
              IF(NTP.GT.NTPMX) THEN
                  WRITE(6,602) NTP,NTPMX
!                 STOP
                  CALL ABEND()
                  ENDIF
              IF(NUSE.GT.0) WRITE(6,604) NUSE,NTP
              IF(NUSE.LE.0) WRITE(6,606) NTP
              IF(IR2.GT.0) WRITE(6,608)
              IF((ILR.GT.1).AND.(DABS(CNN).GT.0.D0))WRITE(6,610)CNN,NCN
!** Read in turning points to be interpolated over
!** RFACT & EFACT are factors required to convert units of input turning
!       points (XI,YI) to Angstroms & cm-1, respectively (may be = 1.d0)
!** Turning points (XI,YI) must be ordered with increasing XI(I)
!** Energy VSHIFT [cm-1] is added to the input potential points to
!   make their absolute energy consistent with VLIM (often VSHIFT=Te).
!-----------------------------------------------------------------------
!             READ(5,*) RFACT, EFACT, VSHIFT
!             READ(5,*) (XI(I), YI(I), I= 1,NTP)
!-----------------------------------------------------------------------
              WRITE(6,612) VSHIFT, RFACT, EFACT
              NROW= (NTP+2)/3
              DO  J= 1,NROW
                  IF(EFACT.LE.10.D0) THEN
                      WRITE(6,614) (XI(I),YI(I),I= J,NTP,NROW)
                    ELSE
                      WRITE(6,616) (XI(I),YI(I),I= J,NTP,NROW)
                    ENDIF
                  ENDDO
                  WRITE(6,624)
              DO  I= 1,NTP
                  YI(I)= YI(I)*EFACT+ VSHIFT
                  XI(I)= XI(I)*RFACT
                  ENDDO
              IF(IR2.GT.0) THEN
                  DO  I= 1,NTP
                      YI(I)= YI(I)*XI(I)**2
                      ENDDO
                  ENDIF
              IF((DABS(YI(NTP)-YI(NTP-1)).LE.0).AND.                    &
     &                              (XI(NTP).LT.RR(NPP))) WRITE(6,618)
              ENDIF
          ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      IF(NTP.GT.0) THEN
          CALL GENINT(LNPT,NPP,RR,VV,NUSE,IR2,NTP,XI,YI,VLIM,ILR,       &
     &                                                        NCN,CNN)
        ELSE
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** If (NTP.le.0) PREPOT uses subroutine POTGEN to generate a fully
!  analytic potential defined by the following read-in parameters.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
!  and for all cases except IPOTL=2, the potential asymptote energy is
!  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
!  potential minimum and  DSCM  the leading (quadratic) potential coeft.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** IPOTL  specifies the type of potential function to be generated.
!** PPAR, QPAR, NSR, NLR & NCMM  integers characterize chosen potential
!** IBOB   specifies whether (if > 0) or not (if .le. 0) atomic mass
!      dependent Born-Oppenheimer breakdown corrections will be included
!- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
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
!      where  PPAR.ge.1  and inputing  RREF.le.0  sets  RREF= REQ
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
!----------------------------------------------------------------------
!++    READ(5,*) IPOTL, PPAR, QPAR, NSR, NLR, NCMM, IBOB
!++    READ(5,*) DSCM, REQ, RREF
!++    IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
!++    IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
!++    IF(IBOB.GT.0) THEN
!++        READ(5,*) MN1R, MN2R, PAD, QAD, NU1, NU2, PNA, NT1, NT2
!++        IF(NU1.GE.0) READ(5,*) U1INF, (U1(I), I=0,NU1)
!++        IF(NU2.GE.0) READ(5,*) U2INF, (U2(I), I=0,NU2)
!++        IF(NT1.GE.0) READ(5,*) T1INF, (T1(I), I=0,NT1)
!++        IF(NT2.GE.0) READ(5,*) T2INF, (T2(I), I=0,NT2)
!++        ENDIF
!++    ENDIF
!-----------------------------------------------------------------------
          NCN= 99
! When debugging, you can print De and the first 3 values of V(R):
!         WRITE(6,*) 'DSCM=',DSCM
!         DO I=1,3
!          WRITE(6,*) RR(I)
!         ENDDO
          WRITE(6,*) ''
          WRITE(6,*) 'Exiting prepot.f'
          WRITE(6,*) 'Entering potgen.f'
          WRITE(6,*) ''
! VV is not yet defined.
          CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,      &
     &  NCN,CNN,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,MMLR,   &
     &  CMM,NCMM,IVSR,IDSTT,RHOAB)
!         CALL POTGEN(LNPT,NPP,IAN1,IAN2,IMN1,IMN2,VLIM,RR,RM2,VV,
!    1                                                        NCN,CNN)
          WRITE(6,*) 'Returned from potgen.f!'
        ENDIF
      IF(LPPOT.NE.0) THEN
!** If desired, on the first pass (i.e. if LNPT > 0) print the potential
          RH= RR(2)-RR(1)
          INPTS= IABS(LPPOT)
          IF(LPPOT.LT.0) THEN
!** Option to write resulting function compactly to channel-8.
              RMIN= RR(1)
! Make sure RMIN is "referenced":
              RR(1)= RMIN
              NLIN= NPP/INPTS+ 1
              WRITE(8,800) NLIN,VLIM
              WRITE(8,802) (RR(I),VV(I),I= 1,NPP,INPTS)
            ELSE
!** Option to print potential & its 1-st three derivatives, the latter
!  calculated by differences, assuming equally spaced RR(I) values.
              DO  I= 1,3
                  RWRB(I)= 0.d0
                  VWRB(I)= 0.d0
                  D1V(I)= 0.d0
                  ENDDO
              WRITE(6,620)
              NLIN= NPP/(2*INPTS)+1
              RH= INPTS*RH
              DO  I= 1,NLIN
                  LWR= 1+ INPTS*(I-1)
                  DO  J= 1,2
                      JWR= LWR+(J-1)*NLIN*INPTS
                      IF(JWR.LE.NPP) THEN
                          RWR(J)= RR(JWR)
                          VWR(J)= VV(JWR)
                          D1V(J)= (VWR(J)-VWRB(J))/(RWR(J)-RWRB(J))
                          VWRB(J)= VWR(J)
                          D2V(J)= (D1V(J)-D1VB(J))/(RWR(J)-RWRB(J))
                          RWRB(J)= RWR(J)
                          D1VB(J)= D1V(J)
                        ELSE
                          RWR(J)= 0.d0
                          VWR(J)= 0.d0
                        ENDIF
                      IF(I.LE.2) THEN
                          D2V(J)= 0.d0
                          IF(I.EQ.1) D1V(J)= 0.d0
                          ENDIF
                      ENDDO
                  WRITE(6,622) (RWR(J),VWR(J),D1V(J),D2V(J),J= 1,2)
                  ENDDO
            ENDIF
          ENDIF
      IF(LNPT.GT.0) WRITE(6,624)
      RETURN
  600 FORMAT(' State has  OMEGA=',i2,'   and energy asymptote:   Y(lim)=&
     &',F12.4,'(cm-1)')
  602 FORMAT(/' **** ERROR in dimensioning of arrays required'          &
     & ,' by GENINT;   No. input points ',I5,' > NTPMX =',I4)
  604 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov&
     &er',I5,' input points' )
  606 FORMAT(' Perform cubic spline interpolation over the',I5,         &
     &  ' input points' )
  608 FORMAT(' Interpolation actually performed over modified input arra&
     &y:   Y(I) * r(I)**2')
  610 FORMAT( ' Beyond read-in points extrapolate to limiting asymptotic&
     & behaviour:'/20x,'Y(r)  =  Y(lim) - (',D16.7,')/r**',I2)
  612 FORMAT(' To make input points Y(i) consistent with  Y(lim),  add' &
     & ,'  Y(shift)=',F12.4/' Scale input points:  (distance)*',        &
     & 1PD16.9,'  &  (energy)*',D16.9/13x,'to get required internal unit&
     &s  [Angstroms & cm-1 for potentials]'/                            &
     &  3('      r(i)         Y(i)  ')/3(3X,11('--')))
  614 FORMAT((3(F13.8,F12.4)))
  616 FORMAT((3(F12.6,F13.8)))
  618 FORMAT(/' !!! CAUTION !!! Last two mesh point  YI  values are equa&
     &l'/17x,'so extrapolation to large  r  will be unreliable !!!'/)
  620 FORMAT(/'  Function and first 2 derivatives by differences'/      &
     &  2('     r       Y(r)     d1Y/dr1    d2Y/dr2')/2(2X,19('--')))
  622 FORMAT(2(0PF8.3,F11.3,1PD11.3,D10.2))
! 622 FORMAT(2(0PF7.2,F12.5,1PD11.3,D10.2))
  624 FORMAT(1x,38('--'))
  800 FORMAT(/I7,' function values with asymptotic value:',F14.6)
  802 FORMAT((1X,1(F12.8,F14.6)))
      END
