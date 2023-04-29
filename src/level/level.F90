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
subroutine LEVEL(RC)
!******  Eigenvalue program  LEVEL 2022 ********************************
!   !!!! Form modified for handling Stolyarov radial variable  !!!!!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Program for calculating eigenvalues and eigenfunctions (and if
!  desired, also various expectation values & matrix elements) of a
!  one-dimensional potential, and/or matrix elements (& Franck-Condon
!  factors) between levels of two different potentials.
!** As with most similar codes, the original version of this program was
!  based on the Franck-Condon intensity program of R.N. Zare, report
!  UCRL-10925(1963), but the present version is massively modified.
!** This program is unique in that it can:  (1) automatically locate &
!      calculate the widths of quasibound levels (orbiting resonances);
!  (2) can calculate diatomic molecule centrifugal distortion constants;
!
!      it will also automatically generate the eigenvalues etc. for all
!      vibrational and/or rotational levels of a given well-behaved
!      single-minimum potential.
!***** Main calling and I/O routines.  Last Updated  28 June 2009 *****

use LEVEL_COMMON, only: ARV, DRDY2, NDIMR, PRV, RVB, SDRDY, YVB, VBZ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Quart, Pi, uToau
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(out) :: RC
! Dimensions for  potential arrays  and  vib. level arrays.
integer(kind=iwp), parameter :: MORDRMX = 20, NTPMX = 1600, RORDR = 7, VIBMX = 400
integer(kind=iwp) :: AFLAG, AUTO1, AUTO2, CHARGE, I, IAN1, IAN2, IBOB, ICOR, IDSTT, IJ(VIBMX), IJD, ILEV1, ILEV2, ILRF, IMN1, &
                     IMN2, INNER, INNOD1, INNOD2, INNR1(0:VIBMX), INNR2(0:VIBMX), IOMEG1, IOMEG2, IPOTL, IQT, IR2F, IRFN, &
                     IV(VIBMX), IV2(VIBMX), IVD, IVS, IVSR, IWR, J, J2DD, J2DL, J2DU, JCT, JDJR, JLEV, JREF, JROT, JROT2, &
                     JWR(VIBMX), KV, KV2, KVIN, LCDC, LPPOT, LPRWF, LRPT, LXPCT, M, MMLR(3), MORDR, NBEG, NBEG2, NCMM, NCN1, NCN2, &
                     NCNF, NEND, NEND2, NFP, NJM, NJMM, NLEV, NLEV1, NLEV2, NLP, NLR, NPP, NRFN, NROW, NSR, NTP, NUMPOT, NUSEF, &
                     PPAR, QPAR, SINNER, VMAX, VMAX1, VMAX2, WARN
real(kind=wp) :: aRVp, BEFF, BFCT, BvWN, BZ, CMM(3), CNN1, CNN2, CNNF, DEJ, DM(0:MORDRMX), DRDY, DREF, DREFP, DSCM, EJ, EJ2, EJP, &
                 EJREF, EO, EO2, EPS, ESLJ(VIBMX), ESOLN(VIBMX), FFAS, GAMA, GB, GBB, GI, GV(0:VIBMX), MASS1, MASS2, MFACTF, &
                 PARM(4), PINV, PMAX1, PMAX2, PW, RCNST(RORDR), REQ, RFACTF, RFLIM, RH, RHOAB, RMIN, RR, RREF, RRp, SL, SOMEG1, &
                 SOMEG2, VD, VDMV, VLIM1, VLIM2, WV, XIF(NTPMX), XX, YH, YH2, YIF(NTPMX), YMAX, YMIN, YMINN, ZK1(0:VIBMX,0:RORDR), &
                 ZK2(0:VIBMX,0:RORDR), ZMU
logical(kind=iwp) :: Skip
character(len=78) :: TITL
character(len=2) :: NAME1, NAME2
real(kind=wp), allocatable :: FAS(:), RFN(:), RM2(:), RM22(:), RRM2(:), RRM22(:), V1(:), V1BZ(:), V2(:), V2BZ(:), VJ(:), WF1(:), &
                              WF2(:)

call mma_allocate(RVB,NDIMR,label='RVB')
call mma_allocate(YVB,NDIMR,label='YVB')
call mma_allocate(DRDY2,NDIMR,label='DRDY2')
call mma_allocate(SDRDY,NDIMR,label='SDRDY')
call mma_allocate(VBZ,NDIMR,label='VBZ')

call mma_allocate(FAS,NDIMR,label='FAS')
call mma_allocate(V1,NDIMR,label='V1')
call mma_allocate(V2,NDIMR,label='V2')
call mma_allocate(VJ,NDIMR,label='VJ')
call mma_allocate(V1BZ,NDIMR,label='V1BZ')
call mma_allocate(V2BZ,NDIMR,label='V2BZ')
call mma_allocate(WF1,NDIMR,label='WF1')
call mma_allocate(WF2,NDIMR,label='WF2')

call mma_allocate(RFN,NDIMR,label='RFN')
call mma_allocate(RRM2,NDIMR,label='RRM2')
call mma_allocate(RM2,NDIMR,label='RM2')
call mma_allocate(RRM22,NDIMR,label='RRM22')
call mma_allocate(RM22,NDIMR,label='RM22')
PINV = One
NLEV2 = -1
AUTO2 = 0
VMAX2 = 0
IOMEG2 = 0
SOMEG2 = Zero
CNN2 = Zero
PMAX2 = Zero
MASS2 = Zero
NCN2 = 0
NEND2 = 0
NBEG2 = 0
INNOD2 = 0
EJREF = Zero
IV2(:) = 0
RCNST(:) = Zero
! Default (Q-branch) defining J-increments for matrix element calcn.
J2DL = 0
J2DU = 0
J2DD = 1
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Begin by reading in the (integer) atomic numbers and mass numbers
!  defining the effective reduced mass of the system considered.
!** IAN1 & IAM2, and IMN1 & IMN2 are, respectively, the atomic numbers
!    and the mass numbers identifying the atoms forming the molecule.
!    Their masses are extracted from data subroutine MASSES and used
!    to generate the the reduced mass ZMU.
!** If  IMN1  or  IMN2  lie outside the range of mass numbers for normal
!  stable isotopes of that species, subroutine MASSES returns the
!  average atomic mass based on the natural isotope abundance.
!** If the read-in value of IAN1 and/or IAN2 is <= 0, then instead of
!  using the MASS table, read an actual particle mass for it/them.
!** CHARGE (integer) is the charge on the molecule (=0 for neutral). If
!   (CHARGE /= 0)  generate & use Watson's charge-adjusted reduced mass.
!** Parameter NUMPOT specifies whether to calculate eigenvalues etc. for
!  a single potential (when NUMPOT <= 1), or to generate two independent
!  potentials & calculate matrix elements coupling levels of one to
!  levels of the other (for NUMPOT >= 2).
!-----------------------------------------------------------------------

! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING:
!write(u6,*) IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS
!write(u6,*) NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB
!write(u6,*) DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1
! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING (ANOTHER WAY):
!write(u6,*) AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
!write(u6,*) 'level has the following after CALL LEVEL_RDINP:'
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
!write(u6,*) 'IOMEG1 = ',IOMEG1
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
!read(u5,*,iostat=istatus)
outer: do
  !if (istatus < 0) then
  !  exit outer
  !else if (istatus > 0) then
  !  call abend()
  !end if
  call LEVEL_RDINP(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS,NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM, &
                   REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)
  !read(u5,*,iostat=istatus) IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT
  !if (istatus < 0) then
  !  exit outer
  !else if (istatus > 0) then
  !  call abend()
  !end if
  !---------------------------------------------------------------------
  ! Subroutine MASSES returns the names of the atoms NAMEi,ground
  ! electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
  ! mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
  if ((IAN1 > 0) .and. (IAN1 <= 109)) then
    call MASSES(IAN1,IMN1,NAME1,MASS1)
  else
    ! If particle-i is not a normal atomic isotope, read a 2-character
    ! name (enclosed between '', as in 'mu') and its actual mass.
    !-------------------------------------------------------------------
    !read(u5,*) NAME1,MASS1
    !-------------------------------------------------------------------
  end if
  if ((IAN2 > 0) .and. (IAN2 <= 109)) then
    call MASSES(IAN2,IMN2,NAME2,MASS2)
  else
    !-------------------------------------------------------------------
    !read(u5,*) NAME2,MASS2
    !-------------------------------------------------------------------
  end if
  ZMU = MASS1*MASS2/(MASS1+MASS2-CHARGE/uToau)
  !=====================================================================
  ! TITL is a title or output header of up to 78 characters, read on a
  !   single line enclosed between single quotes: e.g.  'title of problem'
  !=====================================================================
  !read(u5,*) TITL
  TITL = 'Beginning execution of LEVEL:'
  !---------------------------------------------------------------------
  !** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
  !  wavelength of proton & proton mass (u) from 2002 physical constants.
  BZ = ZMU/16.85762920_wp
  write(u6,605) TITL,ZMU,BZ,MASS1,MASS2
  BvWN = One/BZ
  if (CHARGE /= 0) write(u6,624) CHARGE,CHARGE
  EJ = Zero
  EJ2 = Zero
  LRPT = 1
  ! Lower limit (RMIN) and increment (RH) of integration in (Angstroms).
  ! Upper limit of the reduced variable integration range automatically
  ! set at  YMAX= 1.0 , which corresponds to  RMAX= infinity !!.
  ! A hard wall boundary condition may be imposed at a smaller distance
  ! using an appropriate choice of the read-in level parameter IV (below)
  !!! The radial integration variable is  yp(r;Reff)  with   p= PRV
  ! EPS (cm-1) is the desired eigenvalue convergence criterion
  !---------------------------------------------------------------------
  !read(u5,*) RH,RMIN,PRV,ARV,EPS
  !---------------------------------------------------------------------
  ! NPP = no. of points in potential and wavefunction array.
  !!! First ... calculate new AS radial mesh YH implied but the given RH
  I = int(0.5e7_wp*(PRV/ARV)*RH)
  !.... give YH a rounded-off value (to 8 digits)
  YH = real(I,kind=wp)*1.0e-7_wp
  aRVp = ARV**PRV
  RRp = RMIN**PRV
  YMIN = (RRp-aRVp)/(RRp+aRVp)
  YMAX = One
  ! NPP = no. of points in potential and wavefunction array.
  NPP = int(((YMAX-YMIN)/YH+1.00001_wp))
  if (NDIMR < NPP) then
    write(u6,6604) NDIMR,YH,real(NPP,kind=wp)/real(NDIMR,kind=wp)
    NPP = NDIMR
  end if
  !... reset YMIN slightly to precisely span range
  YMIN = YMAX-(NPP-1)*YH
  YH2 = YH*YH
  BFCT = BZ*YH2
  YMINN = YMIN-YH
  write(u6,604) YMIN,YMAX,YH,PRV,ARV,RMIN,RH,ARV,NAME1,IMN1,NAME2,IMN2
  PINV = One/PRV
  FFAS = YH2*(PINV**2-One)/(Four*aRVp)**2
  do I=2,NPP-1
    YVB(I) = YMINN+I*YH
    RRp = (One+YVB(I))/(One-YVB(I))
    RR = RRp**PINV
    RVB(I) = ARV*RR
    RRM2(I) = One/RVB(I)**2
    RRM22(I) = RRM2(I)
    RRp = RRp*aRVp
    DRDY = RVB(I)*(RRp+aRVp)**2/(Two*pRV*RRp*aRVp)
    DRDY2(I) = DRDY**2
    SDRDY(I) = sqrt(DRDY)
    FAS(I) = FFAS*((RRp+aRVp)**2/RRp)**2
  end do
  ! OPIONALLY WRITE SOME VARIABLES IF DEBUGGING:
  !write(u6,*) 'DRDY=',DRDY
  !write(u6,*) 'RRp=',RRp
  !write(u6,*) 'aRVp=',aRVp
  !write(u6,*) 'pRV=',pRV
  !write(u6,*) 'RRp=',RRp
  !write(u6,*) 'aRVp=',aRVp
  !write(u6,*) 'DRDY2(1)=',DRDY2(1)
  !write(u6,*) 'DRDY2(2)=',DRDY2(2)
  !write(u6,*) 'RVB(1)=',RVB(1)
  !write(u6,*) 'RVB(2)=',RVB(2)
  YVB(1) = YMIN
  RVB(1) = RMIN
  RRM2(1) = RRM2(2)
  DRDY2(1) = DRDY2(2)
  SDRDY(1) = SDRDY(2)
  if (RMIN > Zero) then
    RRM2(1) = One/RMIN**2
    RRp = RMIN**PRV
    DRDY = RVB(1)*(RRp+aRVp)**2/(Two*PRV*RRp*aRVp)
    DRDY2(1) = DRDY**2
    SDRDY(1) = sqrt(DRDY)
  end if
  RRM22(1) = RRM2(1)
  YVB(NPP) = YMAX
  !... 'fake' RMAX value to ensure last 1/R**2 point is stable.
  RVB(NPP) = RVB(NPP-1)+RVB(NPP-1)-RVB(NPP-2)
  RRp = RVB(NPP)**PRV
  DRDY = RVB(NPP)*(RRp+aRVp)**2/(Two*PRV*RRp*aRVp)
  DRDY2(NPP) = DRDY**2
  SDRDY(NPP) = sqrt(DRDY)
  RRM2(NPP) = One/RVB(NPP)
  RRM22(NPP) = RRM2(NPP)
  !For debugging purposes, you can print the first 3 V(R) values:
  !do I=1,3
  !  write(u6,*) RVB(I)
  !end do

  !++ Begin reading appropriate parameters & preparing potential(s)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !+  Subroutine "PREPOT" prepares (and if desired, writes) the potential
  !+  array V(i) (cm-1)  at the NPP distances RVB(i) (Angst).
  ! NPP = no. of points in potential and wavefunction array.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !* If NTP > 0 :  define potential by interpolation over & extrapolation
  !        beyond the NTP read-in turning points using subroutine GENINT.
  !   If NTP <= 0 : generate a (fully analytic) potential in POTGEN.
  !* If LPPOT > 0 : at every |LPPOT|-th point, print potential and
  !      derivatives-by-differences. ***  If  LPPOT < 0  write potential
  !      at every |LPPOT|-th point to channel-8 in a compact format **
  !* OMEGA is the electronic contribution to the angular momentum such
  !  that the reduced centrifugal potential is:  (J*(J+1)-OMEGA**2)/R**2
  !* Set (OMEGA >= 99) if wish to use centrifugal factor for rotation
  !  in two dimensions:   (J**2 - 1/4)/R**2  .
  !* VLIM (cm-1) is the energy associated with the potential asymptote.
  !---------------------------------------------------------------------
  !read(u5,*) NTP,LPPOT,OMEGA,VLIM
  !---------------------------------------------------------------------
  !** For pointwise potentials, PREPOT uses subroutine GENINT to read
  !  points and conditions and interpolate (using subroutines NTRPSR,
  !  SPLINE & SPLINT) and extrapolate to get full potential.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** For a pointwise potential (NTP > 0), now read points & parameters
  !  controlling how the interpolation/extrapolation is to be done.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** NTP (read above) is number of turning points (XI,YI) to be read in.
  !** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
  !    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE <= 0)
  !    interpolate with cubic spline instead of local polynomials.
  !** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
  !   This may help if interpolation has trouble on steep repulsive wall.
  !** ILR specifies how to extrapolate beyond largest input distance XI(i)
  !  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
  !  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
  !  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
  !** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
  !  inverse-power terms beginning with  1/R**NCN}. *** If CNN /= 0 ,
  !  leading coefficient fixed at  CNN ; otherwise get it from points too.
  !* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
  !* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
  !* If ILR > 3 : successive higher power terms differ by factor  1/R
  !
  !** RFACT & EFACT are factors required to convert units of input turning
  !       points (XI,YI) to Angstroms & cm-1, respectively (often = 1.0)
  !** Turning points (XI,YI) must be ordered with increasing XI(I)
  !** Energy VSHIFT (cm-1) is added to the input potential points to
  !   make their absolute energy consistent with VLIM (often VSHIFT=Te).
  !---------------------------------------------------------------------
  !read(u5,*) NUSE,IR2,ILR,NCN,CNN
  !read(u5,*) RFACT,EFACT,VSHIFT
  !read(u5,*) (XI(I),YI(I),I=1,NTP)
  !---------------------------------------------------------------------
  ! NCN1 (returned by PREPOT) is the power of the asymptotically-
  ! dominant inverse-power long range potential term.
  ! VLIM1, VLIM1, V1, NCN1 and CNN1 are not defined yet, but are input parameters for
  ! PREPOT
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  write(u6,*) 'Exiting level'
  write(u6,*) 'Entering prepot'
  write(u6,*) ''
  call PREPOT(LRPT,NPP,IOMEG1,RVB,RRM2,VLIM1,V1,CNN1,NCN1,PPAR,QPAR,NSR,NLR,DSCM,REQ,RREF,PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
  !call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,V1,CNN1,NCN1)
  write(u6,*) 'Successfully made it through Prepot!'
  !OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
  !do I=1,3
  !  write(u6,*) 'V(',I,')=',V1(I)
  !end do
  !write(u6,*) 'V(                 20000)=',V1(20000)
  !write(u6,*) 'V(',NPP,')=',V1(NPP)
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** If (NTP <= 0) PREPOT uses subroutine POTGEN to generate a fully
  !  analytic potential defined by the following read-in parameters.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
  !  and for all cases except IPOTL=2, the potential asymptote energy is
  !  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
  !  potential minimum and  DSCM  the leading (quadratic) potential coeft.
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  !** IPOTL specifies the type of potential function to be generated.
  !** MPAR, NSR & NCMM are integers characterizing the chosen potential
  !** NVARB is number of (real*8) potential parameters read in.
  !** IBOB specifies whether (if > 0) or not (if <= 0) atomic mass
  !      dependent Born-Oppenheimer breakdown corrections will be included
  !- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  !** If IPOTL=1  generate an L.J.(MPAR,NCN) potential.
  !** If IPOTL=2  use Seto's modification of Surkus' GPEF expansion in
  !       z = [R**MPAR - Re**MPAR]/[a*R**MPAR + b*Re**MPAR] where
  !       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
  !       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
  !       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
  !       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
  !  * For Dunham case:  NCN=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
  !  * For SPF case:  NCN=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
  !  * For Ogilvie-Tipping:  NCN=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
  !  * NOTE that for Surkus MPAR < 0 case:  z(MPAR,a,b)= z(|MPAR|,-b,-a)
  !      Generate & return the  D_e  value implied by these coefficients.
  !** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
  !      with exponent factor "beta" defined as a power series of order
  !      (NVARB-1) in  y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR)
  !      with NVARB coefficients PARM(i).    [!! MPAR >= 1 !!]
  !    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
  !*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
  !      modified Morse function with  b= PARM(1)  and C= PARM(2).
  !*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
  !      Morse Oscillator" potential with exponent expansion in (R-Re)]
  ! ...  otherwise, set  MPAR >= 0
  !** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
  !      If MPAR > 0  exponent parameter defined in terms of a polynomial
  !           of order (NVARB-2) with the NVARB coefficients  PARM(j).
  !           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR),
  !           and long-range defined by NCN inverse-power terms
  !      If MPAR = 0  exponent polynomial variable is  2*y_{1}= y{O-T}
  !      If MPAR < 0  exponent polynomial variable is  y_{|MPAR|}
  !      If MPAR <= 0  exponent polynomial connected to limiting inverse-
  !           power potential exponent by exponential switching function
  !           with parameters  Asw= PARM(NVARB-1)  and  RSW= PARM(NVARB).
  !** If IPOTL=5  generate a Double-Exponential Long-Range (DELR)
  !       potential [JCP 119, 7398 (2003)] with additive long-range part
  !       defined by a sum of NCMM damped inverse-power terms, & exponent
  !       polynomial radial variable defined by parameter MPAR (=p)
  !** If IPOTL=6  generate generalized HFD({m_i},i=1,NCMM) potential.
  !       PARM(1-3) are the parameters defining the HFD damping function
  !       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
  !       PARM(4) the quadratic coefficient in the exponent, and
  !       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
  !              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
  !** If IPOTL=7  use Tiemann-type polynomial potential attached to an
  !     inverse-power long-range tail and an 1/R^{12} (or exponential)
  !     inner wall.
  !---------------------------------------------------------------------
  !read(u5,*) IPOTL,MPAR,NSR,NCMM,NVARB,IBOB,DSCM,REQ
  !if (IPOTL >= 4) read(u5,*) (MMLR(I),CMM(I),I=1,NCMM)
  !if ((IPOTL == 4) .or. (IPOT == 7)) read(u5,*) (MMLR(I),CMM(I),I=1,MPAR)
  !if (NVARB > 0) read(u5,*) (PARM(I),I=1,NVARB)
  !if (IBOB > 0) then
  !  read(u5,*) MN1R,MN2R,PAD,MAD,NU1,NU2,PNA,NT1,NT2
  !  if (PAD > 0) then
  !    IF (NU1 >= 0) read(u5,*) U1INF,(U1(I),I=0,NU1)
  !    IF (NU2 >= 0) read(u5,*) U2INF,(U2(I),I=0,NU2)
  !    IF (NT1 >= 0) read(u5,*) T1INF,(T1(I),I=0,NT1)
  !    IF (NT2 >= 0) read(u5,*) T2INF,(T2(I),I=0,NT2)
  !   else
  !    IF (NU1 >= 0) read(u5,*) U1INF,(U1(I),I=0,NU1),Aad1,Rad1
  !    IF (NU2 >= 0) read(u5,*) U2INF,(U2(I),I=0,NU2),Aad2,Rad2
  !    IF (NT1 >= 0) read(u5,*) T1INF,(T1(I),I=0,NT1),Ana1,Rna1
  !    IF (NT2 >= 0) read(u5,*) T2INF,(T2(I),I=0,NT2),Ana2,Rna2
  !  end if
  !end if
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  PW = Two
  if ((NCN1 > 0) .and. (NCN1 /= 2)) PW = Two*NCN1/(NCN1-Two)
  if (real(NCN1,kind=wp) < (Two*PRV+1.9999999_wp)) write(u6,629) (Two*PRV+Two),NCN1
  ! Convert potential in [cm-1] to form appropriate for SCHRQas
  V1BZ(1:NPP) = V1(1:NPP)*BFCT
  V1(1:NPP) = V1BZ(1:NPP)*DRDY2(1:NPP)+FAS(1:NPP)
  V2(1:NPP) = V1(1:NPP)
  RM2(1:NPP) = RRM2(1:NPP)*DRDY2(1:NPP)
  VLIM2 = VLIM1
  ! OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
  !write(u6,*) 'V(R) after converting into form for schrq:'
  !do I=1,3
  !  !write(u6,*) 'V(',I,')=',V1(I)
  !  !write(u6,*) 'FAS(',I,')=',FAS(I)
  !  write(u6,*) 'DRDY2(',I,')=',DRDY2(I)
  !  !write(u6,*) 'RRM2(',I,')=',RRM2(I)
  !end do
  !write(u6,*) 'V(                 20000)=',V1(20000)
  !write(u6,*) 'V(',NPP,')=',V1(NPP)
  if (NUMPOT <= 1) then
    write(u6,636)
    IOMEG2 = IOMEG1
  else
    !write(u6,635)
    write(u6,*) ''
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! For 2-potential Franck-Condon factor calculation, get the second
    ! potential in this second call to PREPOT (uses the same parameter
    ! reading sequence so exhaustively described immediately above).
    !call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,VLIM2,V2,CNN2,NCN2)
    call PREPOT(LRPT,NPP,IOMEG2,RVB,RRM22,VLIM2,V2,CNN2,NCN2,PPAR,QPAR,NSR,NLR,DSCM,REQ,RREF,PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    ! Convert potential (in (cm-1)) to form appropriate for SCHRQas
    V2BZ(1:NPP) = V2(1:NPP)*BFCT
    V2(1:NPP) = V2BZ(1:NPP)*DRDY2(1:NPP)+FAS(1:NPP)
    RM22(1:NPP) = RRM22(1:NPP)*DRDY2(1:NPP)
  end if

  !** NLEV1 is the no. of levels {v=IV(i), J=IJ(i)} of potential-1 which
  !   we wish to find.
  !* IF(NLEV1=0) calculate (and print?) potential, and then quit.
  !* If read-in value of  NLEV1 < 0 , and program (attempts to) find all
  !  vibrational levels of potential-1 up to  v = |NLEV1|.  [This case
  !  assumes  AUTO1 > 0.]
  !** If (AUTO1 > 0) read in only (v,J) quantum numbers of NLEV1 desired
  !  level & subroutine ALFas tries to locate them (normal preferred case).
  !   If (AUTO1 <= 0) also read in trial energy for each level.  In this
  !      case, the  NLEV <= 0  option does not work.
  !   If (AUTO1 <= 0) and vib. quant. No. IV < 0, seek level nearest to
  !      given trial energy but for whatever q. number shows up
  !** If(LCDC > 0) calculate centrifugal distortion constants for each
  !  level via the Tellinghuisen implementation of Hutson's method.
  !  If(LCDC >= 2) also write them compactly to channel-9.
  !** IF(LXPCT=0) calculate no expectation values or matrix elements.
  !* IF(LXPCT = -1) only calculate and write compactly to channel-7 the
  !  eigenvalues and level widths.
  !* IF(LXPCT= 1,2 or -2) calculate expectation values, or if |LXPCT| > 2
  !  the off-diagonal matrix elements, of powers of the distance
  !  coordinate or radial function defined by parameters IRFN & DREF.
  !* For   LXPCT > 0  write all these results to channel-6;  otherwise
  !                  supress most such printing to channel-6.
  !* For  |LXPCT| = 2  write eigenvalues and expectation values in
  !                       compact form on channel-7.
  !* For  |LXPCT| > 2  calculate matrix elements coupling each level
  !  to all (up to VIBMX) preceeding levels of the same potential (for
  !  NUMPOT <= 1), or to NLEV2 (see below) vib. levels of potential-2
  !  (for NUMPOT >= 2), and if (|LXPCT| > 3) write the overall
  !  off-diagonal matrix elements on channel-8.
  !* For  |LXPCT| > 4  also write to channel-7 the matrix elements of the
  !  individual powers of the chosen distance coordinate (or radial fx.)
  !* For  |LXPCT| > 5  WRITE(7,xx) only those matrix element components.
  !** IF(NJM > 0), for each vibrational level, calculate all rotational
  !  levels up to  J=NJM  or predissociation, whichever comes first.
  !  Note that  AUTO1 <= 0  forces  NJM= 0
  !** When (NJM > 0) increase J in increments of JDJR.
  !** IF(IWR /= 0) print error & warning descriptions
  !  IF(IWR >= 1) also print final eigenvalues & node count.
  !  IF(IWR >= 2) also show end-of-range wave function amplitudes
  !  IF(IWR >= 3) print also intermediate trial eigenvalues, etc.
  !** IF(LPRWF > 0) print wave function every LPRWF-th  point.
  !** IF(LPRWF < 0) compactly write to channel-7 every |LPRWF|-th
  !  wave function value.  **  A lead "card" identifies the level, gives
  !  the position of 1-st point and radial mesh, & states No. of  points
  !=====================================================================
  !** INNOD1 specified wave fx. initiation at RMIN.  Normal case of
  !  INNOD1 > 0  gives initiation with wave fx. node @ RMIN.
  !  INNOD1 <= 0  give initiation with  zero slope @ RMIN.  This determines
  !    symmetric eigenfunctions for rare special case when input potential
  !    is half of a precisely symmetric potential with mid-point at RMIN.
  !---------------------------------------------------------------------
  ! ... comment this out if use first version of following READ statement
  INNOD1 = 1
  INNOD2 = INNOD1
  !---------------------------------------------------------------------
  !read(u5,*) NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF,INNOD1
  !---------------------------------------------------------------------
  !read(u5,*) NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
  !---------------------------------------------------------------------
  ! SINNER specifies whether wave function matching occurs at outermost
  ! (SINNER <= 0) or innermost well turning point, to facilitate finding
  ! inner vs. outer wells of a double well potential; Normally controlled
  ! automatically,

  if (LPRWF < 0) write(10,605) TITL

  SINNER = 0
  INNER = SINNER
  if (INNOD1 > 0) write(u6,686) 1
  if (INNOD1 <= 0) write(u6,688) 1
  if (JDJR <= 0) JDJR = 1
  write(u6,612) EPS
  NLEV = NLEV1
  if (NLEV1 <= 0) NLEV = 1
  SOMEG1 = IOMEG1**2
  if (IOMEG1 >= 99) then
    write(u6,609)
  else
    if (IOMEG1 >= 0) write(u6,608) 1,IOMEG1,SOMEG1
    if (IOMEG1 < 0) write(u6,6085) 1,IOMEG1,-IOMEG1
  end if
  VMAX1 = 0
  !** Read the vibrational & rotational quantum numbers IV(i) & IJ(i) [and
  !  if AUTO1 <= 0 also trial energy GV(I)] of the NLEV levels to be found
  !** For  IV(i)  values < -10,  SCHRQ  imposes a hard wall boundary
  !  condition (i.e., a node) at mesh point # |-IV(i)| .
  !---------------------------------------------------------------------
  !if (AUTO1 > 0) read(u5,*) (IV(I),IJ(I),I=1,NLEV)
  !if (AUTO1 <= 0) read(u5,*) (IV(I),IJ(I),GV(I),I=1,NLEV)
  !---------------------------------------------------------------------
  ! IJ(i) for each level i should be read from the input file but
  ! otherwise initialize them explicitly:
  IJ(1:NLEV) = 0
  if (NLEV1 > 0) then
    if (AUTO1 > 0) write(u6,607) NLEV,(IV(I),IJ(I),I=1,NLEV)
    if (AUTO1 <= 0) then
      write(u6,6607) NLEV,(IV(I),IJ(I),GV(I),I=1,NLEV)
      do I=1,NLEV1
        if (IV(I) < VIBMX) ZK1(IV(I),0) = GV(I)
      end do
    end if
    do I=1,NLEV
      if (IV(I) < VIBMX) VMAX1 = max(VMAX1,IV(I))
    end do
    JREF = 0
  else
    if (NLEV1 < (1-VIBMX)) NLEV1 = 1-VIBMX
    VMAX1 = -NLEV1
    NLEV = VMAX1+1
    write(u6,625) IJ(1),NLEV
    JREF = IJ(1)
    IJ(1:NLEV) = JREF
    do I=1,NLEV
      IV(I) = I-1
    end do
  end if
  if (NJM > IJ(1)) write(u6,638) JDJR,NJM
  !if (LCDC > 1) write(9,901) TITL
  if (LXPCT == -1) write(7,723) TITL
  !** MORDR is the highest power of the radial function (or distance
  !  coordinate whose expectation values or matrix elements are to be
  !  calculated.  Program currently dimensioned for (MORDR <= 10).  To
  !  calculate only F-C Factors (when LXPCT>2), set  MORDR = -1.
  !** IRFN & DREF specify the definition of the radial function or
  !  distance coordinate  RFN(R), powers of which are averaged over in
  !  expectation value or matrix element calculations.
  !* If(IRFN <= -10) utilize the USER-CHOSEN and CODED radial function
  !                 generated in Lines #500-504 (below)
  !* If(IRFN = -4)  the function is a power series in  R  premultiplying a
  !          first derivative operator acting on the wavefx of Potential-2
  !* If(IRFN = -3)  the function is the inverse power   1/R**3
  !* If(IRFN = -2)  the function is the inverse power   1/R**2
  !* If(IRFN = -1)  the function is the Dunham coordinate  X=(R-DREF)/DREF
  !* If(IRFN =  0)  the function  RFN(R)  is the distance  R  itself.
  !* If(IRFN = 1-9)  use the Surkus-type variable
  !                 X=(R^p - DREF^p)/(R^p + DREF^p)  where  p= IRFN
  !* For  IRFN = -1 or 1-9,  if  DREF > 0  the read-in DREF value is the
  !  reference length used to define the distance coordinate, while
  !  if  DREF <= 0  determine the value of this reference length by
  !  requiring that the expectation value  X**1  of the distance
  !  coordinate for the first level considered be identically zero.
  !* IF(IRFN > 10) define  RFN(R)   by reading in, interpolating over (and
  !  extrapolating beyond) IRFN read-in values of some known radial
  !  (transition moment) function, whose asymptotic value is DREF.  Do
  !  this in using the same read statements and GENINT subroutine calls
  !  used for generating a numerical potential.
  if ((LXPCT /= 0) .and. (LXPCT /= -1)) then
    !-------------------------------------------------------------------
    !read(u5,*) MORDR,IRFN,DREF
    !-------------------------------------------------------------------
    if (MORDR > MORDRMX) MORDR = MORDRMX
    if (abs(LXPCT) == 2) write(7,724) TITL,MORDR
    !if ((abs(LXPCT) == 4) .or. (abs(LXPCT) == 5)) write(8,824) TITL
    if (abs(LXPCT) >= 5) write(7,725) TITL,MORDR
    if (abs(IRFN) >= 10) then
      MORDR = 1
      DM(0) = Zero
      DM(1) = One
    else
      if (MORDR >= 0) then
        ! Overall calculated matrix elements are for a power series in the
        ! radial function  RFN(i)  (specified by IRFN & DREF), so must read
        ! coefficients  DM(J)  of this power series.
        !---------------------------------------------------------------
        !read(u5,*) (DM(J),J=0,MORDR)
        !---------------------------------------------------------------
      else
        RFN(1:NPP) = One
        if (MORDR < 0) write(u6,617)
      end if
    end if
    ! Define radial function (distance coordinate) operator  RFN(R)  for
    ! expectation values or matrix elements.
    ! First ... for matrix elements of an operator consisting of a power
    ! series in  R  premultiplying the radial derivative of the wavefx.
    if (IRFN == -4) write(u6,650) MORDR
    if (MORDR > 0) then
      ! If  RFN(R)  is the distance itself ...
      if (IRFN == 0) then
        write(u6,614)
        RFN(1:NPP) = RVB(1:NPP)
      end if
      if ((IRFN == -2) .or. (IRFN == -3)) then
        if ((IRFN == 0) .or. (IRFN == -2) .or. (IRFN == -3)) DREF = Zero
        ! If  RFN(R)  is   1/(distance)**|IRFN|  ....
        J = -IRFN
        write(u6,616)-IRFN
        RFN(1:NPP) = One/RVB(1:NPP)**J
      end if
      ! Any other user-defined matrix element argument radial function
      ! may be introduced to the code here, and invoked by:  IRFN= -4
      ! Note that the existing  RVB(i)  array is the radial distances  R .
      if (IRFN <= -10) then
        ! Illustrative user-defined analysis RFN(R) function
        !write(u6,*) 'Print description of function introduced'
        !write(u6,*) 'Use Freedman Pade DMF for CO'
        !do I=1,NPP
        !  !------------------------------------------------------------
        !  RFN(I) = {calculate users chosen radial function}
        !  ! Freedman's DMF for CO  ------------------------------------
        !  RFN(I) = {calculate users chosen radial function}
        !  data coeff_new /-24.6005858_wp,-109.5939637_wp,-524.8233323_wp,4.5194090_wp,19.7954955_wp,6.6011985_wp,19.7206690_wp/
        !  dm = -0.122706_wp*(One+coeff(1)*x+coeff(2)*x*x+coeff(3)*x**3)/(One+coeff(4)*x+coeff(5)*x*x+coeff(6)*x**3+coeff(7)*x**6)
        !  !------------------------------------------------------------
        !  XX = RFN(I)/1.128322714_wp-One
        !  RFN(I)= -0.122706_wp*(One+XX*(-24.6005858_wp+XX*(-109.5939637_wp+XX*(-524.8233323_wp))))/ &
        !          (One+XX*(4.5194090_wp+XX*(19.7954955_wp+XX*(6.6011985_wp+19.7206690_wp*XX**3))))
        !end do
      end if
      if ((IRFN == -1) .or. ((IRFN >= 1) .and. (IRFN <= 9))) then
        ! If  RFN(R)  is the Dunham or Surkus-type distance coordinate
        if (IRFN == -1) write(u6,615)
        if ((IRFN >= 1) .and. (IRFN <= 9)) write(u6,611) IRFN,IRFN,IRFN,IRFN
        if (DREF > Zero) then
          DREFP = DREF**IRFN
          write(u6,613) DREF
          do I=1,NPP
            XX = YMINN+I*YH
            if (IRFN == -1) RFN(I) = (RVB(I)-DREF)/DREF
            if (IRFN >= 1) RFN(I) = (RVB(I)**IRFN-DREFP)/(RVB(I)**IRFN+DREFP)
          end do
        else
          write(u6,610)
        end if
      end if

      !QQQQQQQQQQQ new ... the following not fully tested ...

      ! If  RFN(R)  is defined by interpolating over read-in points, use
      ! potential generating routine to do interpolation/extrapolation.
      if (IRFN >= 10) then
        MORDR = 1
        DM(0) = Zero
        DM(1) = One
        write(u6,603)
        !** If the expectation value/matrix element radial function argument to
        !   be defined by interpolating/extrapolating over read-in points, then
        !   read input analogous to that for a pointwise potential, and then call
        !   interpolation/extrapolation routine GENINT (from PREPOT package)
        !* NRFN is the number of points [XIF(i),YIF(i)] to be read in
        !* RFLIM  is the limiting asymptotic value imposed on the extrapolation
        !* Interpolate with NUSEF-point piecewise polynomials (or splines for
        !    NUSEF <= 0), which are extrapolated to the asymptote as specified by
        !    parameters ILRF, NCNF & CNNF (see read #20).
        !* RFACTF - factor converts read-in distances XIF(i) to angstroms
        !* MFACTF - factor converts read-in moment values YIF(i) to debye.
        !---------------------------------------------------------------
        !read(u5,*) NRFN,RFLIM
        !read(u5,*) NUSEF,ILRF,NCNF,CNNF
        !read(u5,*) RFACTF,MFACTF
        !read(u5,*) (XIF(I),YIF(I),I=1,NRFN)
        ! If you uncomment the above, you better also uncomment the
        ! initialization to 0 below.
        !---------------------------------------------------------------
        MFACTF = 0
        RFACTF = 0
        write(u6,810) NRFN,RFLIM
        if (NUSEF > 0) write(u6,812) NUSEF,NRFN
        if (NUSEF <= 0) write(u6,814) NRFN
        if ((ILRF > 1) .and. (abs(CNNF) > Zero)) write(u6,816) CNNF,NCNF
        write(u6,818) RFACTF,MFACTF
        NROW = (NRFN+2)/3
        do J=1,NROW
          write(u6,820) (XIF(I),YIF(I),I=J,NRFN,NROW)
        end do
        XIF(1:NRFN) = XIF(1:NRFN)*RFACTF
        YIF(1:NRFN) = YIF(1:NRFN)*MFACTF
        IR2F = 0
        call GENINT(LRPT,NPP,RVB,RFN,NUSEF,IR2F,NRFN,XIF,YIF,RFLIM,ILRF,NCNF,CNNF)
      end if
    end if
    if ((MORDR >= 0) .and. (abs(IRFN) <= 9)) write(u6,602) (DM(J),J=0,MORDR)
  end if
  !** For matrix element calculation, couple each level of potential-1 to
  !  up to (see below) NLEV2 other vibrational levels, subject to
  !  rotational selection rules:  DELTA(J)= J2DL to J2DU with increment
  !  J2DD (e.g., -1,+1,+2 for P- and R-branches).
  !** If (AUTO2 > 0) read in only (v,J) quantum numbers of desired levels
  !  and trust subroutine ALFas to locate them (normal preferred case).
  !   If (AUTO2 <= 0) also read in a trial pure vib energy for each level.
  !* For the one-potential case (NUMPOT <= 1), automatically truncate to
  !  avoid redundancy and consider only emission into these NLEV2 levels.
  !* Trial level energies are generated internally.
  !**  IV2(i) are the vibrational quantum numbers of the Potential-2
  !  levels for which matrix elements are desired.
  !**  ZK(IV2(i),0) are the associated pure vibrational trial energies
  !  (which are only read in if AUTO2 <= 0!)
  !=====================================================================
  !** INNOD2 specified wave fx. initiation at RMIN.  Normal case of
  !  INNOD2 > 0  gives initiation with wave fx. node @ RMIN.
  !  INNOD2 <= 0  give initiation with  zero slope @ RMIN.  This determines
  !    symmetric eigenfunctions for rare special case when input potential
  !    is half of a precisely symmetric potential with mid-point at RMIN.
  !=====================================================================
  if (abs(LXPCT) >= 3) then
    !-------------------------------------------------------------------
    !read(u5,*) NLEV2,AUTO2,J2DL,J2DU,J2DD,INNOD2
    !-------------------------------------------------------------------
    !read(u5,*) NLEV2,AUTO2,J2DL,J2DU,J2DD
    !-------------------------------------------------------------------
    if (NLEV2 > VIBMX) NLEV2 = VIBMX
    if (NLEV2 <= 0) then
      write(u6,644) NLEV2
      !stop
      call abend()
    end if
    !-------------------------------------------------------------------
    !if (AUTO2 > 0) read(u5,*) (IV2(I),I=1,NLEV2)
    if (AUTO2 <= 0) then
      !read(u5,*) (IV2(I),ZK2(I,1),I=1,NLEV2)
      !-----------------------------------------------------------------
      !** Give potential-2 trial energy the correct vibrational label
      do I=1,NLEV2
        ZK2(IV2(I),0) = ZK2(I,1)
      end do
    end if
    if (NUMPOT > 1) then
      if (INNOD2 > 0) write(u6,686) 2
      if (INNOD2 <= 0) write(u6,688) 2
    end if
    VMAX2 = 0
    do ILEV2=1,NLEV2
      VMAX2 = max(VMAX2,IV2(ILEV2))
    end do
    if (MORDR < 0) DM(1) = One
    SOMEG2 = IOMEG2**2
    if (J2DD == 0) J2DD = 1
    if (AUTO2 > 0) write(u6,634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),I=1,NLEV2)
    if (AUTO2 <= 0) write(u6,6634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),ZK2(IV2(I),0),I=1,NLEV2)
    if (NUMPOT >= 2) then
      if (IOMEG2 >= 99) then
        write(u6,609)
      else
        if (IOMEG2 >= 0) write(u6,608) 2,IOMEG2,SOMEG2
        if (IOMEG2 < 0) write(u6,6085) 2,IOMEG2,-IOMEG2
      end if
    end if
    write(u6,632)
  end if

  if (AUTO1 > 0) then
    ! If using automatic search for desired levels, subroutine ALFas gets
    ! eigenvalues ZK1(v,0) for desired vibrational levels of Potential-1,
    ! centrifugally-distorted to J=JREF.
    EJREF = JREF*(JREF+1)*YH**2
    ! Replace  [J(J+1)] by  [J(J+1) + |IOMEG1|]  for Li2(A) and like cases.
    if (IOMEG1 < 0) EJREF = EJREF-real(IOMEG1,kind=wp)*YH**2
    VBZ(1:NPP) = V1BZ(1:NPP)+EJREF*RRM2(1:NPP)
    VJ(1:NPP) = V1(1:NPP)+EJREF*RM2(1:NPP)
    ! OPTIONALLY WRITE SOME VALUES WHEN DEBUGGING:
    !write(u6,*) 'EJREF=',EJREF
    !do I=1,3
    !  write(u6,*) 'VJ=',VJ(I)
    !  write(u6,*) 'V1',V1(I)
    !  write(u6,*) 'RM2=',RM2(I)
    !end do
    if ((NLEV1 == 1) .and. (IV(1) > 998)) then
      !** Option to search for very highest level (within 0.0001 cm-1 of Disoc)
      EO = VLIM1-1.0e-4_wp
      KV = IV(1)
      write(u6,*) ''
      write(u6,*) 'Exiting level (1)'
      write(u6,*) 'Entering schrqas (1)'
      write(u6,*) ''
      call SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
      IV(1) = KV
      if (KV >= 0) then
        write(u6,622) IJ(1),KV,VLIM1-EO
        GV(KV) = EO
        VMAX1 = KV
      else
        write(u6,626) J,1.0e-3_wp
        cycle outer
      end if
    else
      VMAX = VMAX1
      AFLAG = JREF
      if ((abs(LXPCT) > 2) .and. (NUMPOT == 1)) VMAX = max(VMAX1,VMAX2)
      write(u6,*) ''
      write(u6,*) 'Exiting level'
      write(u6,*) 'Entering alfas'
      write(u6,*) ''
      !write(u6,*) 'NDP=',NPP
      !write(u6,*) 'YMIN=',YMIN
      !write(u6,*) 'YH=',YH
      !write(u6,*) 'NCN=',NCN1
      !do I=1,3
      !  write(u6,*) 'V=',VJ(I)
      !  write(u6,*) 'SWF=',WF1(I)
      !  write(u6,*) 'GV=',GV(I)
      !  write(u6,*) 'INNR=',INNR1(I)
      !end do
      !write(u6,*) 'VLIM=',VLIM1
      !write(u6,*) 'KVMAX=',VMAX
      !write(u6,*) 'AFLAG=',AFLAG
      !write(u6,*) 'ZMU=',ZMU
      !write(u6,*) 'EPS=',EPS
      !write(u6,*) 'BFCT=',BFCT
      !write(u6,*) 'INNODE=',INNOD1
      !write(u6,*) 'IWR=',IWR
      !write(u6,*) ''
      !deallocate(RVB)
      call ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,VMAX,AFLAG,EPS,GV,BFCT,INNOD1,INNR1,IWR)
      VMAX1 = VMAX
    end if
    !** Get band constants for v=0-VMAX1 for generating trial eigenvalues
    WARN = 0
    do ILEV1=0,VMAX
      KV = ILEV1
      EO = GV(KV)
      INNER = INNR1(KV)
      write(u6,*) ''
      write(u6,*) 'Exiting level (2)'
      write(u6,*) 'Entering schrqas (2)'
      write(u6,*) ''
      write(u6,*) 'Getting band constants for v=',KV
      call SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,WARN,LPRWF)
      write(u6,*) ''
      write(u6,*) 'Exiting level'
      write(u6,*) 'Entering cdjoelas'
      write(u6,*) ''
      ! OPTIONALLY WRITE THE INPUT PARAMETERS FOR DEBUGGING:
      !write(u6,*) 'EO=',EO
      !write(u6,*) 'NBEG=',NBEG
      !write(u6,*) 'NEND=',NEND
      !write(u6,*) 'BvWN=',BvWN
      !write(u6,*) 'YH=',YH
      !write(u6,*) 'WARN=',WARN
      !write(u6,*) 'VJ(1)=',VJ(1)
      !write(u6,*) 'WF1(1)=',WF1(1)
      !write(u6,*) 'RM2(1)=',RM2(1)
      !write(u6,*) 'RCNST(1)=',RCNST(1)
      ! For a-state test case, there's a memory error before CALL CDJOELas for v=2
      !call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
      if (NLEV1 < 0) then
        IV(ILEV1+1) = KV
        IJ(ILEV1+1) = JREF
      end if
      ZK1(ILEV1,0) = GV(ILEV1)
      ZK1(ILEV1,1:7) = RCNST(1:7)
    end do
  end if
  if (abs(LXPCT) > 2) then
    if (AUTO2 > 0) then
      ! If using automatic location for levels of potential-2 (AUTO2 > 0)
      ! for matrix element calculation, also need Potential-2 band constants
      ! (rotational energy derivatives) ... again, calculate them at J=JREF
      if (NUMPOT > 1) then
        AFLAG = JREF
        VBZ(1:NPP) = V2BZ(1:NPP)+EJREF*RRM22(1:NPP)
        VJ(1:NPP) = V2(1:NPP)+EJREF*RM22(1:NPP)
        call ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,VMAX2,AFLAG,EPS,GV,BFCT,INNOD2,INNR2,IWR)
      end if
    end if
    do ILEV2=1,NLEV2
      if (NUMPOT == 1) then
        ! For matrix elements within a single potl., copy above band constants
        ZK2(IV2(ILEV2),0:7) = ZK1(IV2(ILEV2),0:7)
      else
        ! ... otherwise, generate them (as above) with SCHRQ & CDJOEL
        KV = IV2(ILEV2)
        if (AUTO2 > 0) EO = GV(KV)
        if (AUTO2 <= 0) EO = ZK2(KV,0)
        INNER = INNR2(KV)
        write(u6,*) ''
        write(u6,*) 'Exiting level (3)'
        write(u6,*) 'Entering schrqas (3)'
        write(u6,*) ''
        call SCHRQas(KV,JREF,EO,GAMA,PMAX2,VLIM2,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD2,INNER,WARN,LPRWF)
        call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
        ZK2(IV2(ILEV2),0) = EO
        ZK2(IV2(ILEV2),1:7) = RCNST(1:7)
      end if
    end do
  end if
  WARN = 1
  EJREF = EJREF/YH**2
  if (NLEV1 <= 0) NLEV = VMAX1+1

  !===== Begin Actual Potential-1 Eigenvalue Calculation Loop Here =====
  ! Loop to compute eigenvalues ... etc. for NLEV levels of Potential-1
  do ILEV1=1,NLEV
    KV = IV(ILEV1)
    if (KV < 0) exit
    NJMM = max(NJM,IJ(ILEV1))
    JROT = IJ(ILEV1)-JDJR
    IQT = 0
    JCT = 0
    ! If NJM > IJ(ILEV1) loop over range of rotational levels too
    jloop: do JLEV=IJ(ILEV1),NJMM,JDJR
      JROT = JROT+JDJR
      EJ = JROT*(JROT+1)-SOMEG1
      if (IOMEG1 >= 99) EJ = JROT*JROT-Quart
      ! If   IOMEG < 0   centrifugal term is  [J(J+1) + |IOMEG|]
      if (IOMEG1 < 0) EJ = JROT*(JROT+1)-real(IOMEG1,kind=wp)
      ! If appropriate (AUTO1>0) use ALFas results to generate trial eigenvalue
      if (AUTO1 > 0) then
        EO = ZK1(KV,0)
        DEJ = EJ-EJREF
        EJP = One
        do M=1,7
          EJP = EJP*DEJ
          EO = EO+EJP*ZK1(KV,M)
        end do
      else
        !... otherwise - use read-in trial energy
        if (IV(ILEV1) < VIBMX) EO = ZK1(IV(ILEV1),0)
        if (IV(ILEV1) >= VIBMX) EO = GV(ILEV1)
      end if
      Skip = .false.
      if ((AUTO1 <= 0) .and. (abs(ZK1(IV(ILEV1),0)) <= Zero)) then
        call SCATTLEN(JROT,SL,VLIM1,V1,WF1,BFCT,YMIN,YH,NPP,CNN1,NCN1,IWR,LPRWF)
        if (NUMPOT == 1) cycle outer
      else
        ! ... or if JLEV > IJ(ILEV1) ... use local Beff to estimate next level
        if (JLEV > IJ(ILEV1)) then
          BEFF = sum(WF1(NBEG:NEND)**2*RM2(NBEG:NEND))*YH*BvWN
          EO = ESLJ(JCT)+(2*JLEV+1-JDJR)*JDJR*BEFF
        end if
        ! Now add centrifugal term to get effective (radial) potential
        EJ = EJ*YH**2
        VBZ(1:NPP) = V1BZ(1:NPP)+EJ*RRM2(1:NPP)
        VJ(1:NPP) = V1(1:NPP)+EJ*RM2(1:NPP)
        ! Set wall outer boundary condition, if specified by input IV(ILEV1)
        if (KV < -10) then
          WF1(-IV(ILEV1)) = Zero
          WF1(-IV(ILEV1)-1) = -One
        end if
        KVIN = KV
        if (AUTO1 > 0) INNER = INNR1(KV)
        if (SINNER /= 0) INNER = SINNER
        ! Call SCHRQ to find Potential-1 eigenvalue EO and eigenfn. WF1(i)
        write(u6,*) ''
        write(u6,*) 'Exiting level (4)'
        write(u6,*) 'Entering schrqas (4)'
        write(u6,*) ''
        ! The next CALL SCHRQas for a-state failes during v=0, but EO's from when
        ! Entering schrqas (2) were already accurate, so the next step is not necessary.
        do
          call SCHRQas(KV,JROT,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
          if (KV < 0) then
            ! SCHRQ  error condition is  (KV < 0) .
            if (NJM > IJ(ILEV1)) then
              ! ... in automatic search for ever-higher J levels
              if (IQT <= 0) then
                ! ... try one more time with E(trial) slightly below barrier maximum
                IQT = 1
                EO = PMAX1-0.1_wp
                cycle
              else
                KV = KVIN
                cycle jloop
              end if
            end if
            Skip = .true.
            exit
          end if
          !write(u6,*) 'Hello! (for debugging) 1'
          if ((KV /= KVIN) .and. (AUTO1 > 0)) then
            !** If got wrong vib level, do a brute force ALFas calculation to find it.
            KV = KVIN
            AFLAG = JROT
            call ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,KV,AFLAG,EPS,GV,BFCT,INNOD1,INNR1,IWR)
            if (KV == KVIN) then
              EO = GV(KVIN)
            else
              write(u6,618) KVIN,JROT,KV
              KV = KVIN
              cycle jloop
            end if
          else
            exit
          end if
        end do
        if (.not. Skip) then
          !write(u6,*) 'Hello! (for debugging) 2'
          if (KV /= IV(ILEV1)) IV(ILEV1) = KV
          ! If desired, calculate rotational & centrifugal distortion constants
          if (LCDC > 0) then
            if ((IOMEG1 > 0) .and. (JROT == 0)) then
              ! Calculate 'true' rotational constants for rotationless IOMEG>0 case
              write(u6,*) ''
              write(u6,*) 'Exiting level (5)'
              write(u6,*) 'Entering schrqas (5)'
              write(u6,*) ''
              call SCHRQas(KV,0,EO,GAMA,PMAX1,VLIM1,V1,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V1,WF1,RM2,RCNST)
            else
              ! Calculate rotational constants for actual (v,J) level of interest.
              call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
            end if
            if (abs(EO) > One) then
              write(u6,606) KV,JROT,EO,(RCNST(M),M=1,7)
            else
              write(u6,6065) KV,JROT,EO,(RCNST(M),M=1,7)
            end if
            if (LCDC > 1) then
              if (abs(EO) > One) then
                !write(9,902) KV,JROT,EO,(RCNST(M),M=1,7)
              else
                !write(9,904) KV,JROT,EO,(RCNST(M),M=1,7)
              end if
            end if
          end if
          !write(u6,*) 'Hello! (for debugging) 3'
          if (LXPCT == -1) write(7,703) KV,JROT,EO,GAMA
          if (((LXPCT == 1) .or. (abs(LXPCT) == 2)) .or. &
              ((abs(LXPCT) > 2) .and. ((IRFN == -1) .or. (IRFN >= 1) .and. (IRFN >= 9)) .and. (DREF <= Zero))) then
            ! Calculate various expectation values in LEVXPC
            call LEVXPC(KV,JROT,EO,GAMA,NPP,WF1,RFN,VBZ,VLIM1,YH,DREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
            if ((LXPCT > 0) .and. (MORDR > 0)) write(u6,632)
          end if
        end if
      end if
      !write(u6,*) 'Hello! (for debugging) 4'
      if ((.not. Skip) .and. (abs(LXPCT) > 2) .and. (NLEV2 > 0)) then
        !** If desired, now calculate off-diagonal matrix elements, either
        !  between levels of different potentials, IF(NUMPOT >= 2), or between
        !  levels of a single potential, for (NUMPOT <= 1).
        !** First prepare centrifugally distorted potential, trial energy, etc.,
        !  and calculate second wave function and matrix element(s)
        do ILEV2=1,NLEV2
          ! For case of a single potential, avoid redundancy by considering
          ! only emission
          if ((NUMPOT <= 1) .and. (IV2(ILEV2) > KV)) cycle
          ! Loop over J2's allowed by given selection rule.
          iloop: do IJD=J2DL,J2DU,J2DD
            KV2 = IV2(ILEV2)
            KVIN = KV2
            JROT2 = JROT+IJD
            if (JROT2 < 0) cycle
            if ((NUMPOT <= 1) .and. (IV2(ILEV2) == KV) .and. (JROT2 > JROT)) cycle
            EJ2 = JROT2*(JROT2+1)-SOMEG2
            if (IOMEG2 >= 99) EJ2 = JROT2**2-Quart
            !... allow for weird Li2(A) and Li2(c) potential cases
            if (IOMEG2 < 0) EJ2 = JROT2*(JROT2+1)-real(IOMEG2,kind=wp)
            EO2 = ZK2(KV2,0)
            DEJ = EJ2-EJREF
            EJP = One
            ! Use calculated state-2 CDC's to predict trial eigenvalue
            do M=1,7
              EJP = EJP*DEJ
              EO2 = EO2+EJP*ZK2(KV2,M)
            end do
            ! Now ... update to appropriate centrifugally distorted potential
            EJ2 = EJ2*YH*YH
            VBZ(1:NPP) = V2BZ(1:NPP)+EJ2*RRM22(1:NPP)
            VJ(1:NPP) = V2(1:NPP)+EJ2*RM22(1:NPP)
            INNER = INNR2(KV2)
            if (SINNER /= 0) INNER = SINNER
            ICOR = 0
            write(u6,*) ''
            write(u6,*) 'Exiting level (6)'
            write(u6,*) 'Entering schrqas (6)'
            write(u6,*) ''
            do
              call SCHRQas(KV2,JROT2,EO2,GAMA,PMAX2,VLIM2,VJ,WF2,BFCT,EPS,YMIN,YH,NPP,NBEG2,NEND2,INNOD2,INNER,IWR,LPRWF)
              if (KV2 /= KVIN) then
                if (KV2 >= 0) then
                  ! Using CDC's to estimate trial eigenvalue failed:
                  ICOR = ICOR+1
                  if (ICOR <= 2) then
                    ! ... first correction attempt ... use semiclassical dv/dE to improve
                    GB = -One
                    GI = -One
                    ! ... hey RJ!!  shouldn't you update this using SCECOR ?
                    WV = Zero
                    XX = EO2*BFCT
                    do I=NBEG2,NEND2
                      GBB = GB
                      GB = GI
                      GI = XX-VJ(I)
                      if ((GBB > Zero) .and. (GI > Zero)) WV = WV+One/sqrt(GB)
                    end do
                    WV = Two*Pi/(BFCT*WV)
                    EO2 = EO2+WV*(KVIN-KV2)
                    cycle
                  end if
                  write(u6,633) IV2(ILEV2),JROT2,KV2
                end if
                ! ... if that fails, do a brute force ALFas calculation to find it.
                KV2 = KVIN
                AFLAG = JROT2
                call ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,KV2,AFLAG,EPS,GV,BFCT,INNOD2,INNR2,IWR)
                if (KV2 == KVIN) then
                  EO2 = GV(KV2)
                  INNER = INNR2(KV2)
                else
                  write(u6,618) KVIN,JROT,KV2
                  cycle iloop
                end if
              end if
            end do
            if (NBEG > NBEG2) NBEG2 = NBEG
            if (NEND < NEND2) NEND2 = NEND
            !if ((NUMPOT <= 1) .and. (EO2 > (EO+EPS))) exit
            call MATXEL(KV,JROT,IOMEG1,EO,KV2,JROT2,IOMEG2,IRFN,EO2,NBEG2,NEND2,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
          end do iloop
          ! End of Potential-2 rotational selection level loop
        end do
        !++++ End of Potential-2 vibrational level matrix element loop +
      end if

      JCT = JCT+1
      ! ... check to avoid array overflow
      if (JCT > VIBMX) then
        write(u6,637) VIBMX
        !stop
        call abend()
      end if
      JWR(JCT) = JROT
      ESLJ(JCT) = EO
    end do jloop
    !write(u6,*) 'Hello! (for debugging) 5'
    !++ End of Potential-1 loop over NJM-specified J-sublevels
    if (NJM > IJ(ILEV1)) then
      ! Print rotational sublevels generated for vibrational level  ILEV
      NROW = (JCT+4)/5
      write(u6,627) KV
      do J=1,NROW
        write(u6,628) (JWR(I),ESLJ(I),I=J,JCT,NROW)
      end do
      write(u6,641)
    end if
    ESOLN(ILEV1) = ESLJ(1)
  end do
  !++ End of loop over the NLEV Potential-1 input levels
  write(u6,*) ''
  write(u6,*) 'SUMMARY (ALL ENERGIES IN CM-1):'
  if (NLEV1 < 0) then
    NROW = (NLEV+3)/4
    write(u6,623) NLEV,IJ(1)
    do J=1,NROW
      write(u6,630) (IV(I),ESOLN(I),I=J,NLEV,NROW)
    end do
    if ((NLEV > 1) .and. (IJ(1) == 0) .and. (NCN1 > 0) .and. (ESOLN(NLEV) < VLIM1)) then
      ! In (NLEV1 < 0) option, estimate vD using the N-D theory result that:
      ! (vD - v) {is proportional to} (binding energy)**((NCN-2)/(2*NCN))
      VDMV = One/(((VLIM1-ESOLN(NLEV-1))/(VLIM1-ESOLN(NLEV)))**(One/PW)-One)
      !call ADD_INFO('LEVEL_CHKSUM',[CHKSUM],1,2)
      call ADD_INFO('Vibrational energy:',ESOLN,NLEV,2)
      ! IF YOU WANT TO ONLY VERIFY THE ENERGY OF THE HIGHEST LEVEL:
      !call ADD_INFO('Last vibrational energy:',ESOLN(NLEV),1,2)
      ! IF YOU WANT TO ONLY VERIFY THE NUMBER OF LEVELS:
      ! All values must be reals, so if you want to check an integer, e.g. the
      ! number of iterations, then you must convert this to a real:
      !call ADD_INFO('Numer of vibrational levels:',NLEV,1,2)
      ! Use empirical N-D Expression to predict number and (if there are
      ! any) energies of missing levels
      VD = IV(NLEV)+VDMV
      IVD = int(VD)
      if (IVD >= VIBMX) IVD = VIBMX-1
      IVS = IV(NLEV)+1
      write(u6,620) NCN1,VD
      if ((IVD >= IVS) .and. (real(IV(NLEV),kind=wp)/VD > 0.9_wp)) then
        NFP = NLEV+1
        do I=IVS,IVD
          NLEV = NLEV+1
          IV(NLEV) = IV(NLEV-1)+1
          ESOLN(NLEV) = VLIM1-(VLIM1-ESOLN(NLEV-1))*(One-One/VDMV)**PW
          VDMV = VDMV-One
        end do
        NLP = NLEV-NFP+1
        NROW = (NLP+3)/4
        write(u6,621) NLP
        do J=1,NROW
          write(u6,630) (IV(I),ESOLN(I),I=NFP+J-1,NLEV,NROW)
        end do
      end if
    end if
  end if
  if ((NJM <= 0) .and. (NLEV1 >= 0)) then
    NROW = (NLEV+2)/3
    write(u6,619) NLEV
    do J=1,NROW
      write(u6,631) (IV(I),IJ(I),ESOLN(I),I=J,NLEV,NROW)
    end do
  end if
  write(u6,601)
  ! Comment out the following line to read input file again if you want to find
  ! the levels for a different potential. Currently LEVEL_RDINP doesn't
  ! allow us to read the input file a second time though.
  if (.true.) exit outer
end do outer

!-----------------------------------------------------------------------
call mma_deallocate(RVB)
call mma_deallocate(YVB)
call mma_deallocate(DRDY2)
call mma_deallocate(FAS)
call mma_deallocate(SDRDY)
call mma_deallocate(VBZ)

call mma_deallocate(V1)
call mma_deallocate(V2)
call mma_deallocate(VJ)
call mma_deallocate(V1BZ)
call mma_deallocate(V2BZ)
call mma_deallocate(WF1)
call mma_deallocate(WF2)

call mma_deallocate(RFN)
call mma_deallocate(RRM2)
call mma_deallocate(RM2)
call mma_deallocate(RRM22)
call mma_deallocate(RM22)

RC = 0

601 format(1x,79('=')////)
602 format(' Coefficients of expansion for radial matrix element/expectation value argument:'/(5X,5(1PD14.6)))
603 format(/' Expectation value/matrix element arguments are powers of a radial function'/5x, &
           'defined by interpolating over read-in points'//' Transition moment function:'/1x,9('==='))
6604 format(/' !!! NOTE:  array dimension limit   NDIMR=',i5/'    prevents preliminary mesh   YH=',f9.6, &
            '  from spanning range [YMIN,YMAX]'/'    so increase mesh by factor',f7.4/)
604 format(' Integrate from   Ymin=',f11.7,'   to   Ymax=',f7.4,'  with mesh  YH=',f9.7/5x, &
           'based on radial variable  yp(r;a)= (r^p - a^p)/(r^p + a^p)'/5x,'for   p=',f6.3,'    and   a=',F9.6/ &
           ' Range corresponds to   Rmin=',f7.3,' [Angst]    to    Rmax= infinity (!!)'/5x,'and   RH=',F10.7,' [Angst]   at   R=', &
           f11.8//' Potential #1 for ',A2,'(',I3,')-',A2,'(',I3,')'/1x,32('='))
605 format(/A78/40('=='):/' Generate   ZMU=',F15.11,'(u)','   &   BZ=',1PD16.9,'((1/cm-1)(1/Ang**2))'/10x,'from atomic masses:', &
           0Pf16.11,'  & ',F16.11,'(u)')
606 format(' E(v=',i3,', J=',i3,')=',f10.3,'   Bv=',F11.7,'  -Dv=',1PD12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,'   Mv=',D12.4, &
           '   Nv=',D12.4,'   Ov=',D12.4)
6065 format(' E(v=',i3,', J=',i3,')=',f11.7,'  Bv=',1PD11.4,'  -Dv=',D12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,'   Mv=',D12.4, &
           '   Nv=',D12.4,'   Ov=',D12.4)
607 format(/' Solve for the',i4,' vibration-rotation levels of Potential-1:'/'   (v,J) =',6('  (',i3,',',i3,')':)/(10x,6('  (',i3, &
           ',',i3,')':)))
6607 format(/' Solve for',i4,' vibration-rotation levels of Potential-1 using Trial energies:'/(3x,3('   E(',I3,',',I3,')=', &
            F11.2:)))
608 format(/' State-',I1,' electronic angular momentum  OMEGA=',I2/9x,'yields centrifugal potential  [J*(J+1) -',F5.2,']/r**2')
6085 format(/' State-',I1,' electronic angular momentum  OMEGA=',I2/9x,'yields centrifugal potential  [J*(J+1) +',I2,']/r**2')
609 format('  Use centrifugal potential for rotation in two dimensions:   (J**2 - 1/4)/r**2')
610 format(5X,'where DREF defined by requiring  <X**1> = 0  for first level considered')
611 format(/' Matrix element argument expansion vble is   X = ((r^',i1,' - DREF^',i1,')/(r^',i1,' + DREF^',i1,'))')
612 format(/' Eigenvalue convergence criterion is   EPS=',1PD8.1,'(cm-1)'/ &
           ' Airy function at 3-rd turning point is quasibound outer boundary condition')
613 format(5X,'where reference length is held fixed at   DREF =',F13.10,'(Angstroms)')
614 format(/' Matrix element arguments are powers of the distance  r (in Angstroms)')
615 format(/' Matrix element argument expansion variable is:    X = (r - DREF)/DREF')
616 format(/' Matrix element arguments are powers of the squared inverse distance  X = 1/r**',i1)
617 format(/' Matrix element argument is fixed as a constant = 1')
618 format(' *** PROBLEM *** Searching for  v=',i3,' , J=',i3,'  ALFas only found to  v=',i3)
619 format(/' Find the',i4,' vibration-rotation levels:'/3('     v   J      E(v)   ')/3(2x,7('---')))
620 format(/' An  n=',I2,'  N-D theory extrapolation from last 2 levels implies   vD =',F8.3)
621 format(5X,'with the',I4,' missing level(s) predicted to be:'/4('     v     E(v)   ')/4(4x,7('--')))
622 format(/' Search for highest bound  J=',i3,'  level finds  E(v=',i3,') = VLIM -',1PD12.5/)
623 format(/' Find',I4,' Potential-1 vibrational levels with  J=',i3/4('     v     E(v)   ')/4(4x,7('--')))
624 format(4x,'Since the molecule is an ion with charge',SP,I3/6x, &
           "use Watson's charge-adjusted reduced mass   mu = M1*M2/[M1 + M2 - (",i2,')*me]')
625 format(' For  J=',i3,', try to find the first',i4,' vibrational levels of Potential-1')
626 format(/' *** FAIL to find highest bound J=',i3,'  level from trial   E = VLIM -',1PD11.4)
627 format(/' For vibrational level  v =',I3,'   of Potential-1'/1X,5('  J',6X,'E',7X)/1X,5(7('--'),2X))
628 format((1X,5(I3,F11.3,2X)))
629 format(/' *** Note that Radial variable power \alpha optimal for NLR=',f5.2,' > NCN=',i2)
630 format((4(I6,F12.4:)))
631 format((3(I6,I4,F13.5:)))
632 format(1X,79('-'))
633 format(' **** Caution: Search for   v=',I3,'   J=',i3,'  on potential-2 actually found   v=',I3)
634 format(/' Using the rotational selection rule:  delta(J)=',i3,' to',i2,' with increment',i2/ &
           '   calculate matrix elements for coupling to the',I4,' vibrational levels of'/'   Potential-2:   v =',14I4:/(21x,14i4:))
6634 format(/' Using the rotational selection rule:  delta(J)=',i3,' to',i2,' with increment',i2/ &
            '   calculate matrix elements for coupling to the',I4,' vibrational levels of'/'   Potential-2 using trial energies:', &
            2('   E(',I3,')=',F9.2:)/4('   E(',I3,')=',F9.2:))
!635 format(/' Get matrix elements between levels of Potential-1 (above) & Potential-2 (below)'/1X,39('--')/' For Potential #2:'/ &
!           1x,17('='))
636 format(/' Calculate properties of the single potential described above')
637 format(/' *** Array Dimension OVERFLOW ***   (Number of J sublevels) > VIBMX=',i4)
638 format('   and automatically increment  J  in steps of',i3,' to a maximum value of',i4)
641 format(1X,39('++'))
644 format(/' *** Input data ERROR *** matrix element calculation needs  NLEV2=',i3,' > 0')
650 format(/' Matrix element argument is radial first derivative operator premultiplied by'/5x,'a power series in  r  of order',i3)
686 format(' Potential-',i1,' uses inner boundary condition of  zero value  at  RMIN')
688 format(' Potential-',i1,' uses symmetric-well inner boundary condition  of zero slope at RMIN')
!703 format(1X,I4,I5,F13.4,G13.5)
703 format(1X,I4,I5,1PD20.11,D13.5)
723 format(/A78/1x,'Output values of:  v, J, E (Level Width)')
724 format(//A78//'   v   J    E(v,J)     Width       <KE>',6x,'<M(r)>  &  <XI**k>  for k=1 to',i3/2x,38('=='))
725 format(//A78//"   v'  J'",'  v"  J"     FREQ',"    <v',J'| XI**k",' |v",J">  for  k=0  to  MORDR=',i2/2x,37('=='))
810 format(' Transition moment function defined by interpolating over',I4,' read-in points'/5x, &
           'and approaching the asymptotic value',f12.6)
812 format(' Perform',I3,'-point piecewise polynomial interpolation over',I5,' input points')
814 format(' Perform cubic spline interpolation over the',I5,' input points')
816 format('- Beyond read-in points extrapolate to limiting asymptoticbehaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
818 format(' Scale input points:  (distance)*',1PD16.9,'     (moment)*',D16.9/4x,'to get required units  [Angstroms & debye]'/ &
           3('      R(i)         Y(i)  ')/3(3X,11('--')))
820 format((3(F12.6,F13.6)))
!824 format(//A78/30('==')/" Note that (v',J') &",' (v",J") strictly label the upper and lower levels, resp.,'/6x, &
!           'and  E(lower)=E"'/' but  E(2)-E(1)  is:  (energy of State-2 level) - (energy of State-1 level)'//12x,'Band'/ &
!           ' dJ(J")',4x,7hv'   v",'  E(lower)  E(2)-E(1)  A(Einstein)   F-C Factor  ',13h<v'j'|M|v"j"> /1x,3('--'),( &
!           '   -------'),'  --------',3x,4('--'),3x,11('-'),3x,11('-'),3x,11('-') )
!811 format(//A78/30('==')//12x,"   v'","  J'",'    v"','  J"   position    E(upper)    E(lower)',16h   <v'j'|M|v"j">/1x,68('-') )
!901 format(//A78/1x,62('==')/'   v    J',7x,'E',10x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
!902 format(I4,I5,f25.15,f14.10,6(1PD15.7))
!904 format(I4,I5,f25.15,1PD14.7,6(D15.7))

end subroutine LEVEL
