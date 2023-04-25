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

use stdalloc, only: mma_allocate, mma_deallocate
use LEVEL_COMMON

implicit none
integer, intent(OUT) :: RC
!** Dimensions for  potential arrays  and  vib. level arrays.
integer VIBMX, MORDRMX, RORDR, NTPMX
parameter(VIBMX=400,RORDR=7,MORDRMX=20,NTPMX=1600)
integer NDIMR
!parameter (NDIMR=200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive. fmax-stack-var-size=2^20
!parameter (NDIMR=131074) ! MMA_ALLOCATE won't allow PARAMETERs
!real*8 PRV, ARV, RVB(NDIMR), YVB(NDIMR), DRDY2(NDIMR), FAS(NDIMR), SDRDY(NDIMR), VBZ(NDIMR), aRVp
real*8 PRV, ARV, aRVp
!real*8, ALLOCATABLE :: RVB(:), YVB(:), DRDY2(:), FAS(:), SDRDY(:), VBZ(:)
common/BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
integer I, J, M, III, IJD, ILEV1, ILEV2, IOMEG1, IOMEG2, INNOD1, INNOD2, INNER, SINNER, IQT, IWR, IRFN, IVD, IVS, IAN1, IAN2, &
        IMN1, IMN2, GEL1, GEL2, GNS1, GNS2, JDJR, JCT, J2DL, J2DU, J2DD, JROT, JROT2, JLEV, JREF, ICOR, CHARGE, KV, KV2, KVIN, &
        LCDC, LPRWF, LRPT, LXPCT, MORDR, NUSEF, ILRF, IR2F, NUMPOT, NBEG, NBEG2, NEND, NEND2, NPP, NCN1, NCN2, NCNF, NLEV, NLEV1, &
        NLEV2, NJM, NJMM, NFP, NLP, NRFN, NROW, WARN, VMAX, VMAX1, VMAX2, AFLAG, AUTO1, AUTO2, IV(VIBMX), IJ(VIBMX), IV2(VIBMX), &
        JWR(VIBMX), INNR1(0:VIBMX), INNR2(0:VIBMX), NTP, LPPOT, IPOTL, PPAR, QPAR, NSR, NLR, IBOB, NCMM, IVSR, IDSTT, MMLR(3)
!real*8 ZK1(0:VIBMX,0:RORDR), ZK2(0:VIBMX,0:RORDR), RCNST(RORDR), V1(NDIMR), V2(NDIMR), VJ(NDIMR), V1BZ(NDIMR), V2BZ(NDIMR), &
!       WF1(NDIMR), WF2(NDIMR), CMM(3), PARM(4)
real*8 ZK1(0:VIBMX,0:RORDR), ZK2(0:VIBMX,0:RORDR), RCNST(RORDR), CMM(3), PARM(4)
real*8, allocatable :: V1(:), V2(:), VJ(:), V1BZ(:), V2BZ(:), WF1(:), WF2(:)
!real*8 RFN(NDIMR), RRM2(NDIMR), RM2(NDIMR), RRM22(NDIMR), RM22(NDIMR), GV(0:VIBMX), ESOLN(VIBMX), ESLJ(VIBMX), XIF(NTPMX), &
!       YIF(NTPMX), ABUND1, ABUND2, MASS1, MASS2, DM(0:MORDRMX)
real*8 GV(0:VIBMX), ESOLN(VIBMX), ESLJ(VIBMX), XIF(NTPMX), YIF(NTPMX), ABUND1, ABUND2, MASS1, MASS2, DM(0:MORDRMX)
real*8, allocatable :: RFN(:), RRM2(:), RM2(:), RRM22(:), RM22(:)
real*8 BZ, BvWN, BFCT, BEFF, DEJ, EPS, EO, EO2, EJ, EJ2, EJP, EJREF, GAMA, MEL, PMAX1, PMAX2, PW, RH, RMIN, RR, RRp, PINV, DRDY, &
       YH, YH2, YMIN, YMINN, YMAX, DREF, DREFP, CNN1, CNN2, RFLIM, CNNF, RFACTF, MFACTF, SOMEG1, SOMEG2, VLIM1, VLIM2, VD, VDMV, &
       XX, ZMU, GI, GB, GBB, WV, FFAS, SL, DSCM, REQ, RREF, RHOAB
character*78 TITL
character*2 NAME1, NAME2
data MEL/5.4857990945d-4/,YMAX/1.d+00/
!** Default (Q-branch) defining J-increments for matrix element calcn.
data J2DL,J2DU,J2DD/0,0,1/

NDIMR = 131074
call mma_allocate(RVB,NDIMR,label='RVB')
call mma_allocate(YVB,NDIMR,label='YVB')
call mma_allocate(DRDY2,NDIMR,label='DRDY2')
call mma_allocate(FAS,NDIMR,label='FAS')
call mma_allocate(SDRDY,NDIMR,label='SDRDY')
call mma_allocate(VBZ,NDIMR,label='VBZ')

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
PINV = 1.d0
NLEV2 = -1
AUTO2 = 0
VMAX2 = 0
IOMEG2 = 0
SOMEG2 = 0
CNN2 = 0
PMAX2 = 0
ABUND2 = 0
MASS2 = 0
NCN2 = 0
NEND2 = 0
NBEG2 = 0
GNS2 = 0
GEL2 = 0
INNOD2 = 0
EJREF = 0
do I=1,VIBMX
  IV2(I) = 0
end do
do I=1,RORDR
  RCNST(I) = 0
end do
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
!----------------------------------------------------------------------
2 continue
call LEVEL_RDINP(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS,NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM, &
                 REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)
! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING:
!write(6,*) IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS
!write(6,*) NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB
!write(6,*) DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1
! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING (ANOTHER WAY):
!write(6,*) AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
!write(6,*) 'level.f has the following after CALL LEVEL_RDINP:'
!write(6,*) 'IAN1 = ',IAN1
!write(6,*) 'IMN1 = ',IMN1
!write(6,*) 'IAN2 = ',IAN2
!write(6,*) 'IMN2 = ',IMN2
!write(6,*) 'CHARGE = ',CHARGE
!write(6,*) 'NUMPOT = ',NUMPOT
!write(6,*) 'RH = ',RH
!write(6,*) 'RMIN = ',RMIN
!write(6,*) 'PRV = ',PRV
!write(6,*) 'ARV = ',ARV
!write(6,*) 'EPS = ',EPS
!write(6,*) 'NTP = ',NTP
!write(6,*) 'LPPOT = ',LPPOT
!write(6,*) 'IOMEG1 = ',IOMEG1
!write(6,*) 'VLIM = ',VLIM
!write(6,*) 'IPOTL = ',IPOTL
!write(6,*) 'PPAR = ',PPAR
!write(6,*) 'QPAR = ',QPAR
!write(6,*) 'NSR = ',NSR
!write(6,*) 'NLR = ',NLR
!write(6,*) 'IBOB = ',IBOB
!write(6,*) 'DSCM = ',DSCM
!write(6,*) 'REQ = ',REQ
!write(6,*) 'RREF = ',RREF
!write(6,*) 'NCMM = ',NCMM
!write(6,*) 'IVSR = ',IVSR
!write(6,*) 'IDSTT = ',IDSTT
!write(6,*) 'RHOAB = ',RHOAB
!write(6,*) 'MMLR = ',MMLR
!write(6,*) 'CMM = ',CMM
!write(6,*) 'PARM = ',PARM
!write(6,*) 'NLEV1 = ',NLEV1
!write(6,*) 'AUTO1 = ',AUTO1
!write(6,*) 'LCDC = ',LCDC
!write(6,*) 'LXPCT = ',LXPCT
!write(6,*) 'NJM = ',NJM
!write(6,*) 'JDJR = ',JDJR
!write(6,*) 'IWF = ',IWF
!write(6,*) 'LPRWF = ',LPRWF
!read(5,*,END=999)
!2 continue
!read(5,*,END=999) IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT
!----------------------------------------------------------------------
! Subroutine MASSES returns the names of the atoms NAMEi,ground
! electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
! mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
if ((IAN1 > 0) .and. (IAN1 <= 109)) then
  call MASSES(IAN1,IMN1,NAME1,GEL1,GNS1,MASS1,ABUND1)
else
  ! If particle-i is not a normal atomic isotope, read a 2-character
  ! name (enclosed between '', as in 'mu') and its actual mass.
  !----------------------------------------------------------------------
  !read(5,*) NAME1,MASS1
  !----------------------------------------------------------------------
end if
if ((IAN2 > 0) .and. (IAN2 <= 109)) then
  call MASSES(IAN2,IMN2,NAME2,GEL2,GNS2,MASS2,ABUND2)
else
!----------------------------------------------------------------------
!read(5,*) NAME2,MASS2
!----------------------------------------------------------------------
end if
ZMU = MASS1*MASS2/(MASS1+MASS2-CHARGE*MEL)
!=======================================================================
! TITL is a title or output header of up to 78 characters, read on a
!   single line enclosed between single quotes: e.g.  'title of problem'
!=======================================================================
!read(5,*) TITL
TITL = 'Beginning execution of LEVEL:'
!----------------------------------------------------------------------
!** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
!  wavelength of proton & proton mass (u) from 2002 physical constants.
BZ = ZMU/16.85762920d0
write(6,605) TITL,ZMU,BZ,MASS1,MASS2
BvWN = 1.d0/BZ
if (CHARGE /= 0) write(6,624) CHARGE,CHARGE
EJ = 0.d0
EJ2 = 0.d0
LRPT = 1
! Lower limit (RMIN) and increment (RH) of integration in (Angstroms).
! Upper limit of the reduced variable integration range automatically
! set at  YMAX= 1.0 , which corresponds to  RMAX= infinity !!.
! A hard wall boundary condition may be imposed at a smaller distance
! using an appropriate choice of the read-in level parameter IV (below)
!!! The radial integration variable is  yp(r;Reff)  with   p= PRV
! EPS (cm-1) is the desired eigenvalue convergence criterion
!---------------------------------------------------------------------
!read(5,*) RH,RMIN,PRV,ARV,EPS
!---------------------------------------------------------------------
! NPP = no. of points in potential and wavefunction array.
!!! First ... calculate new AS radial mesh YH implied but the given RH
I = int(0.5d7*(PRV/ARV)*RH)
!.... give YH a rounded-off value (to 8 digits)
YH = dble(I)*1.d-07
aRVp = ARV**PRV
RRp = RMIN**PRV
YMIN = (RRp-aRVp)/(RRp+aRVp)
YMAX = 1.d0
! NPP = no. of points in potential and wavefunction array.
NPP = int(((YMAX-YMIN)/YH+1.00001))
if (NDIMR < NPP) then
  !write(6,6604) NDIMR,YH,dfloat(NPP)/dfloat(NDIMR)
  write(6,6604) NDIMR,YH,dble(NPP)/dble(NDIMR)
  NPP = NDIMR
end if
!... reset YMIN slightly to precisely span range
YMIN = YMAX-(NPP-1)*YH
YH2 = YH*YH
BFCT = BZ*YH2
YMINN = YMIN-YH
write(6,604) YMIN,YMAX,YH,PRV,ARV,RMIN,RH,ARV,NAME1,IMN1,NAME2,IMN2
PINV = 1.d0/PRV
FFAS = YH2*(PINV**2-1.d0)/(4.d0*aRVp)**2
do I=2,NPP-1
  YVB(I) = YMINN+I*YH
  RRp = (1.d0+YVB(I))/(1.d0-YVB(I))
  RR = RRp**PINV
  RVB(I) = ARV*RR
  RRM2(I) = 1.d0/RVB(I)**2
  RRM22(I) = RRM2(I)
  RRp = RRp*aRVp
  DRDY = RVB(I)*(RRp+aRVp)**2/(2.d0*pRV*RRp*aRVp)
  DRDY2(I) = DRDY**2
  SDRDY(I) = dsqrt(DRDY)
  FAS(I) = FFAS*((RRp+aRVp)**2/RRp)**2
end do
! OPIONALLY WRITE SOME VARIABLES IF DEBUGGING:
!write(6,*) 'DRDY=',DRDY
!write(6,*) 'RRp=',RRp
!write(6,*) 'aRVp=',aRVp
!write(6,*) 'pRV=',pRV
!write(6,*) 'RRp=',RRp
!write(6,*) 'aRVp=',aRVp
!write(6,*) 'DRDY2(1)=',DRDY2(1)
!write(6,*) 'DRDY2(2)=',DRDY2(2)
!write(6,*) 'RVB(1)=',RVB(1)
!write(6,*) 'RVB(2)=',RVB(2)
YVB(1) = YMIN
RVB(1) = RMIN
RRM2(1) = RRM2(2)
DRDY2(1) = DRDY2(2)
SDRDY(1) = SDRDY(2)
if (RMIN > 0.d0) then
  RRM2(1) = 1.d0/RMIN**2
  RRp = RMIN**PRV
  DRDY = RVB(1)*(RRp+aRVp)**2/(2.d0*PRV*RRp*aRVp)
  DRDY2(1) = DRDY**2
  SDRDY(1) = dsqrt(DRDY)
end if
RRM22(1) = RRM2(1)
YVB(NPP) = YMAX
!... 'fake' RMAX value to ensure last 1/R**2 point is stable.
RVB(NPP) = RVB(NPP-1)+RVB(NPP-1)-RVB(NPP-2)
RRp = RVB(NPP)**PRV
DRDY = RVB(NPP)*(RRp+aRVp)**2/(2.d0*PRV*RRp*aRVp)
DRDY2(NPP) = DRDY**2
SDRDY(NPP) = dsqrt(DRDY)
RRM2(NPP) = 1.d0/RVB(NPP)
RRM22(NPP) = RRM2(NPP)
!For debugging purposes, you can print the first 3 V(R) values:
!do I=1,3
!  write(6,*) RVB(I)
!end do

!++ Begin reading appropriate parameters & preparing potential(s)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+  Subroutine "PREPOT" prepares (and if desired, writes) the potential
!+  array V(i) (cm-1)  at the NPP distances RVB(i) (Angst).
! NPP = no. of points in potential and wavefunction array.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!-----------------------------------------------------------------------
!++   READ(5,*) NTP, LPPOT, OMEGA, VLIM
!----------------------------------------------------------------------
!** For pointwise potentials, PREPOT uses subroutine GENINT to read
!  points and conditions and interpolate (using subroutines NTRPSR,
!  SPLINE & SPLINT) and extrapolate to get full potential.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** For a pointwise potential (NTP > 0), now read points & parameters
!  controlling how the interpolation/extrapolation is to be done.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!       points (XI,YI) to Angstroms & cm-1, respectively (often = 1.d0)
!** Turning points (XI,YI) must be ordered with increasing XI(I)
!** Energy VSHIFT (cm-1) is added to the input potential points to
!   make their absolute energy consistent with VLIM (often VSHIFT=Te).
!-----------------------------------------------------------------------
!read(5,*) NUSE,IR2,ILR,NCN,CNN
!read(5,*) RFACT,EFACT,VSHIFT
!read(5,*) (XI(I),YI(I),I=1,NTP)
!-----------------------------------------------------------------------
! NCN1 (returned by PREPOT) is the power of the asymptotically-
! dominant inverse-power long range potential term.
! VLIM1, VLIM1, V1, NCN1 and CNN1 are not defined yet, but are input parameters for
! PREPOT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
write(6,*) 'Exiting level.f'
write(6,*) 'Entering prepot.f'
write(6,*) ''
call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,V1,CNN1,NCN1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,MMLR, &
            CMM,NCMM,IVSR,IDSTT,RHOAB)
!call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,V1,CNN1,NCN1)
write(6,*) 'Successfully made it through Prepot.f!'
!OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
!do I=1,3
!  write(6,*) 'V(',I,')=',V1(I)
!end do
!write(6,*) 'V(                 20000)=',V1(20000)
!write(6,*) 'V(',NPP,')=',V1(NPP)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** If (NTP <= 0) PREPOT uses subroutine POTGEN to generate a fully
!  analytic potential defined by the following read-in parameters.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
!  and for all cases except IPOTL=2, the potential asymptote energy is
!  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
!  potential minimum and  DSCM  the leading (quadratic) potential coeft.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!----------------------------------------------------------------------
!read(5,*) IPOTL,MPAR,NSR,NCMM,NVARB,IBOB,DSCM,REQ
!if (IPOTL >= 4) read(5,*) (MMLR(I),CMM(I),I=1,NCMM)
!if ((IPOTL == 4) .or. (IPOT == 7)) read(5,*) (MMLR(I),CMM(I),I=1,MPAR)
!if (NVARB > 0) read(5,*) (PARM(I),I=1,NVARB)
!if (IBOB > 0) then
!  read(5,*) MN1R,MN2R,PAD,MAD,NU1,NU2,PNA,NT1,NT2
!  if (PAD > 0) then
!    IF (NU1 >= 0) read(5,*) U1INF,(U1(I),I=0,NU1)
!    IF (NU2 >= 0) read(5,*) U2INF,(U2(I),I=0,NU2)
!    IF (NT1 >= 0) read(5,*) T1INF,(T1(I),I=0,NT1)
!    IF (NT2 >= 0) read(5,*) T2INF,(T2(I),I=0,NT2)
!   else
!    IF (NU1 >= 0) read(5,*) U1INF,(U1(I),I=0,NU1),Aad1,Rad1
!    IF (NU2 >= 0) read(5,*) U2INF,(U2(I),I=0,NU2),Aad2,Rad2
!    IF (NT1 >= 0) read(5,*) T1INF,(T1(I),I=0,NT1),Ana1,Rna1
!    IF (NT2 >= 0) read(5,*) T2INF,(T2(I),I=0,NT2),Ana2,Rna2
!  end if
!end if
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
PW = 2.d0
if ((NCN1 > 0) .and. (NCN1 /= 2)) PW = 2.d0*NCN1/(NCN1-2.d0)
!if (dfloat(NCN1) < (2.d0*PRV+1.9999999d0)) then
if (dble(NCN1) < (2.d0*PRV+1.9999999d0)) then
  write(6,629) (2.d0*PRV+2.),NCN1
  629 format(/' *** Note that Radial variable power \alpha optimal for NLR=',f5.2,' > NCN=',i2)
end if
! Convert potential in [cm-1] to form appropriate for SCHRQas
do I=1,NPP
  V1BZ(I) = V1(I)*BFCT
  V1(I) = V1BZ(I)*DRDY2(I)+FAS(I)
  V2(I) = V1(I)
  RM2(I) = RRM2(I)*DRDY2(I)
end do
VLIM2 = VLIM1
! OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
!write(6,*) 'V(R) after converting into form for schrq.f:'
!do I=1,3
!  !write(6,*) 'V(',I,')=',V1(I)
!  !write(6,*) 'FAS(',I,')=',FAS(I)
!  write(6,*) 'DRDY2(',I,')=',DRDY2(I)
!  !write(6,*) 'RRM2(',I,')=',RRM2(I)
!end do
!write(6,*) 'V(                 20000)=',V1(20000)
!write(6,*) 'V(',NPP,')=',V1(NPP)
if (NUMPOT <= 1) then
  write(6,636)
  IOMEG2 = IOMEG1
else
  !write(6,635)
  write(6,*) ''
  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! For 2-potential Franck-Condon factor calculation, get the second
  ! potential in this second call to PREPOT (uses the same parameter
  ! reading sequence so exhaustively described immediately above).
  !call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,VLIM2,V2,CNN2,NCN2)
  call PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,VLIM2,V2,CNN2,NCN2,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM, &
              MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)

  !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  ! Convert potential (in (cm-1)) to form appropriate for SCHRQas
  do I=1,NPP
    V2BZ(I) = V2(I)*BFCT
    V2(I) = V2BZ(I)*DRDY2(I)+FAS(I)
    RM22(I) = RRM22(I)*DRDY2(I)
  end do
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
!=======================================================================
!** INNOD1 specified wave fx. initiation at RMIN.  Normal case of
!  INNOD1 > 0  gives initiation with wave fx. node @ RMIN.
!  INNOD1 <= 0  give initiation with  zero slope @ RMIN.  This determines
!    symmetric eigenfunctions for rare special case when input potential
!    is half of a precisely symmetric potential with mid-point at RMIN.
!-----------------------------------------------------------------------
! ... comment this out if use first version of following READ statement
INNOD1 = 1
INNOD2 = INNOD1
!-----------------------------------------------------------------------
!read(5,*) NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF,INNOD1
!-----------------------------------------------------------------------
!read(5,*) NLEV1,AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
!-----------------------------------------------------------------------
! SINNER specifies whether wave function matching occurs at outermost
! (SINNER <= 0) or innermost well turning point, to facilitate finding
! inner vs. outer wells of a double well potential; Normally controlled
! automatically,

if (LPRWF < 0) write(10,605) TITL

SINNER = 0
INNER = SINNER
if (INNOD1 > 0) write(6,686) 1
if (INNOD1 <= 0) write(6,688) 1
if (JDJR <= 0) JDJR = 1
write(6,612) EPS
NLEV = NLEV1
if (NLEV1 <= 0) NLEV = 1
SOMEG1 = IOMEG1**2
if (IOMEG1 >= 99) then
  write(6,609)
else
  if (IOMEG1 >= 0) write(6,608) 1,IOMEG1,SOMEG1
  if (IOMEG1 < 0) write(6,6085) 1,IOMEG1,-IOMEG1
end if
VMAX1 = 0
!** Read the vibrational & rotational quantum numbers IV(i) & IJ(i) [and
!  if AUTO1 <= 0 also trial energy GV(I)] of the NLEV levels to be found
!** For  IV(i)  values < -10,  SCHRQ  imposes a hard wall boundary
!  condition (i.e., a node) at mesh point # |-IV(i)| .
!-----------------------------------------------------------------------
!if (AUTO1 > 0) read(5,*) (IV(I),IJ(I),I=1,NLEV)
!if (AUTO1 <= 0) read(5,*) (IV(I),IJ(I),GV(I),I=1,NLEV)
!-----------------------------------------------------------------------
! IJ(i) for each level i should be read from the input file but
! otherwise initialize them explicitly:
do I=1,NLEV
  IJ(I) = 0
end do
if (NLEV1 > 0) then
  if (AUTO1 > 0) write(6,607) NLEV,(IV(I),IJ(I),I=1,NLEV)
  if (AUTO1 <= 0) then
    write(6,6607) NLEV,(IV(I),IJ(I),GV(I),I=1,NLEV)
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
  write(6,625) IJ(1),NLEV
  JREF = IJ(1)
  do I=1,NLEV
    IV(I) = I-1
    IJ(I) = JREF
  end do
end if
if (NJM > IJ(1)) write(6,638) JDJR,NJM
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
  !---------------------------------------------------------------------
  !read(5,*) MORDR,IRFN,DREF
  !---------------------------------------------------------------------
  if (MORDR > MORDRMX) MORDR = MORDRMX
  if (iabs(LXPCT) == 2) write(7,724) TITL,MORDR
  !if ((iabs(LXPCT) == 4) .or. (iabs(LXPCT) == 5)) write(8,824) TITL
  if (iabs(LXPCT) >= 5) write(7,725) TITL,MORDR
  if (iabs(IRFN) >= 10) then
    MORDR = 1
    DM(0) = 0.d0
    DM(1) = 1.d0
  else
    if (MORDR >= 0) then
      ! Overall calculated matrix elements are for a power series in the
      ! radial function  RFN(i)  (specified by IRFN & DREF), so must read
      ! coefficients  DM(J)  of this power series.
      !-------------------------------------------------------------------
      !read(5,*) (DM(J),J=0,MORDR)
      !-------------------------------------------------------------------
    else
      do I=1,NPP
        RFN(I) = 1.d0
      end do
      if (MORDR < 0) write(6,617)
    end if
  end if
  ! Define radial function (distance coordinate) operator  RFN(R)  for
  ! expectation values or matrix elements.
  ! First ... for matrix elements of an operator consisting of a power
  ! series in  R  premultiplying the radial derivative of the wavefx.
  if (IRFN == -4) write(6,650) MORDR
  if (MORDR > 0) then
    ! If  RFN(R)  is the distance itself ...
    if (IRFN == 0) then
      write(6,614)
      do I=1,NPP
        RFN(I) = RVB(I)
      end do
    end if
    if ((IRFN == -2) .or. (IRFN == -3)) then
      if ((IRFN == 0) .or. (IRFN == -2) .or. (IRFN == -3)) DREF = 0.d0
      ! If  RFN(R)  is   1/(distance)**|IRFN|  ....
      J = -IRFN
      write(6,616)-IRFN
      do I=1,NPP
        RFN(I) = 1.d0/RVB(I)**J
      end do
    end if
    ! Any other user-defined matrix element argument radial function
    ! may be introduced to the code here, and invoked by:  IRFN= -4
    ! Note that the existing  RVB(i)  array is the radial distances  R .
    if (IRFN <= -10) then
      ! Illustrative user-defined analysis RFN(R) function
      !write(6,*) 'Print description of function introduced'
      !write(6,*) 'Use Freedman Pade DMF for CO'
      !do I=1,NPP
      !  !--------------------------------------------------------------
      !  RFN(I) = {calculate users chosen radial function}
      !  ! Freedman's DMF for CO  --------------------------------------
      !  RFN(I) = {calculate users chosen radial function}
      !  data coeff_new /-24.6005858d0,-109.5939637d0,-524.8233323d0,4.5194090d0,19.7954955d0,6.6011985d0,19.7206690d0/
      !  dm = -0.122706d0*(1.0+coeff(1)*x+coeff(2)*x*x+coeff(3)*x**3)/(1.0+coeff(4)*x+coeff(5)*x*x+coeff(6)*x**3+coeff(7)*x**6)
      !  !--------------------------------------------------------------
      !  XX = RFN(I)/1.128322714d0-1.d0
      !  RFN(I)= -0.122706d0*(1.d0+XX*(-24.6005858d0+XX*(-109.5939637d0+XX*(-524.8233323d0))))/ &
      !          (1.d0+XX*(4.5194090d0+XX*(19.7954955d0+XX*(6.6011985d0+19.7206690d0*XX**3))))
      !end do
    end if
    if ((IRFN == -1) .or. ((IRFN >= 1) .and. (IRFN <= 9))) then
      ! If  RFN(R)  is the Dunham or Surkus-type distance coordinate
      if (IRFN == -1) write(6,615)
      if ((IRFN >= 1) .and. (IRFN <= 9)) write(6,611) IRFN,IRFN,IRFN,IRFN
      if (DREF > 0.d0) then
        DREFP = DREF**IRFN
        write(6,613) DREF
        do I=1,NPP
          XX = YMINN+I*YH
          if (IRFN == -1) RFN(I) = (RVB(I)-DREF)/DREF
          if (IRFN >= 1) RFN(I) = (RVB(I)**IRFN-DREFP)/(RVB(I)**IRFN+DREFP)
        end do
      else
        write(6,610)
      end if
    end if

    !QQQQQQQQQQQ new ... the following not fully tested ...

    ! If  RFN(R)  is defined by interpolating over read-in points, use
    ! potential generating routine to do interpolation/extrapolation.
    if (IRFN >= 10) then
      MORDR = 1
      DM(0) = 0.d0
      DM(1) = 1.d0
      write(6,603)
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
      !-----------------------------------------------------------------------
      !read(5,*) NRFN,RFLIM
      !read(5,*) NUSEF,ILRF,NCNF,CNNF
      !read(5,*) RFACTF,MFACTF
      !read(5,*) (XIF(I),YIF(I),I=1,NRFN)
      ! If you uncomment the above, you better also uncomment the
      ! initialization to 0 below.
      !-----------------------------------------------------------------------
      MFACTF = 0
      RFACTF = 0
      write(6,810) NRFN,RFLIM
      if (NUSEF > 0) write(6,812) NUSEF,NRFN
      if (NUSEF <= 0) write(6,814) NRFN
      if ((ILRF > 1) .and. (dabs(CNNF) > 0.d0)) write(6,816) CNNF,NCNF
      write(6,818) RFACTF,MFACTF
      NROW = (NRFN+2)/3
      do J=1,NROW
        write(6,820) (XIF(I),YIF(I),I=J,NRFN,NROW)
      end do
      do I=1,NRFN
        XIF(I) = XIF(I)*RFACTF
        YIF(I) = YIF(I)*MFACTF
      end do
      810 format(' Transition moment function defined by interpolating over',I4,' read-in points'/5x, &
                 'and approaching the asymptotic value',f12.6)
      812 format(' Perform',I3,'-point piecewise polynomial interpolation over',I5,' input points')
      814 format(' Perform cubic spline interpolation over the',I5,' input points')
      816 format('- Beyond read-in points extrapolate to limiting asymptoticbehaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
      818 format(' Scale input points:  (distance)*',1PD16.9,'     (moment)*',D16.9/4x, &
                 'to get required units  [Angstroms & debye]'/3('      R(i)         Y(i)  ')/3(3X,11('--')))
      820 format((3(F12.6,F13.6)))
      IR2F = 0
      call GENINT(LRPT,NPP,RVB,RFN,NUSEF,IR2F,NRFN,XIF,YIF,RFLIM,ILRF,NCNF,CNNF)
    end if
  end if
  if ((MORDR >= 0) .and. (iabs(IRFN) <= 9)) write(6,602) (DM(J),J=0,MORDR)
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
!=======================================================================
!** INNOD2 specified wave fx. initiation at RMIN.  Normal case of
!  INNOD2 > 0  gives initiation with wave fx. node @ RMIN.
!  INNOD2 <= 0  give initiation with  zero slope @ RMIN.  This determines
!    symmetric eigenfunctions for rare special case when input potential
!    is half of a precisely symmetric potential with mid-point at RMIN.
!=======================================================================
if (iabs(LXPCT) >= 3) then
  !---------------------------------------------------------------------
  !read(5,*) NLEV2,AUTO2,J2DL,J2DU,J2DD,INNOD2
  !---------------------------------------------------------------------
  !read(5,*) NLEV2,AUTO2,J2DL,J2DU,J2DD
  !---------------------------------------------------------------------
  if (NLEV2 > VIBMX) NLEV2 = VIBMX
  if (NLEV2 <= 0) then
    write(6,644) NLEV2
    go to 997
    !stop
  end if
  !---------------------------------------------------------------------
  !if (AUTO2 > 0) read(5,*) (IV2(I),I=1,NLEV2)
  if (AUTO2 <= 0) then
    !read(5,*) (IV2(I),ZK2(I,1),I=1,NLEV2)
    !-------------------------------------------------------------------
    !** Give potential-2 trial energy the correct vibrational label
    do I=1,NLEV2
      ZK2(IV2(I),0) = ZK2(I,1)
    end do
  end if
  if (NUMPOT > 1) then
    if (INNOD2 > 0) write(6,686) 2
    if (INNOD2 <= 0) write(6,688) 2
  end if
  VMAX2 = 0
  do ILEV2=1,NLEV2
    VMAX2 = max(VMAX2,IV2(ILEV2))
  end do
  if (MORDR < 0) DM(1) = 1.d0
  SOMEG2 = IOMEG2**2
  if (J2DD == 0) J2DD = 1
  if (AUTO2 > 0) write(6,634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),I=1,NLEV2)
  if (AUTO2 <= 0) write(6,6634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),ZK2(IV2(I),0),I=1,NLEV2)
  if (NUMPOT >= 2) then
    if (IOMEG2 >= 99) then
      write(6,609)
    else
      if (IOMEG2 >= 0) write(6,608) 2,IOMEG2,SOMEG2
      if (IOMEG2 < 0) write(6,6085) 2,IOMEG2,-IOMEG2
    end if
  end if
  write(6,632)
end if

if (AUTO1 > 0) then
  ! If using automatic search for desired levels, subroutine ALFas gets
  ! eigenvalues ZK1(v,0) for desired vibrational levels of Potential-1,
  ! centrifugally-distorted to J=JREF.
  EJREF = JREF*(JREF+1)*YH**2
  ! Replace  [J(J+1)] by  [J(J+1) + |IOMEG1|]  for Li2(A) and like cases.
  !if (IOMEG1 < 0) EJREF = EJREF-dfloat(IOMEG1)*YH**2
  if (IOMEG1 < 0) EJREF = EJREF-dble(IOMEG1)*YH**2
  do I=1,NPP
    VBZ(I) = V1BZ(I)+EJREF*RRM2(I)
    VJ(I) = V1(I)+EJREF*RM2(I)
  end do
  ! OPTIONALLY WRITE SOME VALUES WHEN DEBUGGING:
  !write(6,*) 'EJREF=',EJREF
  !do I=1,3
  !  write(6,*) 'VJ=',VJ(I)
  !  write(6,*) 'V1',V1(I)
  !  write(6,*) 'RM2=',RM2(I)
  !end do
  if ((NLEV1 == 1) .and. (IV(1) > 998)) then
    !** Option to search for very highest level (within 0.0001 cm-1 of Disoc)
    EO = VLIM1-0.0001d0
    KV = IV(1)
    write(6,*) ''
    write(6,*) 'Exiting level.f (1)'
    write(6,*) 'Entering schrq.f (1)'
    write(6,*) ''
    call SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
    IV(1) = KV
    if (KV >= 0) then
      write(6,622) IJ(1),KV,VLIM1-EO
      GV(KV) = EO
      VMAX1 = KV
    else
      write(6,626) J,0.001d0
      go to 2
    end if
  else
    VMAX = VMAX1
    AFLAG = JREF
    if ((iabs(LXPCT) > 2) .and. (NUMPOT == 1)) VMAX = max(VMAX1,VMAX2)
    write(6,*) ''
    write(6,*) 'Exiting level.f'
    write(6,*) 'Entering alf.f'
    write(6,*) ''
    !write(6,*) 'NDP=',NPP
    !write(6,*) 'YMIN=',YMIN
    !write(6,*) 'YH=',YH
    !write(6,*) 'NCN=',NCN1
    !do I=1,3
    !  write(6,*) 'V=',VJ(I)
    !  write(6,*) 'SWF=',WF1(I)
    !  write(6,*) 'GV=',GV(I)
    !  write(6,*) 'INNR=',INNR1(I)
    !end do
    !write(6,*) 'VLIM=',VLIM1
    !write(6,*) 'KVMAX=',VMAX
    !write(6,*) 'AFLAG=',AFLAG
    !write(6,*) 'ZMU=',ZMU
    !write(6,*) 'EPS=',EPS
    !write(6,*) 'BFCT=',BFCT
    !write(6,*) 'INNODE=',INNOD1
    !write(6,*) 'IWR=',IWR
    !write(6,*) ''
    !deallocate(RVB)
    call ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,VMAX,AFLAG,ZMU,EPS,GV,BFCT,INNOD1,INNR1,IWR)
    VMAX1 = VMAX
  end if
  !** Get band constants for v=0-VMAX1 for generating trial eigenvalues
  WARN = 0
  do ILEV1=0,VMAX
    KV = ILEV1
    EO = GV(KV)
    INNER = INNR1(KV)
    write(6,*) ''
    write(6,*) 'Exiting level.f (2)'
    write(6,*) 'Entering schrq.f (2)'
    write(6,*) ''
    write(6,*) 'Getting band constants for v=',KV
    call SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,WARN,LPRWF)
    write(6,*) ''
    write(6,*) 'Exiting level.f'
    write(6,*) 'Entering cdjoel.f'
    write(6,*) ''
    ! OPTIONALLY WRITE THE INPUT PARAMETERS FOR DEBUGGING:
    !write(6,*) 'EO=',EO
    !write(6,*) 'NBEG=',NBEG
    !write(6,*) 'NEND=',NEND
    !write(6,*) 'BvWN=',BvWN
    !write(6,*) 'YH=',YH
    !write(6,*) 'WARN=',WARN
    !write(6,*) 'VJ(1)=',VJ(1)
    !write(6,*) 'WF1(1)=',WF1(1)
    !write(6,*) 'RM2(1)=',RM2(1)
    !write(6,*) 'RCNST(1)=',RCNST(1)
    ! For a-state test case, there's a memory error before CALL CDJOELas for v=2
    !call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
    if (NLEV1 < 0) then
      IV(ILEV1+1) = KV
      IJ(ILEV1+1) = JREF
    end if
    ZK1(ILEV1,0) = GV(ILEV1)
    do M=1,7
      ZK1(ILEV1,M) = RCNST(M)
    end do
  end do
end if
if (iabs(LXPCT) > 2) then
  if (AUTO2 > 0) then
    ! If using automatic location for levels of potential-2 (AUTO2 > 0)
    ! for matrix element calculation, also need Potential-2 band constants
    ! (rotational energy derivatives) ... again, calculate them at J=JREF
    if (NUMPOT > 1) then
      AFLAG = JREF
      do I=1,NPP
        VBZ(i) = V2BZ(I)+EJREF*RRM22(I)
        VJ(I) = V2(I)+EJREF*RM22(I)
      end do
      call ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,VMAX2,AFLAG,ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
    end if
  end if
  do ILEV2=1,NLEV2
    if (NUMPOT == 1) then
      ! For matrix elements within a single potl., copy above band constants
      do M=0,7
        ZK2(IV2(ILEV2),M) = ZK1(IV2(ILEV2),M)
      end do
    else
      ! ... otherwise, generate them (as above) with SCHRQ & CDJOEL
      KV = IV2(ILEV2)
      if (AUTO2 > 0) EO = GV(KV)
      if (AUTO2 <= 0) EO = ZK2(KV,0)
      INNER = INNR2(KV)
      write(6,*) ''
      write(6,*) 'Exiting level.f (3)'
      write(6,*) 'Entering schrq.f (3)'
      write(6,*) ''
      call SCHRQas(KV,JREF,EO,GAMA,PMAX2,VLIM2,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD2,INNER,WARN,LPRWF)
      call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
      ZK2(IV2(ILEV2),0) = EO
      do M=1,7
        ZK2(IV2(ILEV2),M) = RCNST(M)
      end do
    end if
  end do
end if
WARN = 1
EJREF = EJREF/YH**2
if (NLEV1 <= 0) NLEV = VMAX1+1

!===== Begin Actual Potential-1 Eigenvalue Calculation Loop Here =======
! Loop to compute eigenvalues ... etc. for NLEV levels of Potential-1
do ILEV1=1,NLEV
  KV = IV(ILEV1)
  if (KV < 0) exit
  NJMM = max(NJM,IJ(ILEV1))
  JROT = IJ(ILEV1)-JDJR
  IQT = 0
  JCT = 0
  ! If NJM > IJ(ILEV1) loop over range of rotational levels too
  do JLEV=IJ(ILEV1),NJMM,JDJR
    JROT = JROT+JDJR
    EJ = JROT*(JROT+1)-SOMEG1
    if (IOMEG1 >= 99) EJ = JROT*JROT-0.25d0
    ! If   IOMEG < 0   centrifugal term is  [J(J+1) + |IOMEG|]
    !if (IOMEG1 < 0) EJ=JROT*(JROT+1)-dfloat(IOMEG1)
    if (IOMEG1 < 0) EJ = JROT*(JROT+1)-dble(IOMEG1)
    ! If appropriate (AUTO1>0) use ALFas results to generate trial eigenvalue
    if (AUTO1 > 0) then
      EO = ZK1(KV,0)
      DEJ = EJ-EJREF
      EJP = 1.d0
      do M=1,7
        EJP = EJP*DEJ
        EO = EO+EJP*ZK1(KV,M)
      end do
    else
      !... otherwise - use read-in trial energy
      if (IV(ILEV1) < VIBMX) EO = ZK1(IV(ILEV1),0)
      if (IV(ILEV1) >= VIBMX) EO = GV(ILEV1)
    end if
    if ((AUTO1 <= 0) .and. (dabs(ZK1(IV(ILEV1),0)) <= 0.d0)) then
      call SCATTLEN(JROT,SL,VLIM1,V1,WF1,BFCT,YMIN,YH,NPP,CNN1,NCN1,IWR,LPRWF)
      if (NUMPOT == 1) goto 2
      goto 104
    end if
    ! ... or if JLEV > IJ(ILEV1) ... use local Beff to estimate next level
    if (JLEV > IJ(ILEV1)) then
      BEFF = 0.d0
      do I=NBEG,NEND
        BEFF = BEFF+WF1(I)**2*RM2(I)
      end do
      BEFF = BEFF*YH*BvWN
      EO = ESLJ(JCT)+(2*JLEV+1-JDJR)*JDJR*BEFF
    end if
    !** Now add centrifugal term to get effective (radial) potential
    EJ = EJ*YH**2
    do J=1,NPP
      VBZ(J) = V1BZ(J)+EJ*RRM2(J)
      VJ(J) = V1(J)+EJ*RM2(J)
    end do
    ! Set wall outer boundary condition, if specified by input IV(ILEV1)
    if (KV < -10) then
      WF1(-IV(ILEV1)) = 0.d0
      WF1(-IV(ILEV1)-1) = -1.d0
    end if
    KVIN = KV
    if (AUTO1 > 0) INNER = INNR1(KV)
    if (SINNER /= 0) INNER = SINNER
    ! Call SCHRQ to find Potential-1 eigenvalue EO and eigenfn. WF1(i)
    write(6,*) ''
    write(6,*) 'Exiting level.f (4)'
    write(6,*) 'Entering schrq.f (4)'
    write(6,*) ''
    ! The next CALL SCHRQas for a-state failes during v=0, but EO's from when
    ! Entering schrq.f (2) were already accurate, so the next step is not necessary.
    100 continue
    call SCHRQas(KV,JROT,EO,GAMA,PMAX1,VLIM1,VJ,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
    if (KV < 0) then
      ! SCHRQ  error condition is  (KV < 0) .
      if (NJM > IJ(ILEV1)) then
        ! ... in automatic search for ever-higher J levels
        if (IQT <= 0) then
          ! ... try one more time with E(trial) slightly below barrier maximum
          IQT = 1
          EO = PMAX1-0.1d0
          go to 100
        else
          KV = KVIN
          go to 130
        end if
      end if
      go to 122
    end if
    !write(6,*) 'Hello! (for debugging) 1'
    if ((KV /= KVIN) .and. (AUTO1 > 0)) then
      !** If got wrong vib level, do a brute force ALFas calculation to find it.
      KV = KVIN
      AFLAG = JROT
      call ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,KV,AFLAG,ZMU,EPS,GV,BFCT,INNOD1,INNR1,IWR)
      if (KV == KVIN) then
        EO = GV(KVIN)
        go to 100
      else
        write(6,618) KVIN,JROT,KV
        KV = KVIN
        go to 130
      end if
    end if
    !write(6,*) 'Hello! (for debugging) 2'
    if (KV /= IV(ILEV1)) IV(ILEV1) = KV
    ! If desired, calculate rotational & centrifugal distortion constants
    if (LCDC > 0) then
      if ((IOMEG1 > 0) .and. (JROT == 0)) then
        ! Calculate 'true' rotational constants for rotationless IOMEG>0 case
        write(6,*) ''
        write(6,*) 'Exiting level.f (5)'
        write(6,*) 'Entering schrq.f (5)'
        write(6,*) ''
        call SCHRQas(KV,0,EO,GAMA,PMAX1,VLIM1,V1,WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
        call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V1,WF1,RM2,RCNST)
      else
        ! Calculate rotational constants for actual (v,J) level of interest.
        call CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,RCNST)
      end if
      if (dabs(EO) > 1.d0) then
        write(6,606) KV,JROT,EO,(RCNST(M),M=1,7)
      else
        write(6,6065) KV,JROT,EO,(RCNST(M),M=1,7)
      end if
      if (LCDC > 1) then
        if (dabs(EO) > 1.d0) then
          !write(9,902) KV,JROT,EO,(RCNST(M),M=1,7)
        else
          !write(9,904) KV,JROT,EO,(RCNST(M),M=1,7)
        end if
      end if
    end if
    !write(6,*) 'Hello! (for debugging) 3'
    if (LXPCT == -1) write(7,703) KV,JROT,EO,GAMA
    if (((LXPCT == 1) .or. (iabs(LXPCT) == 2)) .or. &
        ((iabs(LXPCT) > 2) .and. ((IRFN == -1) .or. (IRFN >= 1) .and. (IRFN >= 9)) .and. (DREF <= 0.d0))) then
      ! Calculate various expectation values in LEVXPC
      call LEVXPC(KV,JROT,EO,GAMA,NPP,WF1,RFN,VBZ,VLIM1,YH,DREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
      if ((LXPCT > 0) .and. (MORDR > 0)) write(6,632)
    end if
    !write(6,*) 'Hello! (for debugging) 4'
    104 continue
    if ((iabs(LXPCT) <= 2) .or. (NLEV2 <= 0)) go to 122
    !** If desired, now calculate off-diagonal matrix elements, either
    !  between levels of different potentials, IF(NUMPOT >= 2), or between
    !  levels of a single potential, for (NUMPOT <= 1).
    !** First prepare centrifugally distorted potential, trial energy, etc.,
    !  and calculate second wave function and matrix element(s)
    do ILEV2=1,NLEV2
      ! For case of a single potential, avoid redundancy by considering
      ! only emission
      if ((NUMPOT <= 1) .and. (IV2(ILEV2) > KV)) go to 120
      ! Loop over J2's allowed by given selection rule.
      do IJD=J2DL,J2DU,J2DD
        KV2 = IV2(ILEV2)
        KVIN = KV2
        JROT2 = JROT+IJD
        if (JROT2 < 0) go to 116
        if ((NUMPOT <= 1) .and. (IV2(ILEV2) == KV) .and. (JROT2 > JROT)) go to 116
        EJ2 = JROT2*(JROT2+1)-SOMEG2
        if (IOMEG2 >= 99) EJ2 = JROT2**2-0.25d0
        !... allow for weird Li2(A) and Li2(c) potential cases
        !if (IOMEG2 < 0) EJ2=JROT2*(JROT2+1)-dfloat(IOMEG2)
        if (IOMEG2 < 0) EJ2 = JROT2*(JROT2+1)-dble(IOMEG2)
        EO2 = ZK2(KV2,0)
        DEJ = EJ2-EJREF
        EJP = 1.d0
        ! Use calculated state-2 CDC's to predict trial eigenvalue
        do M=1,7
          EJP = EJP*DEJ
          EO2 = EO2+EJP*ZK2(KV2,M)
        end do
        ! Now ... update to appropriate centrifugally distorted potential
        EJ2 = EJ2*YH*YH
        do I=1,NPP
          VBZ(I) = V2BZ(I)+EJ2*RRM22(I)
          VJ(I) = V2(I)+EJ2*RM22(I)
        end do
        INNER = INNR2(KV2)
        if (SINNER /= 0) INNER = SINNER
        ICOR = 0
        write(6,*) ''
        write(6,*) 'Exiting level.f (6)'
        write(6,*) 'Entering schrq.f (6)'
        write(6,*) ''
        110 continue
        call SCHRQas(KV2,JROT2,EO2,GAMA,PMAX2,VLIM2,VJ,WF2,BFCT,EPS,YMIN,YH,NPP,NBEG2,NEND2,INNOD2,INNER,IWR,LPRWF)
        if (KV2 /= KVIN) then
          if (KV2 < 0) go to 114
          ! Using CDC's to estimate trial eigenvalue failed:
          ICOR = ICOR+1
          if (ICOR <= 2) then
            ! ... first correction attempt ... use semiclassical dv/dE to improve
            GB = -1.d0
            GI = -1.d0
            ! ... hey RJ!!  shouldn't you update this using SCECOR ?
            WV = 0.d0
            XX = EO2*BFCT
            do I=NBEG2,NEND2
              GBB = GB
              GB = GI
              GI = XX-VJ(I)
              if ((GBB > 0.d0) .and. (GI > 0.d0)) WV = WV+1.d0/dsqrt(GB)
            end do
            WV = 6.2832d0/(BFCT*WV)
            EO2 = EO2+WV*(KVIN-KV2)
            go to 110
          end if
          write(6,633) IV2(ILEV2),JROT2,KV2
          ! ... if that fails, do a brute force ALFas calculation to find it.
          114 continue
          KV2 = KVIN
          AFLAG = JROT2
          call ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,KV2,AFLAG,ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
          if (KV2 == KVIN) then
            EO2 = GV(KV2)
            INNER = INNR2(KV2)
            go to 110
          else
            write(6,618) KVIN,JROT,KV2
            go to 116
          end if
        end if
        if (NBEG > NBEG2) NBEG2 = NBEG
        if (NEND < NEND2) NEND2 = NEND
        !if ((NUMPOT <= 1) .and. (EO2 > (EO+EPS))) go to 120
        call MATXEL(KV,JROT,IOMEG1,EO,KV2,JROT2,IOMEG2,IRFN,EO2,NBEG2,NEND2,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
        116 continue
      end do
      ! End of Potential-2 rotational selection level loop
      120 continue
    end do
    !++++ End of Potential-2 vibrational level matrix element loop +++++++++

    122 continue
    JCT = JCT+1
    ! ... check to avoid array overflow
    if (JCT > VIBMX) then
      write(6,637) VIBMX
      go to 997
      !stop
    end if
    JWR(JCT) = JROT
    ESLJ(JCT) = EO
  end do
  !write(6,*) 'Hello! (for debugging) 5'
  !++ End of Potential-1 loop over NJM-specified J-sublevels
  130 continue
  if (NJM > IJ(ILEV1)) then
    ! Print rotational sublevels generated for vibrational level  ILEV
    NROW = (JCT+4)/5
    write(6,627) KV
    do J=1,NROW
      write(6,628) (JWR(I),ESLJ(I),I=J,JCT,NROW)
    end do
    write(6,641)
  end if
  ESOLN(ILEV1) = ESLJ(1)
end do
!++ End of loop over the NLEV Potential-1 input levels
write(6,*) ''
write(6,*) 'SUMMARY (ALL ENERGIES IN CM-1):'
if (NLEV1 < 0) then
  NROW = (NLEV+3)/4
  write(6,623) NLEV,IJ(1)
  do J=1,NROW
    write(6,630) (IV(I),ESOLN(I),I=J,NLEV,NROW)
  end do
  if ((NLEV > 1) .and. (IJ(1) == 0) .and. (NCN1 > 0) .and. (ESOLN(NLEV) < VLIM1)) then
    ! In (NLEV1 < 0) option, estimate vD using the N-D theory result that:
    ! (vD - v) {is proportional to} (binding energy)**((NCN-2)/(2*NCN))
    VDMV = 1.d0/(((VLIM1-ESOLN(NLEV-1))/(VLIM1-ESOLN(NLEV)))**(1.d0/PW)-1.d0)
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
    write(6,620) NCN1,VD
    !if ((IVD > = IVS) .and. (dfloat(IV(NLEV))/VD > 0.9d0)) then
    if ((IVD >= IVS) .and. (dble(IV(NLEV))/VD > 0.9d0)) then
      NFP = NLEV+1
      do I=IVS,IVD
        NLEV = NLEV+1
        IV(NLEV) = IV(NLEV-1)+1
        ESOLN(NLEV) = VLIM1-(VLIM1-ESOLN(NLEV-1))*(1.d0-1.d0/VDMV)**PW
        VDMV = VDMV-1.d0
      end do
      NLP = NLEV-NFP+1
      NROW = (NLP+3)/4
      write(6,621) NLP
      do J=1,NROW
        III = NFP+J-1
        write(6,630) (IV(I),ESOLN(I),I=III,NLEV,NROW)
      end do
    end if
  end if
end if
if ((NJM <= 0) .and. (NLEV1 >= 0)) then
  NROW = (NLEV+2)/3
  write(6,619) NLEV
  do J=1,NROW
    write(6,631) (IV(I),IJ(I),ESOLN(I),I=J,NLEV,NROW)
  end do
end if
write(6,601)
! The following two lines are to read input file again if you want to find
! the levels for a different potential. Currently LEVEL_RDINP.F90 doesn't
! allow us to read the input file a second time though.
!go to 2
!999 stop
!-------------------------------------------------------------------
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

997 continue
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
!824 format(//A78/30('==')/" Note that (v',J') &",' (v",J") strictly label the upper and lower levels, resp.,'/6x, &
!           'and  E(lower)=E"'/' but  E(2)-E(1)  is:  (energy of State-2 level) - (energy of State-1 level)'//12x,'Band'/ &
!           ' dJ(J")',4x,7hv'   v",'  E(lower)  E(2)-E(1)  A(Einstein)   F-C Factor  ',13h<v'j'|M|v"j"> /1x,3('--'),( &
!           '   -------'),'  --------',3x,4('--'),3x,11('-'),3x,11('-'),3x,11('-') )
!811 format(//A78/30('==')//12x,"   v'","  J'",'    v"','  J"   position    E(upper)    E(lower)',16h   <v'j'|M|v"j">/1x,68('-') )
!901 format(//A78/1x,62('==')/'   v    J',7x,'E',10x,'Bv',11x,'-Dv',13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
!902 format(I4,I5,f25.15,f14.10,6(1PD15.7))
!904 format(I4,I5,f25.15,1PD14.7,6(D15.7))

end subroutine LEVEL
