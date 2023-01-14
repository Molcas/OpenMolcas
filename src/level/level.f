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
c***********************************************************************
c Please inform me of any bugs at nike@hpqc.org or ndattani@uwaterloo.ca
c***********************************************************************
c******  Eigenvalue program  LEVEL 2022 ********************************
c   !!!! Form modified for handling Stolyarov radial variable  !!!!!
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Program for calculating eigenvalues and eigenfunctions (and if
c  desired, also various expectation values & matrix elements) of a
c  one-dimensional potential, and/or matrix elements (& Franck-Condon
c  factors) between levels of two different potentials.
c** As with most similar codes, the original version of this program was
c  based on the Franck-Condon intensity program of R.N. Zare, report
c  UCRL-10925(1963), but the present version is massively modified.
c** This program is unique in that it can:  (1) automatically locate &
c      calculate the widths of quasibound levels (orbiting resonances);
c  (2) can calculate diatomic molecule centrifugal distortion constants;

c      it will also automatically generate the eigenvalues etc. for all
c      vibrational and/or rotational levels of a given well-behaved
c      single-minimum potential.
c***** Main calling and I/O routines.  Last Updated  28 June 2009 *****
      SUBROUTINE LEVEL(RC)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
      USE LEVEL_COMMON
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: RC
c** Dimensions for  potential arrays  and  vib. level arrays.
      INTEGER VIBMX,MORDRMX,RORDR,NTPMX
      PARAMETER (VIBMX=400,RORDR=7,MORDRMX=20,NTPMX= 1600)
c!!---------------------------------------------------------------------
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive. fmax-stack-var-size=2^20
!     PARAMETER (NDIMR= 131074) ! MMA_ALLOCATE won't allow PARAMETERs
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                   SDRDY(NDIMR),VBZ(NDIMR),aRVp
      REAL*8 PRV,ARV,aRVp
!     REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),SDRDY(:),
!    1 VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
c!!---------------------------------------------------------------------
      INTEGER I,J,M,III,IJD,ILEV1,ILEV2,IOMEG1,IOMEG2,INNOD1,INNOD2,
     1 INNER,SINNER,IQT,IWR,IRFN,IVD,IVS,IAN1,IAN2,IMN1,IMN2,GEL1,GEL2,
     2 GNS1,GNS2,JDJR,JCT,J2DL,J2DU,J2DD,JROT,JROT2,JLEV,JREF, ICOR,
     3 CHARGE, KV,KV2,KVIN,LCDC,LPRWF,LRPT,LXPCT,MORDR,NUSEF,ILRF,IR2F,
     4 NUMPOT,NBEG,NBEG2,NEND,NEND2,NPP,NCN1,NCN2,NCNF,NLEV,NLEV1,
     5 NLEV2,NJM,NJMM,NFP,NLP,NRFN,NROW,WARN,VMAX,VMAX1,VMAX2,AFLAG,
     6 AUTO1,AUTO2, IV(VIBMX),IJ(VIBMX),IV2(VIBMX),JWR(VIBMX),
     7 INNR1(0:VIBMX),INNR2(0:VIBMX),NTP,LPPOT,IPOTL,PPAR,QPAR,NSR,NLR,
     8 IBOB,NCMM,IVSR,IDSTT,MMLR(3)
c
!     REAL*8 ZK1(0:VIBMX,0:RORDR),ZK2(0:VIBMX,0:RORDR),RCNST(RORDR),
!    1 V1(NDIMR),V2(NDIMR),VJ(NDIMR),V1BZ(NDIMR),V2BZ(NDIMR),
!    2 WF1(NDIMR),WF2(NDIMR),CMM(3),PARM(4)
      REAL*8 ZK1(0:VIBMX,0:RORDR),ZK2(0:VIBMX,0:RORDR),RCNST(RORDR),
     1 CMM(3),PARM(4)
      REAL*8, ALLOCATABLE ::  V1(:),V2(:),VJ(:),V1BZ(:),V2BZ(:),
     1 WF1(:),WF2(:)
c
!     REAL*8  RFN(NDIMR),RRM2(NDIMR),RM2(NDIMR),RRM22(NDIMR),
!    2  RM22(NDIMR),GV(0:VIBMX),ESOLN(VIBMX),ESLJ(VIBMX),XIF(NTPMX),
!    4  YIF(NTPMX),ABUND1,ABUND2,MASS1,MASS2,DM(0:MORDRMX)
      REAL*8 GV(0:VIBMX),ESOLN(VIBMX),ESLJ(VIBMX),XIF(NTPMX),
     1  YIF(NTPMX),ABUND1,ABUND2,MASS1,MASS2,DM(0:MORDRMX)
      REAL*8, ALLOCATABLE :: RFN(:),RRM2(:),RM2(:),RRM22(:),RM22(:)
      REAL*8 BZ,BvWN,BFCT,BEFF,DEJ,EPS,EO,EO2,EJ,EJ2,EJP,EJREF,GAMA,
     1 MEL,PMAX1,PMAX2,PW,RH,RMIN,RR,RRp,PINV,DRDY,YH,YH2,YMIN,YMINN,
     2 YMAX,DREF,DREFP,CNN1,CNN2,RFLIM,CNNF,RFACTF,MFACTF,SOMEG1,
     3 SOMEG2,VLIM1,VLIM2,VD,VDMV,XX,ZMU,GI,GB,GBB,WV,FFAS,SL,DSCM,
     4 REQ,RREF,RHOAB
c
      CHARACTER*78 TITL
      CHARACTER*2 NAME1,NAME2
c
      DATA MEL/5.4857990945d-4/,YMAX/1.d+00/
c** Default (Q-branch) defining J-increments for matrix element calcn.
      DATA J2DL,J2DU,J2DD/0,0,1/
      NDIMR = 131074
      CALL MMA_ALLOCATE(RVB,NDIMR,label='RVB')
      CALL MMA_ALLOCATE(YVB,NDIMR,label='YVB')
      CALL MMA_ALLOCATE(DRDY2,NDIMR,label='DRDY2')
      CALL MMA_ALLOCATE(FAS,NDIMR,label='FAS')
      CALL MMA_ALLOCATE(SDRDY,NDIMR,label='SDRDY')
      CALL MMA_ALLOCATE(VBZ,NDIMR,label='VBZ')
!
      CALL MMA_ALLOCATE(V1,NDIMR,label='V1')
      CALL MMA_ALLOCATE(V2,NDIMR,label='V2')
      CALL MMA_ALLOCATE(VJ,NDIMR,label='VJ')
      CALL MMA_ALLOCATE(V1BZ,NDIMR,label='V1BZ')
      CALL MMA_ALLOCATE(V2BZ,NDIMR,label='V2BZ')
      CALL MMA_ALLOCATE(WF1,NDIMR,label='WF1')
      CALL MMA_ALLOCATE(WF2,NDIMR,label='WF2')
!
      CALL MMA_ALLOCATE(RFN,NDIMR,label='RFN')
      CALL MMA_ALLOCATE(RRM2,NDIMR,label='RRM2')
      CALL MMA_ALLOCATE(RM2,NDIMR,label='RM2')
      CALL MMA_ALLOCATE(RRM22,NDIMR,label='RRM22')
      CALL MMA_ALLOCATE(RM22,NDIMR,label='RM22')
      PINV=1.d0
      NLEV2=-1
      AUTO2=0
      VMAX2=0
      IOMEG2=0
      SOMEG2=0
      CNN2=0
      PMAX2=0
      ABUND2=0
      MASS2=0
      NCN2=0
      NEND2=0
      NBEG2=0
      GNS2=0
      GEL2=0
      INNOD2=0
      EJREF=0
      DO I=1,VIBMX
       IV2(I)=0
      ENDDO
      DO I=1,RORDR
       RCNST(I)=0
      ENDDO
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Begin by reading in the (integer) atomic numbers and mass numbers
c  defining the effective reduced mass of the system considered.
c** IAN1 & IAM2, and IMN1 & IMN2 are, respectively, the atomic numbers
c    and the mass numbers identifying the atoms forming the molecule.
c    Their masses are extracted from data subroutine MASSES and used
c    to generate the the reduced mass ZMU.
c** If  IMN1  or  IMN2  lie outside the range of mass numbers for normal
c  stable isotopes of that species, subroutine MASSES returns the
c  average atomic mass based on the natural isotope abundance.
c** If the read-in value of IAN1 and/or IAN2 is .LE.0, then instead of
c  using the MASS table, read an actual particle mass for it/them.
c** CHARGE (integer) is the charge on the molecule (=0 for neutral). If
c   (CHARGE.ne.0)  generate & use Watson's charge-adjusted reduced mass.
c** Parameter NUMPOT specifies whether to calculate eigenvalues etc. for
c  a single potential (when NUMPOT.LE.1), or to generate two independent
c  potentials & calculate matrix elements coupling levels of one to
c  levels of the other (for NUMPOT.GE.2).
c----------------------------------------------------------------------
    2 CALL LEVEL_RDINP(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,
     1 ARV,EPS,NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,
     2 DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1,AUTO1,
     3 LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)
! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING:
!     WRITE(6,*) IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,ARV,EPS
!     WRITE(6,*) NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB
!    WRITE(6,*) DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1
! OPTIONALLY WRITE THE INPUT KEYWORDS WHEN DEBUGGING (ANOTHER WAY):
!     WRITE(6,*) AUTO1,LCDC,LXPCT,NJM,JDJR,IWR,LPRWF
!     WRITE(6,*) 'level.f has the following after CALL LEVEL_RDINP:'
!     WRITE(6,*) 'IAN1 = ',IAN1
!     WRITE(6,*) 'IMN1 = ',IMN1
!     WRITE(6,*) 'IAN2 = ',IAN2
!     WRITE(6,*) 'IMN2 = ',IMN2
!     WRITE(6,*) 'CHARGE = ',CHARGE
!     WRITE(6,*) 'NUMPOT = ',NUMPOT
!     WRITE(6,*) 'RH = ',RH
!     WRITE(6,*) 'RMIN = ',RMIN
!     WRITE(6,*) 'PRV = ',PRV
!     WRITE(6,*) 'ARV = ',ARV
!     WRITE(6,*) 'EPS = ',EPS
!     WRITE(6,*) 'NTP = ',NTP
!     WRITE(6,*) 'LPPOT = ',LPPOT
!     WRITE(6,*) 'IOMEG1 = ',IOMEG1
!     WRITE(6,*) 'VLIM = ',VLIM
!     WRITE(6,*) 'IPOTL = ',IPOTL
!     WRITE(6,*) 'PPAR = ',PPAR
!     WRITE(6,*) 'QPAR = ',QPAR
!     WRITE(6,*) 'NSR = ',NSR
!     WRITE(6,*) 'NLR = ',NLR
!     WRITE(6,*) 'IBOB = ',IBOB
!     WRITE(6,*) 'DSCM = ',DSCM
!     WRITE(6,*) 'REQ = ',REQ
!     WRITE(6,*) 'RREF = ',RREF
!     WRITE(6,*) 'NCMM = ',NCMM
!     WRITE(6,*) 'IVSR = ',IVSR
!     WRITE(6,*) 'IDSTT = ',IDSTT
!     WRITE(6,*) 'RHOAB = ',RHOAB
!     WRITE(6,*) 'MMLR = ',MMLR
!     WRITE(6,*) 'CMM = ',CMM
!     WRITE(6,*) 'PARM = ',PARM
!     WRITE(6,*) 'NLEV1 = ',NLEV1
!     WRITE(6,*) 'AUTO1 = ',AUTO1
!     WRITE(6,*) 'LCDC = ',LCDC
!     WRITE(6,*) 'LXPCT = ',LXPCT
!     WRITE(6,*) 'NJM = ',NJM
!     WRITE(6,*) 'JDJR = ',JDJR
!     WRITE(6,*) 'IWF = ',IWF
!     WRITE(6,*) 'LPRWF = ',LPRWF
!     READ(5,*,END=999)
!   2 READ(5,*,END=999) IAN1, IMN1, IAN2, IMN2, CHARGE, NUMPOT
c----------------------------------------------------------------------
c** Subroutine MASSES returns the names of the atoms NAMEi,ground
c  electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
c  mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
      IF((IAN1.GT.0).AND.(IAN1.LE.109)) THEN
          CALL MASSES(IAN1,IMN1,NAME1,GEL1,GNS1,MASS1,ABUND1)
      ELSE
c** If particle-i is not a normal atomic isotope, read a 2-character
c   name (enclosed between '', as in 'mu') and its actual mass.
c----------------------------------------------------------------------
!         READ(5,*) NAME1, MASS1
c----------------------------------------------------------------------
      ENDIF
      IF((IAN2.GT.0).AND.(IAN2.LE.109)) THEN
          CALL MASSES(IAN2,IMN2,NAME2,GEL2,GNS2,MASS2,ABUND2)
      ELSE
c----------------------------------------------------------------------
!         READ(5,*) NAME2, MASS2
c----------------------------------------------------------------------
      ENDIF
      ZMU= MASS1*MASS2/(MASS1+MASS2- CHARGE*MEL)
c=======================================================================
c TITL is a title or output header of up to 78 characters, read on a
c   single line enclosed between single quotes: e.g.  'title of problem'
c=======================================================================
!     READ(5,*) TITL
      TITL = 'Beginning execution of LEVEL:'
c----------------------------------------------------------------------
c** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
c  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920D0
      WRITE(6,605) TITL,ZMU,BZ,MASS1,MASS2
      BvWN= 1.D0/BZ
      IF(CHARGE.NE.0) WRITE(6,624) CHARGE,CHARGE
      EJ= 0.D0
      EJ2= 0.D0
      LRPT= 1
c** Lower limit (RMIN) and increment (RH) of integration in (Angstroms).
c** Upper limit of the reduced variable integration range automatically
c  set at  YMAX= 1.0 , which corresponds to  RMAX= infinity !!.
c* A hard wall boundary condition may be imposed at a smaller distance
c  using an appropriate choice of the read-in level parameter IV (below)
c!! The radial integration variable is  yp(r;Reff)  with   p= PRV
c** EPS (cm-1) is the desired eigenvalue convergence criterion
c---------------------------------------------------------------------
!     READ(5,*) RH, RMIN, PRV, ARV, EPS
c---------------------------------------------------------------------
c** NPP = no. of points in potential and wavefunction array.
c!! First ... calculate new AS radial mesh YH implied but the given RH
      I= INT(0.5d7*(PRV/ARV)*RH)
c.... give YH a rounded-off value (to 8 digits)
      YH= DBLE(I)*1.d-07
      aRVp= ARV**PRV
      RRp= RMIN**PRV
      YMIN= (RRp - aRVp)/(RRp + aRVp)
      YMAX= 1.d0
c** NPP = no. of points in potential and wavefunction array.
      NPP= INT(((YMAX-YMIN)/YH+ 1.00001))
      IF(NDIMR.LT.NPP) THEN
!         WRITE(6,6604)  NDIMR,YH,DFLOAT(NPP)/DFLOAT(NDIMR)
          WRITE(6,6604)  NDIMR,YH,DBLE(NPP)/DBLE(NDIMR)
          NPP= NDIMR
      ENDIF
c... reset YMIN slightly to precisely span range
      YMIN= YMAX - (NPP-1)*YH
      YH2= YH*YH
      BFCT= BZ*YH2
      YMINN= YMIN-YH
      WRITE(6,604) YMIN,YMAX,YH,PRV,ARV,RMIN,RH,ARV,
     1                                           NAME1,IMN1,NAME2,IMN2
      PINV= 1.d0/PRV
      FFAS= YH2*(PINV**2 - 1.d0)/(4.d0*aRVp)**2
      DO  I= 2,NPP-1
          YVB(I)= YMINN + I*YH
          RRp= (1.d0 + YVB(I))/(1.d0 - YVB(I))
          RR= RRp**PINV
          RVB(I)= ARV*RR
          RRM2(I)= 1.D0/RVB(I)**2
          RRM22(I)= RRM2(I)
          RRp= RRp*aRVp
          DRDY= RVB(I)*(RRp + aRVp)**2/(2.d0*pRV*RRp*aRVp)
          DRDY2(I)= DRDY**2
          SDRDY(I)= DSQRT(DRDY)
          FAS(I)= FFAS*((RRp + aRVp)**2/RRp)**2
      ENDDO
! OPIONALLY WRITE SOME VARIABLES IF DEBUGGING:
!     WRITE(6,*) 'DRDY=',DRDY
!     WRITE(6,*) 'RRp=',RRp
!     WRITE(6,*) 'aRVp=',aRVp
!     WRITE(6,*) 'pRV=',pRV
!     WRITE(6,*) 'RRp=',RRp
!     WRITE(6,*) 'aRVp=',aRVp
!     WRITE(6,*) 'DRDY2(1)=',DRDY2(1)
!     WRITE(6,*) 'DRDY2(2)=',DRDY2(2)
!     WRITE(6,*) 'RVB(1)=',RVB(1)
!     WRITE(6,*) 'RVB(2)=',RVB(2)
      YVB(1)= YMIN
      RVB(1)= RMIN
      RRM2(1)= RRM2(2)
      DRDY2(1)= DRDY2(2)
      SDRDY(1)= SDRDY(2)
      IF(RMIN.GT.0.D0) THEN
          RRM2(1)= 1.D0/RMIN**2
          RRp= RMIN**PRV
          DRDY= RVB(1)*(RRp + aRVp)**2/(2.d0*PRV*RRp*aRVp)
          DRDY2(1)= DRDY**2
          SDRDY(1)= DSQRT(DRDY)
      ENDIF
      RRM22(1)= RRM2(1)
      YVB(NPP)= YMAX
c... 'fake' RMAX value to ensure last 1/R**2 point is stable.
      RVB(NPP)= RVB(NPP-1) + RVB(NPP-1) - RVB(NPP-2)
      RRp= RVB(NPP)**PRV
      DRDY= RVB(NPP)*(RRp + aRVp)**2/(2.d0*PRV*RRp*aRVp)
      DRDY2(NPP)= DRDY**2
      SDRDY(NPP)= DSQRT(DRDY)
      RRM2(NPP)= 1.d0/RVB(NPP)
      RRM22(NPP)= RRM2(NPP)
!     For debugging purposes, you can print the first 3 V(R) values:
!     DO I=1,3
!          WRITE(6,*) RVB(I)
!     ENDDO

c
c++ Begin reading appropriate parameters & preparing potential(s)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+  Subroutine "PREPOT" prepares (and if desired, writes) the potential
c+  array V(i) (cm-1)  at the NPP distances RVB(i) (Angst).
c** NPP = no. of points in potential and wavefunction array.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c* If NTP > 0 :  define potential by interpolation over & extrapolation
c        beyond the NTP read-in turning points using subroutine GENINT.
c   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
c* If LPPOT > 0 : at every |LPPOT|-th point, print potential and
c      derivatives-by-differences. ***  If  LPPOT < 0  write potential
c      at every |LPPOT|-th point to channel-8 in a compact format **
c* OMEGA is the electronic contribution to the angular momentum such
c  that the reduced centrifugal potential is:  (J*(J+1)-OMEGA**2)/R**2
c* Set (OMEGA.GE.99) if wish to use centrifugal factor for rotation
c  in two dimensions:   (J**2 - 1/4)/R**2  .
c* VLIM (cm-1) is the energy associated with the potential asymptote.
c-----------------------------------------------------------------------
c++   READ(5,*) NTP, LPPOT, OMEGA, VLIM
c----------------------------------------------------------------------
c** For pointwise potentials, PREPOT uses subroutine GENINT to read
c  points and conditions and interpolate (using subroutines NTRPSR,
c  SPLINE & SPLINT) and extrapolate to get full potential.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For a pointwise potential (NTP > 0), now read points & parameters
c  controlling how the interpolation/extrapolation is to be done.
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** NTP (read above) is number of turning points (XI,YI) to be read in.
c** If NUSE > 0  interpolate with NUSE-point piecewise polynomials
c    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
c    interpolate with cubic spline instead of local polynomials.
c** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
c   This may help if interpolation has trouble on steep repulsive wall.
c** ILR specifies how to extrapolate beyond largest input distance XI(i)
c  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
c  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
c  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
c** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
c  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
c  leading coefficient fixed at  CNN ; otherwise get it from points too.
c* Assume read-in CNN value has units:  [(cm-1)(Angstroms)**'NCN'].
c* If ILR = 2 or 3 , successive higher power terms differ by  1/R**2
c* If ILR > 3 : successive higher power terms differ by factor  1/R
c
c** RFACT & EFACT are factors required to convert units of input turning
c       points (XI,YI) to Angstroms & cm-1, respectively (often = 1.d0)
c** Turning points (XI,YI) must be ordered with increasing XI(I)
c** Energy VSHIFT (cm-1) is added to the input potential points to
c   make their absolute energy consistent with VLIM (often VSHIFT=Te).
c-----------------------------------------------------------------------
c++   READ(5,*) NUSE, IR2, ILR, NCN, CNN
c++   READ(5,*) RFACT, EFACT, VSHIFT
c++   READ(5,*) (XI(I), YI(I), I= 1,NTP)
c-----------------------------------------------------------------------
c** NCN1 (returned by PREPOT) is the power of the asymptotically-
c  dominant inverse-power long range potential term.
c   VLIM1, VLIM1, V1, NCN1 and CNN1 are not defined yet, but are input parameters for
c  PREPOT
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(6,*) 'Exiting level.f'
      WRITE(6,*) 'Entering prepot.f'
      WRITE(6,*) ''
      CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,
     1  V1,CNN1,NCN1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,
     2  MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
!     CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,
!    1                                                   V1,CNN1,NCN1)
      WRITE(6,*) 'Successfully made it through Prepot.f!'
!OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
!     DO I=1,3
!      WRITE(6,*) 'V(',I,')=',V1(I)
!     ENDDO
!     WRITE(6,*) 'V(                 20000)=',V1(20000)
!     WRITE(6,*) 'V(',NPP,')=',V1(NPP)
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** If (NTP.le.0) PREPOT uses subroutine POTGEN to generate a fully
c  analytic potential defined by the following read-in parameters.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c* Potentials generated in cm-1 with equilibrium distance REQ [Angst.],
c  and for all cases except IPOTL=2, the potential asymptote energy is
c  VLIM and well depth is DSCM.  For IPOTL=2, VLIM is the energy at the
c  potential minimum and  DSCM  the leading (quadratic) potential coeft.
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** IPOTL specifies the type of potential function to be generated.
c** MPAR, NSR & NCMM are integers characterizing the chosen potential
c** NVARB is number of (real*8) potential parameters read in.
c** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass
c      dependent Born-Oppenheimer breakdown corrections will be included
c- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
c** If IPOTL=1  generate an L.J.(MPAR,NCN) potential.
c** If IPOTL=2  use Seto's modification of Surkus' GPEF expansion in
c       z = [R**MPAR - Re**MPAR]/[a*R**MPAR + b*Re**MPAR] where
c       a=PARM(NVARB-1) & b=PARM(NVARB), which incorporates Dunham, SPF,
c       O-T and other forms: V(z) = c_0 z^2 [1 + c_1 z + c_2 z^2 + ...]
c       where  c_0 [cm-1] is read in as DSCM, and the first (NVARB-2)
c       PARM(i)'s are the  c_i  (i > 0).  [MPAR is dummy parameter here]
c  * For Dunham case:  NCN=1, PARM(NVARB-1)= 0.0, PARM(NVARB)= 1.0
c  * For SPF case:  NCN=1, PARM(NVARB-1)= 1.0, PARM(NVARB)= 0.0
c  * For Ogilvie-Tipping:  NCN=1, PARM(NVARB-1)= 0.5 = PARM(NVARB)
c  * NOTE that for Surkus MPAR < 0 case:  z(MPAR,a,b)= z(|MPAR|,-b,-a)
c      Generate & return the  D_e  value implied by these coefficients.
c** If IPOTL=3  generate a Morse or Extended Morse Oscillator potential
c      with exponent factor "beta" defined as a power series of order
c      (NVARB-1) in  y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR)
c      with NVARB coefficients PARM(i).    [!! MPAR .ge.1 !!]
c    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
c*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
c      modified Morse function with  b= PARM(1)  and C= PARM(2).
c*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
c      Morse Oscillator" potential with exponent expansion in (R-Re)]
c ...  otherwise, set  MPAR.ge.0
c** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
c      If MPAR > 0  exponent parameter defined in terms of a polynomial
c           of order (NVARB-2) with the NVARB coefficients  PARM(j).
c           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR),
c           and long-range defined by NCN inverse-power terms
c      If MPAR = 0  exponent polynomial variable is  2*y_{1}= y{O-T}
c      If MPAR < 0  exponent polynomial variable is  y_{|MPAR|}
c      If MPAR.le.0  exponent polynomial connected to limiting inverse-
c           power potential exponent by exponential switching function
c           with parameters  Asw= PARM(NVARB-1)  and  RSW= PARM(NVARB).
c** If IPOTL=5  generate a Double-Exponential Long-Range (DELR)
c       potential [JCP 119, 7398 (2003)] with additive long-range part
c       defined by a sum of NCMM damped inverse-power terms, & exponent
c       polynomial radial variable defined by parameter MPAR (=p)
c** If IPOTL=6  generate generalized HFD({m_i},i=1,NCMM) potential.
c       PARM(1-3) are the parameters defining the HFD damping function
c       D(x)=exp[-pparm(1)*(PARM(2)/x - 1)**PARM(3)] {for x < PARM(2)}
c       PARM(4) the quadratic coefficient in the exponent, and
c       PARM(5) is the power of  x=R/Req  multiplying the repulsive term
c              AREP*x**PARM(5) *exp[-beta*x - PARM(4)*x**2] ;
c** If IPOTL=7  use Tiemann-type polynomial potential attached to an
c     inverse-power long-range tail and an 1/R^{12} (or exponential)
c     inner wall.
c----------------------------------------------------------------------
c++     READ(5,*) IPOTL, MPAR, NSR, NCMM, NVARB, IBOB, DSCM, REQ
c++     IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
c++     IF((IPOTL.EQ.4).OR.(IPOT.EQ.7)) READ(5,*) (MMLR(I),CMM(I),I= 1,MPAR)
c++     IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
c++     IF(IBOB.GT.0) THEN
c++         READ(5,*) MN1R, MN2R, PAD, MAD, NU1, NU2, PNA, NT1, NT2
c++         IF(PAD.GT.0) THEN
c++              IF(NU1.GE.0) READ(5,*) U1INF, (U1(I), I=0,NU1)
c++              IF(NU2.GE.0) READ(5,*) U2INF, (U2(I), I=0,NU2)
c++              IF(NT1.GE.0) READ(5,*) T1INF, (T1(I), I=0,NT1)
c++              IF(NT2.GE.0) READ(5,*) T2INF, (T2(I), I=0,NT2)
c++           ELSE
c++              IF(NU1.GE.0) READ(5,*) U1INF,(U1(I), I=0,NU1), Aad1, Rad1
c++              IF(NU2.GE.0) READ(5,*) U2INF,(U2(I), I=0,NU2), Aad2, Rad2
c++              IF(NT1.GE.0) READ(5,*) T1INF,(T1(I), I=0,NT1), Ana1, Rna1
c++              IF(NT2.GE.0) READ(5,*) T2INF,(T2(I), I=0,NT2), Ana2, Rna2
c++           ENDIF
c++         ENDIF
c++     ENDIF
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PW= 2.D0
      IF((NCN1.GT.0).AND.(NCN1.NE.2)) PW= 2.D0*NCN1/(NCN1-2.D0)
!     IF(DFLOAT(NCN1).LT.(2.d0*PRV + 1.9999999d0)) THEN
      IF(DBLE(NCN1).LT.(2.d0*PRV + 1.9999999d0)) THEN
          WRITE(6,629) (2.d0*PRV + 2.),NCN1
  629 FORMAT(/  ' *** Note that Radial variable power \alpha optimal for
     1   NLR=',f5.2,' > NCN=', i2)
      ENDIF
c** Convert potential in [cm-1] to form appropriate for SCHRQas
      DO  I= 1,NPP
          V1BZ(I)= V1(I)*BFCT
          V1(I)= V1BZ(I)*DRDY2(I) + FAS(I)
          V2(I)= V1(I)
          RM2(I)= RRM2(I)*DRDY2(I)
      ENDDO
      VLIM2= VLIM1
!OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
!     WRITE(6,*) 'V(R) after converting into form for schrq.f:'
!     DO I=1,3
!     !WRITE(6,*) 'V(',I,')=',V1(I)
!     !WRITE(6,*) 'FAS(',I,')=',FAS(I)
!      WRITE(6,*) 'DRDY2(',I,')=',DRDY2(I)
!     !WRITE(6,*) 'RRM2(',I,')=',RRM2(I)
!     ENDDO
!     WRITE(6,*) 'V(                 20000)=',V1(20000)
!     WRITE(6,*) 'V(',NPP,')=',V1(NPP)
      IF(NUMPOT.LE.1) THEN
          WRITE(6,636)
          IOMEG2= IOMEG1
      ELSE
!         WRITE(6,635)
          WRITE(6,*) ''
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** For 2-potential Franck-Condon factor calculation, get the second
c  potential in this second call to PREPOT (uses the same parameter
c  reading sequence so exhaustively described immediately above).
!          CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,
!     1                                             VLIM2,V2,CNN2,NCN2)
           CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,
     1  VLIM2,V2,CNN2,NCN2,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,
     2  PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)

c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Convert potential (in (cm-1)) to form appropriate for SCHRQas
          DO  I=1,NPP
              V2BZ(I)= V2(I)*BFCT
              V2(I)= V2BZ(I)*DRDY2(I) + FAS(I)
              RM22(I)= RRM22(I)*DRDY2(I)
          ENDDO
      ENDIF
c
c** NLEV1 is the no. of levels {v=IV(i), J=IJ(i)} of potential-1 which
c   we wish to find.
c* IF(NLEV1=0) calculate (and print?) potential, and then quit.
c* If read-in value of  NLEV1 < 0 , and program (attempts to) find all
c  vibrational levels of potential-1 up to  v = |NLEV1|.  [This case
c  assumes  AUTO1 > 0.]
c** If (AUTO1.gt.0) read in only (v,J) quantum numbers of NLEV1 desired
c  level & subroutine ALFas tries to locate them (normal preferred case).
c   If (AUTO1.le.0) also read in trial energy for each level.  In this
c      case, the  NLEV.le.0  option does not work.
c   If (AUTO1.le.0) and vib. quant. No. IV < 0, seek level nearest to
c      given trial energy but for whatever q. number shows up
c** If(LCDC.gt.0) calculate centrifugal distortion constants for each
c  level via the Tellinghuisen implementation of Hutson's method.
c  If(LCDC.GE.2) also write them compactly to channel-9.
c** IF(LXPCT=0) calculate no expectation values or matrix elements.
c* IF(LXPCT = -1) only calculate and write compactly to channel-7 the
c  eigenvalues and level widths.
c* IF(LXPCT= 1,2 or -2) calculate expectation values, or if |LXPCT| > 2
c  the off-diagonal matrix elements, of powers of the distance
c  coordinate or radial function defined by parameters IRFN & DREF.
c* For   LXPCT > 0  write all these results to channel-6;  otherwise
c                  supress most such printing to channel-6.
c* For  |LXPCT| = 2  write eigenvalues and expectation values in
c                       compact form on channel-7.
c* For  |LXPCT| > 2  calculate matrix elements coupling each level
c  to all (up to VIBMX) preceeding levels of the same potential (for
c  NUMPOT.le.1), or to NLEV2 (see below) vib. levels of potential-2
c  (for NUMPOT.ge.2), and if (|LXPCT| > 3) write the overall
c  off-diagonal matrix elements on channel-8.
c* For  |LXPCT| > 4  also write to channel-7 the matrix elements of the
c  individual powers of the chosen distance coordinate (or radial fx.)
c* For  |LXPCT| > 5  WRITE(7,xx) only those matrix element components.
c** IF(NJM > 0), for each vibrational level, calculate all rotational
c  levels up to  J=NJM  or predissociation, whichever comes first.
c  Note that  AUTO1.le.0  forces  NJM= 0
c** When (NJM.GT.0) increase J in increments of JDJR.
c** IF(IWR.NE.0) print error & warning descriptions
c  IF(IWR.GE.1) also print final eigenvalues & node count.
c  IF(IWR.GE.2) also show end-of-range wave function amplitudes
c  IF(IWR.GE.3) print also intermediate trial eigenvalues, etc.
c** IF(LPRWF.GT.0) print wave function every LPRWF-th  point.
c** IF(LPRWF.LT.0) compactly write to channel-7 every |LPRWF|-th
c  wave function value.  **  A lead "card" identifies the level, gives
c  the position of 1-st point and radial mesh, & states No. of  points
c=======================================================================
c** INNOD1 specified wave fx. initiation at RMIN.  Normal case of
c  INNOD1 > 0  gives initiation with wave fx. node @ RMIN.
c  INNOD1.le.0  give initiation with  zero slope @ RMIN.  This determines
c    symmetric eigenfunctions for rare special case when input potential
c    is half of a precisely symmetric potential with mid-point at RMIN.
c-----------------------------------------------------------------------
c ... comment this out if use first version of following READ statement
      INNOD1= 1
      INNOD2= INNOD1
c-----------------------------------------------------------------------
c     READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF, INNOD1
c-----------------------------------------------------------------------
!     READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF
c-----------------------------------------------------------------------
c** SINNER specifies whether wave function matching occurs at outermost
c  (SINNER.le.0) or innermost well turning point, to facilitate finding
c  inner vs. outer wells of a double well potential; Normally controlled
c  automatically,
c!!
      IF(LPRWF.LT.0) WRITE(10,605) TITL
c!!
      SINNER= 0
      INNER= SINNER
      IF(INNOD1.GT.0) WRITE(6,686) 1
      IF(INNOD1.LE.0) WRITE(6,688) 1
      IF(JDJR.LE.0) JDJR=1
      WRITE(6,612) EPS
      NLEV= NLEV1
      IF(NLEV1.LE.0) NLEV= 1
      SOMEG1= IOMEG1**2
      IF(IOMEG1.GE.99) THEN
          WRITE(6,609)
        ELSE
          IF(IOMEG1.GE.0) WRITE(6,608) 1,IOMEG1,SOMEG1
          IF(IOMEG1.LT.0) WRITE(6,6085) 1,IOMEG1,-IOMEG1
        ENDIF
      VMAX1= 0
c** Read the vibrational & rotational quantum numbers IV(i) & IJ(i) [and
c  if AUTO1.le.0 also trial energy GV(I)] of the NLEV levels to be found
c** For  IV(i)  values < -10,  SCHRQ  imposes a hard wall boundary
c  condition (i.e., a node) at mesh point # |-IV(i)| .
c-----------------------------------------------------------------------
!     IF(AUTO1.GT.0) READ(5,*) (IV(I), IJ(I), I= 1,NLEV)
!     IF(AUTO1.LE.0) READ(5,*) (IV(I), IJ(I), GV(I), I= 1,NLEV)
c-----------------------------------------------------------------------
! IJ(i) for each level i should be read from the input file but
! otherwise initialize them explicitly:
      DO I= 1,NLEV
       IJ(I)=0
      ENDDO
      IF(NLEV1.GT.0) THEN
          IF(AUTO1.GT.0) WRITE(6,607) NLEV,(IV(I),IJ(I),I=1,NLEV)
          IF(AUTO1.LE.0) THEN
              WRITE(6,6607) NLEV,(IV(I),IJ(I),GV(I),I=1,NLEV)
              DO  I= 1,NLEV1
                  IF(IV(I).LT.VIBMX) ZK1(IV(I),0)= GV(I)
                  ENDDO
              ENDIF
          DO I= 1,NLEV
              IF(IV(I).LT.VIBMX) VMAX1= MAX(VMAX1,IV(I))
              ENDDO
          JREF= 0
        ELSE
          IF(NLEV1.LT.(1-VIBMX)) NLEV1= 1- VIBMX
          VMAX1= -NLEV1
          NLEV= VMAX1+ 1
          WRITE(6,625) IJ(1),NLEV
          JREF= IJ(1)
          DO  I= 1,NLEV
              IV(I)= I-1
              IJ(I)= JREF
              ENDDO
        ENDIF
      IF(NJM.GT.IJ(1)) WRITE(6,638) JDJR,NJM
!     IF(LCDC.GT.1)  WRITE(9,901) TITL
      IF(LXPCT.EQ.-1) WRITE(7,723) TITL
c** MORDR is the highest power of the radial function (or distance
c  coordinate whose expectation values or matrix elements are to be
c  calculated.  Program currently dimensioned for (MORDR.LE.10).  To
c  calculate only F-C Factors (when LXPCT>2), set  MORDR = -1.
c** IRFN & DREF specify the definition of the radial function or
c  distance coordinate  RFN(R), powers of which are averaged over in
c  expectation value or matrix element calculations.
c* If(IRFN .le. -10) utilize the USER-CHOSEN and CODED radial function
c                 generated in Lines #500-504 (below)
c* If(IRFN = -4)  the function is a power series in  R  premultiplying a
c          first derivative operator acting on the wavefx of Potential-2
c* If(IRFN = -3)  the function is the inverse power   1/R**3
c* If(IRFN = -2)  the function is the inverse power   1/R**2
c* If(IRFN = -1)  the function is the Dunham coordinate  X=(R-DREF)/DREF
c* If(IRFN =  0)  the function  RFN(R)  is the distance  R  itself.
c* If(IRFN = 1-9)  use the Surkus-type variable
c                 X=(R^p - DREF^p)/(R^p + DREF^p)  where  p= IRFN
c* For  IRFN = -1 or 1-9,  if  DREF.gt.0  the read-in DREF value is the
c  reference length used to define the distance coordinate, while
c  if  DREF.le.0  determine the value of this reference length by
c  requiring that the expectation value  X**1  of the distance
c  coordinate for the first level considered be identically zero.
c* IF(IRFN > 10) define  RFN(R)   by reading in, interpolating over (and
c  extrapolating beyond) IRFN read-in values of some known radial
c  (transition moment) function, whose asymptotic value is DREF.  Do
c  this in using the same read statements and GENINT subroutine calls
c  used for generating a numerical potential.
      IF((LXPCT.NE.0).AND.(LXPCT.NE.-1)) THEN
c-----------------------------------------------------------------------
!         READ(5,*) MORDR, IRFN, DREF
c-----------------------------------------------------------------------
          IF(MORDR.GT.MORDRMX) MORDR= MORDRMX
          IF(IABS(LXPCT).EQ.2) WRITE(7,724) TITL,MORDR
!         IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5)) WRITE(8,824) TITL
          IF(IABS(LXPCT).GE.5) WRITE(7,725) TITL,MORDR
          IF(IABS(IRFN).GE.10) THEN
              MORDR= 1
              DM(0)= 0.d0
              DM(1)= 1.d0
            ELSE
              IF(MORDR.GE.0) THEN
c** Overall calculated matrix elements are for a power series in the
c  radial function  RFN(i)  (specified by IRFN & DREF), so must read
c  coefficients  DM(J)  of this power series.
c-----------------------------------------------------------------------
!                 READ(5,*) (DM(J), J= 0,MORDR)
c-----------------------------------------------------------------------
                ELSE
                  DO  I= 1,NPP
                      RFN(I)= 1.D0
                      ENDDO
                  IF(MORDR.LT.0) WRITE(6,617)
                ENDIF
            ENDIF
c** Define radial function (distance coordinate) operator  RFN(R)  for
c  expectation values or matrix elements.
c** First ... for matrix elements of an operator consisting of a power
c    series in  R  premultiplying the radial derivative of the wavefx.
          IF(IRFN.EQ.-4) WRITE(6,650) MORDR
          IF(MORDR.GT.0) THEN
c** If  RFN(R)  is the distance itself ...
              IF(IRFN.EQ.0) THEN
                  WRITE(6,614)
                  DO  I= 1,NPP
                      RFN(I)= RVB(I)
                      ENDDO
                  ENDIF
              IF((IRFN.EQ.-2).OR.(IRFN.EQ.-3)) THEN
                  IF((IRFN.EQ.0).OR.(IRFN.EQ.-2).OR.(IRFN.EQ.-3))
     1                                                      DREF= 0.D0
c** If  RFN(R)  is   1/(distance)**|IRFN|  ....
                  J= -IRFN
                  WRITE(6,616) -IRFN
                  DO  I= 1,NPP
                      RFN(I)= 1.d0/RVB(I)**J
                      ENDDO
                  ENDIF
c%% Any other user-defined matrix element argument radial function
c   may be introduced to the code here, and invoked by:  IRFN= -4
c  Note that the existing  RVB(i)  array is the radial distances  R .
              IF(IRFN.LE.-10) THEN
c&& Illustrative user-defined analysis RFN(R) function
c&&               WRITE(6,*) 'Print description of function introduced'
c&&               WRITE(6,*) 'Use Freedman Pade DMF for CO'
c&&               DO  I= 1,NPP
c%%---------------------------------------------------------------------
c%%                  RFN(I)= {calculate users chosen radial function}
c%% Freedman's DMF for CO  ---------------------------------------------
c%%                  RFN(I)= {calculate users chosen radial function}
c&&   data coeff_new /-24.6005858d0,-109.5939637d0,-524.8233323d0,
c&&  +                 4.5194090d0,19.7954955d0,
c&&  +                 6.6011985d0,19.7206690d0/
c&&   dm = -0.122706d0*(1.+coeff(1)*x+coeff(2)*x*x+coeff(3)*x**3)/
c&&  +                   (1.+coeff(4)*x+coeff(5)*x*x+coeff(6)*x**3
c&&  +                    + coeff(7)*x**6)
c&&---------------------------------------------------------------------
c&&                   XX= RFN(I)/1.128322714d0 - 1.d0
c&&                   RFN(I)= -0.122706d0*(1.d0+ XX*(-24.6005858d0
c&&  1                  + XX*(-109.5939637d0 + XX*(-524.8233323d0))))/
c&&  2                    (1.d0 + XX*(4.5194090d0 + XX*(19.7954955d0 +
c&&  3                        XX*(6.6011985d0 + 19.7206690d0*XX**3))))
c&&                   ENDDO
                  ENDIF
              IF((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9))) THEN
c** If  RFN(R)  is the Dunham or Surkus-type distance coordinate
                  IF(IRFN.EQ.-1) WRITE(6,615)
                  IF((IRFN.GE.1).AND.(IRFN.LE.9)) WRITE(6,611)
     1                                             IRFN,IRFN,IRFN,IRFN
                  IF(DREF.GT.0.D0) THEN
                      DREFP= DREF**IRFN
                      WRITE(6,613) DREF
                      DO  I=1,NPP
                          XX= YMINN+I*YH
                          IF(IRFN.EQ.-1) RFN(I)= (RVB(I)- DREF)/DREF
                          IF(IRFN.GE.1) RFN(I)= (RVB(I)**IRFN- DREFP)
     1                                        /(RVB(I)**IRFN + DREFP)
                          ENDDO
                    ELSE
                      WRITE(6,610)
                    ENDIF
                  ENDIF
c
cQQQQQQQQQQQ new ... the following not fully tested ...
c
c** If  RFN(R)  is defined by interpolating over read-in points, use
c  potential generating routine to do interpolation/extrapolation.
              IF(IRFN.GE.10) THEN
                  MORDR= 1
                  DM(0)= 0.d0
                  DM(1)= 1.d0
                  WRITE(6,603)
c** If the expectation value/matrix element radial function argument to
c   be defined by interpolating/extrapolating over read-in points, then
c   read input analogous to that for a pointwise potential, and then call
c   interpolation/extrapolation routine GENINT (from PREPOT package)
c* NRFN is the number of points [XIF(i),YIF(i)] to be read in
c* RFLIM  is the limiting asymptotic value imposed on the extrapolation
c* Interpolate with NUSEF-point piecewise polynomials (or splines for
c    NUSEF.le.0), which are extrapolated to the asymptote as specified by
c    parameters ILRF, NCNF & CNNF (see read #20).
c* RFACTF - factor converts read-in distances XIF(i) to angstroms
c* MFACTF - factor converts read-in moment values YIF(i) to debye.
c-----------------------------------------------------------------------
!                 READ(5,*) NRFN, RFLIM
!                 READ(5,*) NUSEF, ILRF, NCNF, CNNF
!                 READ(5,*) RFACTF, MFACTF
!                 READ(5,*) (XIF(I), YIF(I), I= 1,NRFN)
! If you uncomment the above, you better also uncomment the
! initialization to 0 below.
c-----------------------------------------------------------------------
                  MFACTF=0
                  RFACTF=0
                  WRITE(6,810) NRFN, RFLIM
                  IF(NUSEF.GT.0) WRITE(6,812) NUSEF, NRFN
                  IF(NUSEF.LE.0) WRITE(6,814) NRFN
                  IF((ILRF.GT.1).AND.(DABS(CNNF).GT.0.D0))
     1                                         WRITE(6,816) CNNF, NCNF
                  WRITE(6,818) RFACTF, MFACTF
                  NROW= (NRFN+ 2)/3
                  DO  J= 1,NROW
                      WRITE(6,820) (XIF(I), YIF(I), I=J, NRFN, NROW)
                  ENDDO
                  DO  I= 1,NRFN
                      XIF(I)= XIF(I)*RFACTF
                      YIF(I)= YIF(I)*MFACTF
              ENDDO
  810 FORMAT(' Transition moment function defined by interpolating over'
     1  ,I4,' read-in points'/5x,'and approaching the asymptotic value',
     2  f12.6)
  812 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov
     1er',I5,' input points' )
  814 FORMAT(' Perform cubic spline interpolation over the',I5,' input p
     1oints' )
  816 FORMAT('- Beyond read-in points extrapolate to limiting asymptotic
     1 behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  818 FORMAT(' Scale input points:  (distance)*',1PD16.9,'   &  (moment)
     1*',D16.9/4x,'to get required units  [Angstroms & debye]'/
     3  3('      R(i)         Y(i)  ')/3(3X,11('--')))
  820 FORMAT((3(F12.6,F13.6)))
                  IR2F= 0
                  CALL GENINT(LRPT,NPP,RVB,RFN,NUSEF,IR2F,NRFN,XIF,YIF,
     1                                           RFLIM,ILRF,NCNF,CNNF)
                  ENDIF
              ENDIF
          IF((MORDR.GE.0).AND.(IABS(IRFN).LE.9))
     1                                  WRITE(6,602) (DM(J),J=0,MORDR)
          ENDIF
c** For matrix element calculation, couple each level of potential-1 to
c  up to (see below) NLEV2 other vibrational levels, subject to
c  rotational selection rules:  DELTA(J)= J2DL to J2DU with increment
c  J2DD (e.g., -1,+1,+2 for P- and R-branches).
c** If (AUTO2.gt.0) read in only (v,J) quantum numbers of desired levels
c  and trust subroutine ALFas to locate them (normal preferred case).
c   If (AUTO2.le.0) also read in a trial pure vib energy for each level.
c* For the one-potential case (NUMPOT.LE.1), automatically truncate to
c  avoid redundancy and consider only emission into these NLEV2 levels.
c* Trial level energies are generated internally.
c**  IV2(i) are the vibrational quantum numbers of the Potential-2
c  levels for which matrix elements are desired.
c**  ZK(IV2(i),0) are the associated pure vibrational trial energies
c  (which are only read in if AUTO2.le.0!)
c=======================================================================
c** INNOD2 specified wave fx. initiation at RMIN.  Normal case of
c  INNOD2 > 0  gives initiation with wave fx. node @ RMIN.
c  INNOD2.le.0  give initiation with  zero slope @ RMIN.  This determines
c    symmetric eigenfunctions for rare special case when input potential
c    is half of a precisely symmetric potential with mid-point at RMIN.
c=======================================================================
      IF(IABS(LXPCT).GE.3) THEN
c-----------------------------------------------------------------------
ccc       READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD, INNOD2
c-----------------------------------------------------------------------
!         READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD
c-----------------------------------------------------------------------
          IF(NLEV2.GT.VIBMX) NLEV2= VIBMX
          IF(NLEV2.LE.0) THEN
              WRITE(6,644) NLEV2
              GO TO 997
!             STOP
          ENDIF
c----------------------------------------------------------------------
!         IF(AUTO2.GT.0) READ(5,*) (IV2(I), I= 1,NLEV2)
          IF(AUTO2.LE.0) THEN
!             READ(5,*) (IV2(I), ZK2(I,1), I= 1,NLEV2)
c----------------------------------------------------------------------
c** Give potential-2 trial energy the correct vibrational label
              DO  I= 1,NLEV2
                  ZK2(IV2(I),0)= ZK2(I,1)
              ENDDO
          ENDIF
          IF(NUMPOT.GT.1) THEN
              IF(INNOD2.GT.0) WRITE(6,686) 2
              IF(INNOD2.LE.0) WRITE(6,688) 2
          ENDIF
          VMAX2= 0
          DO  ILEV2= 1,NLEV2
              VMAX2= MAX(VMAX2,IV2(ILEV2))
          ENDDO
          IF(MORDR.LT.0) DM(1)= 1.d0
          SOMEG2= IOMEG2**2
          IF(J2DD.EQ.0) J2DD= 1
          IF(AUTO2.GT.0) WRITE(6,634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),
     1                                                     I= 1,NLEV2)
          IF(AUTO2.LE.0) WRITE(6,6634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),
     1                                      ZK2(IV2(I),0), I= 1,NLEV2)
          IF(NUMPOT.GE.2) THEN
              IF(IOMEG2.GE.99) THEN
                  WRITE(6,609)
              ELSE
                  IF(IOMEG2.GE.0) WRITE(6,608) 2,IOMEG2,SOMEG2
                  IF(IOMEG2.LT.0) WRITE(6,6085) 2,IOMEG2,-IOMEG2
              ENDIF
          ENDIF
          WRITE(6,632)
      ENDIF
c
      IF(AUTO1.GT.0) THEN
c** If using automatic search for desired levels, subroutine ALFas gets
c  eigenvalues ZK1(v,0) for desired vibrational levels of Potential-1,
c  centrifugally-distorted to J=JREF.
          EJREF= JREF*(JREF+1)*YH**2
c** Replace  [J(J+1)] by  [J(J+1) + |IOMEG1|]  for Li2(A) and like cases.
!         IF(IOMEG1.LT.0) EJREF= EJREF - DFLOAT(IOMEG1)*YH**2
          IF(IOMEG1.LT.0) EJREF= EJREF - DBLE(IOMEG1)*YH**2
          DO  I= 1,NPP
              VBZ(I)= V1BZ(I) + EJREF*RRM2(I)
              VJ(I)= V1(I) + EJREF*RM2(I)
          ENDDO
! OPTIONALLY WRITE SOME VALUES WHEN DEBUGGING:
!         WRITE(6,*) 'EJREF=',EJREF
!         DO I=1,3
!             WRITE(6,*) 'VJ=',VJ(I)
!             WRITE(6,*) 'V1',V1(I)
!             WRITE(6,*) 'RM2=',RM2(I)
!         ENDDO
          IF((NLEV1.EQ.1).AND.(IV(1).gt.998)) THEN
c** Option to search for very highest level (within 0.0001 cm-1 of Disoc)
              EO= VLIM1- 0.0001d0
              KV= IV(1)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f (1)'
              WRITE(6,*) 'Entering schrq.f (1)'
              WRITE(6,*) ''
              CALL SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,
     1      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              IV(1)= KV
              IF(KV.GE.0) THEN
                  WRITE(6,622) IJ(1),KV,VLIM1-EO
                  GV(KV)= EO
                  VMAX1= KV
              ELSE
                  WRITE(6,626) J, 0.001d0
                  GO TO 2
              ENDIF
          ELSE
              VMAX= VMAX1
              AFLAG= JREF
              IF((IABS(LXPCT).GT.2).AND.(NUMPOT.EQ.1)) VMAX=
     1                                                MAX(VMAX1,VMAX2)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f'
              WRITE(6,*) 'Entering alf.f'
              WRITE(6,*) ''
!             WRITE(6,*) 'NDP=',NPP
!             WRITE(6,*) 'YMIN=',YMIN
!             WRITE(6,*) 'YH=',YH
!             WRITE(6,*) 'NCN=',NCN1
!             DO I=1,3
!              WRITE(6,*) 'V=',VJ(I)
!              WRITE(6,*) 'SWF=',WF1(I)
!              WRITE(6,*) 'GV=',GV(I)
!             WRITE(6,*) 'INNR=',INNR1(I)
!             ENDDO
!             WRITE(6,*) 'VLIM=',VLIM1
!             WRITE(6,*) 'KVMAX=',VMAX
!             WRITE(6,*) 'AFLAG=',AFLAG
!             WRITE(6,*) 'ZMU=',ZMU
!             WRITE(6,*) 'EPS=',EPS
!             WRITE(6,*) 'BFCT=',BFCT
!             WRITE(6,*) 'INNODE=',INNOD1
!             WRITE(6,*) 'IWR=',IWR
!             WRITE(6,*) ''
!             DEALLOCATE(RVB)
              CALL ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,VMAX,AFLAG,ZMU,
     1                               EPS,GV,BFCT,INNOD1,INNR1,IWR)
              VMAX1= VMAX
          ENDIF
c** Get band constants for v=0-VMAX1 for generating trial eigenvalues
          WARN=  0
          DO  ILEV1= 0,VMAX
              KV= ILEV1
              EO= GV(KV)
              INNER= INNR1(KV)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f (2)'
              WRITE(6,*) 'Entering schrq.f (2)'
              WRITE(6,*) ''
              WRITE(6,*) 'Getting band constants for v=',KV
              CALL SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,
     1     WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,WARN,LPRWF)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f'
              WRITE(6,*) 'Entering cdjoel.f'
              WRITE(6,*) ''
! OPTIONALLY WRITE THE INPUT PARAMETERS FOR DEBUGGING:
!             WRITE(6,*) 'EO=',EO
!             WRITE(6,*) 'NBEG=',NBEG
!             WRITE(6,*) 'NEND=',NEND
!             WRITE(6,*) 'BvWN=',BvWN
!             WRITE(6,*) 'YH=',YH
!             WRITE(6,*) 'WARN=',WARN
!             WRITE(6,*) 'VJ(1)=',VJ(1)
!             WRITE(6,*) 'WF1(1)=',WF1(1)
!             WRITE(6,*) 'RM2(1)=',RM2(1)
!             WRITE(6,*) 'RCNST(1)=',RCNST(1)
! For a-state test case, there's a memory error before CALL CDJOELas for v=2
!             CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,RM2,
!    1                                                          RCNST)
              IF(NLEV1.LT.0) THEN
                  IV(ILEV1+1)= KV
                  IJ(ILEV1+1)= JREF
                  ENDIF
              ZK1(ILEV1,0)= GV(ILEV1)
              DO  M= 1,7
                  ZK1(ILEV1,M)= RCNST(M)
                  ENDDO
              ENDDO
          ENDIF
      IF(IABS(LXPCT).GT.2) THEN
          IF(AUTO2.GT.0) THEN
c** If using automatic location for levels of potential-2 (AUTO2 > 0)
c  for matrix element calculation, also need Potential-2 band constants
c  (rotational energy derivatives) ... again, calculate them at J=JREF
              IF(NUMPOT.GT.1) THEN
                  AFLAG= JREF
                  DO  I= 1,NPP
                      VBZ(i)= V2BZ(I) + EJREF*RRM22(I)
                      VJ(I)= V2(I) + EJREF*RM22(I)
                      ENDDO
                  CALL ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,VMAX2,AFLAG,
     1                            ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
                  ENDIF
              ENDIF
          DO  ILEV2= 1,NLEV2
              IF(NUMPOT.EQ.1) THEN
c** For matrix elements within a single potl., copy above band constants
                  DO  M= 0,7
                      ZK2(IV2(ILEV2),M)= ZK1(IV2(ILEV2),M)
                      ENDDO
                ELSE
c ... otherwise, generate them (as above) with SCHRQ & CDJOEL
                  KV= IV2(ILEV2)
                  IF(AUTO2.GT.0) EO= GV(KV)
                  IF(AUTO2.LE.0) EO= ZK2(KV,0)
                  INNER= INNR2(KV)
                  WRITE(6,*) ''
                  WRITE(6,*) 'Exiting level.f (3)'
                  WRITE(6,*) 'Entering schrq.f (3)'
                  WRITE(6,*) ''
                  CALL SCHRQas(KV,JREF,EO,GAMA,PMAX2,VLIM2,VJ,
     1     WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD2,INNER,WARN,LPRWF)
                  CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,
     1                                                      RM2,RCNST)
                  ZK2(IV2(ILEV2),0)= EO
                  DO  M= 1,7
                      ZK2(IV2(ILEV2),M)= RCNST(M)
                      ENDDO
                ENDIF
              ENDDO
          ENDIF
      WARN= 1
      EJREF= EJREF/YH**2
      IF(NLEV1.LE.0) NLEV= VMAX1+1
c
c===== Begin Actual Potential-1 Eigenvalue Calculation Loop Here =======
c** Loop to compute eigenvalues ... etc. for NLEV levels of Potential-1
      DO 190 ILEV1= 1,NLEV
          KV= IV(ILEV1)
          IF(KV.LT.0) EXIT
          NJMM= MAX(NJM,IJ(ILEV1))
          JROT= IJ(ILEV1)- JDJR
          IQT= 0
          JCT= 0
c** If NJM > IJ(ILEV1) loop over range of rotational levels too
          DO  JLEV= IJ(ILEV1),NJMM,JDJR
              JROT= JROT+ JDJR
              EJ= JROT*(JROT+1) - SOMEG1
              IF(IOMEG1.GE.99) EJ= JROT*JROT - 0.25D0
c** If   IOMEG < 0   centrifugal term is  [J(J+1) + |IOMEG|]
!             IF(IOMEG1.LT.0) EJ= JROT*(JROT+1) - DFLOAT(IOMEG1)
              IF(IOMEG1.LT.0) EJ= JROT*(JROT+1) - DBLE(IOMEG1)
c** If appropriate (AUTO1>0) use ALFas results to generate trial eigenvalue
              IF(AUTO1.GT.0) THEN
                  EO= ZK1(KV,0)
                  DEJ= EJ- EJREF
                  EJP= 1.d0
                  DO M= 1,7
                      EJP= EJP*DEJ
                      EO= EO+ EJP*ZK1(KV,M)
                      ENDDO
                ELSE
c... otherwise - use read-in trial energy
                  IF(IV(ILEV1).LT.VIBMX) EO= ZK1(IV(ILEV1),0)
                  IF(IV(ILEV1).GE.VIBMX) EO= GV(ILEV1)
                ENDIF
              IF((AUTO1.LE.0).AND.(DABS(ZK1(IV(ILEV1),0)).LE.0.d0)) THEN
                  CALL SCATTLEN(JROT,SL,VLIM1,V1,WF1,BFCT,YMIN,YH,NPP,
     1                 CNN1,NCN1,IWR,LPRWF)
                  IF(NUMPOT.EQ.1) GOTO 2
                  GOTO 104
                  ENDIF
c ... or if JLEV > IJ(ILEV1) ... use local Beff to estimate next level
              IF(JLEV.GT.IJ(ILEV1)) THEN
                  BEFF= 0.d0
                  DO  I= NBEG,NEND
                      BEFF= BEFF+ WF1(I)**2*RM2(I)
                      ENDDO
                  BEFF= BEFF*YH*BvWN
                  EO=  ESLJ(JCT)+ (2*JLEV+ 1- JDJR)*JDJR*BEFF
                  ENDIF
c** Now add centrifugal term to get effective (radial) potential
              EJ= EJ*YH**2
              DO  J= 1,NPP
                  VBZ(J)= V1BZ(J) + EJ*RRM2(J)
                  VJ(J)= V1(J) + EJ*RM2(J)
                  ENDDO
c** Set wall outer boundary condition, if specified by input IV(ILEV1)
              IF(KV.LT.-10) THEN
                  WF1(-IV(ILEV1))= 0.D0
                  WF1(-IV(ILEV1)-1)= -1.D0
                  ENDIF
              KVIN= KV
              IF(AUTO1.GT.0) INNER= INNR1(KV)
              IF(SINNER.NE.0) INNER= SINNER
c** Call SCHRQ to find Potential-1 eigenvalue EO and eigenfn. WF1(i)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f (4)'
              WRITE(6,*) 'Entering schrq.f (4)'
              WRITE(6,*) ''
! The next CALL SCHRQas for a-state failes during v=0, but EO's from when
! Entering schrq.f (2) were already accurate, so the next step is not necessary.
  100         CALL SCHRQas(KV,JROT,EO,GAMA,PMAX1,VLIM1,VJ,
     1      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              IF(KV.LT.0) THEN
c** SCHRQ  error condition is  (KV.LT.0) .
                  IF(NJM.GT.IJ(ILEV1)) THEN
c ... in automatic search for ever-higher J levels
                      IF(IQT.LE.0) THEN
c ... try one more time with E(trial) slightly below barrier maximum
                          IQT= 1
                          EO= PMAX1- 0.1d0
                          GO TO 100
                      ELSE
                          KV= KVIN
                          GO TO 130
                      ENDIF
                  ENDIF
                  GO TO 122
              ENDIF
!             WRITE(6,*) 'Hello! (for debugging) 1'
              IF((KV.NE.KVIN).AND.
     1                        ((AUTO1.GT.0))) THEN
c** If got wrong vib level, do a brute force ALFas calculation to find it.
                  KV= KVIN
                  AFLAG= JROT
                  CALL ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,KV,AFLAG,ZMU,
     1                                EPS,GV,BFCT,INNOD1,INNR1,IWR)
                  IF(KV.EQ.KVIN) THEN
                      EO= GV(KVIN)
                      GO TO 100
                  ELSE
                      WRITE(6,618) KVIN,JROT,KV
                      KV= KVIN
                      GO TO 130
                  ENDIF
              ENDIF
!             WRITE(6,*) 'Hello! (for debugging) 2'
              IF(KV.NE.IV(ILEV1)) IV(ILEV1)= KV
c** If desired, calculate rotational & centrifugal distortion constants
              IF(LCDC.GT.0) THEN
                  IF((IOMEG1.GT.0).AND.(JROT.EQ.0)) THEN
c** Calculate 'true' rotational constants for rotationless IOMEG>0 case
                      WRITE(6,*) ''
                      WRITE(6,*) 'Exiting level.f (5)'
                      WRITE(6,*) 'Entering schrq.f (5)'
                      WRITE(6,*) ''
                      CALL SCHRQas(KV,0,EO,GAMA,PMAX1,VLIM1,V1,
     1      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
                      CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V1,
     1                                                  WF1,RM2,RCNST)
                  ELSE
c** Calculate rotational constants for actual (v,J) level of interest.
                      CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,
     1                                                  WF1,RM2,RCNST)
                  ENDIF
                  IF(DABS(EO).GT.1.d0) THEN
                      WRITE(6,606) KV,JROT,EO,(RCNST(M),M=1,7)
                  ELSE
                      WRITE(6,6065) KV,JROT,EO,(RCNST(M),M=1,7)
                  ENDIF
                  IF(LCDC.GT.1) THEN
                      IF(DABS(EO).GT.1.d0) THEN
!                         WRITE(9,902) KV,JROT,EO,(RCNST(M),M=1,7)
                      ELSE
!                         WRITE(9,904) KV,JROT,EO,(RCNST(M),M=1,7)
                      ENDIF
                  ENDIF
              ENDIF
!             WRITE(6,*) 'Hello! (for debugging) 3'
              IF(LXPCT.EQ.-1)  WRITE(7,703) KV,JROT,EO,GAMA
              IF(((LXPCT.EQ.1).OR.(IABS(LXPCT).EQ.2)).OR.
     1                  ((IABS(LXPCT).GT.2).AND.((IRFN.EQ.-1).OR.
     2          (IRFN.GE.1).AND.(IRFN.GE.9)).AND.(DREF.LE.0.d0))) THEN
c** Calculate various expectation values in LEVXPC
                  CALL LEVXPC(KV,JROT,EO,GAMA,NPP,WF1,RFN,VBZ,VLIM1,
     1                    YH,DREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
                  IF((LXPCT.GT.0).AND.(MORDR.GT.0)) WRITE(6,632)
              ENDIF
!             WRITE(6,*) 'Hello! (for debugging) 4'
  104         IF((IABS(LXPCT).LE.2).OR.(NLEV2.LE.0)) GO TO 122
c** If desired, now calculate off-diagonal matrix elements, either
c  between levels of different potentials, IF(NUMPOT.GE.2), or between
c  levels of a single potential, for (NUMPOT.LE.1).
c** First prepare centrifugally distorted potential, trial energy, etc.,
c  and calculate second wave function and matrix element(s)
              DO 120 ILEV2= 1,NLEV2
c** For case of a single potential, avoid redundancy by considering
c  only emission
                  IF((NUMPOT.LE.1).AND.(IV2(ILEV2).GT.KV)) GO TO 120
c** Loop over J2's allowed by given selection rule.
                  DO 116 IJD= J2DL,J2DU,J2DD
                      KV2= IV2(ILEV2)
                      KVIN= KV2
                      JROT2= JROT+IJD
                      IF(JROT2.LT.0) GO TO 116
                      IF((NUMPOT.LE.1).AND.(IV2(ILEV2).EQ.KV).AND.
     1                                      (JROT2.GT.JROT)) GO TO 116
                      EJ2= JROT2*(JROT2+1)- SOMEG2
                      IF(IOMEG2.GE.99) EJ2=JROT2**2-0.25D0
c... allow for weird Li2(A) and Li2(c) potential cases
!                     IF(IOMEG2.LT.0) EJ2=JROT2*(JROT2+1)-DFLOAT(IOMEG2)
                      IF(IOMEG2.LT.0) EJ2=JROT2*(JROT2+1)-DBLE(IOMEG2)
                      EO2= ZK2(KV2,0)
                      DEJ= EJ2- EJREF
                      EJP= 1.d0
c** Use calculated state-2 CDC's to predict trial eigenvalue
                      DO  M= 1,7
                          EJP= EJP*DEJ
                          EO2= EO2+ EJP*ZK2(KV2,M)
                      ENDDO
c** Now ... update to appropriate centrifugally distorted potential
                      EJ2= EJ2*YH*YH
                      DO  I=1,NPP
                          VBZ(I)= V2BZ(I)+ EJ2*RRM22(I)
                          VJ(I)= V2(I)+ EJ2*RM22(I)
                      ENDDO
                      INNER= INNR2(KV2)
                      IF(SINNER.NE.0) INNER= SINNER
                      ICOR= 0
                      WRITE(6,*) ''
                      WRITE(6,*) 'Exiting level.f (6)'
                      WRITE(6,*) 'Entering schrq.f (6)'
                      WRITE(6,*) ''
  110                 CALL SCHRQas(KV2,JROT2,EO2,GAMA,PMAX2,VLIM2,VJ,
     1    WF2,BFCT,EPS,YMIN,YH,NPP,NBEG2,NEND2,INNOD2,INNER,IWR,LPRWF)
                      IF(KV2.NE.KVIN) THEN
                          IF(KV2.LT.0) GO TO 114
c** Using CDC's to estimate trial eigenvalue failed:
                          ICOR= ICOR+1
                          IF(ICOR.LE.2) THEN
c ... first correction attempt ... use semiclassical dv/dE to improve
                              GB= -1.d0
                              GI= -1.d0
c ... hey RJ!!  shouldn't you update this using SCECOR ?
                              WV= 0.d0
                              XX= EO2*BFCT
                              DO  I= NBEG2,NEND2
                                  GBB= GB
                                  GB= GI
                                  GI= XX-VJ(I)
                                  IF((GBB.GT.0.d0).AND.(GI.GT.0.d0))
     1                                          WV= WV+ 1.d0/DSQRT(GB)
                              ENDDO
                              WV= 6.2832d0/(BFCT*WV)
                              EO2= EO2+ WV*(KVIN- KV2)
                              GO TO 110
                          ENDIF
                          WRITE(6,633) IV2(ILEV2),JROT2,KV2
c ... if that fails, do a brute force ALFas calculation to find it.
  114                     KV2= KVIN
                          AFLAG= JROT2
                          CALL ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,KV2,
     1                      AFLAG,ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
                          IF(KV2.EQ.KVIN) THEN
                              EO2= GV(KV2)
                              INNER= INNR2(KV2)
                              GO TO 110
                          ELSE
                              WRITE(6,618) KVIN,JROT,KV2
                              GO TO 116
                          ENDIF
                      ENDIF
                      IF(NBEG.GT.NBEG2) NBEG2= NBEG
                      IF(NEND.LT.NEND2) NEND2= NEND
c                     IF((NUMPOT.LE.1).AND.(EO2.GT.(EO+EPS))) GO TO 120
                      CALL MATXEL(KV,JROT,IOMEG1,EO,KV2,JROT2,IOMEG2,
     1 IRFN,EO2,NBEG2,NEND2,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
  116                 CONTINUE
c** End of Potential-2 rotational selection level loop
  120             CONTINUE
c++++ End of Potential-2 vibrational level matrix element loop +++++++++
c
  122         CONTINUE
              JCT= JCT+1
c ... check to avoid array overflow
              IF(JCT.GT.VIBMX) THEN
                  WRITE(6,637)  VIBMX
                  GO TO 997
!                 STOP
                  ENDIF
              JWR(JCT)= JROT
              ESLJ(JCT)= EO
              ENDDO
!             WRITE(6,*) 'Hello! (for debugging) 5'
c++ End of Potential-1 loop over NJM-specified J-sublevels
  130     IF(NJM.GT.IJ(ILEV1)) THEN
c** Print rotational sublevels generated for vibrational level  ILEV
              NROW=(JCT+4)/5
              WRITE(6,627) KV
              DO  J=1,NROW
                  WRITE(6,628) (JWR(I),ESLJ(I),I=J,JCT,NROW)
                  ENDDO
              WRITE(6,641)
              ENDIF
          ESOLN(ILEV1)= ESLJ(1)
  190     CONTINUE
c++ End of loop over the NLEV Potential-1 input levels
      WRITE(6,*) ''
      WRITE(6,*) 'SUMMARY (ALL ENERGIES IN CM-1):'
      IF(NLEV1.LT.0) THEN
          NROW=(NLEV+3)/4
          WRITE(6,623) NLEV,IJ(1)
          DO  J=1,NROW
              WRITE(6,630) (IV(I),ESOLN(I),I=J,NLEV,NROW)
          ENDDO
          IF((NLEV.GT.1).AND.(IJ(1).EQ.0).AND.(NCN1.GT.0)
     1                               .AND.(ESOLN(NLEV).LT.VLIM1)) THEN
c** In (NLEV1 < 0) option, estimate vD using the N-D theory result that:
c     (vD - v) {is proportional to} (binding energy)**((NCN-2)/(2*NCN))
              VDMV=1.D0/(((VLIM1-ESOLN(NLEV-1))/
     1                         (VLIM1-ESOLN(NLEV)))**(1.D0/PW) - 1.D0)
!          CALL ADD_INFO('LEVEL_CHKSUM',[CHKSUM],1,2)
           CALL ADD_INFO('Vibrational energy:',ESOLN,NLEV,2)
! IF YOU WANT TO ONLY VERIFY THE ENERGY OF THE HIGHEST LEVEL:
!          CALL ADD_INFO('Last vibrational energy:',ESOLN(NLEV),1,2)
! IF YOU WANT TO ONLY VERIFY THE NUMBER OF LEVELS:
! All values must be reals, so if you want to check an integer, e.g. the
! number of iterations, then you must convert this to a real:
!          CALL ADD_INFO('Numer of vibrational levels:',NLEV,1,2)
c** Use empirical N-D Expression to predict number and (if there are
c  any) energies of missing levels
              VD= IV(NLEV)+VDMV
              IVD= INT(VD)
              IF(IVD.GE.VIBMX) IVD= VIBMX-1
              IVS= IV(NLEV)+1
              WRITE(6,620) NCN1,VD
!             IF((IVD.GE.IVS).AND.(DFLOAT(IV(NLEV))/VD.GT.0.9d0)) THEN
              IF((IVD.GE.IVS).AND.(DBLE(IV(NLEV))/VD.GT.0.9d0)) THEN
                  NFP= NLEV+1
                  DO  I= IVS,IVD
                      NLEV= NLEV+1
                      IV(NLEV)= IV(NLEV-1)+1
                      ESOLN(NLEV)= VLIM1 - (VLIM1 - ESOLN(NLEV-1))*
     1                                          (1.D0 - 1.D0/VDMV)**PW
                      VDMV= VDMV-1.D0
                  ENDDO
                  NLP= NLEV-NFP+1
                  NROW= (NLP+3)/4
                  WRITE(6,621) NLP
                  DO  J= 1,NROW
                      III= NFP+J-1
                      WRITE(6,630) (IV(I),ESOLN(I),I= III,NLEV,NROW)
                  ENDDO
              ENDIF
          ENDIF
      ENDIF
      IF((NJM.LE.0).AND.(NLEV1.GE.0)) THEN
          NROW=  (NLEV+2)/3
          WRITE(6,619) NLEV
          DO  J= 1,NROW
              WRITE(6,631) (IV(I),IJ(I),ESOLN(I),I= J,NLEV,NROW)
          ENDDO
      ENDIF
      WRITE(6,601)
! The following two lines are to read input file again if you want to find
! the levels for a different potential. Currently LEVEL_RDINP.F90 doesn't
! allow us to read the input file a second time though.
!     GO TO 2
! 999 STOP
c-------------------------------------------------------------------
      CALL MMA_DEALLOCATE(RVB)
      CALL MMA_DEALLOCATE(YVB)
      CALL MMA_DEALLOCATE(DRDY2)
      CALL MMA_DEALLOCATE(FAS)
      CALL MMA_DEALLOCATE(SDRDY)
      CALL MMA_DEALLOCATE(VBZ)
!
      CALL MMA_DEALLOCATE(V1)
      CALL MMA_DEALLOCATE(V2)
      CALL MMA_DEALLOCATE(VJ)
      CALL MMA_DEALLOCATE(V1BZ)
      CALL MMA_DEALLOCATE(V2BZ)
      CALL MMA_DEALLOCATE(WF1)
      CALL MMA_DEALLOCATE(WF2)
!
      CALL MMA_DEALLOCATE(RFN)
      CALL MMA_DEALLOCATE(RRM2)
      CALL MMA_DEALLOCATE(RM2)
      CALL MMA_DEALLOCATE(RRM22)
      CALL MMA_DEALLOCATE(RM22)
  601 FORMAT(1x,79('=')////)
  602 FORMAT( ' Coefficients of expansion for radial matrix element/expe
     1ctation value argument:'/(5X,5(1PD14.6)))
  603 FORMAT(/' Expectation value/matrix element arguments are powers of
     1 a radial function'/5x,'defined by interpolating over read-in poin
     2ts'//' Transition moment function:'/1x,9('==='))
 6604 FORMAT(/' !!! NOTE:  array dimension limit   NDIMR=',i5/'    preve
     1nts preliminary mesh   YH=',f9.6,'  from spanning range [YMIN,YMAX
     2]'/'    so increase mesh by factor',f7.4/)
  604 FORMAT(' Integrate from   Ymin=',f11.7,'   to   Ymax=',f7.4,
     1 '  with mesh  YH=',f9.7/5x,'based on radial variable  yp(r;a)= (r
     2^p - a^p)/(r^p + a^p)'/5x,'for   p=',f6.3,'    and   a=',F9.6/
     3 ' Range corresponds to   Rmin=',f7.3,' [Angst]    to    Rmax= inf
     4inity (!!)'/5x,'and   RH=',F10.7,' [Angst]   at   R=',f11.8//
     5 ' Potential #1 for ',A2,'(',I3,')-',A2,'(',I3,')'/1x,32('='))
  605 FORMAT(/A78/40('=='):/' Generate   ZMU=',F15.11,'(u)',
     1  '   &   BZ=',1PD16.9,'((1/cm-1)(1/Ang**2))'/
     2  10x,'from atomic masses:',0Pf16.11,'  & ',F16.11,'(u)')
  606 FORMAT(' E(v=',i3,', J=',i3,')=',f10.3,'   Bv=',F11.7,
     1  '  -Dv=',1PD12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,
     2  '   Mv=',D12.4,'   Nv=',D12.4,'   Ov=',D12.4)
 6065 FORMAT(' E(v=',i3,', J=',i3,')=',f11.7,'  Bv=',1PD11.4,
     1  '  -Dv=',D12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,
     2  '   Mv=',D12.4,'   Nv=',D12.4,'   Ov=',D12.4)
  607 FORMAT(/' Solve for the',i4,' vibration-rotation levels of Potenti
     1al-1:'/'   (v,J) =',6('  (',i3,',',i3,')':)/(10x,6('  (',i3,',',
     2  i3,')':)))
 6607 FORMAT(/' Solve for',i4,' vibration-rotation levels of Potential-1
     1 using Trial energies:'/(3x,3('   E(',I3,',',I3,')=',
     2 F11.2:)))
  608 FORMAT(/' State-',I1,' electronic angular momentum  OMEGA=',I2/
     1  9x,'yields centrifugal potential  [J*(J+1) -',F5.2,']/r**2' )
 6085 FORMAT(/' State-',I1,' electronic angular momentum  OMEGA=',I2/
     1  9x,'yields centrifugal potential  [J*(J+1) +',I2,']/r**2' )
  609 FORMAT('  Use centrifugal potential for rotation in two dimensions
     1:   (J**2 - 1/4)/r**2')
  610 FORMAT(5X, 'where DREF defined by requiring  <X**1> = 0  for first
     1 level considered')
  611 FORMAT(/' Matrix element argument expansion vble is   X = ((r^',
     1  i1,' - DREF^',i1,')/(r^',i1,' + DREF^',i1,'))')
  612 FORMAT(/' Eigenvalue convergence criterion is   EPS=',1PD8.1,
     1 '(cm-1)'/' Airy function at 3-rd turning point is quasibound oute
     2r boundary condition')
  613 FORMAT(5X,'where reference length is held fixed at   DREF =',
     1 F13.10,'(Angstroms)')
  614 FORMAT(/' Matrix element arguments are powers of the distance  r (
     1in Angstroms)')
  615 FORMAT(/' Matrix element argument expansion variable is:    X = (r
     1 - DREF)/DREF')
  616 FORMAT(/' Matrix element arguments are powers of the squared inver
     1se distance  X = 1/r**',i1)
  617 FORMAT(/' Matrix element argument is fixed as a constant = 1')
  618 FORMAT(' *** PROBLEM *** Searching for  v=',i3,' , J=',i3,
     1 '  ALFas only found to  v=',i3)
  619 FORMAT(/' Find the',i4,' vibration-rotation levels:'/
     1  3('     v   J      E(v)   ')/3(2x,7('---')))
  620 FORMAT(/' An  n=',I2,'  N-D theory extrapolation from last 2 level
     1s implies   vD =',F8.3)
  621 FORMAT(5X,'with the',I4,' missing level(s) predicted to be:'/
     1  4('     v     E(v)   ')/4(4x,7('--')))
  622 FORMAT(/' Search for highest bound  J=',i3,'  level finds  E(v=',
     1  i3,') = VLIM -',1PD12.5/)
  623 FORMAT(/' Find',I4,' Potential-1 vibrational levels with  J=',i3/
     1  4('     v     E(v)   ')/4(4x,7('--')))
  624 FORMAT(4x,'Since the molecule is an ion with charge',SP,I3/6x,"use
     1 Watson's charge-adjusted reduced mass   mu = M1*M2/[M1 + M2 - (",
     2  i2,')*me]')
  625 FORMAT(' For  J=',i3,', try to find the first',i4,' vibrational le
     1vels of Potential-1')
  626 FORMAT(/' *** FAIL to find highest bound J=',i3,'  level from tria
     1l   E = VLIM -',1PD11.4)
  627 FORMAT(/' For vibrational level  v =',I3,'   of Potential-1'/
     1 1X,5('  J',6X,'E',7X)/1X,5(7('--'),2X))
  628 FORMAT((1X,5(I3,F11.3,2X)))
  630 FORMAT((4(I6,F12.4:)))
  631 FORMAT((3(I6,I4,F13.5:)))
  632 FORMAT(1X,79('-'))
  633 FORMAT(' **** Caution: Search for   v=',I3,'   J=',i3,
     1  '  on potential-2 actually found   v=',I3)
  634 FORMAT(/' Using the rotational selection rule:  delta(J)=',
     1 i3,' to',i2,' with increment',i2/'   calculate matrix elements fo
     2r coupling to the',I4,' vibrational levels of'/
     3 '   Potential-2:   v =',14I4:/(21x,14i4:))
 6634 FORMAT(/' Using the rotational selection rule:  delta(J)=',
     1 i3,' to',i2,' with increment',i2/'   calculate matrix elements fo
     2r coupling to the',I4,' vibrational levels of'/'   Potential-2 usi
     3ng trial energies:',2('   E(',I3,')=',F9.2:)/4('   E(',I3,')=',
     4 F9.2:))
! 635 FORMAT(/' Get matrix elements between levels of Potential-1 (above
!    1) & Potential-2 (below)'/1X,39('--')/' For Potential #2:'/
!    2  1x,17('='))
  636 FORMAT(/' Calculate properties of the single potential described a
     1bove')
  637 FORMAT(/' *** Array Dimension OVERFLOW ***   (Number of J sublevel
     1s) > VIBMX=',i4)
  638 FORMAT('   and automatically increment  J  in steps of',i3, ' to a
     1 maximum value of',i4)
  641 FORMAT(1X,39('++'))
  644 FORMAT(/' *** Input data ERROR *** matrix element calculation need
     1s  NLEV2=',i3,' > 0')
  650 FORMAT(/' Matrix element argument is radial first derivative opera
     1tor premultiplied by'/5x,'a power series in  r  of order',i3)
  686 FORMAT(' Potential-',i1,' uses inner boundary condition of  zero v
     1alue  at  RMIN')
  688 FORMAT(' Potential-',i1,' uses symmetric-well inner boundary condi
     1tion  of zero slope at RMIN')
cc703 FORMAT(1X,I4,I5,F13.4,G13.5)
  703 FORMAT(1X,I4,I5,1PD20.11,D13.5)
  723 FORMAT(/A78/1x,'Output values of:  v, J, E & (Level Width)')
  724 FORMAT(//A78//'   v   J    E(v,J)     Width       <KE>',
     1  6x,'<M(r)>  &  <XI**k>  for k=1 to',i3/2x,38('=='))
  725 FORMAT(//A78//"   v'  J'",'  v"  J"     FREQ',"    <v',J'| XI**k",
     1  ' |v",J">  for  k=0  to  MORDR=',i2/2x,37('=='))
! 824 FORMAT(//A78/30('==')/" Note that (v',J') &",' (v",J") strictly la
!    1bel the upper and lower levels, resp.,'/6x,'and  E(lower)=E"'/
!    2 ' but  E(2)-E(1)  is:  (energy of State-2 level) - (energy of Sta
!    3te-1 level)'//12x,'Band'/' dJ(J")',4x,7hv'   v",'  E(lower)  E(2)-
!    4E(1)  A(Einstein)   F-C Factor  ',13h<v'j'|M|v"j"> /
!    5 1x,3('--'),('   -------'),'  --------',3x,
!    6 4('--'),3x,11('-'),3x,11('-'),3x,11('-') )
cc811 FORMAT(//A78/30('==')//12x,"   v'","  J'",'    v"','  J"',
cc   1 '   position    E(upper)    E(lower)',16h   <v'j'|M|v"j">/
cc   2 1x,68('-') )
! 901 FORMAT(//A78/1x,62('==')/'   v    J',7x,'E',10x,'Bv',11x,'-Dv',
!    1  13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
! 902 FORMAT(I4,I5,f25.15,f14.10,6(1PD15.7))
! 904 FORMAT(I4,I5,f25.15,1PD14.7,6(D15.7))
!     END
  997 RC = 0
      END SUBROUTINE LEVEL
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE LEVXPC(KV,JR,EPR,GAMA,NPP,WF,RFN,V,VLIM,YH,DREF,
     1                             NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Calculates expectation values of the kinetic energy and of X**IP
c  (IP=1,MORDR), denoted XPTKE and XPCTR(IP), respectively, for level
c  v=KV, J=JR, E=EPR(cm-1), using wave function WF(i), (i=NBEG,NEND).
c** Assumes units of length are (Angstroms) .
c** Division by BFCT converts potential V(I) to units (cm-1).
c** If (|LXPCT| = 2  or  4), "punch" (WRITE(7,XXX)) results to channel-7
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c!!!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
      REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),SDRDY(:),
     1 VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
c!!!!
      INTEGER I,K,IRFN,IPNCH,ITRY,JR,KV,LXPCT,NPP,NBEG,NEND,MORDR
      REAL*8  WF(NPP),RFN(NPP),V(NPP),XPCTR(0:11),DM(0:20)
      REAL*8 BFCT,DS,DRT,DMR,DER,EPR,EINN,GAMA,YH,DREF,
     1  RR,RXPCT,SS2,SF2,VLIM,XPTKE,PINV
c
      NDIMR= 131074
      CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
      CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
      CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
      CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
      CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
      CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
      EINN= BFCT*EPR
      IPNCH=0
      IF((IABS(LXPCT).EQ.2).OR.(IABS(LXPCT).GE.4)) IPNCH=1
c** MORDR is the highest-power expectation value considered.
      IF(MORDR.GT.11) MORDR=11
      ITRY=20
      IF(((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9)))
     1                                    .AND. (DREF.LE.0.D0)) ITRY=0
c** Start by calculating contributions at end points
    2 SS2=WF(NBEG)**2 *DRDY2(NBEG)
      SF2=WF(NEND)**2 *DRDY2(NEND)
      XPTKE= 0.5D0*(SS2*(EINN-V(NBEG)) + SF2*(EINN-V(NEND)))
      IF(MORDR.GT.0) THEN
          XPCTR(0)= 1.d0/YH
          DO  K=1,MORDR
              SS2=SS2*RFN(NBEG)
              SF2=SF2*RFN(NEND)
              XPCTR(K)=0.5D0*(SS2+SF2)
              ENDDO
          ENDIF
      IF(IRFN.GT.-4) THEN
c** For regular expectation values of a radial function ...
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2 *DRDY2(I)
              XPTKE= XPTKE+ DS*(EINN-V(I))
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
        ELSE
c** For expectation values involving partial derivative operator ...
          DO  K= 0,MORDR
              XPCTR(K)= 0.d0
              ENDDO
          DO  I=NBEG+1,NEND-1
              DS=WF(I)**2 *DRDY2(I)
              XPTKE= XPTKE+ DS*(EINN-V(I))
              DS= WF(I)*(WF(I+1)- WF(I-1))*DRDY2(I)
              IF(MORDR.GT.0) THEN
                  RR= RFN(I)
                  DO  K=1,MORDR
                      DS=DS*RR
                      XPCTR(K)=XPCTR(K)+DS
                      ENDDO
                  ENDIF
              ENDDO
          DO  K= 0,MORDR
              XPCTR(K)= XPCTR(K)/(2.d0*YH)
              ENDDO
        ENDIF
      XPTKE= XPTKE*YH/BFCT
      IF(MORDR.LT.0) GO TO 99
      DMR= 0.d0
      DO  K=0,MORDR
          XPCTR(K)=XPCTR(K)*YH
          DMR= DMR+ DM(K)*XPCTR(K)
          ENDDO
      IF((LXPCT.EQ.1).OR.(IABS(LXPCT).EQ.2)) THEN
          IF(EPR.LE.VLIM) WRITE(6,600) KV,JR,EPR,DMR,XPTKE
          IF(EPR.GT.VLIM) WRITE(6,602) KV,JR,EPR,DMR,XPTKE,GAMA
          IF(IABS(IRFN).LE.9) WRITE(6,604) (K,XPCTR(K),K=1,MORDR)
          IF(IPNCH.GE.1) WRITE(7,701) KV,JR,EPR,GAMA,XPTKE,DMR,
     1                                        (XPCTR(K),K=1,MORDR)
          ENDIF
      IF(ITRY.GT.19) GO TO 99
c** If appropriate, iteratively correct DREF value till distance
c  coordinate expectation value is identically zero.
      IF(IRFN.EQ.-1) THEN
c** For Dunham expansion parameter, define revised function here
          DREF=XPCTR(1)
          DRT=DREF
          WRITE(6,603) ITRY,DRT,DREF
          DO  I= 1,NPP
              RVB(I)= RVB(I)/DREF - 1.D0
              ENDDO
          ITRY=99
          GO TO 2
          ENDIF
c** For Surkus-type expansion parameter, define revised function
      ITRY=ITRY+1
      IF(ITRY.EQ.1) THEN
          RXPCT= XPCTR(1)
          DREF= 0.D0
          DRT= RXPCT
        ELSE
          DER= -IRFN/(2.d0*DREF)
          DRT= -XPCTR(1)/DER
        ENDIF
      DREF=DREF+DRT
      WRITE(6,603) ITRY,DRT,DREF
c** Redefine Surkus-type distance variable RFN using new DREF
      PINV= 1.d0/PRV
      WRITE(6,*) PINV ! Make sure it's "referneced" in THIS subroutine.
      DO  I= 1,NPP
          RFN(I)= (RVB(I)**IRFN - DREF**IRFN)/(RVB(I)**IRFN+ DREF**IRFN)
          ENDDO
      IF(DABS(DRT/DREF).GE.1.D-12) GO TO 2
      CALL MMA_DEALLOCATE(RVB)
      CALL MMA_DEALLOCATE(YVB)
      CALL MMA_DEALLOCATE(DRDY2)
      CALL MMA_DEALLOCATE(FAS)
      CALL MMA_DEALLOCATE(SDRDY)
      CALL MMA_DEALLOCATE(VBZ)
   99 RETURN
  600 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,
     1  '   <KE>=',F11.3)
  602 FORMAT(' E(v=',i3,', J=',i3,')=',f11.3,'   <M(r)>=',G18.10,
     1 '   <KE>=',F11.3/'   Tunneling predissociation  Width(FWHM)=',
     2 G13.6,'    <X**',I2,'>=',F13.8)
  604 FORMAT((8x,3('   <X**',I2,'>=',F13.8:)))
  603 FORMAT(' On iteration #',I2,'  change DREF by',1PD10.2,
     1  '  to   DREF=',0PF13.10,' [Angstroms]')
  701 FORMAT(2I4,F11.3,G11.4,F11.3,3(F12.7)/(5X,6F12.7))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE MATXEL(KV1,JROT1,IOMEG1,EO1,KV2,JROT2,IOMEG2,IRFN,EO2,
     1  NBEG,NEND,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
!     USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
c** Subroutine to calculate matrix elements of powers of the distance
c  coordinate between vib. eigenfunction WF1(i) for v=KV1, J=JROT1 of
c  potential-1 & WF2(I), corresponding to KV2 & JROT2 of potentl.-2
      INTEGER I,J,IOMEG1,IOMEG2,IOMUP,IOMLW,IRFN,JROT1,JROT2,JUP,JLW,
     1 KV1,KV2,KVUP,KVLW,LXPCT,NBEG,NEND,MORDR,NDIMR
      REAL*8  ZMAT(0:20),WF1(NEND),WF2(NEND),RFN(NEND),DM(0:MORDR),
     1  DRDY2(NDIMR)
      REAL*8  AEINST,DEG,DME,DSM,EO1,EO2,ELW,FCF,FREQ,OMUP,RH,RI,SJ,
     1    ZJUP
      CHARACTER*1  DJ(-3:3)
      DATA DJ/'N','O','P','Q','R','S','T'/
      ZMAT(0)= 0.D0
      IF(MORDR.GE.1) THEN
          DO  J= 1,MORDR
              ZMAT(J)= 0.D0
              ENDDO
          ENDIF
      IF(IRFN.NE.-4) THEN
c** For regular power series or function matrix elements ...
          DO  I=NBEG,NEND
              DSM= WF2(I)*WF1(I) * DRDY2(I)
              ZMAT(0)= ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J= 1,MORDR
                      DSM= DSM*RI
                      ZMAT(J)= ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
        ELSE
c** For partial derivative matrix elements ...
          DO  I= NBEG+1,NEND-1
              DSM= WF1(I)*(WF2(I+1)- WF2(I-1)) * DRDY2(I)
              ZMAT(0)= ZMAT(0)+DSM
              RI= RFN(I)
              IF(MORDR.GE.1) THEN
                  DO  J= 1,MORDR
                      DSM= DSM*RI
                      ZMAT(J)= ZMAT(J)+DSM
                      ENDDO
                  ENDIF
              ENDDO
          DO  J= 0,MORDR
              ZMAT(J)= ZMAT(J)/(2.d0*RH)
              ENDDO
        ENDIF
      DME= 0.D0
      FCF= (ZMAT(0)*RH)**2
      IF(MORDR.GE.0) THEN
          DO  J= 0,MORDR
              ZMAT(J)= ZMAT(J)*RH
              DME= DME+DM(J)*ZMAT(J)
              ENDDO
          ENDIF
      FREQ= EO2-EO1
      ELW= DMIN1(EO1,EO2)
c** Now calculate the Honl-London Factor for the particular transition
c   Factors updated as per Hansson & Watson JMS (2005).
      SJ= 0.D0
      KVUP= KV1
      KVLW= KV2
      JUP= JROT1
      JLW= JROT2
      IOMUP= MAX(IOMEG1,0)
      IOMLW= MAX(IOMEG2,0)
      IF(EO2.GT.EO1) THEN
          KVUP= KV2
          KVLW= KV1
          JUP= JROT2
          JLW= JROT1
          IOMUP= MAX(IOMEG2,0)
          IOMLW= MAX(IOMEG1,0)
          ENDIF
      ZJUP= JUP
      OMUP= IOMUP
      DEG= 2*JUP+ 1
      IF((JLW.LT.IOMLW).OR.(JUP.LT.IOMUP)) GO TO 50
      IF(IOMUP.EQ.IOMLW) THEN
c** Factors for  DELTA(LAMBDA) = 0  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+ OMUP)*(JUP- IOMUP)/ZJUP
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))
     1                       SJ= DEG*OMUP**2/(ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1)) SJ= (ZJUP+1.D0+OMUP)*(JUP+1-IOMUP)/
     1                                                     (ZJUP+1.D0)
          ENDIF
      IF(IOMUP.EQ.(IOMLW+1)) THEN
c** Factors for  DELTA(LAMBDA) = +1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (ZJUP+OMUP)*(JUP-1+IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))
     1       SJ= (ZJUP+OMUP)*(JUP+1-IOMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1))
     1       SJ= (JUP+1-IOMUP)*(ZJUP+2.D0-OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF
      IF(IOMUP.LT.IOMLW) THEN
c** Factors for  DELTA(LAMBDA) = -1  transitions of spin singlets
          IF(JUP.EQ.(JLW+1)) SJ= (JUP-IOMUP)*(JUP-1-IOMUP)/(2.D0*ZJUP)
          IF((JUP.EQ.JLW).AND.(JUP.GT.0))
     1      SJ= (JUP-IOMUP)*(ZJUP+1.D0+OMUP)*DEG/(2.D0*ZJUP*(ZJUP+1.D0))
          IF(JUP.EQ.(JLW-1))
     1           SJ= (ZJUP+1.D0+OMUP)*(ZJUP+2.D0+OMUP)/(2.D0*ZJUP+2.D0)
          ENDIF
c... finally, include Hansson-Watson  w0/w1  term in Honl-London factor
      IF((MIN(IOMUP,IOMLW).EQ.0).and.(IOMUP.NE.IOMLW)) SJ= SJ+SJ
c
c** For FREQ in  cm-1  and dipole moment in  debye , AEINST is the
c  absolute Einstein radiative emission rate (s-1) , using the
c  rotational intensity factors for sigma-sigma transitions.
   50 CONTINUE
      AEINST = DABS(3.1361891D-7 *DABS(FREQ)**3*DME**2 * SJ/DEG)
      IF(LXPCT.GT.0) THEN
          WRITE(6,600) KV1,JROT1,EO1,KV2,JROT2,EO2
          IF(IABS(IRFN).LE.9) WRITE(6,602) (J,ZMAT(J),J= 0,MORDR)
          WRITE(6,604) FCF,DME,FREQ,AEINST
          WRITE(6,606)
          ENDIF
      IF((IABS(LXPCT).EQ.4).OR.(IABS(LXPCT).EQ.5).AND.(SJ.GT.0.D0)) THEN
          IF(IABS(JUP-JLW).LE.3) WRITE(8,801) DJ(JUP-JLW),JLW,KVUP,
     1                                    KVLW,ELW,FREQ,AEINST,FCF,DME
c... Special printout for Hui/LeRoy N2 Quadrupole paper [JCP 1XX (2007)]
ccc       E00= 1175.7693d0
cc        WRITE(11,811) -FREQ,KVUP,JUP,KVLW,JLW,-FREQ,ELW-FREQ-E00,
cc   1                                      ELW-E00,DME**2
cc811 FORMAT(F12.4,2I4,I6,I4,3f12.4,1PD15.6)
          IF(IABS(JUP-JLW).GT.3) WRITE(8,802) JUP-JLW,JLW,KVUP,
     1                                    KVLW,ELW,FREQ,AEINST,FCF,DME
          ENDIF
      IF(IABS(LXPCT).GE.5)
c    1         WRITE(7,701) KVUP,JUP,KVLW,JLW,(ZMAT(J),J=0,MORDR)
     1         WRITE(7,701) KVUP,JUP,KVLW,JLW,FREQ,(ZMAT(J),J=0,MORDR)
      RETURN
  600 FORMAT(' Coupling   E(v=',I3,', J=',I3,')=',F12.4,'   to   E(v=',
     1 I3,', J=',I3,')=',F12.4)
  602 FORMAT(5x,'Moment matrix elements:',2('   <X**',I2,'>=',F14.10:),
     1  1x/(3x,3('   <X**',I2,'>=',F14.10:),1x))
  604 FORMAT(' FCF=',1PD11.4,'   <M>=',D12.5,'   d(E)=',0PF10.2,
     1  '   A(Einst)=',1PD11.4,' s-1')
  606 FORMAT(1X,79('+'))
  701 FORMAT(4I4,F12.4,4F12.8:/(4X,6F12.8))
c 701 FORMAT(4I4,6F12.8:/(16X,6F12.8))
  801 FORMAT(1x,A1,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
  802 FORMAT(i2,'(',I3,')  ',I3,' -',I3,F10.2,F11.2,3(1PD14.5))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
c***********************************************************************
      SUBROUTINE CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V,WF0,RM2,RCNST)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  Subroutine solving the linear inhomogeneous differential equations
c  formulated by J.M. Hutson [J.Phys.B14, 851 (1982)] for treating
c  centrifugal distortion as a perturbation, to determine centrifugal
c  distortion constants of a diatomic molecule.  Uses the algorithm of
c  J. Tellinghuisen [J.Mol.Spectrosc. 122, 455 (1987)].  The current
c  version calculates Bv, Dv, Hv, Lv, Mv, Nv and Ov and writes them out,
c  but does not return values to the calling program.
c
c** On entry:   EO    is the eigenvalue (in units [cm-1])
c               NBEG & NEND  the mesh point range over which the input
c               NDIMR  is dimension of arrays  V(i), WF0(i), ... etc.
c wavefunction  WF0  (in units 1/sqrt(Ang))  has non-negligible values
c               BvWn  is the numerical factor (hbar^2/2mu) [cm-1 Ang^2]
c               YH    is the integration stepsize
c               WARN  is an integer flag: > 0 print internal warnings,
c               V(i)  is the effective potential (including centrifugal
c                     term if calculation performed at  J > 0) in
c                     'internal' units, including the factor  YH**2/BvWN
c               RM2(i) is the array  (r')^2/(distance**2)
c** On exit:    RCNST(i)  is the set of 7 rotational constants: Bv, -Dv,
c                       Hv, Lv, Mv, Nv & Ov
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c** Dimension:  potential arrays  and  vib. level arrays.
c!!
      INTEGER NDIMR
!     PARAMETER (NDIMR= 200001)
! A limit set by the -fmax-stack-var-size in OpenMolcas is making arrays
! of the above size too large. If we can't get that increased, we could
! use an ALLOCATABLE array or use -frecursive.
!     PARAMETER (NDIMR= 131074)
!     REAL*8 PRV,ARV,RVB(NDIMR),YVB(NDIMR),DRDY2(NDIMR),FAS(NDIMR),
!    1                                         SDRDY(NDIMR),VBZ(NDIMR)
      REAL*8 PRV,ARV
      REAL*8, ALLOCATABLE :: RVB(:),YVB(:),DRDY2(:),FAS(:),
     1                                         SDRDY(:),VBZ(:)
      COMMON /BLKAS/PRV,ARV!,RVB,YVB,DRDY2,SDRDY,FAS,VBZ
c!!
      INTEGER I,M,IPASS,M1,M2,NBEG,NEND,WARN
!      REAL*8 V(NEND),WF0(NEND),RM2(NEND),P(NDIMR),WF1(NDIMR),
!     1                                            WF2(NDIMR),RCNST(7)
      REAL*8 V(NEND),WF0(NEND),RM2(NEND),RCNST(7)
      REAL*8, ALLOCATABLE :: P(:),WF1(:),WF2(:)
      REAL*8 BvWN,DV,DVV,HVV,HV2,LVV,LV2,MVV,MV2,NVV,OVV,EO,E,YH,YH2,
     1  ZTW,AR,R2IN,G2,G3,P0,P1,P2,P3,PI,PIF,PRS,PRT,V1,V2,V3,Y1,Y2,Y3,
     2  TSTHv,TSTLv,TSTMv,AMB,AMB1,AMB2,
     3  OV,OV01,OV02,OV03,OV11,OV12,OV13,OV22,OV23,OV33,
     4  PER01,PER02,PER03,PER11,PER12,PER13,PER22,PER23,PER33,R2XX
c
      NDIMR = 131074
      CALL MMA_ALLOCATE(RVB,NDIMR,LABEL='RVB')
      CALL MMA_ALLOCATE(YVB,NDIMR,LABEL='YVB')
      CALL MMA_ALLOCATE(DRDY2,NDIMR,LABEL='DRDY2')
      CALL MMA_ALLOCATE(FAS,NDIMR,LABEL='FAS')
      CALL MMA_ALLOCATE(SDRDY,NDIMR,LABEL='SDRDY')
      CALL MMA_ALLOCATE(VBZ,NDIMR,LABEL='VBZ')
!
      CALL MMA_ALLOCATE(P,NDIMR,LABEL='P')
      CALL MMA_ALLOCATE(WF1,NDIMR,LABEL='WF1')
      CALL MMA_ALLOCATE(WF2,NDMIR,LABEL='WF2')
      P0=0
      MV2=0
      LV2=0
      G3=0
      IF(NEND.GT.NDIMR) THEN
          WRITE(6,602) NEND,NDIMR
          RETURN
          ENDIF
      ZTW= 1.D0/12.d0
      YH2 = YH*YH
      DV = YH2*ZTW
      E= EO*YH2/BvWN
      IPASS = 1
      OV01 = 0.D0
      OV02 = 0.D0
      OV03 = 0.D0
      OV11 = 0.D0
      OV22 = 0.D0
      OV12 = 0.D0
      OV33 = 0.D0
      OV23 = 0.D0
      OV13 = 0.D0
      PER01 = 0.D0
      PER02 = 0.D0
      PER03 = 0.D0
      PER11 = 0.D0
      PER12 = 0.D0
      PER13 = 0.D0
      PER22 = 0.D0
      PER23 = 0.D0
      PER33 = 0.D0
c** First, calculate the expectation value of  1/r**2  and hence Bv
      R2IN= 0.5D0*(RM2(NBEG)*WF0(NBEG)**2 + RM2(NEND)*WF0(NEND)**2)
      DO   I= NBEG+1, NEND-1
         R2IN= R2IN+ RM2(I)*WF0(I)**2
         ENDDO
      R2IN = R2IN*YH
      RCNST(1)= R2IN*BvWN
c
c** On First pass  IPASS=1  and calculate first-order wavefx., Dv & Hv
c  On second pass  IPASS=2  and calculate second-order wavefx., Lv & Mv
c  On third pass   IPASS=3  and calculate third-order wavefx., Nv & Ov
c
   10 P1= 0.D0
      P2= 0.D0
c
c     P1= WF0(NEND)
c     P2= WF0(NEND-1)
c
      P(NEND) = P1
      P(NEND-1) = P2
      V1 = V(NEND) - E*DRDY2(NEND)
      V2 = V(NEND-1) - E*DRDY2(NEND-1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1)
     1                   - DV*(RM2(NEND) - R2IN*DRDY2(NEND))*WF0(NEND)
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF0(NEND-1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)
     1                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF1(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF1(NEND-1)
     1                                 - DVV*WF0(NEND-1)*DRDY2(NEND-1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF1(NEND) - HVV*WF0(NEND)
     1                     - (RM2(NEND) - R2IN*DRDY2(NEND))*WF2(NEND))
          G2 = (RM2(NEND-1) - R2IN*DRDY2(NEND-1))*WF2(NEND-1)
     1             - (DVV*WF1(NEND-1) + HVV*WF0(NEND-1))*DRDY2(NEND-1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      M= NEND-1
c** Now - integrate inward from outer end of range
      DO  I = NBEG+2,NEND
          M = M-1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          R2XX= R2IN*DRDY2(M)
          IF(IPASS.EQ.1) G3 = (RM2(M)- R2XX)*WF0(M)
          IF(IPASS.EQ.2) G3 = (RM2(M)-R2XX)*WF1(M) - DVV*WF0(M)*DRDY2(M)
          IF(IPASS.EQ.3) G3 = (RM2(M)- R2XX)*WF2(M)
     1                            - (DVV*WF1(M) + HVV*WF0(M))*DRDY2(M)
          V3 = V(M) - E*DRDY2(M)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          IF(V3.LT.0.D0)  GO TO 32
          P(M) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          ENDDO
      GO TO 90
c** Escaped loop at outer turning point:  initialize outward integration
   32 PRS = P3
      PRT = P(M+1)
      P1 = 0.D0
      P2 = 0.D0
c
c     P1 = WF0(NBEG)
c     P2 = WF0(NBEG+1)
c
      P(NBEG) = P1
      P(NBEG+1) = P2
      V1 = V(NBEG) - E*DRDY2(NBEG)
      V2 = V(NBEG+1) - E*DRDY2(NBEG+1)
      IF(IPASS.EQ.1) THEN
          Y1 = P1*(1.D0 - ZTW*V1)
     1                   - DV*(RM2(NBEG) - R2IN*DRDY2(NBEG))*WF0(NBEG)
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF0(NBEG+1)
        ELSEIF(IPASS.EQ.2) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*(DVV*WF0(NEND)*DRDY2(NBEG)
     1                     - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF1(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF1(NBEG+1)
     1                                 - DVV*WF0(NBEG+1)*DRDY2(NBEG+1)
        ELSEIF(IPASS.EQ.3) THEN
          Y1 = P1*(1.D0 - ZTW*V1) + DV*((DVV*WF1(NEND) + HVV*WF0(NEND))
     1       *DRDY2(NBEG)  - (RM2(NBEG) - R2IN*DRDY2(NBEG))*WF2(NBEG))
          G2 = (RM2(NBEG+1) - R2IN*DRDY2(NBEG+1))*WF2(NBEG+1)
     1             - (DVV*WF1(NBEG+1) + HVV*WF0(NBEG+1))*DRDY2(NBEG+1)
        ENDIF
      Y2 = P2*(1.D0 - ZTW*V2) - DV*G2
      AR = 0.D0
      M1 = M+1
c** Now ... integrate outward from inner end of range
      DO  I = NBEG+2,M1
          Y3 = Y2 + Y2 - Y1 + YH2*G2 + V2*P2
          P0 = WF0(I)
          R2XX= R2IN*DRDY2(I)
          IF(IPASS.EQ.1) G3 = (RM2(I)-R2XX)*P0
          IF(IPASS.EQ.2) G3 = (RM2(I)-R2XX)*WF1(I) - DVV*P0*DRDY2(I)
          IF(IPASS.EQ.3) G3 = (RM2(I)-R2XX)*WF2(I)
     1                                - (DVV*WF1(I) + HVV*P0)*DRDY2(I)
          V3 = V(I) - E*DRDY2(I)
          P3 = (Y3 + DV*G3)/(1.D0 - ZTW*V3)
          P(I) = P3
          Y1 = Y2
          Y2 = Y3
          V2 = V3
          P2 = P3
          G2 = G3
          AR = AR + P0*P3*DRDY2(I)
          ENDDO
c** Average for 2 adjacent mesh points to get Joel's "(a-b)"
      AMB2 = (P3-PRT)/P0
      AMB1 = (P(M)-PRS)/WF0(M)
      AMB = (AMB1+AMB2)*0.5D0
      M2 = M+2
c** Find the rest of the overlap with zero-th order solution ...
      DO  I = M2,NEND
          P0 = WF0(I)
          PI = P(I) + AMB*P0
          P(I) = PI
          AR = AR + PI*P0*DRDY2(I)
          ENDDO
      OV = AR*YH
      DO  I = NBEG,NEND
          P0 = WF0(I)
c ... and project out contribution of zero'th-order part of solution
          PI = P(I) - OV*P0
          PIF = PI*RM2(I)
          IF(IPASS.EQ.1) THEN
c** Now - on first pass accumulate integrals for Dv and Hv
              WF1(I) = PI
              OV01 = OV01 + PI*P0 * drdy2(i)
              OV11 = OV11 + PI*PI * drdy2(i)
              PER01 = PER01 + PIF*P0
              PER11 = PER11 + PI*PIF
            ELSEIF(IPASS.EQ.2) THEN
c ... and on next pass, accumulate integrals for Lv and Mv
              WF2(I) = PI
              P1 = WF1(I)
              OV02 = OV02 + PI*P0 * drdy2(i)
              OV12 = OV12 + PI*P1 * drdy2(i)
              OV22 = OV22 + PI*PI * drdy2(i)
              PER02 = PER02 + PIF*P0
              PER12 = PER12 + PIF*P1
              PER22 = PER22 + PI*PIF
            ELSEIF(IPASS.EQ.3) THEN
c ... and on next pass, accumulate integrals for Nv and Ov
              P1 = WF1(I)
              P2 = WF2(I)
              OV03 = OV03 + PI*P0 * drdy2(i)
              OV13 = OV13 + PI*P1 * drdy2(i)
              OV23 = OV23 + PI*P2 * drdy2(i)
              OV33 = OV33 + PI*PI * drdy2(i)
              PER03 = PER03 + PIF*P0
              PER13 = PER13 + PIF*P1
              PER23 = PER23 + PIF*P2
              PER33 = PER33 + PIF*PI
            ENDIF
          ENDDO
      IF(IPASS.EQ.1) THEN
          DVV = YH*PER01
          HVV = YH*(PER11 - R2IN*OV11)
          IPASS = 2
          RCNST(2) = DVV*BvWN
          RCNST(3) = HVV*BvWn
          GO TO 10
        ELSEIF(IPASS.EQ.2) THEN
          HV2 = YH*PER02*BvWN
          LVV = YH*(PER12 - R2IN*OV12 - DVV*OV11)
          MVV = YH*(PER22 - R2IN*OV22 - 2.D0*DVV*OV12 - HVV*OV11)
          IPASS = 3
          RCNST(4) = LVV*BvWN
          RCNST(5) = MVV*BvWN
          GO TO 10
        ELSEIF(IPASS.EQ.3) THEN
          LV2 = YH*PER03*BvWN
          MV2 = YH*(PER13 - R2IN*OV13 - DVV*OV12 - HVV*OV11)*BvWN
          NVV = YH*(PER23 - R2IN*OV23 - DVV*(OV13 + OV22)
     1                                     - 2.D0*HVV*OV12 - LVV*OV11)
          OVV = YH*(PER33 - R2IN*OV33 - 2.D0*DVV*OV23
     1             - HVV*(2.D0*OV13+ OV22) - 2.D0*LVV*OV12 - MVV*OV11)
          RCNST(6) = NVV*BvWN
          RCNST(7) = OVV*BvWN
        ENDIF
      IF(WARN.GT.0) THEN
          IF(DMAX1(DABS(OV01),DABS(OV02),DABS(OV01)).GT.1.D-9)
     1                                     WRITE(6,604) OV01,OV02,OV03
          TSTHV= dabs(RCNST(3)/HV2-1.D0)
          TSTLV= dabs(RCNST(4)/LV2-1.D0)
          TSTMV= dabs(RCNST(5)/MV2-1.D0)
          IF(DMAX1(TSTHV,TSTLV,TSTMV).GT.1.d-5)
     1                                  WRITE(6,603) TSTHV,TSTLV,TSTMV
          ENDIF
      CALL MMA_DEALLOCATE(RVB)
      CALL MMA_DEALLOCATE(YVB)
      CALL MMA_DEALLOCATE(DRDY2)
      CALL MMA_DEALLOCATE(FAS)
      CALL MMA_DEALLOCATE(SDRDY)
      CALL MMA_DEALLOCATE(VBZ)
!
      CALL MMA_DEALLOCATE(P)
      CALL MMA_DEALLOCATE(WF1)
      CALL MMA_DEALLOCATE(WF2)
      RETURN
   90 WRITE(6,601) EO
      RETURN
  601 FORMAT(' *** ERROR in CDJOEL *** for input energy  E =',f12.4,
     1   '  never reach outer turning point')
  602 FORMAT(/' *** Dimensioning PROBLEM in CDJOEL ***   NEND=',i6,
     1  ' > NDIMR=',i6)
  603 FORMAT(' ** CAUTION ** Comparison tests for Hv, Lv & Mv give:',
     1 3(1Pd9.1))
  604 FORMAT(' ** CAUTION ** CDJOEL orthogonality tests OV01,OV02 & OV03
     1:',3(1Pd9.1))
      END
c23456789 123456789 123456789 123456789 123456789 123456789 123456789 12
