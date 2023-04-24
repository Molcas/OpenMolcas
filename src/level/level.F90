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

!      it will also automatically generate the eigenvalues etc. for all
!      vibrational and/or rotational levels of a given well-behaved
!      single-minimum potential.
!***** Main calling and I/O routines.  Last Updated  28 June 2009 *****
      SUBROUTINE LEVEL(RC)
      USE STDALLOC, ONLY: MMA_ALLOCATE, MMA_DEALLOCATE
      USE LEVEL_COMMON
      IMPLICIT NONE
      INTEGER, INTENT(OUT) :: RC
!** Dimensions for  potential arrays  and  vib. level arrays.
      INTEGER VIBMX,MORDRMX,RORDR,NTPMX
      PARAMETER (VIBMX=400,RORDR=7,MORDRMX=20,NTPMX= 1600)
!!!---------------------------------------------------------------------
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
!!!---------------------------------------------------------------------
      INTEGER I,J,M,III,IJD,ILEV1,ILEV2,IOMEG1,IOMEG2,INNOD1,INNOD2,    &
     & INNER,SINNER,IQT,IWR,IRFN,IVD,IVS,IAN1,IAN2,IMN1,IMN2,GEL1,GEL2, &
     & GNS1,GNS2,JDJR,JCT,J2DL,J2DU,J2DD,JROT,JROT2,JLEV,JREF, ICOR,    &
     & CHARGE, KV,KV2,KVIN,LCDC,LPRWF,LRPT,LXPCT,MORDR,NUSEF,ILRF,IR2F, &
     & NUMPOT,NBEG,NBEG2,NEND,NEND2,NPP,NCN1,NCN2,NCNF,NLEV,NLEV1,      &
     & NLEV2,NJM,NJMM,NFP,NLP,NRFN,NROW,WARN,VMAX,VMAX1,VMAX2,AFLAG,    &
     & AUTO1,AUTO2, IV(VIBMX),IJ(VIBMX),IV2(VIBMX),JWR(VIBMX),          &
     & INNR1(0:VIBMX),INNR2(0:VIBMX),NTP,LPPOT,IPOTL,PPAR,QPAR,NSR,NLR, &
     & IBOB,NCMM,IVSR,IDSTT,MMLR(3)
!
!     REAL*8 ZK1(0:VIBMX,0:RORDR),ZK2(0:VIBMX,0:RORDR),RCNST(RORDR),
!    1 V1(NDIMR),V2(NDIMR),VJ(NDIMR),V1BZ(NDIMR),V2BZ(NDIMR),
!    2 WF1(NDIMR),WF2(NDIMR),CMM(3),PARM(4)
      REAL*8 ZK1(0:VIBMX,0:RORDR),ZK2(0:VIBMX,0:RORDR),RCNST(RORDR),    &
     & CMM(3),PARM(4)
      REAL*8, ALLOCATABLE ::  V1(:),V2(:),VJ(:),V1BZ(:),V2BZ(:),        &
     & WF1(:),WF2(:)
!
!     REAL*8  RFN(NDIMR),RRM2(NDIMR),RM2(NDIMR),RRM22(NDIMR),
!    2  RM22(NDIMR),GV(0:VIBMX),ESOLN(VIBMX),ESLJ(VIBMX),XIF(NTPMX),
!    4  YIF(NTPMX),ABUND1,ABUND2,MASS1,MASS2,DM(0:MORDRMX)
      REAL*8 GV(0:VIBMX),ESOLN(VIBMX),ESLJ(VIBMX),XIF(NTPMX),           &
     &  YIF(NTPMX),ABUND1,ABUND2,MASS1,MASS2,DM(0:MORDRMX)
      REAL*8, ALLOCATABLE :: RFN(:),RRM2(:),RM2(:),RRM22(:),RM22(:)
      REAL*8 BZ,BvWN,BFCT,BEFF,DEJ,EPS,EO,EO2,EJ,EJ2,EJP,EJREF,GAMA,    &
     & MEL,PMAX1,PMAX2,PW,RH,RMIN,RR,RRp,PINV,DRDY,YH,YH2,YMIN,YMINN,   &
     & YMAX,DREF,DREFP,CNN1,CNN2,RFLIM,CNNF,RFACTF,MFACTF,SOMEG1,       &
     & SOMEG2,VLIM1,VLIM2,VD,VDMV,XX,ZMU,GI,GB,GBB,WV,FFAS,SL,DSCM,     &
     & REQ,RREF,RHOAB
!
      CHARACTER*78 TITL
      CHARACTER*2 NAME1,NAME2
!
      DATA MEL/5.4857990945d-4/,YMAX/1.d+00/
!** Default (Q-branch) defining J-increments for matrix element calcn.
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
!** If the read-in value of IAN1 and/or IAN2 is .LE.0, then instead of
!  using the MASS table, read an actual particle mass for it/them.
!** CHARGE (integer) is the charge on the molecule (=0 for neutral). If
!   (CHARGE.ne.0)  generate & use Watson's charge-adjusted reduced mass.
!** Parameter NUMPOT specifies whether to calculate eigenvalues etc. for
!  a single potential (when NUMPOT.LE.1), or to generate two independent
!  potentials & calculate matrix elements coupling levels of one to
!  levels of the other (for NUMPOT.GE.2).
!----------------------------------------------------------------------
    2 CALL LEVEL_RDINP(IAN1,IMN1,IAN2,IMN2,CHARGE,NUMPOT,RH,RMIN,PRV,   &
     & ARV,EPS,NTP,LPPOT,IOMEG1,VLIM1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,     &
     & DSCM,REQ,RREF,NCMM,IVSR,IDSTT,RHOAB,MMLR,CMM,PARM,NLEV1,AUTO1,   &
     & LCDC,LXPCT,NJM,JDJR,IWR,LPRWF)
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
!----------------------------------------------------------------------
!** Subroutine MASSES returns the names of the atoms NAMEi,ground
!  electronic state degeneracy GELi, nuclear spin degeneracy GNSi,
!  mass MASSi, and isotopic abundance ABUNDi for a given atomic isotope.
      IF((IAN1.GT.0).AND.(IAN1.LE.109)) THEN
          CALL MASSES(IAN1,IMN1,NAME1,GEL1,GNS1,MASS1,ABUND1)
      ELSE
!** If particle-i is not a normal atomic isotope, read a 2-character
!   name (enclosed between '', as in 'mu') and its actual mass.
!----------------------------------------------------------------------
!         READ(5,*) NAME1, MASS1
!----------------------------------------------------------------------
      ENDIF
      IF((IAN2.GT.0).AND.(IAN2.LE.109)) THEN
          CALL MASSES(IAN2,IMN2,NAME2,GEL2,GNS2,MASS2,ABUND2)
      ELSE
!----------------------------------------------------------------------
!         READ(5,*) NAME2, MASS2
!----------------------------------------------------------------------
      ENDIF
      ZMU= MASS1*MASS2/(MASS1+MASS2- CHARGE*MEL)
!=======================================================================
! TITL is a title or output header of up to 78 characters, read on a
!   single line enclosed between single quotes: e.g.  'title of problem'
!=======================================================================
!     READ(5,*) TITL
      TITL = 'Beginning execution of LEVEL:'
!----------------------------------------------------------------------
!** Numerical factor  16.85762920 (+/- 0.00000011) based on Compton
!  wavelength of proton & proton mass (u) from 2002 physical constants.
      BZ= ZMU/16.85762920D0
      WRITE(6,605) TITL,ZMU,BZ,MASS1,MASS2
      BvWN= 1.D0/BZ
      IF(CHARGE.NE.0) WRITE(6,624) CHARGE,CHARGE
      EJ= 0.D0
      EJ2= 0.D0
      LRPT= 1
!** Lower limit (RMIN) and increment (RH) of integration in (Angstroms).
!** Upper limit of the reduced variable integration range automatically
!  set at  YMAX= 1.0 , which corresponds to  RMAX= infinity !!.
!* A hard wall boundary condition may be imposed at a smaller distance
!  using an appropriate choice of the read-in level parameter IV (below)
!!! The radial integration variable is  yp(r;Reff)  with   p= PRV
!** EPS (cm-1) is the desired eigenvalue convergence criterion
!---------------------------------------------------------------------
!     READ(5,*) RH, RMIN, PRV, ARV, EPS
!---------------------------------------------------------------------
!** NPP = no. of points in potential and wavefunction array.
!!! First ... calculate new AS radial mesh YH implied but the given RH
      I= INT(0.5d7*(PRV/ARV)*RH)
!.... give YH a rounded-off value (to 8 digits)
      YH= DBLE(I)*1.d-07
      aRVp= ARV**PRV
      RRp= RMIN**PRV
      YMIN= (RRp - aRVp)/(RRp + aRVp)
      YMAX= 1.d0
!** NPP = no. of points in potential and wavefunction array.
      NPP= INT(((YMAX-YMIN)/YH+ 1.00001))
      IF(NDIMR.LT.NPP) THEN
!         WRITE(6,6604)  NDIMR,YH,DFLOAT(NPP)/DFLOAT(NDIMR)
          WRITE(6,6604)  NDIMR,YH,DBLE(NPP)/DBLE(NDIMR)
          NPP= NDIMR
      ENDIF
!... reset YMIN slightly to precisely span range
      YMIN= YMAX - (NPP-1)*YH
      YH2= YH*YH
      BFCT= BZ*YH2
      YMINN= YMIN-YH
      WRITE(6,604) YMIN,YMAX,YH,PRV,ARV,RMIN,RH,ARV,                    &
     &                                           NAME1,IMN1,NAME2,IMN2
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
!... 'fake' RMAX value to ensure last 1/R**2 point is stable.
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

!
!++ Begin reading appropriate parameters & preparing potential(s)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!+  Subroutine "PREPOT" prepares (and if desired, writes) the potential
!+  array V(i) (cm-1)  at the NPP distances RVB(i) (Angst).
!** NPP = no. of points in potential and wavefunction array.
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!* If NTP > 0 :  define potential by interpolation over & extrapolation
!        beyond the NTP read-in turning points using subroutine GENINT.
!   If NTP.le.0 : generate a (fully analytic) potential in POTGEN.
!* If LPPOT > 0 : at every |LPPOT|-th point, print potential and
!      derivatives-by-differences. ***  If  LPPOT < 0  write potential
!      at every |LPPOT|-th point to channel-8 in a compact format **
!* OMEGA is the electronic contribution to the angular momentum such
!  that the reduced centrifugal potential is:  (J*(J+1)-OMEGA**2)/R**2
!* Set (OMEGA.GE.99) if wish to use centrifugal factor for rotation
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
!    (usually choose NUSE even, say, = 6, 8 or 10). ***  If(NUSE.LE.0)
!    interpolate with cubic spline instead of local polynomials.
!** If IR2 > 0 , interpolate over  YI*XI**2 ; otherwise on  YI  itself
!   This may help if interpolation has trouble on steep repulsive wall.
!** ILR specifies how to extrapolate beyond largest input distance XI(i)
!  If ILR < 0 , fit last 3 points to:  VLIM - A*exp(-b*(R-R0)**2)
!  If ILR = 0 , fit last 3 points to:  VLIM - A*R**p *exp(-b*R)
!  If ILR = 1 : fit last two points to:  VLIM - A/R**B .
!** If(ILR > 1) fit last turning points to:  VLIM - sum{of ILR
!  inverse-power terms beginning with  1/R**NCN}. *** If CNN.ne.0 ,
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
!++   READ(5,*) NUSE, IR2, ILR, NCN, CNN
!++   READ(5,*) RFACT, EFACT, VSHIFT
!++   READ(5,*) (XI(I), YI(I), I= 1,NTP)
!-----------------------------------------------------------------------
!** NCN1 (returned by PREPOT) is the power of the asymptotically-
!  dominant inverse-power long range potential term.
!   VLIM1, VLIM1, V1, NCN1 and CNN1 are not defined yet, but are input parameters for
!  PREPOT
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      WRITE(6,*) 'Exiting level.f'
      WRITE(6,*) 'Entering prepot.f'
      WRITE(6,*) ''
      CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,   &
     &  V1,CNN1,NCN1,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,PARM,   &
     &  MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)
!     CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG1,RVB,RRM2,VLIM1,
!    1                                                   V1,CNN1,NCN1)
      WRITE(6,*) 'Successfully made it through Prepot.f!'
!OPTIONALLY WRITE FIRST FEW v(r) VALUES, THE LAST ONE AND A MIDDLE ONE
!     DO I=1,3
!      WRITE(6,*) 'V(',I,')=',V1(I)
!     ENDDO
!     WRITE(6,*) 'V(                 20000)=',V1(20000)
!     WRITE(6,*) 'V(',NPP,')=',V1(NPP)
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** If (NTP.le.0) PREPOT uses subroutine POTGEN to generate a fully
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
!** IBOB specifies whether (if > 0) or not (if .le. 0) atomic mass
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
!      with NVARB coefficients PARM(i).    [!! MPAR .ge.1 !!]
!    * For conventional "simple" Morse potential,  NVARB=1 & MPAR dummy
!*  Special option #1: set  MPAR= -1  to produce Wei Hua's 4-parameter
!      modified Morse function with  b= PARM(1)  and C= PARM(2).
!*  Special option #2: set  MPAR= -2  to produce Coxon's "Generalized
!      Morse Oscillator" potential with exponent expansion in (R-Re)]
! ...  otherwise, set  MPAR.ge.0
!** If IPOTL=4  generate an MLR potential [Mol.Phys. 105, 691 (2007)]
!      If MPAR > 0  exponent parameter defined in terms of a polynomial
!           of order (NVARB-2) with the NVARB coefficients  PARM(j).
!           in y_{MPAR}= (R**MPAR - Re**MPAR)/(R**MPAR + Re**MPAR),
!           and long-range defined by NCN inverse-power terms
!      If MPAR = 0  exponent polynomial variable is  2*y_{1}= y{O-T}
!      If MPAR < 0  exponent polynomial variable is  y_{|MPAR|}
!      If MPAR.le.0  exponent polynomial connected to limiting inverse-
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
!++     READ(5,*) IPOTL, MPAR, NSR, NCMM, NVARB, IBOB, DSCM, REQ
!++     IF(IPOTL.GE.4) READ(5,*) (MMLR(I), CMM(I),I= 1,NCMM)
!++     IF((IPOTL.EQ.4).OR.(IPOT.EQ.7)) READ(5,*) (MMLR(I),CMM(I),I= 1,MPAR)
!++     IF(NVARB.GT.0)  READ(5,*) (PARM(I), I=1,NVARB)
!++     IF(IBOB.GT.0) THEN
!++         READ(5,*) MN1R, MN2R, PAD, MAD, NU1, NU2, PNA, NT1, NT2
!++         IF(PAD.GT.0) THEN
!++              IF(NU1.GE.0) READ(5,*) U1INF, (U1(I), I=0,NU1)
!++              IF(NU2.GE.0) READ(5,*) U2INF, (U2(I), I=0,NU2)
!++              IF(NT1.GE.0) READ(5,*) T1INF, (T1(I), I=0,NT1)
!++              IF(NT2.GE.0) READ(5,*) T2INF, (T2(I), I=0,NT2)
!++           ELSE
!++              IF(NU1.GE.0) READ(5,*) U1INF,(U1(I), I=0,NU1), Aad1, Rad1
!++              IF(NU2.GE.0) READ(5,*) U2INF,(U2(I), I=0,NU2), Aad2, Rad2
!++              IF(NT1.GE.0) READ(5,*) T1INF,(T1(I), I=0,NT1), Ana1, Rna1
!++              IF(NT2.GE.0) READ(5,*) T2INF,(T2(I), I=0,NT2), Ana2, Rna2
!++           ENDIF
!++         ENDIF
!++     ENDIF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      PW= 2.D0
      IF((NCN1.GT.0).AND.(NCN1.NE.2)) PW= 2.D0*NCN1/(NCN1-2.D0)
!     IF(DFLOAT(NCN1).LT.(2.d0*PRV + 1.9999999d0)) THEN
      IF(DBLE(NCN1).LT.(2.d0*PRV + 1.9999999d0)) THEN
          WRITE(6,629) (2.d0*PRV + 2.),NCN1
  629 FORMAT(/  ' *** Note that Radial variable power \alpha optimal for&
     &   NLR=',f5.2,' > NCN=', i2)
      ENDIF
!** Convert potential in [cm-1] to form appropriate for SCHRQas
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
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** For 2-potential Franck-Condon factor calculation, get the second
!  potential in this second call to PREPOT (uses the same parameter
!  reading sequence so exhaustively described immediately above).
!          CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,
!     1                                             VLIM2,V2,CNN2,NCN2)
           CALL PREPOT(LRPT,IAN1,IAN2,IMN1,IMN2,NPP,IOMEG2,RVB,RRM22,   &
     &  VLIM2,V2,CNN2,NCN2,IPOTL,PPAR,QPAR,NSR,NLR,IBOB,DSCM,REQ,RREF,  &
     &  PARM,MMLR,CMM,NCMM,IVSR,IDSTT,RHOAB)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!** Convert potential (in (cm-1)) to form appropriate for SCHRQas
          DO  I=1,NPP
              V2BZ(I)= V2(I)*BFCT
              V2(I)= V2BZ(I)*DRDY2(I) + FAS(I)
              RM22(I)= RRM22(I)*DRDY2(I)
          ENDDO
      ENDIF
!
!** NLEV1 is the no. of levels {v=IV(i), J=IJ(i)} of potential-1 which
!   we wish to find.
!* IF(NLEV1=0) calculate (and print?) potential, and then quit.
!* If read-in value of  NLEV1 < 0 , and program (attempts to) find all
!  vibrational levels of potential-1 up to  v = |NLEV1|.  [This case
!  assumes  AUTO1 > 0.]
!** If (AUTO1.gt.0) read in only (v,J) quantum numbers of NLEV1 desired
!  level & subroutine ALFas tries to locate them (normal preferred case).
!   If (AUTO1.le.0) also read in trial energy for each level.  In this
!      case, the  NLEV.le.0  option does not work.
!   If (AUTO1.le.0) and vib. quant. No. IV < 0, seek level nearest to
!      given trial energy but for whatever q. number shows up
!** If(LCDC.gt.0) calculate centrifugal distortion constants for each
!  level via the Tellinghuisen implementation of Hutson's method.
!  If(LCDC.GE.2) also write them compactly to channel-9.
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
!  NUMPOT.le.1), or to NLEV2 (see below) vib. levels of potential-2
!  (for NUMPOT.ge.2), and if (|LXPCT| > 3) write the overall
!  off-diagonal matrix elements on channel-8.
!* For  |LXPCT| > 4  also write to channel-7 the matrix elements of the
!  individual powers of the chosen distance coordinate (or radial fx.)
!* For  |LXPCT| > 5  WRITE(7,xx) only those matrix element components.
!** IF(NJM > 0), for each vibrational level, calculate all rotational
!  levels up to  J=NJM  or predissociation, whichever comes first.
!  Note that  AUTO1.le.0  forces  NJM= 0
!** When (NJM.GT.0) increase J in increments of JDJR.
!** IF(IWR.NE.0) print error & warning descriptions
!  IF(IWR.GE.1) also print final eigenvalues & node count.
!  IF(IWR.GE.2) also show end-of-range wave function amplitudes
!  IF(IWR.GE.3) print also intermediate trial eigenvalues, etc.
!** IF(LPRWF.GT.0) print wave function every LPRWF-th  point.
!** IF(LPRWF.LT.0) compactly write to channel-7 every |LPRWF|-th
!  wave function value.  **  A lead "card" identifies the level, gives
!  the position of 1-st point and radial mesh, & states No. of  points
!=======================================================================
!** INNOD1 specified wave fx. initiation at RMIN.  Normal case of
!  INNOD1 > 0  gives initiation with wave fx. node @ RMIN.
!  INNOD1.le.0  give initiation with  zero slope @ RMIN.  This determines
!    symmetric eigenfunctions for rare special case when input potential
!    is half of a precisely symmetric potential with mid-point at RMIN.
!-----------------------------------------------------------------------
! ... comment this out if use first version of following READ statement
      INNOD1= 1
      INNOD2= INNOD1
!-----------------------------------------------------------------------
!     READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF, INNOD1
!-----------------------------------------------------------------------
!     READ(5,*) NLEV1, AUTO1, LCDC, LXPCT, NJM, JDJR, IWR, LPRWF
!-----------------------------------------------------------------------
!** SINNER specifies whether wave function matching occurs at outermost
!  (SINNER.le.0) or innermost well turning point, to facilitate finding
!  inner vs. outer wells of a double well potential; Normally controlled
!  automatically,
!!!
      IF(LPRWF.LT.0) WRITE(10,605) TITL
!!!
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
!** Read the vibrational & rotational quantum numbers IV(i) & IJ(i) [and
!  if AUTO1.le.0 also trial energy GV(I)] of the NLEV levels to be found
!** For  IV(i)  values < -10,  SCHRQ  imposes a hard wall boundary
!  condition (i.e., a node) at mesh point # |-IV(i)| .
!-----------------------------------------------------------------------
!     IF(AUTO1.GT.0) READ(5,*) (IV(I), IJ(I), I= 1,NLEV)
!     IF(AUTO1.LE.0) READ(5,*) (IV(I), IJ(I), GV(I), I= 1,NLEV)
!-----------------------------------------------------------------------
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
!** MORDR is the highest power of the radial function (or distance
!  coordinate whose expectation values or matrix elements are to be
!  calculated.  Program currently dimensioned for (MORDR.LE.10).  To
!  calculate only F-C Factors (when LXPCT>2), set  MORDR = -1.
!** IRFN & DREF specify the definition of the radial function or
!  distance coordinate  RFN(R), powers of which are averaged over in
!  expectation value or matrix element calculations.
!* If(IRFN .le. -10) utilize the USER-CHOSEN and CODED radial function
!                 generated in Lines #500-504 (below)
!* If(IRFN = -4)  the function is a power series in  R  premultiplying a
!          first derivative operator acting on the wavefx of Potential-2
!* If(IRFN = -3)  the function is the inverse power   1/R**3
!* If(IRFN = -2)  the function is the inverse power   1/R**2
!* If(IRFN = -1)  the function is the Dunham coordinate  X=(R-DREF)/DREF
!* If(IRFN =  0)  the function  RFN(R)  is the distance  R  itself.
!* If(IRFN = 1-9)  use the Surkus-type variable
!                 X=(R^p - DREF^p)/(R^p + DREF^p)  where  p= IRFN
!* For  IRFN = -1 or 1-9,  if  DREF.gt.0  the read-in DREF value is the
!  reference length used to define the distance coordinate, while
!  if  DREF.le.0  determine the value of this reference length by
!  requiring that the expectation value  X**1  of the distance
!  coordinate for the first level considered be identically zero.
!* IF(IRFN > 10) define  RFN(R)   by reading in, interpolating over (and
!  extrapolating beyond) IRFN read-in values of some known radial
!  (transition moment) function, whose asymptotic value is DREF.  Do
!  this in using the same read statements and GENINT subroutine calls
!  used for generating a numerical potential.
      IF((LXPCT.NE.0).AND.(LXPCT.NE.-1)) THEN
!-----------------------------------------------------------------------
!         READ(5,*) MORDR, IRFN, DREF
!-----------------------------------------------------------------------
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
!** Overall calculated matrix elements are for a power series in the
!  radial function  RFN(i)  (specified by IRFN & DREF), so must read
!  coefficients  DM(J)  of this power series.
!-----------------------------------------------------------------------
!                 READ(5,*) (DM(J), J= 0,MORDR)
!-----------------------------------------------------------------------
                ELSE
                  DO  I= 1,NPP
                      RFN(I)= 1.D0
                      ENDDO
                  IF(MORDR.LT.0) WRITE(6,617)
                ENDIF
            ENDIF
!** Define radial function (distance coordinate) operator  RFN(R)  for
!  expectation values or matrix elements.
!** First ... for matrix elements of an operator consisting of a power
!    series in  R  premultiplying the radial derivative of the wavefx.
          IF(IRFN.EQ.-4) WRITE(6,650) MORDR
          IF(MORDR.GT.0) THEN
!** If  RFN(R)  is the distance itself ...
              IF(IRFN.EQ.0) THEN
                  WRITE(6,614)
                  DO  I= 1,NPP
                      RFN(I)= RVB(I)
                      ENDDO
                  ENDIF
              IF((IRFN.EQ.-2).OR.(IRFN.EQ.-3)) THEN
                  IF((IRFN.EQ.0).OR.(IRFN.EQ.-2).OR.(IRFN.EQ.-3))       &
     &                                                      DREF= 0.D0
!** If  RFN(R)  is   1/(distance)**|IRFN|  ....
                  J= -IRFN
                  WRITE(6,616) -IRFN
                  DO  I= 1,NPP
                      RFN(I)= 1.d0/RVB(I)**J
                      ENDDO
                  ENDIF
!%% Any other user-defined matrix element argument radial function
!   may be introduced to the code here, and invoked by:  IRFN= -4
!  Note that the existing  RVB(i)  array is the radial distances  R .
              IF(IRFN.LE.-10) THEN
!&& Illustrative user-defined analysis RFN(R) function
!&&               WRITE(6,*) 'Print description of function introduced'
!&&               WRITE(6,*) 'Use Freedman Pade DMF for CO'
!&&               DO  I= 1,NPP
!%%---------------------------------------------------------------------
!%%                  RFN(I)= {calculate users chosen radial function}
!%% Freedman's DMF for CO  ---------------------------------------------
!%%                  RFN(I)= {calculate users chosen radial function}
!&&   data coeff_new /-24.6005858d0,-109.5939637d0,-524.8233323d0,
!&&  +                 4.5194090d0,19.7954955d0,
!&&  +                 6.6011985d0,19.7206690d0/
!&&   dm = -0.122706d0*(1.+coeff(1)*x+coeff(2)*x*x+coeff(3)*x**3)/
!&&  +                   (1.+coeff(4)*x+coeff(5)*x*x+coeff(6)*x**3
!&&  +                    + coeff(7)*x**6)
!&&---------------------------------------------------------------------
!&&                   XX= RFN(I)/1.128322714d0 - 1.d0
!&&                   RFN(I)= -0.122706d0*(1.d0+ XX*(-24.6005858d0
!&&  1                  + XX*(-109.5939637d0 + XX*(-524.8233323d0))))/
!&&  2                    (1.d0 + XX*(4.5194090d0 + XX*(19.7954955d0 +
!&&  3                        XX*(6.6011985d0 + 19.7206690d0*XX**3))))
!&&                   ENDDO
                  ENDIF
              IF((IRFN.EQ.-1).OR.((IRFN.GE.1).AND.(IRFN.LE.9))) THEN
!** If  RFN(R)  is the Dunham or Surkus-type distance coordinate
                  IF(IRFN.EQ.-1) WRITE(6,615)
                  IF((IRFN.GE.1).AND.(IRFN.LE.9)) WRITE(6,611)          &
     &                                             IRFN,IRFN,IRFN,IRFN
                  IF(DREF.GT.0.D0) THEN
                      DREFP= DREF**IRFN
                      WRITE(6,613) DREF
                      DO  I=1,NPP
                          XX= YMINN+I*YH
                          IF(IRFN.EQ.-1) RFN(I)= (RVB(I)- DREF)/DREF
                          IF(IRFN.GE.1) RFN(I)= (RVB(I)**IRFN- DREFP)   &
     &                                        /(RVB(I)**IRFN + DREFP)
                          ENDDO
                    ELSE
                      WRITE(6,610)
                    ENDIF
                  ENDIF
!
!QQQQQQQQQQQ new ... the following not fully tested ...
!
!** If  RFN(R)  is defined by interpolating over read-in points, use
!  potential generating routine to do interpolation/extrapolation.
              IF(IRFN.GE.10) THEN
                  MORDR= 1
                  DM(0)= 0.d0
                  DM(1)= 1.d0
                  WRITE(6,603)
!** If the expectation value/matrix element radial function argument to
!   be defined by interpolating/extrapolating over read-in points, then
!   read input analogous to that for a pointwise potential, and then call
!   interpolation/extrapolation routine GENINT (from PREPOT package)
!* NRFN is the number of points [XIF(i),YIF(i)] to be read in
!* RFLIM  is the limiting asymptotic value imposed on the extrapolation
!* Interpolate with NUSEF-point piecewise polynomials (or splines for
!    NUSEF.le.0), which are extrapolated to the asymptote as specified by
!    parameters ILRF, NCNF & CNNF (see read #20).
!* RFACTF - factor converts read-in distances XIF(i) to angstroms
!* MFACTF - factor converts read-in moment values YIF(i) to debye.
!-----------------------------------------------------------------------
!                 READ(5,*) NRFN, RFLIM
!                 READ(5,*) NUSEF, ILRF, NCNF, CNNF
!                 READ(5,*) RFACTF, MFACTF
!                 READ(5,*) (XIF(I), YIF(I), I= 1,NRFN)
! If you uncomment the above, you better also uncomment the
! initialization to 0 below.
!-----------------------------------------------------------------------
                  MFACTF=0
                  RFACTF=0
                  WRITE(6,810) NRFN, RFLIM
                  IF(NUSEF.GT.0) WRITE(6,812) NUSEF, NRFN
                  IF(NUSEF.LE.0) WRITE(6,814) NRFN
                  IF((ILRF.GT.1).AND.(DABS(CNNF).GT.0.D0))              &
     &                                         WRITE(6,816) CNNF, NCNF
                  WRITE(6,818) RFACTF, MFACTF
                  NROW= (NRFN+ 2)/3
                  DO  J= 1,NROW
                      WRITE(6,820) (XIF(I), YIF(I), I=J, NRFN, NROW)
                  ENDDO
                  DO  I= 1,NRFN
                      XIF(I)= XIF(I)*RFACTF
                      YIF(I)= YIF(I)*MFACTF
              ENDDO
  810 FORMAT(' Transition moment function defined by interpolating over'&
     &  ,I4,' read-in points'/5x,'and approaching the asymptotic value',&
     &  f12.6)
  812 FORMAT(' Perform',I3,'-point piecewise polynomial interpolation ov&
     &er',I5,' input points' )
  814 FORMAT(' Perform cubic spline interpolation over the',I5,' input p&
     &oints' )
  816 FORMAT('- Beyond read-in points extrapolate to limiting asymptotic&
     & behaviour:'/20x,'Y(R)  =  Y(lim) - (',D16.7,')/R**',I2)
  818 FORMAT(' Scale input points:  (distance)*',1PD16.9,'   &  (moment)&
     &*',D16.9/4x,'to get required units  [Angstroms & debye]'/         &
     &  3('      R(i)         Y(i)  ')/3(3X,11('--')))
  820 FORMAT((3(F12.6,F13.6)))
                  IR2F= 0
                  CALL GENINT(LRPT,NPP,RVB,RFN,NUSEF,IR2F,NRFN,XIF,YIF, &
     &                                           RFLIM,ILRF,NCNF,CNNF)
                  ENDIF
              ENDIF
          IF((MORDR.GE.0).AND.(IABS(IRFN).LE.9))                        &
     &                                  WRITE(6,602) (DM(J),J=0,MORDR)
          ENDIF
!** For matrix element calculation, couple each level of potential-1 to
!  up to (see below) NLEV2 other vibrational levels, subject to
!  rotational selection rules:  DELTA(J)= J2DL to J2DU with increment
!  J2DD (e.g., -1,+1,+2 for P- and R-branches).
!** If (AUTO2.gt.0) read in only (v,J) quantum numbers of desired levels
!  and trust subroutine ALFas to locate them (normal preferred case).
!   If (AUTO2.le.0) also read in a trial pure vib energy for each level.
!* For the one-potential case (NUMPOT.LE.1), automatically truncate to
!  avoid redundancy and consider only emission into these NLEV2 levels.
!* Trial level energies are generated internally.
!**  IV2(i) are the vibrational quantum numbers of the Potential-2
!  levels for which matrix elements are desired.
!**  ZK(IV2(i),0) are the associated pure vibrational trial energies
!  (which are only read in if AUTO2.le.0!)
!=======================================================================
!** INNOD2 specified wave fx. initiation at RMIN.  Normal case of
!  INNOD2 > 0  gives initiation with wave fx. node @ RMIN.
!  INNOD2.le.0  give initiation with  zero slope @ RMIN.  This determines
!    symmetric eigenfunctions for rare special case when input potential
!    is half of a precisely symmetric potential with mid-point at RMIN.
!=======================================================================
      IF(IABS(LXPCT).GE.3) THEN
!-----------------------------------------------------------------------
!cc       READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD, INNOD2
!-----------------------------------------------------------------------
!         READ(5,*) NLEV2, AUTO2, J2DL, J2DU, J2DD
!-----------------------------------------------------------------------
          IF(NLEV2.GT.VIBMX) NLEV2= VIBMX
          IF(NLEV2.LE.0) THEN
              WRITE(6,644) NLEV2
              GO TO 997
!             STOP
          ENDIF
!----------------------------------------------------------------------
!         IF(AUTO2.GT.0) READ(5,*) (IV2(I), I= 1,NLEV2)
          IF(AUTO2.LE.0) THEN
!             READ(5,*) (IV2(I), ZK2(I,1), I= 1,NLEV2)
!----------------------------------------------------------------------
!** Give potential-2 trial energy the correct vibrational label
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
          IF(AUTO2.GT.0) WRITE(6,634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),     &
     &                                                     I= 1,NLEV2)
          IF(AUTO2.LE.0) WRITE(6,6634) J2DL,J2DU,J2DD,NLEV2,(IV2(I),    &
     &                                      ZK2(IV2(I),0), I= 1,NLEV2)
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
!
      IF(AUTO1.GT.0) THEN
!** If using automatic search for desired levels, subroutine ALFas gets
!  eigenvalues ZK1(v,0) for desired vibrational levels of Potential-1,
!  centrifugally-distorted to J=JREF.
          EJREF= JREF*(JREF+1)*YH**2
!** Replace  [J(J+1)] by  [J(J+1) + |IOMEG1|]  for Li2(A) and like cases.
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
!** Option to search for very highest level (within 0.0001 cm-1 of Disoc)
              EO= VLIM1- 0.0001d0
              KV= IV(1)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f (1)'
              WRITE(6,*) 'Entering schrq.f (1)'
              WRITE(6,*) ''
              CALL SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,              &
     &      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
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
              IF((IABS(LXPCT).GT.2).AND.(NUMPOT.EQ.1)) VMAX=            &
     &                                                MAX(VMAX1,VMAX2)
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
              CALL ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,VMAX,AFLAG,ZMU,  &
     &                               EPS,GV,BFCT,INNOD1,INNR1,IWR)
              VMAX1= VMAX
          ENDIF
!** Get band constants for v=0-VMAX1 for generating trial eigenvalues
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
              CALL SCHRQas(KV,JREF,EO,GAMA,PMAX1,VLIM1,VJ,              &
     &     WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,WARN,LPRWF)
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
!** If using automatic location for levels of potential-2 (AUTO2 > 0)
!  for matrix element calculation, also need Potential-2 band constants
!  (rotational energy derivatives) ... again, calculate them at J=JREF
              IF(NUMPOT.GT.1) THEN
                  AFLAG= JREF
                  DO  I= 1,NPP
                      VBZ(i)= V2BZ(I) + EJREF*RRM22(I)
                      VJ(I)= V2(I) + EJREF*RM22(I)
                      ENDDO
                  CALL ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,VMAX2,AFLAG, &
     &                            ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
                  ENDIF
              ENDIF
          DO  ILEV2= 1,NLEV2
              IF(NUMPOT.EQ.1) THEN
!** For matrix elements within a single potl., copy above band constants
                  DO  M= 0,7
                      ZK2(IV2(ILEV2),M)= ZK1(IV2(ILEV2),M)
                      ENDDO
                ELSE
! ... otherwise, generate them (as above) with SCHRQ & CDJOEL
                  KV= IV2(ILEV2)
                  IF(AUTO2.GT.0) EO= GV(KV)
                  IF(AUTO2.LE.0) EO= ZK2(KV,0)
                  INNER= INNR2(KV)
                  WRITE(6,*) ''
                  WRITE(6,*) 'Exiting level.f (3)'
                  WRITE(6,*) 'Entering schrq.f (3)'
                  WRITE(6,*) ''
                  CALL SCHRQas(KV,JREF,EO,GAMA,PMAX2,VLIM2,VJ,          &
     &     WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD2,INNER,WARN,LPRWF)
                  CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,WF1,       &
     &                                                      RM2,RCNST)
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
!
!===== Begin Actual Potential-1 Eigenvalue Calculation Loop Here =======
!** Loop to compute eigenvalues ... etc. for NLEV levels of Potential-1
      DO 190 ILEV1= 1,NLEV
          KV= IV(ILEV1)
          IF(KV.LT.0) EXIT
          NJMM= MAX(NJM,IJ(ILEV1))
          JROT= IJ(ILEV1)- JDJR
          IQT= 0
          JCT= 0
!** If NJM > IJ(ILEV1) loop over range of rotational levels too
          DO  JLEV= IJ(ILEV1),NJMM,JDJR
              JROT= JROT+ JDJR
              EJ= JROT*(JROT+1) - SOMEG1
              IF(IOMEG1.GE.99) EJ= JROT*JROT - 0.25D0
!** If   IOMEG < 0   centrifugal term is  [J(J+1) + |IOMEG|]
!             IF(IOMEG1.LT.0) EJ= JROT*(JROT+1) - DFLOAT(IOMEG1)
              IF(IOMEG1.LT.0) EJ= JROT*(JROT+1) - DBLE(IOMEG1)
!** If appropriate (AUTO1>0) use ALFas results to generate trial eigenvalue
              IF(AUTO1.GT.0) THEN
                  EO= ZK1(KV,0)
                  DEJ= EJ- EJREF
                  EJP= 1.d0
                  DO M= 1,7
                      EJP= EJP*DEJ
                      EO= EO+ EJP*ZK1(KV,M)
                      ENDDO
                ELSE
!... otherwise - use read-in trial energy
                  IF(IV(ILEV1).LT.VIBMX) EO= ZK1(IV(ILEV1),0)
                  IF(IV(ILEV1).GE.VIBMX) EO= GV(ILEV1)
                ENDIF
              IF((AUTO1.LE.0).AND.(DABS(ZK1(IV(ILEV1),0)).LE.0.d0)) THEN
                  CALL SCATTLEN(JROT,SL,VLIM1,V1,WF1,BFCT,YMIN,YH,NPP,  &
     &                 CNN1,NCN1,IWR,LPRWF)
                  IF(NUMPOT.EQ.1) GOTO 2
                  GOTO 104
                  ENDIF
! ... or if JLEV > IJ(ILEV1) ... use local Beff to estimate next level
              IF(JLEV.GT.IJ(ILEV1)) THEN
                  BEFF= 0.d0
                  DO  I= NBEG,NEND
                      BEFF= BEFF+ WF1(I)**2*RM2(I)
                      ENDDO
                  BEFF= BEFF*YH*BvWN
                  EO=  ESLJ(JCT)+ (2*JLEV+ 1- JDJR)*JDJR*BEFF
                  ENDIF
!** Now add centrifugal term to get effective (radial) potential
              EJ= EJ*YH**2
              DO  J= 1,NPP
                  VBZ(J)= V1BZ(J) + EJ*RRM2(J)
                  VJ(J)= V1(J) + EJ*RM2(J)
                  ENDDO
!** Set wall outer boundary condition, if specified by input IV(ILEV1)
              IF(KV.LT.-10) THEN
                  WF1(-IV(ILEV1))= 0.D0
                  WF1(-IV(ILEV1)-1)= -1.D0
                  ENDIF
              KVIN= KV
              IF(AUTO1.GT.0) INNER= INNR1(KV)
              IF(SINNER.NE.0) INNER= SINNER
!** Call SCHRQ to find Potential-1 eigenvalue EO and eigenfn. WF1(i)
              WRITE(6,*) ''
              WRITE(6,*) 'Exiting level.f (4)'
              WRITE(6,*) 'Entering schrq.f (4)'
              WRITE(6,*) ''
! The next CALL SCHRQas for a-state failes during v=0, but EO's from when
! Entering schrq.f (2) were already accurate, so the next step is not necessary.
  100         CALL SCHRQas(KV,JROT,EO,GAMA,PMAX1,VLIM1,VJ,              &
     &      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
              IF(KV.LT.0) THEN
!** SCHRQ  error condition is  (KV.LT.0) .
                  IF(NJM.GT.IJ(ILEV1)) THEN
! ... in automatic search for ever-higher J levels
                      IF(IQT.LE.0) THEN
! ... try one more time with E(trial) slightly below barrier maximum
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
              IF((KV.NE.KVIN).AND.                                      &
     &                        ((AUTO1.GT.0))) THEN
!** If got wrong vib level, do a brute force ALFas calculation to find it.
                  KV= KVIN
                  AFLAG= JROT
                  CALL ALFas(NPP,YMIN,YH,NCN1,VJ,WF1,VLIM1,KV,AFLAG,ZMU,&
     &                                EPS,GV,BFCT,INNOD1,INNR1,IWR)
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
!** If desired, calculate rotational & centrifugal distortion constants
              IF(LCDC.GT.0) THEN
                  IF((IOMEG1.GT.0).AND.(JROT.EQ.0)) THEN
!** Calculate 'true' rotational constants for rotationless IOMEG>0 case
                      WRITE(6,*) ''
                      WRITE(6,*) 'Exiting level.f (5)'
                      WRITE(6,*) 'Entering schrq.f (5)'
                      WRITE(6,*) ''
                      CALL SCHRQas(KV,0,EO,GAMA,PMAX1,VLIM1,V1,         &
     &      WF1,BFCT,EPS,YMIN,YH,NPP,NBEG,NEND,INNOD1,INNER,IWR,LPRWF)
                      CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,V1,       &
     &                                                  WF1,RM2,RCNST)
                  ELSE
!** Calculate rotational constants for actual (v,J) level of interest.
                      CALL CDJOELas(EO,NBEG,NEND,BvWN,YH,WARN,VJ,       &
     &                                                  WF1,RM2,RCNST)
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
              IF(((LXPCT.EQ.1).OR.(IABS(LXPCT).EQ.2)).OR.               &
     &                  ((IABS(LXPCT).GT.2).AND.((IRFN.EQ.-1).OR.       &
     &          (IRFN.GE.1).AND.(IRFN.GE.9)).AND.(DREF.LE.0.d0))) THEN
!** Calculate various expectation values in LEVXPC
                  CALL LEVXPC(KV,JROT,EO,GAMA,NPP,WF1,RFN,VBZ,VLIM1,    &
     &                    YH,DREF,NBEG,NEND,LXPCT,MORDR,DM,IRFN,BFCT)
                  IF((LXPCT.GT.0).AND.(MORDR.GT.0)) WRITE(6,632)
              ENDIF
!             WRITE(6,*) 'Hello! (for debugging) 4'
  104         IF((IABS(LXPCT).LE.2).OR.(NLEV2.LE.0)) GO TO 122
!** If desired, now calculate off-diagonal matrix elements, either
!  between levels of different potentials, IF(NUMPOT.GE.2), or between
!  levels of a single potential, for (NUMPOT.LE.1).
!** First prepare centrifugally distorted potential, trial energy, etc.,
!  and calculate second wave function and matrix element(s)
              DO 120 ILEV2= 1,NLEV2
!** For case of a single potential, avoid redundancy by considering
!  only emission
                  IF((NUMPOT.LE.1).AND.(IV2(ILEV2).GT.KV)) GO TO 120
!** Loop over J2's allowed by given selection rule.
                  DO 116 IJD= J2DL,J2DU,J2DD
                      KV2= IV2(ILEV2)
                      KVIN= KV2
                      JROT2= JROT+IJD
                      IF(JROT2.LT.0) GO TO 116
                      IF((NUMPOT.LE.1).AND.(IV2(ILEV2).EQ.KV).AND.      &
     &                                      (JROT2.GT.JROT)) GO TO 116
                      EJ2= JROT2*(JROT2+1)- SOMEG2
                      IF(IOMEG2.GE.99) EJ2=JROT2**2-0.25D0
!... allow for weird Li2(A) and Li2(c) potential cases
!                     IF(IOMEG2.LT.0) EJ2=JROT2*(JROT2+1)-DFLOAT(IOMEG2)
                      IF(IOMEG2.LT.0) EJ2=JROT2*(JROT2+1)-DBLE(IOMEG2)
                      EO2= ZK2(KV2,0)
                      DEJ= EJ2- EJREF
                      EJP= 1.d0
!** Use calculated state-2 CDC's to predict trial eigenvalue
                      DO  M= 1,7
                          EJP= EJP*DEJ
                          EO2= EO2+ EJP*ZK2(KV2,M)
                      ENDDO
!** Now ... update to appropriate centrifugally distorted potential
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
  110                 CALL SCHRQas(KV2,JROT2,EO2,GAMA,PMAX2,VLIM2,VJ,   &
     &    WF2,BFCT,EPS,YMIN,YH,NPP,NBEG2,NEND2,INNOD2,INNER,IWR,LPRWF)
                      IF(KV2.NE.KVIN) THEN
                          IF(KV2.LT.0) GO TO 114
!** Using CDC's to estimate trial eigenvalue failed:
                          ICOR= ICOR+1
                          IF(ICOR.LE.2) THEN
! ... first correction attempt ... use semiclassical dv/dE to improve
                              GB= -1.d0
                              GI= -1.d0
! ... hey RJ!!  shouldn't you update this using SCECOR ?
                              WV= 0.d0
                              XX= EO2*BFCT
                              DO  I= NBEG2,NEND2
                                  GBB= GB
                                  GB= GI
                                  GI= XX-VJ(I)
                                  IF((GBB.GT.0.d0).AND.(GI.GT.0.d0))    &
     &                                          WV= WV+ 1.d0/DSQRT(GB)
                              ENDDO
                              WV= 6.2832d0/(BFCT*WV)
                              EO2= EO2+ WV*(KVIN- KV2)
                              GO TO 110
                          ENDIF
                          WRITE(6,633) IV2(ILEV2),JROT2,KV2
! ... if that fails, do a brute force ALFas calculation to find it.
  114                     KV2= KVIN
                          AFLAG= JROT2
                          CALL ALFas(NPP,YMIN,YH,NCN2,VJ,WF2,VLIM2,KV2, &
     &                      AFLAG,ZMU,EPS,GV,BFCT,INNOD2,INNR2,IWR)
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
!                     IF((NUMPOT.LE.1).AND.(EO2.GT.(EO+EPS))) GO TO 120
                      CALL MATXEL(KV,JROT,IOMEG1,EO,KV2,JROT2,IOMEG2,   &
     & IRFN,EO2,NBEG2,NEND2,LXPCT,MORDR,DM,RH,NDIMR,DRDY2,WF1,WF2,RFN)
  116                 CONTINUE
!** End of Potential-2 rotational selection level loop
  120             CONTINUE
!++++ End of Potential-2 vibrational level matrix element loop +++++++++
!
  122         CONTINUE
              JCT= JCT+1
! ... check to avoid array overflow
              IF(JCT.GT.VIBMX) THEN
                  WRITE(6,637)  VIBMX
                  GO TO 997
!                 STOP
                  ENDIF
              JWR(JCT)= JROT
              ESLJ(JCT)= EO
              ENDDO
!             WRITE(6,*) 'Hello! (for debugging) 5'
!++ End of Potential-1 loop over NJM-specified J-sublevels
  130     IF(NJM.GT.IJ(ILEV1)) THEN
!** Print rotational sublevels generated for vibrational level  ILEV
              NROW=(JCT+4)/5
              WRITE(6,627) KV
              DO  J=1,NROW
                  WRITE(6,628) (JWR(I),ESLJ(I),I=J,JCT,NROW)
                  ENDDO
              WRITE(6,641)
              ENDIF
          ESOLN(ILEV1)= ESLJ(1)
  190     CONTINUE
!++ End of loop over the NLEV Potential-1 input levels
      WRITE(6,*) ''
      WRITE(6,*) 'SUMMARY (ALL ENERGIES IN CM-1):'
      IF(NLEV1.LT.0) THEN
          NROW=(NLEV+3)/4
          WRITE(6,623) NLEV,IJ(1)
          DO  J=1,NROW
              WRITE(6,630) (IV(I),ESOLN(I),I=J,NLEV,NROW)
          ENDDO
          IF((NLEV.GT.1).AND.(IJ(1).EQ.0).AND.(NCN1.GT.0)               &
     &                               .AND.(ESOLN(NLEV).LT.VLIM1)) THEN
!** In (NLEV1 < 0) option, estimate vD using the N-D theory result that:
!     (vD - v) {is proportional to} (binding energy)**((NCN-2)/(2*NCN))
              VDMV=1.D0/(((VLIM1-ESOLN(NLEV-1))/                        &
     &                         (VLIM1-ESOLN(NLEV)))**(1.D0/PW) - 1.D0)
!          CALL ADD_INFO('LEVEL_CHKSUM',[CHKSUM],1,2)
           CALL ADD_INFO('Vibrational energy:',ESOLN,NLEV,2)
! IF YOU WANT TO ONLY VERIFY THE ENERGY OF THE HIGHEST LEVEL:
!          CALL ADD_INFO('Last vibrational energy:',ESOLN(NLEV),1,2)
! IF YOU WANT TO ONLY VERIFY THE NUMBER OF LEVELS:
! All values must be reals, so if you want to check an integer, e.g. the
! number of iterations, then you must convert this to a real:
!          CALL ADD_INFO('Numer of vibrational levels:',NLEV,1,2)
!** Use empirical N-D Expression to predict number and (if there are
!  any) energies of missing levels
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
                      ESOLN(NLEV)= VLIM1 - (VLIM1 - ESOLN(NLEV-1))*     &
     &                                          (1.D0 - 1.D0/VDMV)**PW
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
!-------------------------------------------------------------------
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
  602 FORMAT( ' Coefficients of expansion for radial matrix element/expe&
     &ctation value argument:'/(5X,5(1PD14.6)))
  603 FORMAT(/' Expectation value/matrix element arguments are powers of&
     & a radial function'/5x,'defined by interpolating over read-in poin&
     &ts'//' Transition moment function:'/1x,9('==='))
 6604 FORMAT(/' !!! NOTE:  array dimension limit   NDIMR=',i5/'    preve&
     &nts preliminary mesh   YH=',f9.6,'  from spanning range [YMIN,YMAX&
     &]'/'    so increase mesh by factor',f7.4/)
  604 FORMAT(' Integrate from   Ymin=',f11.7,'   to   Ymax=',f7.4,      &
     & '  with mesh  YH=',f9.7/5x,'based on radial variable  yp(r;a)= (r&
     &^p - a^p)/(r^p + a^p)'/5x,'for   p=',f6.3,'    and   a=',F9.6/    &
     & ' Range corresponds to   Rmin=',f7.3,' [Angst]    to    Rmax= inf&
     &inity (!!)'/5x,'and   RH=',F10.7,' [Angst]   at   R=',f11.8//     &
     & ' Potential #1 for ',A2,'(',I3,')-',A2,'(',I3,')'/1x,32('='))
  605 FORMAT(/A78/40('=='):/' Generate   ZMU=',F15.11,'(u)',            &
     &  '   &   BZ=',1PD16.9,'((1/cm-1)(1/Ang**2))'/                    &
     &  10x,'from atomic masses:',0Pf16.11,'  & ',F16.11,'(u)')
  606 FORMAT(' E(v=',i3,', J=',i3,')=',f10.3,'   Bv=',F11.7,            &
     &  '  -Dv=',1PD12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,              &
     &  '   Mv=',D12.4,'   Nv=',D12.4,'   Ov=',D12.4)
 6065 FORMAT(' E(v=',i3,', J=',i3,')=',f11.7,'  Bv=',1PD11.4,           &
     &  '  -Dv=',D12.4,'   Hv=',D12.4/8x,'   Lv=',D12.4,                &
     &  '   Mv=',D12.4,'   Nv=',D12.4,'   Ov=',D12.4)
  607 FORMAT(/' Solve for the',i4,' vibration-rotation levels of Potenti&
     &al-1:'/'   (v,J) =',6('  (',i3,',',i3,')':)/(10x,6('  (',i3,',',  &
     &  i3,')':)))
 6607 FORMAT(/' Solve for',i4,' vibration-rotation levels of Potential-1&
     & using Trial energies:'/(3x,3('   E(',I3,',',I3,')=',             &
     & F11.2:)))
  608 FORMAT(/' State-',I1,' electronic angular momentum  OMEGA=',I2/   &
     &  9x,'yields centrifugal potential  [J*(J+1) -',F5.2,']/r**2' )
 6085 FORMAT(/' State-',I1,' electronic angular momentum  OMEGA=',I2/   &
     &  9x,'yields centrifugal potential  [J*(J+1) +',I2,']/r**2' )
  609 FORMAT('  Use centrifugal potential for rotation in two dimensions&
     &:   (J**2 - 1/4)/r**2')
  610 FORMAT(5X, 'where DREF defined by requiring  <X**1> = 0  for first&
     & level considered')
  611 FORMAT(/' Matrix element argument expansion vble is   X = ((r^',  &
     &  i1,' - DREF^',i1,')/(r^',i1,' + DREF^',i1,'))')
  612 FORMAT(/' Eigenvalue convergence criterion is   EPS=',1PD8.1,     &
     & '(cm-1)'/' Airy function at 3-rd turning point is quasibound oute&
     &r boundary condition')
  613 FORMAT(5X,'where reference length is held fixed at   DREF =',     &
     & F13.10,'(Angstroms)')
  614 FORMAT(/' Matrix element arguments are powers of the distance  r (&
     &in Angstroms)')
  615 FORMAT(/' Matrix element argument expansion variable is:    X = (r&
     & - DREF)/DREF')
  616 FORMAT(/' Matrix element arguments are powers of the squared inver&
     &se distance  X = 1/r**',i1)
  617 FORMAT(/' Matrix element argument is fixed as a constant = 1')
  618 FORMAT(' *** PROBLEM *** Searching for  v=',i3,' , J=',i3,        &
     & '  ALFas only found to  v=',i3)
  619 FORMAT(/' Find the',i4,' vibration-rotation levels:'/             &
     &  3('     v   J      E(v)   ')/3(2x,7('---')))
  620 FORMAT(/' An  n=',I2,'  N-D theory extrapolation from last 2 level&
     &s implies   vD =',F8.3)
  621 FORMAT(5X,'with the',I4,' missing level(s) predicted to be:'/     &
     &  4('     v     E(v)   ')/4(4x,7('--')))
  622 FORMAT(/' Search for highest bound  J=',i3,'  level finds  E(v=', &
     &  i3,') = VLIM -',1PD12.5/)
  623 FORMAT(/' Find',I4,' Potential-1 vibrational levels with  J=',i3/ &
     &  4('     v     E(v)   ')/4(4x,7('--')))
  624 FORMAT(4x,'Since the molecule is an ion with charge',SP,I3/6x,"use&
     & Watson's charge-adjusted reduced mass   mu = M1*M2/[M1 + M2 - (",&
     &  i2,')*me]')
  625 FORMAT(' For  J=',i3,', try to find the first',i4,' vibrational le&
     &vels of Potential-1')
  626 FORMAT(/' *** FAIL to find highest bound J=',i3,'  level from tria&
     &l   E = VLIM -',1PD11.4)
  627 FORMAT(/' For vibrational level  v =',I3,'   of Potential-1'/     &
     & 1X,5('  J',6X,'E',7X)/1X,5(7('--'),2X))
  628 FORMAT((1X,5(I3,F11.3,2X)))
  630 FORMAT((4(I6,F12.4:)))
  631 FORMAT((3(I6,I4,F13.5:)))
  632 FORMAT(1X,79('-'))
  633 FORMAT(' **** Caution: Search for   v=',I3,'   J=',i3,            &
     &  '  on potential-2 actually found   v=',I3)
  634 FORMAT(/' Using the rotational selection rule:  delta(J)=',       &
     & i3,' to',i2,' with increment',i2/'   calculate matrix elements fo&
     &r coupling to the',I4,' vibrational levels of'/                   &
     & '   Potential-2:   v =',14I4:/(21x,14i4:))
 6634 FORMAT(/' Using the rotational selection rule:  delta(J)=',       &
     & i3,' to',i2,' with increment',i2/'   calculate matrix elements fo&
     &r coupling to the',I4,' vibrational levels of'/'   Potential-2 usi&
     &ng trial energies:',2('   E(',I3,')=',F9.2:)/4('   E(',I3,')=',   &
     & F9.2:))
! 635 FORMAT(/' Get matrix elements between levels of Potential-1 (above
!    1) & Potential-2 (below)'/1X,39('--')/' For Potential #2:'/
!    2  1x,17('='))
  636 FORMAT(/' Calculate properties of the single potential described a&
     &bove')
  637 FORMAT(/' *** Array Dimension OVERFLOW ***   (Number of J sublevel&
     &s) > VIBMX=',i4)
  638 FORMAT('   and automatically increment  J  in steps of',i3, ' to a&
     & maximum value of',i4)
  641 FORMAT(1X,39('++'))
  644 FORMAT(/' *** Input data ERROR *** matrix element calculation need&
     &s  NLEV2=',i3,' > 0')
  650 FORMAT(/' Matrix element argument is radial first derivative opera&
     &tor premultiplied by'/5x,'a power series in  r  of order',i3)
  686 FORMAT(' Potential-',i1,' uses inner boundary condition of  zero v&
     &alue  at  RMIN')
  688 FORMAT(' Potential-',i1,' uses symmetric-well inner boundary condi&
     &tion  of zero slope at RMIN')
!c703 FORMAT(1X,I4,I5,F13.4,G13.5)
  703 FORMAT(1X,I4,I5,1PD20.11,D13.5)
  723 FORMAT(/A78/1x,'Output values of:  v, J, E & (Level Width)')
  724 FORMAT(//A78//'   v   J    E(v,J)     Width       <KE>',          &
     &  6x,'<M(r)>  &  <XI**k>  for k=1 to',i3/2x,38('=='))
  725 FORMAT(//A78//"   v'  J'",'  v"  J"     FREQ',"    <v',J'| XI**k",&
     &  ' |v",J">  for  k=0  to  MORDR=',i2/2x,37('=='))
! 824 FORMAT(//A78/30('==')/" Note that (v',J') &",' (v",J") strictly la
!    1bel the upper and lower levels, resp.,'/6x,'and  E(lower)=E"'/
!    2 ' but  E(2)-E(1)  is:  (energy of State-2 level) - (energy of Sta
!    3te-1 level)'//12x,'Band'/' dJ(J")',4x,7hv'   v",'  E(lower)  E(2)-
!    4E(1)  A(Einstein)   F-C Factor  ',13h<v'j'|M|v"j"> /
!    5 1x,3('--'),('   -------'),'  --------',3x,
!    6 4('--'),3x,11('-'),3x,11('-'),3x,11('-') )
!c811 FORMAT(//A78/30('==')//12x,"   v'","  J'",'    v"','  J"',
!c   1 '   position    E(upper)    E(lower)',16h   <v'j'|M|v"j">/
!c   2 1x,68('-') )
! 901 FORMAT(//A78/1x,62('==')/'   v    J',7x,'E',10x,'Bv',11x,'-Dv',
!    1  13x,'Hv',13x,'Lv',13x,'Mv',13x,'Nv',13x,'Ov'/1x,62('=='))
! 902 FORMAT(I4,I5,f25.15,f14.10,6(1PD15.7))
! 904 FORMAT(I4,I5,f25.15,1PD14.7,6(D15.7))
!     END
  997 RC = 0
      END SUBROUTINE LEVEL
