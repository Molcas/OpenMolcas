!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2018,2019, Roland Lindh                                *
!               2019, Mickael G. Delcey                                *
!               2019, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,JBNUM,EigVec)

use Index_Functions, only: nTri_Elem
use RASSI_AUX, only: iDisk_TDM
use kVectors, only: e_Vector, k_Vector, nk_Vector
use do_grid, only: Do_Lebedev
use nq_Grid, only: Pax
use Cntrl, only: DIPR, Do_Pol, Do_SK, ICOMP, IRREP, L_Eff, LOOPDIVIDE, LOOPMAX, lSym1, lSym2, LUTDM, MLTPLT, NPROP, nQuad, NSTATE, &
                 OSThr_Dipr, OSThr_QIPR, PNAME, PrRaw, PrWeight, PTYPE, QIPR, REDUCELOOP, RSPR, RSThr, TMGR_Thrs
use Symmetry_Info, only: MUL, nIrrep
use rassi_data, only: NBASF, NBST, NTDMZZ
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use RASSIWfn, only: wfn_SOS_TM
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Four, Half, Pi, auTofs, c_in_au, Debye, gElectron
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: PROP(NSTATE,NSTATE,NPROP)
integer(kind=iwp), intent(in) :: NSS, JBNUM(NSTATE)
real(kind=wp), intent(in) :: USOR(NSS,NSS), USOI(NSS,NSS), ENSOR(NSS), EigVec(NSTATE,NSTATE)
integer(kind=iwp) :: I, iCar, IDISK, IEMPTY, IEND, iEnd_, IGO, iGrp, IJ_, ijSO, IJSS(4), IMSS, IOFF(8), iOpt, iPrint, IPROP, iPrP, &
                     IPRTMOM(14), iQuad, iQuad_, ISF, ISM, ISO, ISS, ISSM, iStart_, iState, iSy12, ITYPE, iVec, iVec_, J, JEnd, &
                     jEnd_, jGrp, JMSS, Job, Job1, Job2, JSF, JSM, JSO, JSS, JSSM, JSTart, jStart_, K, KP, lRaw, lRaw_, MASK, &
                     MaxGrp1, MaxGrp2, NDIFF, NGROUP1, NGROUP2, NIP, NMAX2, nQuad_, NSCR, nTmp, NVEC
real(kind=wp) :: A, ANG, CST, Dummy(1), EDiff, EDiff_, EnSOR1, EnSOR2, F, F_Check, F_temp, kPhase(2), OSThr, R, R_Check, R_temp, &
                 RefEne, RKNorm, RNG, RNorm, Tau, TCPU1, TCPU2, TEMP, Thrs, ThrsParse, TM1, TM2, TM3, TM_2, TM_C(3), TM_I(3), &
                 TM_R(3), TWALL1, TWALL2, UK(3), wavevector(3), Weight
logical(kind=iwp) :: TMOgroup
character(len=52) :: STLNE2
character(len=8) :: LABEL
integer(kind=iwp), allocatable :: iMask(:), ISS_INDEX(:), iSSMask(:), jMask(:), jSSMask(:), TMOgrp1(:), TMOgrp2(:)
real(kind=wp), allocatable :: Aux(:,:), DXI(:,:,:), DXIM(:,:,:), DXR(:,:,:), DXRM(:,:,:), IP(:), OscStr(:,:), pol_Vector(:,:), &
                              RAW(:), Rquad(:,:), SCR(:,:), TDMZZ(:), TMI(:,:,:), TMP(:,:), TMR(:,:,:), VSOI(:,:), VSOR(:,:), &
                              WDMZZ(:)
#ifdef _HDF5_
integer(kind=iwp) :: ijSO_, ip_kVector, ip_TMi, ip_TMr, ip_W, nData, nij
real(kind=wp), allocatable :: Storage(:,:,:,:)
#endif
real(kind=wp), parameter :: AFACTOR = Two/c_in_au**3/(auTofs*1.0e-15_wp), AU2REDR = 200.0_wp*Debye, THRSH = 1.0e-10_wp
real(kind=wp), external :: DDot_

#define _TIME_TMOM_
#ifdef _TIME_TMOM_
call CWTime(TCpu1,TWall1)
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Before we start we need to backtransform the coefficients of the
! SO states from the basis of the SF states which diagonalize the
! SF Hamiltonian to the basis of the original (input) SF states.
! This since all transition moments, whether retrieved from disk or
! recomputed, are in the basis of the original SF states.

call mma_allocate(VSOR,NSS,NSS,Label='VSOR')
call mma_allocate(VSOI,NSS,NSS,Label='VSOI')

call USOTRANS(USOR,USOI,NSS,EigVec,NSTATE,VSOR,VSOI)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _TIME_TMOM_
call CWTime(TCpu1,TWall1)
#endif
! Compute transition strengths for spin-orbit states:

! Initial setup for exact operator

! printing threshold
OSTHR = 1.0e-5_wp
if (DIPR) OSTHR = OSTHR_DIPR
if (DIPR) write(u6,30) 'Dipole printing threshold changed to ',OSTHR
! Again to avoid total negative transition strengths
if (QIPR) then
  OSTHR = OSTHR_QIPR
  write(u6,49) 'Printing threshold changed to ',OSTHR,' since quadrupole threshold is given'
end if
! Rotatory strength threshold
if (RSPR) then
  write(u6,30) 'Rotatory strength printing threshold changed to ',RSTHR
else
  RSTHR = 1.0e-07_wp !Default
end if

! Reducing the loop over states - good for X-rays
! At the moment memory is not reduced

JEND = NSS
if (REDUCELOOP) then
  IEND = LOOPDIVIDE
  JSTART = LOOPDIVIDE+1
  if (LOOPMAX > 0) JEND = min(NSS,LOOPDIVIDE+LOOPMAX)
else
  IEND = NSS
  JSTART = 1
end if

!***********************************************************************
!                                                                      *
!     Start of section for transition moments                          *
!                                                                      *
!     This section has two parts. (1) for matrix elements computed by  *
!     Seward, i.e. for a specific wavevector, (2) for the computation  *
!     of the isotropic oscillator strength.                            *
!                                                                      *
!***********************************************************************

! Find the section of transition moments in the property list.

! The operator is split in 4 different component, each with three
! elements corresponding to differentiation in the x, y, and z
! direction. The four parts are labels as:
! TMOM  RS: The symmetric part of the real comp. of the op.
! TMOM  RA: The asymmetric part of the real comp. of the op.
! TMOM  IS: The symmetric part of the imaginary comp. of the op.
! TMOM  IA: The asymmetric part of the imaginary comp. of the op.

!***********************************************************************
!                                                                      *
!     Section (1): Computation of k specific oscillator strength.      *
!     Section (2): Computation of the isotropic oscillator strength.   *
!                                                                      *
!***********************************************************************

! Here we will use a Lebedev grid to integrate over all possible
! directions of the wave vector, k. The property integrals will be
! computed on the fly and traced with the density to generate the
! corresponding values in the PROP matrix.

! Find the slot on the one-electron file where we will store the
! on-the-fly generated property integrals.

IPRTMOM(:) = -1
do IPROP=1,NPROP
  if (PNAME(IPROP) == 'TMOM  RS') then
    if (IPRTMOM(0+ICOMP(IPROP)) == -1) IPRTMOM(0+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  IS') then
    if (IPRTMOM(3+ICOMP(IPROP)) == -1) IPRTMOM(3+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  RA') then
    if (IPRTMOM(6+ICOMP(IPROP)) == -1) IPRTMOM(6+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM  IA') then
    if (IPRTMOM(9+ICOMP(IPROP)) == -1) IPRTMOM(9+ICOMP(IPROP)) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM0  R') then
    if (IPRTMOM(13) == -1) IPRTMOM(13) = IPROP
  end if
  if (PNAME(IPROP) == 'TMOM0  I') then
    if (IPRTMOM(14) == -1) IPRTMOM(14) = IPROP
  end if
end do
if (any(IPRTMOM == -1)) then
  call mma_deallocate(VSOR)
  call mma_deallocate(VSOI)
end if

! Initiate the Seward environment

nDiff = 0
call IniSew(.false.,nDiff)

! Generate the quadrature points.

! In the spin-coupled case, wave functions are complex and there is
! not a simple relation between oscillator and rotatory strengths for
! k and -k vectors, but there is between the integrals computed in
! the spin-free basis, so we compute them only once and obtain separately
! the results for k and -k, given by kPhase

kPhase = [One,-One]
if (Do_SK) then
  nQuad = 1
  nVec = nk_Vector
  call mma_allocate(Rquad,4,nQuad,label='SK')
  if (.not. (PRRAW .or. PRWEIGHT)) kPhase(2) = Zero
else
  call unitmat(Pax,3)
  call Do_Lebedev(L_Eff,nQuad,Rquad,4)
  Rquad(4,:) = Half*Rquad(4,:)
  nVec = 1
end if
if (Do_Pol) call mma_allocate(pol_Vector,3,nVec*nQuad,Label='POL')

! Initialize for density matrices.

call mma_Allocate(TDMZZ,nTDMZZ,Label='TDMZZ')
call mma_Allocate(WDMZZ,nTDMZZ,Label='WDMZZ')
nSCR = nTri_Elem(NBST)
call mma_allocate(SCR,nSCR,4,LABEL='SCR')

! Allocate some temporary arrays for handling the
! properties of the spin-orbit states.

call mma_allocate(DXR,NSS,NSS,3,Label='DXR')
call mma_allocate(DXI,NSS,NSS,3,Label='DXI')
call mma_allocate(DXRM,NSS,NSS,3,Label='DXRM')
call mma_allocate(DXIM,NSS,NSS,3,Label='DXIM')
call mma_Allocate(TMP,NSS,NSS,Label='TMP')
call mma_allocate(TMR,NSS,NSS,3,Label='TMR')
call mma_allocate(TMI,NSS,NSS,3,Label='TMI')

! ALLOCATE A BUFFER FOR READING ONE-ELECTRON INTEGRALS
NIP = 4+nTri_Elem(NBST)
call mma_allocate(IP,NIP,Label='IP')
#ifdef _HDF5_

! Allocate vector to store all individual transition moments.
! We do this for
! all unique pairs ISO-JSO, iSO=/=JSO nTri_Elem(NSS-1)
!     all k-vectors (2*nQuad or 2*nVec)
!         we store:
!             the weight (1)
!             the k-vector (3)
!             we projected transition vector (real and imaginary parts) (2*3)

nIJ = nTri_Elem(nSS-1)
ip_w = 1
ip_kvector = ip_w+1
ip_TMR = ip_kvector+3
ip_TMI = ip_TMR+3
nData = ip_TMI+3-1
call mma_allocate(Storage,nData,2*nQuad,nIJ,nVec,label='Storage')
Storage(:,:,:,:) = Zero
#endif
!MGD group transitions
TMOgroup = .false.
ngroup1 = IEND
ngroup2 = JEND-JSTART+1
nmax2 = 1
if (REDUCELOOP .and. (TMGr_thrs >= Zero)) then
  TMOgroup = .true.
  THRS = TMGr_thrs
  i = IEND
  RefEne = 0
  TAU = -1
  ngroup2 = 1
  do j=JSTART,JEND
    if (ENSOR(J)-Refene > TAU) then
      NGROUP2 = NGROUP2+1
      Refene = ENSOR(J)
      ediff = Refene-ENSOR(I)
      TAU = ediff*THRS
    end if
  end do
  call mma_Allocate(TMOgrp2,NGROUP2,Label='TMOgrp2')
  ngroup2 = 0
  TAU = -1
  RefEne = 0
  do j=JSTART,JEND
    if (ENSOR(J)-Refene > TAU) then
      NGROUP2 = NGROUP2+1
      TMOgrp2(NGROUP2) = J
      Refene = ENSOR(J)
      ediff = Refene-ENSOR(I)
      TAU = ediff*THRS
    end if
  end do
  TMOgrp2(ngroup2+1) = NSS+1

  j = JSTART
  Refene = ENSOR(j)
  TAU = -1
  ngroup1 = 1
  do i=IEND,1,-1
    if (Refene-ENSOR(i) > TAU) then
      ngroup1 = ngroup1+1
      Refene = ENSOR(i)
      ediff = ENSOR(j)-Refene
      Tau = ediff*THRS
    end if
  end do
  call mma_Allocate(TMOgrp1,NGROUP1,Label='TMOgrp1')
  Ntmp = Ngroup1
  Ngroup1 = Ngroup1-1
  Refene = ENSOR(j)
  TAU = -1
  do i=IEND,1,-1
    if (Refene-ENSOR(i) > TAU) then
      TMOgrp1(ntmp) = i+1
      ntmp = ntmp-1
      Refene = ENSOR(i)
      ediff = ENSOR(j)-Refene
      Tau = ediff*THRS
    end if
  end do
  TMOgrp1(1) = 1
  maxgrp1 = maxval(TMOgrp1(2:ngroup1+1)-TMOgrp1(1:ngroup1))
  maxgrp2 = maxval(TMOgrp2(2:ngroup2+1)-TMOgrp2(1:ngroup2))
  nmax2 = maxgrp1*maxgrp2
end if

! Array for printing contributions from different directions

call mma_allocate(RAW,2*NQUAD*6*nmax2,Label='RAW')
LRAW = 1
call mma_allocate(OSCSTR,2,nmax2,Label='OSCSTR')
call mma_allocate(Aux,8,nmax2,Label='Aux')

do iVec=1,nVec

  if (Do_SK) then
    Rquad(1:3,1) = k_Vector(:,iVec)
    Rquad(4,1) = One   ! Dummy weight
  end if

  iPrint = 0
  IJSO = 0

  ThrSparse = 1.0e-12_wp

  call mma_allocate(iMask,NState,LABEL='iMask')
  call mma_allocate(jMask,NState,LABEL='jMask')
  call mma_allocate(iSSMask,NSS,LABEL='iSSMask')
  call mma_allocate(jSSMask,NSS,LABEL='jSSMask')

  call mma_allocate(ISS_INDEX,NState+1,LABEL='ISS_INDEX')
  ISS_INDEX(1) = 0
  do iState=1,NState
    Job = JBNUM(iState)
    ISS_INDEX(iState+1) = ISS_INDEX(IState)+MLTPLT(Job)
  end do

  do igrp=1,ngroup1

    if (TMOgroup) then
      istart_ = TMOgrp1(igrp)
      iend_ = TMOgrp1(igrp+1)-1
      ENSOR1 = Half*(ENSOR(istart_)+ENSOR(iend_))
    else
      istart_ = igrp
      iend_ = igrp
      ENSOR1 = ENSOR(istart_)
    end if
    ! Screening
    iMask(:) = 0
    iSSMask(:) = 0
    ISM = 0
    ISSM = 0
    ISFLoop: do ISF=1,NState
      do ISS=ISS_INDEX(ISF)+1,ISS_INDEX(ISF+1)
        do ISO=istart_,iend_
          Temp = VSOR(ISS,ISO)**2+VSOI(ISS,ISO)**2
          if (Temp > ThrSparse) then
            ISM = ISM+1
            iMask(ISM) = ISF
            iSSMask(ISSM+1:ISSM+ISS_INDEX(ISF+1)-ISS_INDEX(ISF)) = [(IMSS,IMSS=ISS_INDEX(ISF)+1,ISS_INDEX(ISF+1))]
            ISSM = ISSM+ISS_INDEX(ISF+1)-ISS_INDEX(ISF)
            cycle ISFLoop
          end if
        end do
      end do
    end do ISFLoop

    do jgrp=1,ngroup2

      if (TMOgroup) then
        jstart_ = TMOgrp2(jgrp)
        jend_ = TMOgrp2(jgrp+1)-1
        ENSOR2 = Half*(ENSOR(jstart_)+ENSOR(jend_))
      else
        jstart_ = jgrp+jstart-1
        jend_ = jgrp+jstart-1
        ENSOR2 = ENSOR(jstart_)
      end if
      EDIFF_ = ENSOR2-ENSOR1
      ! Screening
      jMask(:) = 0
      jSSMask(:) = 0
      JSM = 0
      JSSM = 0
      JSFLoop: do JSF=1,NState
        do JSS=ISS_INDEX(JSF)+1,ISS_INDEX(JSF+1)
          do JSO=jstart_,jend_
            Temp = VSOR(JSS,JSO)**2+VSOI(JSS,JSO)**2
            if (Temp > ThrSparse) then
              JSM = JSM+1
              jMask(JSM) = JSF
              jSSMask(JSSM+1:JSSM+ISS_INDEX(JSF+1)-ISS_INDEX(JSF)) = [(JMSS,JMSS=ISS_INDEX(JSF)+1,ISS_INDEX(JSF+1))]
              JSSM = JSSM+ISS_INDEX(JSF+1)-ISS_INDEX(JSF)
              cycle JSFLoop
            end if
          end do
        end do
      end do JSFLoop

      IJSS(1) = istart_
      IJSS(2) = iend_
      IJSS(3) = jstart_
      IJSS(4) = jend_

      write(STLNE2,'(A33,I5,A5,I5)') 'Trans. intensities for SO groups ',igrp,' and ',jgrp
      call StatusLine('RASSI: ',STLNE2)

      if (abs(EDIFF_) <= 1.0e-8_wp) cycle
      if (EDIFF_ < Zero) cycle
      IJSO = IJSO+1

      ! The energy difference is used to define the norm of the
      ! wave vector.

      rkNorm = abs(EDIFF_)/c_in_au

      ! For the case the energy difference is negative we
      ! need to change the sign of the B.s term to get
      ! consistency.

      cst = sign(Half,EDIFF_)*gElectron

      ! Iterate over the quadrature points.

      ! Initialize output arrays

      OscStr(:,:) = Zero
      RAW(:) = Zero
      Aux(:,:) = Zero

      do iQuad=1,nQuad
        iVec_ = (iVec-1)*nQuad+iQuad

        ! Generate the wavevector associated with this quadrature
        ! point and pick up the associated quadrature weight.

        UK(:) = Rquad(1:3,iQuad)
        wavevector(:) = rkNorm*UK(:)

        ! Note that the weights are normalized to integrate to
        ! 4*pi over the solid angles.

        Weight = Rquad(4,iQuad)
        if (.not. Do_SK) Weight = Weight/(Four*Pi)

        ! Generate the polarization vector

        if (Do_Pol) then
          pol_Vector(:,iVec_) = e_Vector-DDot_(3,UK,1,e_Vector,1)*UK
          rNorm = DDot_(3,pol_Vector(:,iVec_),1,pol_Vector(:,iVec_),1)
          if (rNorm > 1.0e-12_wp) then
            pol_Vector(:,iVec_) = pol_Vector(:,iVec_)/sqrt(rNorm)
          else
            pol_Vector(:,iVec_) = Zero
          end if
        end if

        ! Generate the property integrals associated with this
        ! direction of the wave vector k.

        iOpt = 2
        call TMOMInt(wavevector,iOpt)

        !***************************************************************
        !                                                              *
        !      Recompute the needed properties for all the spin-free   *
        !      states.                                                 *
        !                                                              *
        !***************************************************************

        do IPRP=1,14
          PROP(:,:,IPRTMOM(IPRP)) = Zero
        end do

        do ISS=1,ISM

          ! Does this spin-free state contribute to any of the
          ! two spin states? Check the corresponding coefficients.

          i = iMask(ISS)
          do JSS=1,JSM
            j = jMask(JSS)

            ! COMBINED SYMMETRY OF STATES:
            JOB1 = JBNUM(I)
            JOB2 = JBNUM(J)
            LSYM1 = IRREP(JOB1)
            LSYM2 = IRREP(JOB2)
            ISY12 = MUL(LSYM1,LSYM2)
            ! THE SYMMETRY CHECK MASK:
            MASK = 2**(ISY12-1)
            ! FIRST SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS OF
            ! TDMSCR
            call mk_IOFF(IOFF,nIrrep,NBASF,ISY12)

            ! Pick up the transition density between the two
            ! states from disc. Generated in GTDMCTL.

            IDISK = iDisk_TDM(I,J,1)
            iEmpty = iDisk_TDM(I,J,2)
            iOpt = 2
            iGO = ibset(ibset(0,0),2)
            call dens2file(TDMZZ,Dummy,WDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,I,J)
            call MK_TWDM(nIrrep,TDMZZ,WDMZZ,nTDMZZ,SCR,nSCR,IOFF,NBASF,ISY12)

            ! Compute the transition property of the property
            ! integrals between the two states.

            do IPRP=1,14
              IPROP = IPRTMOM(IPRP)
              ITYPE = 0
              if (PTYPE(IPROP) == 'HERMSING') ITYPE = 1
              if (PTYPE(IPROP) == 'ANTISING') ITYPE = 2
              if (PTYPE(IPROP) == 'HERMTRIP') ITYPE = 3
              if (PTYPE(IPROP) == 'ANTITRIP') ITYPE = 4
              LABEL = PNAME(IPROP)
              call MK_PROP(PROP,IPROP,I,J,LABEL,ITYPE,IP,NIP,SCR,nSCR,MASK,ISY12,IOFF)
            end do

          end do ! J
        end do ! I

        ! DXR & DXI hold the component that does not change from k to -k
        ! DXRM & DXIM hold the component that changes sign

        DXR(:,:,:) = Zero
        DXI(:,:,:) = Zero
        DXRM(:,:,:) = Zero
        DXIM(:,:,:) = Zero
        do iCar=1,3

          ! (1) the spin-free part.
          !     Note that the integrals are not computed for the
          !     momentum operator but rather for nabla. That is
          !     as we assemble to transition momentum we have to
          !     remember to put in a factor of -i.

          ! The real part (symmetric and anti-symmetric) becomes imaginary

          call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(0+iCar),0,ISS_INDEX,iMask,ISM,jMask,JSM)
          DXI(:,:,iCar) = DXI(:,:,iCar)-TMP(:,:)
          call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(6+iCar),0,ISS_INDEX,iMask,ISM,jMask,JSM)
          DXI(:,:,iCar) = DXI(:,:,iCar)-TMP(:,:)

          ! The imaginary part (symmetric and anti-symmetric) becomes real

          call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(3+iCar),0,ISS_INDEX,iMask,ISM,jMask,JSM)
          DXRM(:,:,iCar) = DXRM(:,:,iCar)+TMP(:,:)
          call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(9+iCar),0,ISS_INDEX,iMask,ISM,jMask,JSM)
          DXRM(:,:,iCar) = DXRM(:,:,iCar)+TMP(:,:)

          ! (2) the spin-dependent part, magnetic

          !     iCar=1: 1/2(S(+)+S(-))
          !     iCar=2: 1/2i(S(+)-S(-))
          !     iCar=3: Sz

          !     Here we need to be very careful. The operator is
          !     similar to the L.S operator, the spin-orbit coupling.
          !     However, here we have that the operator is B.S. Thus
          !     we can do this in a similar fashion as the spin-orbit
          !     coupling, but with some important difference.

          !     For the spin-orbit coupling the L operator is divided
          !     into x-, y-, and z-components, that is the operator
          !     will be represented by three different integrals.
          !     For the integrals over B it is similar but the x-, y-,
          !     and z-components are constants outside of the integral.
          !     For example, the x-component is expressed as
          !     (k x e_l)_x <0|e^(i k.r)|n>. In this section we will
          !     handle the (k x e_l)_x part outside the loop over the
          !     Cartesian components.

          !     We have to note one further difference, the integrals
          !     are complex in this case.

          !     Let us now compute the contributions T(i), i=x,y,z

          !     In the equations we find that the y-component is imaginary,
          !     see the equations on page 234 in Per Ake's paper. Hence,
          !     the y-component is treated slightly differently.

          !     Actually, here we compute (<0|e^(i k.r)|n> x k), which
          !     will be later dotted with e_l.

          !     i*g/2*(s_y*k_z-s_z*k_y) -> T_x
          !     i*g/2*(s_z*k_x-s_x*k_z) -> T_y
          !     i*g/2*(s_x*k_y-s_y*k_x) -> T_z

          !     Note also that the "I" parts should change sign with kPhase,
          !     but so does wavevector, so the net result is the opposite.

          if (iCar == 1) then
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXIM(:,:,3) = DXIM(:,:,3)+wavevector(2)*cst*TMP(:,:)
            DXIM(:,:,2) = DXIM(:,:,2)-wavevector(3)*cst*TMP(:,:)
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXR(:,:,3) = DXR(:,:,3)-wavevector(2)*cst*TMP(:,:)
            DXR(:,:,2) = DXR(:,:,2)+wavevector(3)*cst*TMP(:,:)
          else if (iCar == 2) then
            ! For the y-component we have to interchange the real and
            ! the imaginary components. The real component gets a
            ! minus sign due to the product i*i=-1
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXI(:,:,1) = DXI(:,:,1)-wavevector(3)*cst*TMP(:,:)
            DXI(:,:,3) = DXI(:,:,3)+wavevector(1)*cst*TMP(:,:)
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXRM(:,:,1) = DXRM(:,:,1)-wavevector(3)*cst*TMP(:,:)
            DXRM(:,:,3) = DXRM(:,:,3)+wavevector(1)*cst*TMP(:,:)
          else if (iCar == 3) then
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(13),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXIM(:,:,2) = DXIM(:,:,2)+wavevector(1)*cst*TMP(:,:)
            DXIM(:,:,1) = DXIM(:,:,1)-wavevector(2)*cst*TMP(:,:)
            call SMMAT_MASKED(PROP,TMP,NSS,IPRTMOM(14),iCar,ISS_INDEX,iMask,ISM,jMask,JSM)
            DXR(:,:,2) = DXR(:,:,2)-wavevector(1)*cst*TMP(:,:)
            DXR(:,:,1) = DXR(:,:,1)+wavevector(2)*cst*TMP(:,:)
          end if
        end do

        ! Now we can compute the transition moments for k and -k

        do kp=1,2

          if (abs(kPhase(kp)) < Half) cycle
          TMR(:,:,:) = DXR(:,:,:)+kPhase(kp)*DXRM(:,:,:)
          TMI(:,:,:) = DXI(:,:,:)+kPhase(kp)*DXIM(:,:,:)
          do iCar=1,3
            call ZTRNSF_MASKED(NSS,VSOR,VSOI,TMR(:,:,iCar),TMI(:,:,iCar),IJSS,iSSMask,ISSM,jSSMask,JSSM)
          end do

          IJ_ = 0
          do ISO=istart_,iend_
            do JSO=jstart_,jend_
              IJ_ = IJ_+1
              EDIFF = ENSOR(JSO)-ENSOR(ISO)

              ! Store the vectors for this direction

              TM_R(:) = TMR(ISO,JSO,:)
              TM_I(:) = TMI(ISO,JSO,:)
#             ifdef _HDF5_
              ! Get proper triangular index
              IJSO_ = nTri_Elem(JSO-2)+ISO
              iQuad_ = 2*(iQuad-1)+kp
              Storage(ip_w,iQuad_,IJSO_,iVec) = Weight
              Storage(ip_kvector:ip_kvector+2,iQuad_,IJSO_,iVec) = kPhase(kp)*Wavevector(:)
              Storage(ip_TMR:ip_TMR+2,iQuad_,IJSO_,iVec) = TM_R(:)
              Storage(ip_TMI:ip_TMI+2,iQuad_,IJSO_,iVec) = TM_I(:)
#             endif

              ! Project out the k direction from the real and imaginary components

              TM_R(:) = TM_R(:)-DDot_(3,TM_R,1,UK,1)*UK(:)
              TM_I(:) = TM_I(:)-DDot_(3,TM_I,1,UK,1)*UK(:)

              ! Implicitly integrate over all directions of the
              ! polarization vector to get the average value.

              TM1 = DDot_(3,TM_R,1,TM_R,1)
              TM2 = DDot_(3,TM_I,1,TM_I,1)
              TM_2 = Half*(TM1+TM2)

              ! Compute maximum and minimum oscillator strengths
              ! and the corresponding polarization vectors

              if (Do_SK) then
                TM3 = DDot_(3,TM_R,1,TM_I,1)
                Rng = sqrt((TM1-TM2)**2+Four*TM3**2)
                Aux(1,ij_) = TM_2+Half*Rng
                Aux(5,ij_) = TM_2-Half*Rng
                ! The direction for the maximum
                Ang = Half*atan2(Two*TM3,TM1-TM2)
                Aux(2:4,ij_) = Aux(2:4,ij_)+cos(Ang)*TM_R(:)+sin(Ang)*TM_I(:)
                ! Normalize and compute the direction for the minimum
                ! as a cross product with k
                rNorm = DDot_(3,Aux(2,ij_),1,Aux(2,ij_),1)
                if (rNorm > 1.0e-12_wp) then
                  Aux(2:4,ij_) = Aux(2:4,ij_)/sqrt(rNorm)
                  Aux(6,ij_) = Aux(3,ij_)*UK(3)-Aux(4,ij_)*UK(2)
                  Aux(7,ij_) = Aux(4,ij_)*UK(1)-Aux(2,ij_)*UK(3)
                  Aux(8,ij_) = Aux(2,ij_)*UK(2)-Aux(3,ij_)*UK(1)
                  rNorm = DDot_(3,Aux(6,ij_),1,Aux(6,ij_),1)
                  Aux(6:8,ij_) = Aux(6:8,ij_)/sqrt(rNorm)
                else
                  Aux(2:4,ij_) = Zero
                  Aux(6:8,ij_) = Zero
                end if
              end if

              ! Oscillator strength for a specific polarization vector

              if (Do_Pol) then
                TM1 = DDot_(3,TM_R,1,pol_Vector(1,iVec_),1)
                TM2 = DDot_(3,TM_I,1,pol_Vector(1,iVec_),1)
                TM_2 = TM1*TM1+TM2*TM2
              end if

              ! Compute the oscillator strength

              F_Temp = Two*TM_2/EDIFF
              if (Do_SK) then
                Aux(1,ij_) = Two*Aux(1,ij_)/EDIFF
                Aux(5,ij_) = Two*Aux(5,ij_)/EDIFF
              end if

              ! Compute the rotatory strength, note that it depends on kPhase

              TM_C(1) = TM_R(2)*TM_I(3)-TM_R(3)*TM_I(2)
              TM_C(2) = TM_R(3)*TM_I(1)-TM_R(1)*TM_I(3)
              TM_C(3) = TM_R(1)*TM_I(2)-TM_R(2)*TM_I(1)
              TM_2 = Two*kPhase(kp)*DDot_(3,TM_C,1,UK,1)
              R_Temp = 0.75_wp*c_in_au/EDIFF**2*TM_2
              R_Temp = R_Temp*AU2REDR

              ! Save the raw oscillator and rotatory strengths in a given direction

              NQUAD_ = 2*NQUAD
              LRAW_ = LRAW+6*NQUAD_*(ij_-1)
              IQUAD_ = 2*(IQUAD-1)+(kp-1)
              RAW(LRAW_+IQUAD_+0*NQUAD_) = F_Temp
              RAW(LRAW_+IQUAD_+1*NQUAD_) = R_Temp

              ! Save the direction and weight too

              RAW(LRAW_+IQUAD_+2*NQUAD_) = UK(1)*kPhase(kp)
              RAW(LRAW_+IQUAD_+3*NQUAD_) = UK(2)*kPhase(kp)
              RAW(LRAW_+IQUAD_+4*NQUAD_) = UK(3)*kPhase(kp)
              RAW(LRAW_+IQUAD_+5*NQUAD_) = Weight

              ! Do not accumulate if not doing an isotropic integration

              if (Do_SK .and. (kp > 1)) cycle

              ! Accumulate to the isotropic oscillator strength

              OscStr(1,IJ_) = OscStr(1,IJ_)+Weight*F_Temp

              ! Accumulate to the isotropic rotatory strength

              OscStr(2,IJ_) = OscStr(2,IJ_)+Weight*R_Temp
            end do
          end do

        end do ! kp

      end do ! iQuad

      IJ_ = 0
      do ISO=istart_,iend_
        do JSO=jstart_,jend_
          IJ_ = IJ_+1
          EDIFF = ENSOR(JSO)-ENSOR(ISO)
          F = OscStr(1,IJ_)
          R = OscStr(2,IJ_)

          call Add_Info('ITMS(SO)',[F],1,6)
          call Add_Info('ROTS(SO)',[R],1,4)

          if (Do_Pol) then
            F_CHECK = abs(Aux(1,ij_))
            R_CHECK = Zero ! dummy assign
          else
            F_CHECK = abs(F)
            R_CHECK = abs(R)
          end if
          if ((F_CHECK < OSTHR) .and. (R_CHECK < RSTHR)) cycle
          A = (AFACTOR*EDIFF**2)*F

          if (iPrint == 0) then
            write(u6,*)
            if (Do_SK) then
              call CollapseOutput(1,'Transition moment strengths (SO states):')
              write(u6,'(3X,A)') '----------------------------------------'
              if (Do_Pol) then
                iVec_ = (iVec-1)*nQuad+1
                write(u6,'(4x,a,3F10.6)') 'Direction of the polarization: ',(pol_vector(k,iVec),k=1,3)
              else
                write(u6,'(4x,a)') 'The oscillator strength is integrated over all directions of the polarization vector'
              end if
              write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',Rquad(1:3,1)
            else
              call CollapseOutput(1,'Isotropic transition moment strengths (SO states):')
              write(u6,'(3X,A)') '--------------------------------------------------'
            end if
            if (OSTHR > Zero) write(u6,45) 'For osc. strength at least',OSTHR,'and red. rot. strength  at least',RSTHR
            write(u6,*)
            if (.not. Do_SK) then
              write(u6,'(4x,a,I4,a)') 'Integrated over ',2*nQuad,' directions of the wave vector'
              write(u6,'(4x,a)') 'The oscillator strength is integrated over all directions of the polarization vector'
              write(u6,*)
            end if
            write(u6,31) 'From','To','Osc. strength','Red. rot. str.','Total A (sec-1)'
            write(u6,32)
            iPrint = 1
          end if

          ! Regular print

          ! Don't print osc. str. if below threshold
          if (F_CHECK < OSTHR) then
            write(u6,46) ISO,JSO,'below threshold',R,A
            ! Don't print rot. str. if below threshold
          else if (R_CHECK < RSTHR) then
            write(u6,47) ISO,JSO,F,'below threshold',A
          else
            write(u6,33) ISO,JSO,F,R,A
          end if

          if (Do_SK) then
            write(u6,50) 'maximum',Aux(1,ij_),'for polarization direction:',Aux(2,ij_),Aux(3,ij_),Aux(4,ij_)
            write(u6,50) 'minimum',Aux(5,ij_),'for polarization direction:',Aux(6,ij_),Aux(7,ij_),Aux(8,ij_)
          end if

          ! Printing raw (unweighted) and direction for every transition

          if (PRRAW) then
            write(u6,*)
            write(u6,*)
            write(u6,34) 'From','To','Raw osc. str.','Red. rot. str.','kx','ky','kz'
            write(u6,35)
            NQUAD_ = 2*NQUAD
            LRAW_ = LRAW+6*NQUAD_*(ij_-1)
            do IQUAD=1,NQUAD
              do kp=1,2
                if (abs(kPhase(kp)) < Half) cycle
                IQUAD_ = 2*(IQUAD-1)+(kp-1)
                write(u6,33) ISO,JSO,RAW(LRAW_+IQUAD_+0*NQUAD_),RAW(LRAW_+IQUAD_+1*NQUAD_),RAW(LRAW_+IQUAD_+2*NQUAD_), &
                             RAW(LRAW_+IQUAD_+3*NQUAD_),RAW(LRAW_+IQUAD_+4*NQUAD_)
              end do
            end do
            write(u6,35)
            write(u6,*)
          end if

          ! Printing weighted and direction for every transition

          if (PRWEIGHT) then
            write(u6,*)
            write(u6,*)
            write(u6,34) 'From','To','Weig. osc. str.','Red. rot. str.','kx','ky','kz'
            write(u6,35)
            NQUAD_ = 2*NQUAD
            LRAW_ = LRAW+6*NQUAD_*(ij_-1)
            do IQUAD=1,NQUAD
              do kp=1,2
                if (abs(kPhase(kp)) < Half) cycle
                IQUAD_ = 2*(IQUAD-1)+(kp-1)
                Weight = RAW(LRAW_+IQUAD_+5*NQUAD_)
                write(u6,33) ISO,JSO,RAW(LRAW_+IQUAD_+0*NQUAD_)*Weight,RAW(LRAW_+IQUAD_+1*NQUAD_)*Weight, &
                             RAW(LRAW_+IQUAD_+2*NQUAD_),RAW(LRAW_+IQUAD_+3*NQUAD_),RAW(LRAW_+IQUAD_+4*NQUAD_)
              end do
            end do
            write(u6,35)
            write(u6,*)
          end if

        end do
      end do

    end do
  end do
  call mma_deallocate(iMask)
  call mma_deallocate(jMask)
  call mma_deallocate(iSSMask)
  call mma_deallocate(jSSMask)
  call mma_deallocate(ISS_INDEX)

  if (iPrint == 1) then
    write(u6,32)
    if (Do_SK) then
      call CollapseOutput(0,'Transition moment strengths (SO states):')
    else
      call CollapseOutput(0,'Isotropic transition moment strengths (SO states):')
    end if
  end if

end do ! iVec
#ifdef _TIME_TMOM_
call CWTime(TCpu2,TWall2)
write(u6,*) 'Time for TMOM : ',TCpu2-TCpu1,TWall2-TWall1
#endif

#ifdef _HDF5_
call mh5_put_dset(wfn_sos_tm,pack(Storage,.true.))
call mma_deallocate(Storage)
#endif

! Do some cleanup

call mma_deallocate(RAW)
call mma_deallocate(IP)
call mma_deallocate(DXR)
call mma_deallocate(DXI)
call mma_deallocate(DXRM)
call mma_deallocate(DXIM)
call mma_deallocate(TMP)
if (TMOgroup) then
  call mma_DeAllocate(TMOgrp1)
  call mma_DeAllocate(TMOgrp2)
end if
if (Do_Pol) call mma_deallocate(pol_Vector)
call mma_deallocate(TMR)
call mma_deallocate(TMI)
call mma_deallocate(OSCSTR)
call mma_deallocate(Aux)

call mma_deallocate(SCR)
call mma_deAllocate(TDMZZ)
call mma_deAllocate(WDMZZ)
call mma_deAllocate(Rquad)
call ClsSew()

call mma_deallocate(VSOR)
call mma_deallocate(VSOI)

30 format(5X,A,1X,ES15.8)
31 format(5X,2(1X,A4),4X,3(1X,A15))
32 format(5X,63('-'))
33 format(5X,2(1X,I4),5X,5(1X,ES15.8))
34 format(5X,2(1X,A4),5X,5(1X,A15))
35 format(5X,95('-'))
45 format(4X,2(A,1X,ES15.8,1X))
46 format(5X,2(1X,I4),5X,(1X,A15),2(1X,ES15.8))
47 format(5X,2(1X,I4),5X,(1X,ES15.8),(1X,A15),(1X,ES15.8))
49 format(5X,A,1X,ES15.8,1X,A)
50 format(10X,A7,3X,1(1X,ES15.8),5X,A27,3(1X,F7.4))

!***********************************************************************
!                                                                      *
!     End of section for transition moments                            *
!                                                                      *
!***********************************************************************

#ifdef _TIME_TMOM_
call CWTime(TCpu2,TWall2)
write(u6,*) 'Time for TMOM(SO) : ',TCpu2-TCpu1,TWall2-TWall1
#endif

end subroutine PRPROP_TM_Exact
