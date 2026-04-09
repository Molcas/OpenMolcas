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
! Copyright (C) 2021, Rulin Feng                                       *
!***********************************************************************
!***************************************************
!         Get TDM in AO basis from SO states
!****************************************************
! This routine is modified from sonatorbm to give the
! full transition density matrix (TDM) in atomic orbital
! basis, which is not a nbtri sized matrix as in sonatorbm
! but a full nbst**2 sized matrix.

subroutine MAKETDMAO(CHARTYPE,USOR,USOI,ASS,BSS,NSS,iOpt,ROTMAT,DENSOUT)

use rassi_aux, only: idisk_TDM
use rassi_global_arrays, only: JBNUM
use rassi_data, only: NBASF, NBST, NTDMZZ
use Cntrl, only: IRREP, LSYM1, LSYM2, LuTDM, MLTPLT, NSTATE
use Symmetry_Info, only: MUL, nIrrep
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
character(len=8) :: CHARTYPE
integer(kind=iwp) :: ASS, BSS, NSS, iOpt
real(kind=wp) :: USOR(NSS,NSS), USOI(NSS,NSS), ROTMAT(3,3), DENSOUT(6,nbst**2)
integer(kind=iwp) :: I, IDIR, IDISK, IEMPTY, IGO, IJ, IOFF(8), ISF, ISS, ISY1, ISY12, ISY12_MA, ISY2, ITD, ITYPE, J, JOB, JOB1, &
                     JOB2, KSF, KSS, LSF, LSS, MPLET, MPLETK, MPLETL, MSPROJ, MSPROJK, MSPROJL, NB1, NB1_F, NB1_I, NB2, NB2_F, &
                     NB2_I, NBSTS
real(kind=wp) :: CG0, CGM, CGP, CGX, CGY, FACT, S1, S2, SM1, SM2, TDM, UIL, UIR, URL, URR
integer, allocatable :: MAPST(:), MAPSP(:), MAPMS(:)
real(kind=wp), allocatable :: SDMXYZR(:,:), SDMXYZI(:,:), SDMXYZR2(:,:), SDMXYZI2(:,:), TMPR(:), TMPI(:), SCR(:), TDMZZ(:)
real(kind=wp), external :: DCLEBS

nbsts = nbst**2

! VV: dummy initialization
CGY = -1
CGX = -1
CG0 = -1
! Initialize
DENSOUT(:,:) = Zero
! Get the proper type of the property
ITYPE = 0
if (CHARTYPE == 'HERMSING') ITYPE = 1
if (CHARTYPE == 'ANTISING') ITYPE = 2
if (CHARTYPE == 'HERMTRIP') ITYPE = 3
if (CHARTYPE == 'ANTITRIP') ITYPE = 4
if (ITYPE == 0) then
  write(u6,*) 'RASSI/SONATORB internal error.'
  write(u6,*) 'Erroneous property type:',CHARTYPE
  call ABEND()
end if

! The following creates an array that is used to
! map a specific spin state to the corresponding
! spin-free state and to its spin
! (see prprop and others)

call mma_allocate(MAPST,NSS,Label='MAPST')
call mma_allocate(MAPSP,NSS,Label='MAPSP')
call mma_allocate(MAPMS,NSS,Label='MAPMS')

ISS = 0
do ISF=1,NSTATE
  JOB = JBNUM(ISF)
  MPLET = MLTPLT(JOB)

  do MSPROJ=-MPLET+1,MPLET-1,2
    ISS = ISS+1
    MAPST(ISS) = ISF
    MAPSP(ISS) = MPLET
    MAPMS(ISS) = MSPROJ
  end do
end do

! Allocate some arrays
! SDMXYZR, etc    DM/TDM for this iteration
! SDMXYZR2, etc   Accumulated DM/TDM
! TMPR,I          Temporary array for U*AU multiplication
! TDMZZ           DM/TDM as read from file
! SCR             Scratch for expansion of TDMZZ
call mma_allocate(SDMXYZR,nbsts,3,Label='SDMXYZR')
call mma_allocate(SDMXYZI,nbsts,3,Label='SDMXYZI')
SDMXYZR(:,:) = Zero
SDMXYZI(:,:) = Zero

call mma_allocate(SDMXYZR2,nbsts,3,Label='SDMXYZR2')
call mma_allocate(SDMXYZI2,nbsts,3,Label='SDMXYZI2')
SDMXYZR2(:,:) = Zero
SDMXYZI2(:,:) = Zero

call mma_allocate(TMPR,nbsts,Label='TMPR')
call mma_allocate(TMPI,nbsts,Label='TMPI')
TMPR(:) = Zero
TMPI(:) = Zero

call mma_allocate(SCR,nbsts,Label='SCR')
! zeroed inside the loop

call mma_allocate(TDMZZ,NTDMZZ,Label='TDMZZ')
TDMZZ(:) = Zero

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MAIN LOOP OVER KSF/LSF
! WRITTEN AS IN PRPROP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CORRESPONDING SPIN-FREE STATES OF THE
! REQUESTED SPIN STATES
!ASF = MAPST(ASS)
!BSF = MAPST(BSS)

do KSS=1,NSS
  KSF = MAPST(KSS)
  MPLETK = MAPSP(KSS)
  MSPROJK = MAPMS(KSS)

  do LSS=1,NSS
    LSF = MAPST(LSS)
    MPLETL = MAPSP(LSS)
    MSPROJL = MAPMS(LSS)

    JOB1 = JBNUM(KSF)
    JOB2 = JBNUM(LSF)
    LSYM1 = IRREP(JOB1)
    LSYM2 = IRREP(JOB2)
    ISY12 = MUL(LSYM1,LSYM2)

    ! SET UP AN OFFSET TABLE FOR SYMMETRY BLOCKS
    call mk_IOFF(IOFF,nIrrep,NBASF,ISY12)

    ! These are going to be zero, so head them off at the pass
    if ((ITYPE <= 2) .and. ((MPLETK /= MPLETL) .or. (MSPROJK /= MSPROJL))) cycle

    !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
    ! Transition density matrices, TDMZZ, in AO basis.
    ! WDMZZ similar, but WE-reduced 'triplet' densities.
    ! TDMZZ will store either, depending on the type
    TDMZZ(:) = Zero

    !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! IDTDM: TOC array for transition 1-matrices
    ! TDMZZ is stored on disk from i = 1, NSTATE j=1, i
    ! so swap if needed
    iEmpty = iDisk_TDM(KSF,LSF,2)
    IDISK = iDisk_TDM(KSF,LSF,1)
    iOpt = 2
    if (ITYPE >= 3) then
      iGo = 4
      call dens2file(TDMZZ,TDMZZ,TDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
      ! NOTE-the TD matrix as read in has an incorrect sign
      TDMZZ(:) = -TDMZZ(:)
    else
      iGo = 1
      call dens2file(TDMZZ,TDMZZ,TDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
    end if

    ! Anti-hermitian properties need a little fixing
    if (((ITYPE == 2) .or. (ITYPE == 4)) .and. (KSF <= LSF)) TDMZZ(:) = -TDMZZ(:)

    ! CALCULATE THE SYMMETRIC AND ANTISYMMETRIC FOLDED TRANS D MATRICES
    ! AND SIMILAR WE-REDUCED SPIN DENSITY MATRICES
    SCR(:) = Zero
    !ccccccccccccc
    ! nbsq and NTDMZZ are the same, but they are the
    ! sizes of the symmetry-adapted matrices
    ! what we need is a NBST**2-sized matrix without symmetry
    ! Thus we expand the final TDM in C1 symmetry
    ! Leave the zero matrix elements as they are
    !cccccccccccccc

    ! Expand into C1
    ITD = 0
    NB1_i = 0
    NB1_f = 0
    do ISY1=1,nIrrep
      NB2_i = 0
      NB2_f = 0
      NB1 = NBASF(ISY1)
      NB1_f = NB1_i+NB1
      do ISY2=1,nIrrep
        ISY12_ma = MUL(ISY1,ISY2)
        NB2 = NBASF(ISY2)
        NB2_f = NB2_i+NB2
        if (ISY12_ma == ISY12) then
          do J=NB2_i+1,NB2_f
            do I=NB1_i+1,NB1_f
              ITD = ITD+1
              TDM = TDMZZ(ITD)
              IJ = I+NBST*(J-1)
              SCR(IJ) = SCR(IJ)+TDM
            end do
          end do
        end if
        NB2_i = NB2_i+NB2
      end do
      NB1_i = NB1_i+NB1
    end do

    ! ie, see how AMFI is processed in soeig
    SDMXYZR(:,:) = Zero
    SDMXYZI(:,:) = Zero

    if (ITYPE >= 3) then
      S1 = Half*real(MPLETK-1,kind=wp)
      SM1 = Half*real(MSPROJK,kind=wp)
      S2 = Half*real(MPLETL-1,kind=wp)
      SM2 = Half*real(MSPROJL,kind=wp)
      FACT = One/sqrt(real(MPLETK,kind=wp))
      if (MPLETK == MPLETL-2) FACT = -FACT

      CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
      CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
      CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
      CGX = sqrt(Half)*(CGM-CGP)
      CGY = sqrt(Half)*(CGM+CGP)
    end if

    if (((ITYPE == 1) .or. (ITYPE == 2)) .and. (MPLETK == MPLETL) .and. (MSPROJK == MSPROJL)) then
      SDMXYZR(:,3) = SDMXYZR(:,3)+SCR(:)
    else if ((ITYPE == 3) .or. (ITYPE == 4)) then
      if (iOpt == 1) then
        SDMXYZR(:,1) = SDMXYZR(:,1)+(CGX*ROTMAT(1,1)+CGY*ROTMAT(2,1)+CG0*ROTMAT(3,1))*SCR(:)

        SDMXYZR(:,2) = SDMXYZR(:,2)+(CGX*ROTMAT(1,2)+CGY*ROTMAT(2,2)+CG0*ROTMAT(3,2))*SCR(:)

        SDMXYZR(:,3) = SDMXYZR(:,3)+(CGX*ROTMAT(1,3)+CGY*ROTMAT(2,3)+CG0*ROTMAT(3,3))*SCR(:)
      else
        SDMXYZR(:,1) = SDMXYZR(:,1)+CGX*SCR(:)
        SDMXYZI(:,2) = SDMXYZI(:,2)+CGY*SCR(:)
        SDMXYZR(:,3) = SDMXYZR(:,3)+CG0*SCR(:)
      end if
    end if

    !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! SPINORBIT
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccc
    ! Sign of the left-hand imaginary part is handled
    ! when doing DAXPY
    URR = USOR(LSS,BSS)
    UIR = USOI(LSS,BSS)
    URL = USOR(KSS,ASS)
    UIL = USOI(KSS,ASS)

    do IDIR=1,3
      TMPR(:) = Zero
      TMPI(:) = Zero

      ! right side
      TMPR(:) = TMPR(:)+URR*SDMXYZR(:,IDIR)-UIR*SDMXYZI(:,IDIR)
      TMPI(:) = TMPI(:)+UIR*SDMXYZR(:,IDIR)+URR*SDMXYZI(:,IDIR)

      ! left side
      SDMXYZR2(:,IDIR) = SDMXYZR2(:,IDIR)+URL*TMPR(:)+UIL*TMPI(:)
      SDMXYZI2(:,IDIR) = SDMXYZI2(:,IDIR)+URL*TMPI(:)-UIL*TMPR(:)
    end do
    !ccccccccccccccccccccc
    ! END SPINORBIT STUFF
    !ccccccccccccccccccccc

    ! END MAIN LOOP OVER STATES (KSS,LSS)
  end do
end do

! Store this density to DENSOUT
if ((ITYPE == 3) .or. (ITYPE == 4)) then
  DENSOUT(1,:) = SDMXYZR2(:,1)
  DENSOUT(2,:) = SDMXYZR2(:,2)
  DENSOUT(3,:) = SDMXYZR2(:,3)
  DENSOUT(4,:) = SDMXYZI2(:,1)
  DENSOUT(5,:) = SDMXYZI2(:,2)
  DENSOUT(6,:) = SDMXYZI2(:,3)
else
  DENSOUT(1,:) = SDMXYZR2(:,3)
  DENSOUT(2,:) = SDMXYZR2(:,3)
  DENSOUT(3,:) = SDMXYZR2(:,3)
  DENSOUT(4,:) = SDMXYZI2(:,3)
  DENSOUT(5,:) = SDMXYZI2(:,3)
  DENSOUT(6,:) = SDMXYZI2(:,3)
end if

! Free memory
call mma_deallocate(SCR)
call mma_deallocate(TDMZZ)

call mma_deallocate(TMPI)
call mma_deallocate(TMPR)

call mma_deallocate(SDMXYZR)
call mma_deallocate(SDMXYZI)

call mma_deallocate(SDMXYZR2)
call mma_deallocate(SDMXYZI2)

call mma_deallocate(MAPST)
call mma_deallocate(MAPSP)
call mma_deallocate(MAPMS)

end subroutine MAKETDMAO
