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
use stdalloc, only: mma_allocate, mma_deallocate
use Cntrl, only: NSTATE, LSYM1, LSYM2, IRREP, MLTPLT
use cntrl, only: LuTDM
use Symmetry_Info, only: nSym => nIrrep, MUL
use rassi_data, only: NBST, NBASF, NTDMZZ
use Constants, only: Zero, One, Half
use Definitions, only: wp, u6

implicit none
character(len=8) CHARTYPE!,LABEL
integer NSS
real*8 USOR(NSS,NSS), USOI(NSS,NSS)
integer ASS, BSS, iOpt
real*8 ROTMAT(3,3)
real*8 DENSOUT(6,nbst**2)
integer IOFF(8)
integer, allocatable :: MAPST(:), MAPSP(:), MAPMS(:)
real*8, allocatable, target :: SDMXR(:), SDMXI(:), SDMYR(:), SDMYI(:), SDMZR(:), SDMZI(:)
real*8, allocatable, target :: SDMXR2(:), SDMXI2(:), SDMYR2(:), SDMYI2(:), SDMZR2(:), SDMZI2(:)
type A2_array
  real*8, pointer :: A2(:)
end type A2_Array
type(A2_array) :: pZMR(3), pZMI(3)
type(A2_array) :: pZMR2(3), pZMI2(3)
real*8, allocatable :: TMPR(:), TMPI(:), SCR(:), TDMZZ(:)
real*8 CGY, CGX, CG0, TDM, S1, S2, SM1, SM2, FACT, CGM, CGP, URR, UIR, URL, UIL, DCLEBS
integer NBSTS, ITYPE, ISS, ISF, JOB, MPLET, MSPROJ, KSS, KSF, MPLETK, MSPROJK, LSS, LSF, MPLETL, MSPROJL, JOB1, JOB2, ISY12, &
        IEMPTY, IDISK, IGO, ITD, NB1_I, NB1_F, ISY1, NB2_I, NB2_F, NB1, ISY2, ISY12_MA, NB2, J, I, IJ, IDIR

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
! SDMXR, etc      DM/TDM for this iteration
! SDMXR2, etc     Accumulated DM/TDM
! TMPR,I          Temporary array for U*AU multiplication
! TDMZZ           DM/TDM as read from file
! SCR             Scratch for expansion of TDMZZ
call mma_allocate(SDMXR,nbsts,Label='SDMXR')
call mma_allocate(SDMXI,nbsts,Label='SDMXI')
call mma_allocate(SDMYR,nbsts,Label='SDMYR')
call mma_allocate(SDMYI,nbsts,Label='SDMYI')
call mma_allocate(SDMZR,nbsts,Label='SDMZR')
call mma_allocate(SDMZI,nbsts,Label='SDMZI')
SDMXR(:) = Zero
SDMXI(:) = Zero
SDMYR(:) = Zero
SDMYI(:) = Zero
SDMZR(:) = Zero
SDMZI(:) = Zero
pZMR(1)%A2 => SDMXR(:)
pZMR(2)%A2 => SDMYR(:)
pZMR(3)%A2 => SDMZR(:)
pZMI(1)%A2 => SDMXI(:)
pZMI(2)%A2 => SDMYI(:)
pZMI(3)%A2 => SDMZI(:)

call mma_allocate(SDMXR2,nbsts,Label='SDMXR2')
call mma_allocate(SDMXI2,nbsts,Label='SDMXI2')
call mma_allocate(SDMYR2,nbsts,Label='SDMYR2')
call mma_allocate(SDMYI2,nbsts,Label='SDMYI2')
call mma_allocate(SDMZR2,nbsts,Label='SDMZR2')
call mma_allocate(SDMZI2,nbsts,Label='SDMZI2')
SDMXR2(:) = Zero
SDMXI2(:) = Zero
SDMYR2(:) = Zero
SDMYI2(:) = Zero
SDMZR2(:) = Zero
SDMZI2(:) = Zero
pZMR2(1)%A2 => SDMXR2(:)
pZMR2(2)%A2 => SDMYR2(:)
pZMR2(3)%A2 => SDMZR2(:)
pZMI2(1)%A2 => SDMXI2(:)
pZMI2(2)%A2 => SDMYI2(:)
pZMI2(3)%A2 => SDMZI2(:)

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
    call mk_IOFF(IOFF,nSYM,NBASF,ISY12)

    ! These are going to be zero, so head them off at the pass
    if ((ITYPE <= 2) .and. ((MPLETK /= MPLETL) .or. (MSPROJK /= MSPROJL))) goto 2200

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
      call DSCAL_(NTDMZZ,-One,TDMZZ,1)
    else
      iGo = 1
      call dens2file(TDMZZ,TDMZZ,TDMZZ,nTDMZZ,LUTDM,IDISK,iEmpty,iOpt,iGo,KSF,LSF)
    end if

    ! Anti-hermitian properties need a little fixing
    if (((ITYPE == 2) .or. (ITYPE == 4)) .and. (KSF <= LSF)) call DSCAL_(NTDMZZ,-One,TDMZZ,1)

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
    do ISY1=1,NSYM
      NB2_i = 0
      NB2_f = 0
      NB1 = NBASF(ISY1)
      NB1_f = NB1_i+NB1
      do ISY2=1,NSYM
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
    SDMXR(:) = Zero
    SDMXI(:) = Zero
    SDMYR(:) = Zero
    SDMYI(:) = Zero
    SDMZR(:) = Zero
    SDMZI(:) = Zero

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
      call DAXPY_(nbsts,One,SCR,1,SDMZR,1)
    else if ((ITYPE == 3) .or. (ITYPE == 4)) then
      if (iOpt == 1) then
        call DAXPY_(nbsts,CGX*ROTMAT(1,1),SCR,1,SDMXR,1)
        call DAXPY_(nbsts,CGY*ROTMAT(2,1),SCR,1,SDMXI,1)
        call DAXPY_(nbsts,CG0*ROTMAT(3,1),SCR,1,SDMXR,1)

        call DAXPY_(nbsts,CGX*ROTMAT(1,2),SCR,1,SDMYR,1)
        call DAXPY_(nbsts,CGY*ROTMAT(2,2),SCR,1,SDMYI,1)
        call DAXPY_(nbsts,CG0*ROTMAT(3,2),SCR,1,SDMYR,1)

        call DAXPY_(nbsts,CGX*ROTMAT(1,3),SCR,1,SDMZR,1)
        call DAXPY_(nbsts,CGY*ROTMAT(2,3),SCR,1,SDMZI,1)
        call DAXPY_(nbsts,CG0*ROTMAT(3,3),SCR,1,SDMZR,1)
      else
        call DAXPY_(nbsts,CGX,SCR,1,SDMXR,1)
        call DAXPY_(nbsts,CGY,SCR,1,SDMYI,1)
        call DAXPY_(nbsts,CG0,SCR,1,SDMZR,1)
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
      call DAXPY_(nbsts,URR,pZMR(IDIR)%A2,1,TMPR,1)
      call DAXPY_(nbsts,-UIR,pZMI(IDIR)%A2,1,TMPR,1)
      call DAXPY_(nbsts,UIR,pZMR(IDIR)%A2,1,TMPI,1)
      call DAXPY_(nbsts,URR,pZMI(IDIR)%A2,1,TMPI,1)

      ! left side
      call DAXPY_(nbsts,URL,TMPR,1,pZMR2(IDIR)%A2,1)
      call DAXPY_(nbsts,UIL,TMPI,1,pZMR2(IDIR)%A2,1)
      call DAXPY_(nbsts,URL,TMPI,1,pZMI2(IDIR)%A2,1)
      call DAXPY_(nbsts,-UIL,TMPR,1,pZMI2(IDIR)%A2,1)
    end do
    !ccccccccccccccccccccc
    ! END SPINORBIT STUFF
    !ccccccccccccccccccccc

    ! END MAIN LOOP OVER STATES (KSS,LSS)
2200 continue

  end do
end do

! Store this density to DENSOUT
if ((ITYPE == 3) .or. (ITYPE == 4)) then
  do I=1,nbsts
    DENSOUT(1,I) = SDMXR2(I)
    DENSOUT(2,I) = SDMYR2(I)
    DENSOUT(3,I) = SDMZR2(I)
    DENSOUT(4,I) = SDMXI2(I)
    DENSOUT(5,I) = SDMYI2(I)
    DENSOUT(6,I) = SDMZI2(I)
  end do
else
  do I=1,nbsts
    DENSOUT(1,I) = SDMZR2(I)
    DENSOUT(2,I) = SDMZR2(I)
    DENSOUT(3,I) = SDMZR2(I)
    DENSOUT(4,I) = SDMZI2(I)
    DENSOUT(5,I) = SDMZI2(I)
    DENSOUT(6,I) = SDMZI2(I)
  end do
end if

! Free memory
call mma_deallocate(SCR)
call mma_deallocate(TDMZZ)

nullify(pZMI(3)%A2,pZMI(2)%A2,pZMI(1)%A2,pZMR(3)%A2,pZMR(2)%A2,pZMR(1)%A2)
call mma_deallocate(SDMZI)
call mma_deallocate(SDMZR)
call mma_deallocate(SDMYI)
call mma_deallocate(SDMYR)
call mma_deallocate(SDMXI)
call mma_deallocate(SDMXR)

call mma_deallocate(TMPI)
call mma_deallocate(TMPR)

nullify(pZMI2(3)%A2,pZMI2(2)%A2,pZMI2(1)%A2,pZMR2(3)%A2,pZMR2(2)%A2,pZMR2(1)%A2)
call mma_deallocate(SDMZI2)
call mma_deallocate(SDMZR2)
call mma_deallocate(SDMYI2)
call mma_deallocate(SDMYR2)
call mma_deallocate(SDMXI2)
call mma_deallocate(SDMXR2)

call mma_deallocate(MAPST)
call mma_deallocate(MAPSP)
call mma_deallocate(MAPMS)

end subroutine MAKETDMAO
