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

subroutine SONATORBM(CHARTYPE,USOR,USOI,ASS,BSS,NSS,iOpt,ROTMAT,DENSOUT)

use rassi_aux, only: idisk_TDM
use rassi_global_arrays, only: JBNUM
use Cntrl, only: IRREP, LSYM1, LSYM2, LuTDM, MLTPLT, NSTATE
use Symmetry_Info, only: nSym => nIrrep, MUL
use rassi_data, only: NBASF, NBTRI, NTDMZZ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
character(len=*) :: CHARTYPE
integer(kind=iwp) :: ASS, BSS, NSS, iOpt
real(kind=wp) :: USOR(NSS,NSS), USOI(NSS,NSS), ROTMAT(3,3), DENSOUT(6,NBTRI)
integer(kind=iwp) :: I, IDIR, IDISK, IEMPTY, IGO, IJ, IOF, IOFF(8), ISF, ISS, ISY, ISY1, ISY12, ISY2, ITD, ITYPE, J, JOB, JOB1, &
                     JOB2, KSF, KSS, LSF, LSS, MPLET, MPLETK, MPLETL, MSPROJ, MSPROJK, MSPROJL, NB, NB1, NB2
real(kind=wp) :: CG0, CGM, CGP, CGX, CGY, DCLEBS, FACT, S1, S2, SM1, SM2, TDM, UIL, UIR, URL, URR
integer(kind=iwp), allocatable :: MAPST(:), MAPSP(:), MAPMS(:)
real(kind=wp), allocatable :: SCR(:), SDMXYZI(:,:), SDMXYZI2(:,:), SDMXYZR(:,:), SDMXYZR2(:,:), TDMZZ(:), TMPI(:), TMPR(:)

! VV: dummy initialization
CGY = -1
CGX = -1
CG0 = -1
! Get the proper type of the property
ITYPE = 0
if (CHARTYPE == 'HERMSING') ITYPE = 1
if (CHARTYPE == 'ANTISING') ITYPE = 2
if (CHARTYPE == 'HERMTRIP') ITYPE = 3
if (CHARTYPE == 'ANTITRIP') ITYPE = 4
if (ITYPE == 0) then
  write(u6,*) 'RASSI/SONATORB internal error.'
  write(u6,*) 'Erroneous property type:',trim(CHARTYPE)
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
! SDMXYSR, etc    DM/TDM for this iteration
! SDMXYZR2, etc   Accumulated DM/TDM
! TMPR,I          Temporary array for U*AU multiplication
! TDMZZ           DM/TDM as read from file
! SCR             Scratch for expansion of TDMZZ

call mma_allocate(SDMXYZR,NBTRI,3,Label='SDMXYZR')
call mma_allocate(SDMXYZI,NBTRI,3,Label='SDMXYZI')
call mma_allocate(SDMXYZR2,NBTRI,3,Label='SDMXYZR2')
call mma_allocate(SDMXYZI2,NBTRI,3,Label='SDMXYZI2')

SDMXYZR2(:,:) = Zero
SDMXYZI2(:,:) = Zero

call mma_allocate(TMPR,NBTRI,Label='TMPR')
call mma_allocate(TMPI,NBTRI,Label='TMPI')
TMPR(:) = Zero
TMPI(:) = Zero

call mma_allocate(SCR,NBTRI,Label='SCR')
! zeroed inside the loop
!SCR(:) = Zero

call mma_allocate(TDMZZ,NTDMZZ,Label='TDMZZ')
TDMZZ(:) = Zero

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! MAIN LOOP OVER KSF/LSF
! WRITTEN AS IN PRPROP
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
! CORRESPONDING SPIN-FREE STATES OF THE
! REQUESTED SPIN STATES
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
    call DCOPY_(NTDMZZ,[Zero],0,TDMZZ,1)

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

    ! This code expands the NTDMZZ-size matrix into
    ! an NBTRI-sized matrix
    ! SPECIAL CASE: DIAGONAL SYMMETRY BLOCKS.
    ! NOTE: During the folding, all off-diagonal values
    !  are doubled. This factor of 2 is made up when
    !  multiplying just the triangle with the AO matrix (DDOT)
    !  rather than the entire matrix

    ! TODO - get rid of if statements in inner loop
    !cccccccccccccccccc

    if (ISY12 == 1) then
      ! DIAGONAL SYMMETRY BLOCKS
      IOF = 0
      ITD = 0
      do ISY=1,NSYM
        NB = NBASF(ISY)
        if (NB == 0) cycle
        do J=1,NB
          do I=1,NB
            ITD = ITD+1
            TDM = TDMZZ(ITD)
            if (I >= J) then
              IJ = IOF+(I*(I-1))/2+J
              if (I > J) then
                if (ITYPE == 2) SCR(IJ) = SCR(IJ)+TDM
                if (ITYPE == 4) SCR(IJ) = SCR(IJ)+TDM
              end if
            else
              IJ = IOF+(J*(J-1))/2+I
              if (ITYPE == 2) SCR(IJ) = SCR(IJ)-TDM
              if (ITYPE == 4) SCR(IJ) = SCR(IJ)-TDM
            end if
            if (ITYPE == 1) SCR(IJ) = SCR(IJ)+TDM
            if (ITYPE == 3) SCR(IJ) = SCR(IJ)+TDM
          end do
        end do
        IOF = IOF+(NB*(NB+1))/2
      end do
    else
      ! GENERAL CASE, NON-DIAGONAL SYMMETRY BLOCKS
      ! THEN LOOP OVER ELEMENTS OF TDMZZ
      ITD = 0
      do ISY1=1,NSYM
        NB1 = NBASF(ISY1)
        if (NB1 == 0) cycle
        ISY2 = MUL(ISY1,ISY12)
        NB2 = NBASF(ISY2)
        if (NB2 == 0) cycle
        if (ISY1 > ISY2) then
          do J=1,NB2
            do I=1,NB1
              ITD = ITD+1
              TDM = TDMZZ(ITD)
              IJ = IOFF(ISY1)+I+NB1*(J-1)
              SCR(IJ) = SCR(IJ)+TDM
            end do
          end do
        else
          do J=1,NB2
            do I=1,NB1
              ITD = ITD+1
              TDM = TDMZZ(ITD)
              IJ = IOFF(ISY2)+J+NB2*(I-1)
              SCR(IJ) = SCR(IJ)-TDM
            end do
          end do
        end if
      end do
    end if

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
      call DAXPY_(NBTRI,One,SCR,1,SDMXYZR(:,3),1)
    else if ((ITYPE == 3) .or. (ITYPE == 4)) then
      if (iOpt == 1) then
        call DAXPY_(NBTRI,CGX*ROTMAT(1,1),SCR,1,SDMXYZR(:,1),1)
        call DAXPY_(NBTRI,CGY*ROTMAT(2,1),SCR,1,SDMXYZI(:,1),1)
        call DAXPY_(NBTRI,CG0*ROTMAT(3,1),SCR,1,SDMXYZR(:,1),1)

        call DAXPY_(NBTRI,CGX*ROTMAT(1,2),SCR,1,SDMXYZR(:,2),1)
        call DAXPY_(NBTRI,CGY*ROTMAT(2,2),SCR,1,SDMXYZI(:,2),1)
        call DAXPY_(NBTRI,CG0*ROTMAT(3,2),SCR,1,SDMXYZR(:,2),1)

        call DAXPY_(NBTRI,CGX*ROTMAT(1,3),SCR,1,SDMXYZR(:,3),1)
        call DAXPY_(NBTRI,CGY*ROTMAT(2,3),SCR,1,SDMXYZI(:,3),1)
        call DAXPY_(NBTRI,CG0*ROTMAT(3,3),SCR,1,SDMXYZR(:,3),1)
      else
        call DAXPY_(NBTRI,CGX,SCR,1,SDMXYZR(:,1),1)
        call DAXPY_(NBTRI,CGY,SCR,1,SDMXYZI(:,2),1)
        call DAXPY_(NBTRI,CG0,SCR,1,SDMXYZR(:,3),1)
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
      call DCOPY_(NBTRI,[Zero],0,TMPR,1)
      call DCOPY_(NBTRI,[Zero],0,TMPI,1)

      ! right side
      call DAXPY_(NBTRI,URR,SDMXYZR(:,IDIR),1,TMPR,1)
      call DAXPY_(NBTRI,-UIR,SDMXYZI(:,IDIR),1,TMPR,1)
      call DAXPY_(NBTRI,UIR,SDMXYZR(:,IDIR),1,TMPI,1)
      call DAXPY_(NBTRI,URR,SDMXYZI(:,IDIR),1,TMPI,1)

      ! left side
      call DAXPY_(NBTRI,URL,TMPR,1,SDMXYZR2(:,IDIR),1)
      call DAXPY_(NBTRI,UIL,TMPI,1,SDMXYZR2(:,IDIR),1)
      call DAXPY_(NBTRI,URL,TMPI,1,SDMXYZI2(:,IDIR),1)
      call DAXPY_(NBTRI,-UIL,TMPR,1,SDMXYZI2(:,IDIR),1)
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
  do I=1,NBTRI
    DENSOUT(1,I) = SDMXYZR2(I,1)
    DENSOUT(2,I) = SDMXYZR2(I,2)
    DENSOUT(3,I) = SDMXYZR2(I,3)
    DENSOUT(4,I) = SDMXYZI2(I,1)
    DENSOUT(5,I) = SDMXYZI2(I,2)
    DENSOUT(6,I) = SDMXYZI2(I,3)
  end do
else
  do I=1,NBTRI
    DENSOUT(1,I) = SDMXYZR2(I,3)
    DENSOUT(2,I) = SDMXYZR2(I,3)
    DENSOUT(3,I) = SDMXYZR2(I,3)
    DENSOUT(4,I) = SDMXYZI2(I,3)
    DENSOUT(5,I) = SDMXYZI2(I,3)
    DENSOUT(6,I) = SDMXYZI2(I,3)
  end do
end if

! Free memory
call mma_deallocate(SCR)
call mma_deallocate(TDMZZ)

call mma_deallocate(SDMXYZR)
call mma_deallocate(SDMXYZI)

call mma_deallocate(TMPI)
call mma_deallocate(TMPR)

call mma_deallocate(SDMXYZR2)
call mma_deallocate(SDMXYZI2)

call mma_deallocate(MAPMS)
call mma_deallocate(MAPSP)
call mma_deallocate(MAPST)

end subroutine SONATORBM
