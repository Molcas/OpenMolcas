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
! Copyright (C) 1987, Per Ake Malmqvist                                *
!               2018, Jesper Norell                                    *
!               2018, Joel Creutzberg                                  *
!               2023, Ignacio Fdez. Galvan                             *
!***********************************************************************

subroutine SODYSORB(NSS,USOR,USOI,DYSAMPS,NZ,SOENE)

use Index_Functions, only: nTri_Elem
use rassi_global_arrays, only: JBNUM, SFDYS, SODYSAMPS, SODYSAMPSI, SODYSAMPSR
use OneDat, only: sNoNuc, sNoOri
use Cntrl, only: DYSEXPORT, DYSEXPSO, MLTPLT, NSTATE
use Symmetry_Info, only: nIrrep
use rassi_data, only: NBASF
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NSS, NZ
real(kind=wp) :: USOR(NSS,NSS), USOI(NSS,NSS), DYSAMPS(NSTATE,NSTATE), SOENE(NSS)
integer(kind=iwp) :: ICMP, IEIG, IOPT, IRC, ISTATE, ISY, ISYLAB, JEIG, JOB1, JSTATE, LUNIT, MPLET1, MSPROJ, NB, NDUM, NOFF, NSZZ, &
                     ORBNUM, SFI, SFJ, SFTOT, SODYSCIND, SOTOT, ZI, ZJ
real(kind=wp) :: AMPLITUDE, CII, CIMAG, CIR, CJI, CJR, CREAL, MSPROJI, MSPROJJ
character(len=80) :: TITLE
character(len=30) :: Filename
character(len=8) :: Label
integer(kind=iwp), allocatable :: SO2SF(:)
real(kind=wp), allocatable :: AMPS(:), DYSEN(:), MSPROJS(:), SODYSCMOI(:), SODYSCMOR(:), SODYSCOFSI(:), SODYSCOFSR(:), SZZ(:), &
                              SZZFULL(:,:)
integer(kind=iwp), external :: IsFreeUnit

! +++ J.Norell 2018

! Calculates spin-orbit Dyson orbitals
! The routine was in some part adapted from DO_SONATORB

! Computes SO Dyson amplitudes by expanding the SF results with
! the SO eigenvectors (of the complex Hamiltonian)
! 1. (Fast): Compute the SO amplitudes directly from the SF amplitudes
!    (approximation) for all states
! 2. (Slower): Compute the full SO Dyson orbitals for the requested
!    initial states and export them to .molden format. SO amplitudes
!    are correctly calculated for these states.

! IFG: Added DysOrb export, not sure it's correct

! ****************************************************************

! Setup SO2SF list which contains the original SF state numbers
! as a function of the SO state number
! And MSPROJS which saves their ms projections for later use

call mma_allocate(SO2SF,NSS,Label='SO2SF')
call mma_allocate(MSPROJS,NSS,Label='MSPROJS')
SOTOT = 0
SFTOT = 0
do ISTATE=1,NSTATE
  JOB1 = JBNUM(ISTATE)
  MPLET1 = MLTPLT(JOB1)
  SFTOT = SFTOT+1

  do MSPROJ=-MPLET1+1,MPLET1-1,2
    SOTOT = SOTOT+1
    SO2SF(SOTOT) = SFTOT
    MSPROJS(SOTOT) = MSPROJ

  end do ! DO MSPROJ1=-MPLET1+1,MPLET1-1,2
end do ! DO ISTATE=1,NSTATE

if (.not. DYSEXPORT) then ! Approximative amplitude calculation

  ! Now construct the SF dysamps in the multiplicity expanded basis
  ! (initially all real, therefore put into SODYSAMPSR)
  SODYSAMPSR(:,:) = Zero
  SODYSAMPSI(:,:) = Zero
  do JSTATE=1,NSS
    do ISTATE=JSTATE+1,NSS
      SFJ = SO2SF(JSTATE)
      SFI = SO2SF(ISTATE)
      SODYSAMPSR(JSTATE,ISTATE) = DYSAMPS(SFJ,SFI)
      SODYSAMPSR(ISTATE,JSTATE) = DYSAMPS(SFJ,SFI)
    end do
  end do

  ! Now perform the transformation from SF dysamps to SO dysamps
  ! by combining the multiplicity expanded SF dysamps with the
  ! SO eigenvector in the ZTRNSF routine.
  call ZTRNSF(NSS,USOR,USOI,SODYSAMPSR,SODYSAMPSI)

  ! Compute the magnitude of the complex amplitudes as an approximation
  SODYSAMPSR(:,:) = SODYSAMPSR*SODYSAMPSR
  SODYSAMPSI(:,:) = SODYSAMPSI*SODYSAMPSI
  SODYSAMPS(:,:) = sqrt(SODYSAMPSR+SODYSAMPSI)

end if ! Approximative amplitude calculation

! ****************************************************************

if (DYSEXPORT) then

  ! Export part of the routine and exact calculation of amplitudes

  ! Read in the atomic overlap matrix, that will be needed below for
  ! for normalization of DOs
  ! (Code from mksxy)
  NSZZ = sum(nTri_Elem(NBASF(1:nIrrep)))
  call mma_allocate(SZZ,NSZZ,Label='SZZ')
  IRC = -1
  IOPT = ibset(ibset(0,sNoOri),sNoNuc)
  ICMP = 1
  ISYLAB = 1
  Label = 'MLTPL  0'
  call RDONE(IRC,IOPT,Label,ICMP,SZZ,ISYLAB)
  if (IRC /= 0) then
    write(u6,*)
    write(u6,*) '      *** ERROR IN SUBROUTINE SODYSORB ***'
    write(u6,*) '     OVERLAP INTEGRALS ARE NOT AVAILABLE'
    write(u6,*)
    call ABEND()
  end if

  call mma_allocate(SZZFULL,NZ,NZ,Label='SZZFULL')

  ! SZZ is originally given in symmetry-blocked triangular form,
  ! lets make it a full matrix for convenience
  SZZFULL(:,:) = Zero
  NDUM = 1
  NOFF = 0
  do ISY=1,nIrrep
    NB = NBASF(ISY)
    do ZJ=1,NB
      do ZI=1,ZJ
        SZZFULL(ZJ+NOFF,ZI+NOFF) = SZZ(NDUM)
        SZZFULL(ZI+NOFF,ZJ+NOFF) = SZZ(NDUM)
        NDUM = NDUM+1
      end do
    end do
    NOFF = NOFF+NB
  end do
  call mma_deallocate(SZZ)

  ! ****************************************************************

  ! Multiply together with the SO eigenvector coefficients with the SF
  ! Dyson orbital coefficients in the atomic basis to obtain the full
  ! SO Dyson orbitals

  ! Multiply together with the SO eigenvector coefficients with the SF
  ! Dyson orbital coefficients in the atomic basis to obtain
  ! SO Dyson orbitals

  call mma_allocate(SODYSCMOR,NZ*NSS,Label='SODYSCMOR')
  call mma_allocate(SODYSCMOI,NZ*NSS,Label='SODYSCMOI')
  call mma_allocate(AMPS,NSS,Label='AMPS')
  call mma_allocate(DYSEN,NSS,Label='DYSEN')
  call mma_allocate(SODYSCOFSR,NZ,Label='SODYSCOFSR')
  call mma_allocate(SODYSCOFSI,NZ,Label='SODYSCOFSI')

  SODYSAMPS(:,:) = Zero
  ! For all requested initial states J and all final states I
  do JSTATE=1,DYSEXPSO

    ! For each initial state JSTATE up to DYSEXPSFSO we will
    ! gather all the obtained Dysorbs
    ! and export to a shared .molden file
    SODYSCIND = 0 ! Orbital coeff. index
    ORBNUM = 0 ! Dysorb index for given JSTATE
    SODYSCMOR(:) = Zero ! Real orbital coefficients
    SODYSCMOI(:) = Zero ! Imaginary orbital coefficients
    DYSEN(:) = Zero ! Orbital energies
    AMPS(:) = Zero ! Transition amplitudes (shown as occupations)

    do ISTATE=JSTATE+1,NSS

      ! Reset values for next state combination
      SODYSCOFSR(:) = Zero
      SODYSCOFSI(:) = Zero

      ! Iterate over the eigenvector components of both states
      do JEIG=1,NSS

        ! Coefficient of first state
        CJR = USOR(JEIG,JSTATE)
        CJI = USOI(JEIG,JSTATE)
        ! Find the corresponding SF states
        SFJ = SO2SF(JEIG)

        do IEIG=1,NSS

          ! Coefficient of second state
          CIR = USOR(IEIG,ISTATE)
          CII = USOI(IEIG,ISTATE)
          ! Find the corresponding SF states
          SFI = SO2SF(IEIG)

          ! Check change in ms projection
          MSPROJJ = MSPROJS(JEIG)
          MSPROJI = MSPROJS(IEIG)
          ! Check |delta ms|=0.5 selection rule
          if (abs(MSPROJJ-MSPROJI) /= 1) cycle

          if (DYSAMPS(SFJ,SFI) > 1.0e-5_wp) then
            ! Multiply together coefficients
            CREAL = CJR*CIR+CJI*CII
            CIMAG = CJR*CII-CJI*CIR
            ! Multiply with the corresponding SF Dyson orbital
            SODYSCOFSR(:) = SODYSCOFSR(:)+CREAL*SFDYS(:,SFJ,SFI)
            SODYSCOFSI(:) = SODYSCOFSI(:)+CIMAG*SFDYS(:,SFJ,SFI)
          end if

        end do ! IEIG
      end do ! JEIG

      ! Normalize the overlap of SODYSCOFS expanded orbitals with the
      ! atomic overlap matrix SZZ to obtain correct amplitudes

      AMPLITUDE = Zero

      do ZI=1,NZ
        AMPLITUDE = AMPLITUDE+sum(SZZFULL(:,ZI)*(SODYSCOFSR(:)*(SODYSCOFSR(ZI)-SODYSCOFSI(ZI))+ &
                                                 SODYSCOFSI(:)*(SODYSCOFSR(ZI)+SODYSCOFSI(ZI))))
      end do ! ZI

      AMPLITUDE = sqrt(AMPLITUDE)
      SODYSAMPS(JSTATE,ISTATE) = AMPLITUDE
      SODYSAMPS(ISTATE,JSTATE) = AMPLITUDE

      ! Export Re and Im part of the coefficients
      if (AMPLITUDE > 1.0e-5_wp) then
        SODYSCMOR(SODYSCIND+1:SODYSCIND+NZ) = SODYSCOFSR(:)
        SODYSCMOI(SODYSCIND+1:SODYSCIND+NZ) = SODYSCOFSI(:)
        SODYSCIND = SODYSCIND+NZ
        ORBNUM = ORBNUM+1
        DYSEN(ORBNUM) = SOENE(ISTATE)-SOENE(JSTATE)
        AMPS(ORBNUM) = AMPLITUDE*AMPLITUDE
      end if

    end do ! ISTATE

    ! If at least one orbital was found, export it/them
    if (ORBNUM > 0) then
      write(filename,'(A,I0,A3)') 'MD_DYS.SO.',JSTATE,'.Re'
      call Molden_DysOrb(filename,DYSEN,AMPS,SODYSCMOR,ORBNUM,NZ)
      write(filename,'(A,I0,A3)') 'MD_DYS.SO.',JSTATE,'.Im'
      call Molden_DysOrb(filename,DYSEN,AMPS,SODYSCMOI,ORBNUM,NZ)

      ! This does not work for SO-Dyson orbitals, because they may
      ! contain contributions from several irreps.
      ! Either that's a bug elsewhere in the code, or the InpOrb
      ! format is not adequate for these orbitals.
      write(filename,'(A,I0,A3)') 'DYSORB.SO.',JSTATE,'.Re'
      LUNIT = IsFreeUnit(50)
      write(TITLE,'(A,I0,A)') '* Spin-orbit Dyson orbitals for state ',JSTATE,' (real part)'
      call WRVEC_DYSON(filename,LUNIT,nIrrep,NBASF,ORBNUM,SODYSCMOR,AMPS,DYSEN,trim(TITLE),NZ)
      write(filename,'(A,I0,A3)') 'DYSORB.SO.',JSTATE,'.Im'
      write(TITLE,'(A,I0,A)') '* Spin-orbit Dyson orbitals for state ',JSTATE,' (imaginary part)'
      call WRVEC_DYSON(filename,LUNIT,nIrrep,NBASF,ORBNUM,SODYSCMOI,AMPS,DYSEN,trim(TITLE),NZ)
      close(LUNIT)
    end if

  end do ! JSTATE

end if

call mma_deallocate(SODYSCMOR)
call mma_deallocate(SODYSCMOI)
call mma_deallocate(AMPS)
call mma_deallocate(DYSEN)
call mma_deallocate(SODYSCOFSR)
call mma_deallocate(SODYSCOFSI)
call mma_deallocate(SZZFULL)

! ****************************************************************

! Free all the allocated memory

call mma_deallocate(SO2SF)
call mma_deallocate(MSPROJS)

end subroutine SODYSORB
