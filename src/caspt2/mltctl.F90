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
! Copyright (C) 1998, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1998  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine MLTCTL(HEFF,EIGVEC,U0)

use Index_Functions, only: nTri_Elem
use caspt2_global, only: iPrGlb
use PrintLevel, only: TERSE, USUAL, VERBOSE
use caspt2_module, only: ENERGY, HZERO, IfChol, IFDW, IFRMS, IFXMS, JMS, MSTATE, NSTATE, PT2Method
use SC_NEVPT2, only: Do_FIC, Do_SC, SC_prop, ENERGY_SC, HEFF_SC, UEFF_SC
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: HEFF(NSTATE,NSTATE), EIGVEC(NSTATE,NSTATE)
real(kind=wp), intent(in) :: U0(Nstate,Nstate)
integer(kind=iwp) :: I, IEND, II0, IJ, ISTA, J, LAXITY, NHTRI, NUMAT, iloop, nloop
real(kind=wp) :: DSHIFT, tmp
character(len=8) :: INLAB
character(len=5) :: variant
real(kind=wp), allocatable :: HTRI(:), UMAT(:,:), Utmp(:,:)
logical(kind=iwp) :: PrRes
integer(kind=iwp), external :: Cho_X_GetTol

if (IPRGLB >= TERSE) then
  call CollapseOutput(1,'Multi-State '//trim(PT2Method)//' section:')
  write(u6,'(A)') repeat('*',80)
  write(u6,*) ' MULTI-STATE '//trim(PT2Method)//' SECTION'
  if (IPRGLB >= USUAL) then
    write(u6,'(A)') repeat('-',80)
    write(u6,*)
  end if
end if

! Write out the effective Hamiltonian, for use in e.g. RASSI:
INLAB = 'HEFF'
if (HZERO /= 'DYALL' .or. .not. SC_prop) then
  call put_darray(INLAB,HEFF,NSTATE**2)
else
  call put_darray(INLAB,HEFF_SC,NSTATE**2)
end if

! Analyze the effective Hamiltonian:
DSHIFT = Zero
if (HEFF(1,1) <= -100.0_wp) DSHIFT = -real(ceiling(HEFF(1,1)),kind=wp)
do I=1,NSTATE
  HEFF(I,I) = HEFF(I,I)-DSHIFT
end do

if ((IPRGLB >= TERSE) .and. (DSHIFT /= Zero)) write(u6,*) ' Output diagonal energies have been shifted. Add ',DSHIFT

nloop = 1
if (HZERO == 'DYALL' .and. Do_SC) nloop = 2

do iloop = 1, nloop
  PrRes = .true.
  if (iloop == 1 .and. .not. Do_FIC) PrRes = .false.

  if (((IPRGLB >= VERBOSE) .or. JMS) .and. PrRes) then
    write(u6,*) ' Effective Hamiltonian matrix (Asymmetric):'
    do ISTA=1,NSTATE,5
      IEND = min(ISTA+4,NSTATE)
      write(u6,*)
      write(u6,'(1x,5I16)') (MSTATE(I),I=ISTA,IEND)
      do J=1,NSTATE
        write(u6,'(1x,I3,3X,5F16.8)') MSTATE(J),(HEFF(J,I),I=ISTA,IEND)
      end do
    end do
  end if

  ! Diagonalize:
  ! Use a symmetrized matrix, in triangular storage:
  NUMAT = NSTATE**2
  NHTRI = (NUMAT+NSTATE)/2
  call mma_allocate(UMAT,NSTATE,NSTATE,LABEL='UMAT')
  call mma_allocate(HTRI,NHTRI,LABEL='HTRI')
  IJ = 0
  do I=1,NSTATE
    HTRI(IJ+1:IJ+I) = Half*(HEFF(I,1:I)+HEFF(1:I,I))
    IJ = IJ+I
  end do
  if (IPRGLB >= USUAL .and. PrRes) then
    write(u6,*)
    write(u6,*) ' Effective Hamiltonian matrix (Symmetric):'
    do ISTA=1,NSTATE,5
      IEND = min(ISTA+4,NSTATE)
      write(u6,*)
      write(u6,'(1x,5I16)') (MSTATE(I),I=ISTA,IEND)
      do I=ISTA,NSTATE
        II0 = nTri_Elem(I-1)
        write(u6,'(1x,I3,3X,5F16.8)') MSTATE(I),(HTRI(II0+J),J=ISTA,min(I,IEND))
      end do
    end do
  end if
  call unitmat(UMAT,NSTATE)
  call NIDiag(HTRI,UMAT,NSTATE,NSTATE)
  call JACORD(HTRI,UMAT,NSTATE,NSTATE)
  do I=1,NSTATE
    ENERGY(I) = DSHIFT+HTRI(nTri_Elem(I))
    EIGVEC(:,I) = UMAT(:,I)
  end do
  call mma_deallocate(UMAT)
  call mma_deallocate(HTRI)

  if (IPRGLB >= TERSE) then
    if (IFRMS) then
      variant = 'RMS'
    else if (IFXMS .and. IFDW) then
      variant = 'XDW'
    else if (IFXMS) then
      variant = 'XMS'
    else if (IFDW) then
      variant = 'DW '
    else
      variant = 'MS '
    end if
    if (HZERO /= 'DYALL') then !! CASPT2
      write(u6,*)
      write(u6,'(6X,A,A)') ' Total ',trim(variant)//'-CASPT2 energies:'
      do I=1,NSTATE
        call PrintResult(u6,'(6x,A,I3,5X,A,F16.8)',trim(variant)//'-CASPT2 Root',I,'Total energy:',ENERGY(I),1)
      end do
    else if (HZERO == 'DYALL') then !! NEVPT2
      if (iloop == 1) variant = 'QD-PC'
      if (iloop == 2) variant = 'QD-SC'
      write(u6,*)
      write(u6,'(6X,A,A)') ' Total ',trim(variant)//'-NEVPT2 energies:'
      do I=1,NSTATE
        call PrintResult(u6,'(6x,A,I3,5X,A,F16.8)',trim(variant)//'-NEVPT2 Root',I,'Total energy:',ENERGY(I),1)
      end do
    end if
  end if

  if (IPRGLB >= USUAL .and. PrRes) then
    write(u6,*)
    write(u6,'(6X,A)') ' Eigenvectors:'
    do ISTA=1,NSTATE,5
      IEND = min(ISTA+4,NSTATE)
      do J=1,NSTATE
        write(u6,'(6x,5F16.8)') (EIGVEC(J,I),I=ISTA,IEND)
      end do
      write(u6,*)
    end do
    if (IFXMS .or. IFRMS) then
      ! Transform eigenvectors into the original input basis
      call mma_allocate(Utmp,Nstate,Nstate,Label='Utmp')
      call dgemm_('N','N',Nstate,Nstate,Nstate,One,U0,Nstate,eigvec,Nstate,Zero,Utmp,Nstate)
      write(u6,'(6X,A)') ' In terms of the input states:'
      do ISTA=1,NSTATE,5
        IEND = min(ISTA+4,NSTATE)
        do J=1,NSTATE
          write(u6,'(6x,5F16.8)') (Utmp(J,I),I=ISTA,IEND)
        end do
        write(u6,*)
      end do
      call mma_deallocate(Utmp)
    end if
    call CollapseOutput(0,'Multi-State CASPT2 section:')
    write(u6,*)
  end if

  ! Restore original effective Hamiltonian
  do I=1,NSTATE
    HEFF(I,I) = HEFF(I,I)+DSHIFT
  end do

  if (nloop == 2) then
    do i = 1, nstate
      do j = 1, nstate
        tmp = HEFF(j,i)
        HEFF(j,i) = HEFF_SC(j,i)
        HEFF_SC(j,i) = tmp
        tmp = EIGVEC(j,i)
        EIGVEC(j,i) = UEFF_SC(j,i)
        UEFF_SC(j,i) = tmp
      end do
    end do
    call dswap_(NSTATE,ENERGY,1,ENERGY_SC,1)
    if (iloop == 1) then
      do I=1,NSTATE
        HEFF(I,I) = HEFF(I,I)-DSHIFT
      end do
    end if
  end if
end do

! In automatic verification calculations, the precision is lower
! in case of Cholesky calculation.
LAXITY = 8
if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
if (HZERO /= 'DYALL' .or. .not. SC_prop) then
  call Add_Info('E_MSPT2',ENERGY,nState,LAXITY)
else
  call Add_Info('E_MSPT2',ENERGY_SC,nState,LAXITY)
  ENERGY(1:NSTATE) = ENERGY_SC(1:NSTATE)
  !! Replace the effective H and eigenvec with those for SC
  HEFF(1:NSTATE,1:NSTATE) = HEFF_SC(1:NSTATE,1:NSTATE)
  EIGVEC(1:NSTATE,1:NSTATE) = UEFF_SC(1:NSTATE,1:NSTATE)
end if

end subroutine MLTCTL
