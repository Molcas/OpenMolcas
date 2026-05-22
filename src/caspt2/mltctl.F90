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

use caspt2_global, only: iPrGlb
use PrintLevel, only: TERSE, USUAL, VERBOSE
use caspt2_module, only: ENERGY, IfChol, IFDW, IFRMS, IFXMS, JMS, MSTATE, NSTATE
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Half
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(inout) :: HEFF(NSTATE,NSTATE)
real(kind=wp), intent(out) :: EIGVEC(NSTATE,NSTATE)
real(kind=wp), intent(in) :: U0(Nstate,Nstate)
integer(kind=iwp) :: I, IEND, II0, IJ, ISTA, J, LAXITY, NHTRI, NUMAT
real(kind=wp) :: DSHIFT
character(len=8) :: INLAB
character(len=3) :: variant
real(kind=wp), allocatable :: HTRI(:), UMAT(:,:), Utmp(:,:)
integer(kind=iwp), external :: Cho_X_GetTol

if (IPRGLB >= TERSE) then
  call CollapseOutput(1,'Multi-State CASPT2 section:')
  write(u6,'(A)') repeat('*',80)
  write(u6,*) ' MULTI-STATE CASPT2 SECTION'
  if (IPRGLB >= USUAL) then
    write(u6,'(A)') repeat('-',80)
    write(u6,*)
  end if
end if

! Write out the effective Hamiltonian, for use in e.g. RASSI:
INLAB = 'HEFF'
call put_darray(INLAB,HEFF,NSTATE**2)

! Analyze the effective Hamiltonian:
DSHIFT = Zero
if (HEFF(1,1) <= -100.0_wp) DSHIFT = -real(ceiling(HEFF(1,1)),kind=wp)
do I=1,NSTATE
  HEFF(I,I) = HEFF(I,I)-DSHIFT
end do

if ((IPRGLB >= TERSE) .and. (DSHIFT /= Zero)) write(u6,*) ' Output diagonal energies have been shifted. Add ',DSHIFT

if ((IPRGLB >= VERBOSE) .or. JMS) then
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
if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,*) ' Effective Hamiltonian matrix (Symmetric):'
  do ISTA=1,NSTATE,5
    IEND = min(ISTA+4,NSTATE)
    write(u6,*)
    write(u6,'(1x,5I16)') (MSTATE(I),I=ISTA,IEND)
    do I=ISTA,NSTATE
      II0 = (I*(I-1))/2
      write(u6,'(1x,I3,3X,5F16.8)') MSTATE(I),(HTRI(II0+J),J=ISTA,min(I,IEND))
    end do
  end do
end if
call unitmat(UMAT,NSTATE)
call NIDiag(HTRI,UMAT,NSTATE,NSTATE)
call JACORD(HTRI,UMAT,NSTATE,NSTATE)
do I=1,NSTATE
  ENERGY(I) = DSHIFT+HTRI((I*(I+1))/2)
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
  write(u6,*)
  write(u6,'(6X,A,A)') ' Total ',trim(variant)//'-CASPT2 energies:'
  do I=1,NSTATE
    call PrintResult(u6,'(6x,A,I3,5X,A,F16.8)',trim(variant)//'-CASPT2 Root',I,'Total energy:',ENERGY(I),1)
  end do
end if

if (IPRGLB >= USUAL) then
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

! In automatic verification calculations, the precision is lower
! in case of Cholesky calculation.
LAXITY = 8
if (IfChol) LAXITY = Cho_X_GetTol(LAXITY)
call Add_Info('E_MSPT2',ENERGY,nState,LAXITY)

end subroutine MLTCTL
