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
! Copyright (C) 2012, Per Ake Malmqvist                                *
!               2019, Stefano Battaglia                                *
!***********************************************************************

subroutine xdwinit(Heff,H0,U0,nState)

use PrintLevel, only: DEBUG, INSANE, USUAL, VERBOSE
use caspt2_global, only: CMO, CMO_Internal, CMOPT2, do_grad, DREF, FIFA, FIMO, iPrGlb, LUONEM, NCMO
use caspt2_module, only: CIThr, iAd1m, iSCF, mState, nAshT, nConf, NoTri, STSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: NState
real(kind=wp), intent(inout) :: Heff(Nstate,Nstate)
real(kind=wp), intent(out) :: H0(Nstate,Nstate), U0(Nstate,Nstate)
integer(kind=iwp) :: I, iDisk, iState, J
real(kind=wp) :: FIJ, wgt
logical(kind=iwp) :: Initiate
real(kind=wp), allocatable :: CI(:), CIRef(:,:), CIXMS(:), DAVE(:), HONE(:)

! Allocate memory for CI array state averaged 1-RDM
call mma_allocATE(CI,NCONF,Label='CI')
call mma_allocate(DAVE,size(DREF),Label='DAVE')
DAVE(:) = Zero

! Set the weight for the density averaging
wgt = One/real(Nstate,kind=wp)

! Loop over all states to compute the state-average density matrix
do Istate=1,Nstate

  if (ISCF /= 0) then
    ! Special case for a single Slater determinant
    CI(1) = One
  else
    ! Get the CI array
    call loadCI(CI,Istate)
  end if

  ! Compute 1-particle active density matrix GAMMA1
  call POLY1(CI,nConf)

  ! Restructure GAMMA1 as DREF array
  call GETDREF(DREF,size(DREF))

  ! Average the density
  DAVE(:) = DAVE(:)+DREF(:)
  !call DAXPY_(SIZE(DREF),wgt,DREF,1,DAVE,1)

end do
DAVE(:) = Wgt*DAVE(:)

if (IPRGLB >= INSANE) then
  write(u6,*) ' State-average 1-RDM'
  do I=1,NASHT
    write(u6,'(1x,14f10.6)') (DAVE((I*(I-1))/2+J),J=1,I)
  end do
  write(u6,*)
end if

! Copy the state-average 1-RDM into DREF and release memory
! for both the CI array and DAVE
DREF(:) = DAVE(:)
call mma_deallocate(CI)
call mma_deallocate(DAVE)

! Load CASSCF MO coefficients
call mma_allocate(CMO_Internal,NCMO,Label='CMO_Internal')
CMO => CMO_Internal
iDisk = IAD1M(1)
call ddafile(LUONEM,2,CMO,NCMO,iDisk)

! Allocate memory for HONE and call for the computation of the
! one- and two-electron in the CMO basis.
call mma_allocate(HONE,NoTri,Label='HONE')
Initiate = .true.

! Build the state-average Fock matrix in MO basis
call MkFock(CMO,nCMO,FIMO,size(FIMO),FIFA,size(FIFA),DREF,size(DREF),HONE,size(HONE),Initiate)

call mma_deallocate(HONE)

! Loop again over all states to compute H0 in the model space
! Loop over ket functions
do J=1,Nstate
  ! Loop over bra functions
  do I=1,Nstate
    ! Compute matrix element <I|F|J> and store it into H0
    FIJ = Zero
    call FOPAB(FIFA,size(FIFA),I,J,FIJ)
    H0(I,J) = FIJ
  end do
end do
! End of loop over states

if (IPRGLB >= USUAL) then
  write(u6,*)
  write(u6,*) ' H0 in the original model space basis:'
  call prettyprint(H0,Nstate,Nstate)
end if

! Diagonalize H0 in the model space
call eigen(H0,U0,Nstate)

! Transform the Fock matrix in the new basis
call transmat(H0,U0,Nstate)
if (IPRGLB >= USUAL) then
  write(u6,*) ' H0 eigenvectors:'
  call prettyprint(U0,Nstate,Nstate)
end if
if (IPRGLB >= DEBUG) then
  write(u6,*) ' H0 in the rotated model space basis:'
  call prettyprint(H0,Nstate,Nstate)
end if

! As well as Heff
call transmat(Heff,U0,Nstate)
if (IPRGLB >= VERBOSE) then
  write(u6,*) ' Heff[1] in the rotated model space basis:'
  call prettyprint(Heff,Nstate,Nstate)
end if

! Mix the CI arrays according to the H0 eigenvectors. Assume we can
! put all the original ones in memory, but put the resulting vectors
! one by one in a buffer.
if (IPRGLB >= VERBOSE) then
  write(u6,'(A)') ' The CASSCF states are now rotated according to the H0 eigenvectors'
  write(u6,*)
end if

call mma_allocate(CIref,nConf,Nstate,Label='CIRef')
! Load the CI arrays into memory
do I=1,Nstate
  call loadCI(CIref(:,I),I)
end do

call mma_allocate(CIXMS,Nconf,Label='CIXMS')
do J=1,Nstate
  ! Transform the states
  call dgemm_('N','N',Nconf,1,Nstate,One,CIREF,Nconf,U0(:,J),Nstate,Zero,CIXMS,Nconf)

  ! Write the rotated CI coefficients back into LUCIEX and REPLACE the
  ! original unrotated CASSCF states. Note that the original states
  ! are still available in the JobIph file
  call writeCI(CIXMS,J)

  if (IPRGLB >= VERBOSE) then
    write(u6,'(1x,a,i3)') ' The CI coefficients of rotated model state nr. ',MSTATE(J)
    call PRWF_CP2(STSYM,NCONF,CIXMS,CITHR)
  end if
end do

if (do_grad) CMOPT2(:) = CMO(:)

! Release all memory
call mma_deallocate(CIRef)
call mma_deallocate(CIXMS)
call mma_deallocate(CMO_Internal)
nullify(CMO)

end subroutine xdwinit

