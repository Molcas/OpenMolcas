!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICEnSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2019, Oskar Weser                                      *
!***********************************************************************

module orthonormalization

use blockdiagonal_matrices, only: t_blockdiagonal, new, delete, from_raw, to_raw, from_symm_raw, blocksizes
use linalg_mod, only: Gram_Schmidt, Lowdin, Canonical
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, u6

implicit none
private

#include "intent.fh"

public :: t_ON_scheme, ON_scheme, ON_scheme_values, orthonormalize

! TODO: Should be changed to default construction in the future.
! As of July 2019 the Sun and PGI compiler have problems.
type :: t_ON_scheme_values
  integer :: no_ON, Gram_Schmidt, Lowdin, Canonical
end type t_ON_scheme_values

! TODO: Dear fellow MOLCAS developer of the future:
! Please replace the following explicit constructur
! with the default constructor
!     instance = type()
! as soon as possible.
! As of July 2019 the Sun compiler requires explicit construction
! for parameter variables. (Which is wrong IMHO.)
type(t_ON_scheme_values), parameter :: ON_scheme_values = t_ON_scheme_values(no_ON=1,Gram_Schmidt=2,Lowdin=3,Canonical=4)

type :: t_ON_scheme
  integer :: val = ON_scheme_values%Gram_Schmidt
end type t_ON_scheme

type(t_ON_scheme) :: ON_scheme = t_ON_scheme()

!>  @brief
!>    Orthonormalize a basis.
!>
!>  @author
!>    Oskar Weser
!>
!>  @details
!>  Reads in the overlap matrix and orthonormalizes accordingly.
!>
!>  @param[in] basis Can be raw memory i.e. a 1D double array
!>    or a t_blockdiagonal matrix.
!>    return type depends on input type.
!>  @param[in] scheme Optional argument. The orthonormalization scheme to use.
!>    The possibilities are Gram_Schmidt, Lowdin, Canonical, or no_ON
!>    (no_orthonormalization)
!>    as given by orthonormalization::ON_scheme_values.
!>    For a detailed explanation see \cite szabo_ostlund (p. 143).
interface orthonormalize
  module procedure :: orthonormalize_raw, orthonormalize_blocks
end interface orthonormalize

contains

subroutine orthonormalize_blocks(basis,scheme,ONB)

  use general_data, only: nSym, nBAs, nDel, nDelt, nSSH, nOrb
  use rasscf_global, only: nSec, nOrbt, nTot3, nTot4

  type(t_blockdiagonal), intent(in) :: basis(:)
  type(t_ON_scheme), intent(in) :: scheme
  type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
  type(t_blockdiagonal), allocatable :: S(:)
  integer :: n_to_ON(nSym), n_new(nSym)

  ! gfortran -O0 warning bug
  allocate(S(0))
  call new(S,blocksizes=blocksizes(basis))
  call read_S(S)

  select case (scheme%val)
    case (ON_scheme_values%no_ON)
      !continue
    case (ON_scheme_values%Lowdin)
      call Lowdin_Blocks(basis,S,ONB)
    case (ON_scheme_values%Canonical)
      n_to_ON(:) = nBas(:nSym)-nDel(:nSym)
      call Canonical_Blocks(basis,S,n_to_ON,ONB,n_new)
      call update_orb_numbers(n_to_ON,n_new,nDel,nSSH,nOrb,nDelt,nSec,nOrbt,nTot3,nTot4)
    case (ON_scheme_values%Gram_Schmidt)
      n_to_ON(:) = nBas(:nSym)-nDel(:nSym)
      call Gram_Schmidt_Blocks(basis,S,n_to_ON,ONB,n_new)
      call update_orb_numbers(n_to_ON,n_new,nDel,nSSH,nOrb,nDelt,nSec,nOrbt,nTot3,nTot4)
  end select

  call delete(S)

end subroutine orthonormalize_blocks

subroutine orthonormalize_raw(CMO,scheme,ONB_v)

  use general_data, only: nBas, nSym

  real(wp), intent(in) :: CMO(:)
  type(t_ON_scheme), intent(in) :: scheme
  real(wp), intent(out) :: ONB_v(:)
  type(t_blockdiagonal), allocatable :: basis(:), ONB(:)

  ! gfortran -O0 warning bug
  allocate(basis(0))
  call new(basis,blocksizes=nBAS(:nSym))
  call new(ONB,blocksizes=nBAS(:nSym))

  call from_raw(CMO,basis)
  call orthonormalize(basis,scheme,ONB)
  call to_raw(ONB,ONB_v)

  call delete(ONB)
  call delete(basis)

end subroutine orthonormalize_raw

! TODO: It would be nice, to use `impure elemental`
! instead of the manual looping.
! As of July 2019 some compilers don't support it.
subroutine Lowdin_Blocks(basis,S,ONB)

  type(t_blockdiagonal), intent(in) :: basis(:), S(:)
  type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
  integer :: i

  do i=1,size(basis)
    call Lowdin(basis(i)%blck,ONB(i)%blck,S(i)%blck)
  end do

end subroutine Lowdin_Blocks

! TODO: It would be nice, to use `impure elemental`
! instead of the manual looping.
! As of July 2019 some compilers don't support it.
subroutine Canonical_Blocks(basis,S,n_to_ON,ONB,n_new)

  type(t_blockdiagonal), intent(in) :: basis(:), S(:)
  integer, intent(in) :: n_to_ON(:)
  type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
  integer, intent(out) :: n_new(:)
  integer :: i

  do i=1,size(basis)
    call Canonical(basis(i)%blck,n_to_ON(i),ONB(i)%blck,n_new(i),S(i)%blck)
  end do

end subroutine Canonical_Blocks

! TODO: It would be nice, to use `impure elemental`
! instead of the manual overloading.
! As of July 2019 some compilers don't support it.
subroutine Gram_Schmidt_Blocks(basis,S,n_to_ON,ONB,n_new)

  type(t_blockdiagonal), intent(in) :: basis(:), S(:)
  integer, intent(in) :: n_to_ON(:)
  type(t_blockdiagonal), intent(_OUT_) :: ONB(:)
  integer, intent(out) :: n_new(:)
  integer :: i

  do i=1,size(basis)
    call Gram_Schmidt(basis(i)%blck,n_to_ON(i),ONB(i)%blck,n_new(i),S(i)%blck)
  end do

end subroutine Gram_Schmidt_Blocks

subroutine update_orb_numbers(n_to_ON,nNew,nDel,nSSH,nOrb,nDelt,nSec,nOrbt,nTot3,nTot4)

  use general_data, only: nSym
  use PrintLevel, only: USUAL
  use output_ras, only: IPRLOC
# include "warnings.h"

  integer, intent(in) :: n_to_ON(:), nNew(:)
  integer, intent(inout) :: nDel(:), nSSH(:), nOrb(:), nDelt, nSec, nOrbt, nTot3, nTot4
  integer :: iSym, remove(nSym), total_remove

  remove = n_to_ON(:nSym)-nNew(:nSym)
  total_remove = sum(remove(:nSym))
  if (total_remove > 0) then
    do iSym=1,nSym
      if (nSSH(iSym) < remove(iSym)) then
        call WarningMessage(2,'Orthonormalization Error')
        write(u6,*) 'Exact or very near linear dependence '
        write(u6,*) 'forces RASSCF to stop execution.'
        write(u6,*) 'Symmetry block:',iSym
        write(u6,*) 'Effective NR of orthonormal orbs:',nNew(iSym)
        write(u6,*) 'Earlier number of deleted orbs:',nDel(iSym)
        write(u6,*) 'Earlier number of secondary orbs:',nSSH(iSym)
        write(u6,*) 'New number of deleted orbs:',nDel(iSym)+remove(iSym)
        write(u6,*) 'New number of secondary orbs:',nSSH(iSym)-remove(iSym)
        call quit(_RC_GENERAL_ERROR_)
      else if (iPRLoc(1) >= USUAL) then
        call WarningMessage(1,'Orthonormalization Warning')
        write(u6,*) 'Exact or very near linear dependence'
        write(u6,*) 'forces RASSCF to delete additional orbitals.'
        write(u6,*) 'Symmetry block:',iSym
        write(u6,*) 'Earlier number of deleted orbs =',nDel(iSym)
        write(u6,*) 'New number of deleted orbs =',nDel(iSym)+remove(iSym)
      end if
    end do
    nDel(:nSym) = nDel(:nSym)+remove(:nSym)
    nSSH(:nSym) = nSSH(:nSym)-remove(:nSym)
    nOrb(:nSym) = nOrb(:nSym)-remove(:nSym)
    nDelt = nDelt+total_remove
    nSec = nSec-total_remove
    nOrbt = nOrbt-total_remove
    nTot3 = sum((nOrb(:nSym)+nOrb(:nSym)**2)/2)
    nTot4 = sum(nOrb(:nSym)**2)
  end if

end subroutine update_orb_numbers

subroutine read_raw_S(S_buffer)

  use OneDat, only: sNoOri

  real(wp), intent(inout) :: S_buffer(:)
  integer :: i_Rc, i_Opt, i_Component, i_SymLbl
  character(len=8) :: Label
# include "warnings.h"

  i_Rc = 0
  i_Opt = ibset(0,sNoOri)
  i_Component = 1
  i_SymLbl = 1
  Label = 'Mltpl  0'
  call RdOne(i_Rc,i_Opt,Label,i_Component,S_buffer,i_SymLbl)
  if (i_rc /= 0) then
    write(u6,*) ' RASSCF is trying to orthonormalize orbitals but'
    write(u6,*) ' could not read overlaps from ONEINT. Something'
    write(u6,*) ' is wrong with the file, or possibly with the'
    write(u6,*) ' program. Please check.'
    call quit(_RC_IO_ERROR_READ_)
  end if

end subroutine read_raw_S

subroutine read_S(S)

  use general_data, only: nBas, nSym, nActEl
  use rasscf_global, only: nFr, nIn, Tot_Nuc_Charge
  use PrintLevel, only: USUAL
  use output_ras, only: IPRLOC

  type(t_blockdiagonal) :: S(nSym)
  integer :: size_S_buffer
  real(wp) :: Mol_Charge
  real(wp), allocatable :: S_buffer(:)
# include "warnings.h"

  size_S_buffer = sum(nBas(:nSym)*(nBas(:nSym)+1)/2)
  call mma_allocate(S_buffer,size_S_buffer+4)
  call read_raw_S(S_buffer)
  Tot_Nuc_Charge = S_buffer(size_S_buffer+4)
  call from_symm_raw(S_buffer,S)
  call mma_deallocate(S_buffer)

  Mol_Charge = Tot_Nuc_Charge-real(2*(nFr+nIn)+nActEl,kind=wp)
  call put_dscalar('Total Charge    ',Mol_Charge)
  if (IPRLOC(1) >= USUAL) then
    write(u6,*)
    write(u6,'(6x,A,f8.2)') 'Total molecular charge',Mol_Charge
  end if

end subroutine read_S

end module orthonormalization
