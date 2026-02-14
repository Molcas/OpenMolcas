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
! Copyright (C) 2024, Matthew R. Hennefarth                            *
!***********************************************************************

module ontop_functional

use Functionals, only: Get_Func_Type
use nq_Info, only: meta_GGA_Type1, meta_GGA_Type2, Other_Type
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
private

integer(kind=iwp), parameter :: Invalid_OTF = -1, Translated = 1, FullyTranslated = 2

type :: OTFNAL_t
  real(kind=wp) :: lambda = Zero
  character(len=80) :: otxc = ''
  character(len=80) :: xc = ''
  contains
  procedure :: is_hybrid
  procedure :: energy_ot
end type

interface OTFNAL_t
  module procedure :: new
end interface

public :: FullyTranslated, get_base, Invalid_OTF, OTFNAL_t, Translated

contains

function new(otxc,lambda)

  type(OTFNAL_t) :: new
  character(len=80), intent(in) :: otxc
  real(kind=wp), intent(in) :: lambda

  new%lambda = lambda
  new%otxc = otxc
  call upcase(new%otxc)

  if (.not. valid_otxc(new%otxc)) then
    call warningmessage(2,'Wrong on-top functional for MC-PDFT')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' Current on-top functionals are:   '
    write(u6,*) ' T :  translated functionals       '
    write(u6,*) ' FT:  fully translated functionals '
    write(u6,*) ' e.g. T:PBE for tPBE functional    '
    write(u6,*) ' **********************************'
    call abend()
  end if

  new%xc = get_base(new%otxc)

  if (.not. supported_functional(new%xc)) then
    call warningmessage(2,'functional type not supported')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' Type-2 meta-GGA (functionals with '
    write(u6,*) ' orbital Laplacians as ingredients)'
    write(u6,*) ' are not supported in MC-PDFT.     '
    write(u6,*) ' **********************************'
    call abend()
  end if

  if (is_hybrid_xc(new%xc)) then
    call warningmessage(2,'Hybrid functionals not supported')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' MC-PDFT does not translate hybrid '
    write(u6,*) ' functionals. If you want to run   '
    write(u6,*) ' hybrid MC-PDFT, use the LAMBda    '
    write(u6,*) ' keyword instead.                  '
    write(u6,*) '                                   '
    write(u6,*) ' EXAMPLE:                          '
    write(u6,*) '  tPBE0 = 75% tPBE + 25% MCSCF.    '
    write(u6,*) ' Usage:                            '
    write(u6,*) '  FUNC=T:PBE                       '
    write(u6,*) '  LAMB=0.25                        '
    write(u6,*) ' **********************************'
    call abend()
  end if

  if (is_ft_meta(new%otxc)) then
    call warningmessage(2,'fully translated meta-GGA not supported')
    write(u6,*) ' ************* ERROR **************'
    write(u6,*) ' meta-GGA functional is currently  '
    write(u6,*) ' not compatible with full transla- '
    write(u6,*) ' tion as the full translation of   '
    write(u6,*) ' kinetic energy is not developed   '
    write(u6,*) ' **********************************'
    call abend()
  end if

end function new

function valid_otxc(otxc)

  logical(kind=iwp) :: valid_otxc
  character(len=80), intent(in) :: otxc

  valid_otxc = (otxc(1:2) == 'T:') .or. (otxc(1:3) == 'FT:')

end function valid_otxc

function is_hybrid_xc(xc_base)

  logical(kind=iwp) :: is_hybrid_xc
  character(len=80), intent(in) :: xc_base
  real(kind=wp), external :: get_exfac

  is_hybrid_xc = (abs(get_exfac(xc_base)) > 1.0e-8_wp)

end function is_hybrid_xc

function is_hybrid(self)

  logical(kind=iwp) :: is_hybrid
  class(OTFNAL_t), intent(in) :: self

  is_hybrid = (self%lambda > Zero)

end function is_hybrid

function get_base(otxc)

  character(len=80) :: get_base
  character(len=80), intent(in) :: otxc

  get_base = otxc(index(otxc,'T:')+2:)

end function get_base

function is_ft_meta(otxc)

  logical(kind=iwp) :: is_ft_meta
  character(len=80), intent(in) :: otxc

  is_ft_meta = (Get_Func_type(get_base(otxc)) == meta_GGA_type1) .and. (get_ontop_type(otxc) == FullyTranslated)

end function is_ft_meta

function supported_functional(base_func)

  logical(kind=iwp) :: supported_functional
  character(len=80), intent(in) :: base_func
  integer(kind=iwp) :: func_type

  func_type = Get_Func_Type(base_func)
  if ((func_type == meta_GGA_Type2) .or. (func_type == Other_Type)) then
    supported_functional = .false.
  else
    supported_functional = .true.
  end if

end function supported_functional

! -1: invalid option
!  1: translated
!  2: fully-translated
!  3: converted (saved for future DCFT)
function get_ontop_type(otxc)

  integer(kind=iwp) :: get_ontop_type
  character(len=80), intent(in) :: otxc

  if (otxc(1:2) == 'T:') then
    get_ontop_type = Translated
  else if (otxc(1:3) == 'FT:') then
    get_ontop_type = FullyTranslated
  else
    get_ontop_type = Invalid_OTF
  end if

end function get_ontop_type

function energy_ot(self,folded_dm1,folded_dm1s,casdm2,charge)

  use rctfld_module, only: lrf
  use rasscf_global, only: dftfock, exfac, nacpr2, noneq, potnuc
  use general_data, only: ispin, nash, nfro, nish, nsym, ntot1
  use stdalloc, only: mma_allocate, mma_deallocate

  real(kind=wp) :: energy_ot
  class(OTFNAL_t), intent(in) :: self
  real(kind=wp), intent(in) :: folded_dm1(ntot1), folded_dm1s(ntot1), casdm2(nacpr2)
  integer(kind=iwp), intent(in) :: charge
  integer(kind=iwp) :: charge_
  logical(kind=iwp) :: first, dff, do_dft
  real(kind=wp), allocatable :: dummy1(:), dummy2(:)

  call put_darray('D1ao',folded_dm1,ntot1)
  call put_darray('D1sao',folded_dm1s,ntot1)
  !call put_darray('D1mo',casdm1,nacpar)
  call put_darray('P2mo',casdm2,nacpr2)

  call mma_allocate(dummy1,ntot1,label='dummy1')
  call mma_allocate(dummy2,ntot1,label='dummy2')
  dummy1(:) = Zero
  dummy2(:) = Zero

  first = .true.
  dff = .false.
  do_dft = .true.

  call put_iarray('nfro',nfro,nsym)
  call put_iarray('nash',nash,nsym)
  call put_iarray('nish',nish,nsym)

  ! Perhaps ideally, we should reword how drvxv (and its children) handles the AO to MO transformation on the grid. It seems like
  ! perhaps we are doing redundant transformations by retransforming AOs
  charge_ = charge
  call drvxv(dummy1,dummy2,folded_dm1,potnuc,ntot1,first,dff,noneq,lrf,self%otxc,exfac,charge_,ispin,dftfock,do_dft)

  call mma_deallocate(dummy1)
  call mma_deallocate(dummy2)

  call get_dscalar('CASDFT energy',energy_ot)

end function energy_ot

end module
