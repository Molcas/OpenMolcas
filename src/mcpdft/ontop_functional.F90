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
  use definitions,only:wp,iwp
  use Functionals, only: Get_Func_Type
  use nq_Info, only: meta_GGA_Type1, meta_GGA_Type2, Other_Type

  implicit none
  private

  integer(iwp):: Invalid_OTF=-1, Translated=1, FullyTranslated=2

  type :: OTFNAL_t
    real(kind=wp) :: lambda = 0.0d0
    character(len=80) :: otxc = ""
    character(len=80) :: xc = ""
  contains
    procedure :: is_hybrid
  endtype

  interface OTFNAL_t
    module procedure :: new
  endinterface

  public :: OTFNAL_t, Invalid_OTF, Translated, FullyTranslated, get_base, func_type

contains

  type(OTFNAL_t) function new(otxc,lambda)
    use definitions,only:u6
    implicit none

    real(kind=wp),intent(in) :: lambda
    character(len=80),intent(in) :: otxc

    new%lambda = lambda
    new%otxc = otxc
    call upcase(new%otxc)

    if(.not. valid_otxc(new%otxc)) then
      call warningmessage(2,"Wrong on-top functional for MC-PDFT")
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' Current on-top functionals are:   '
      write(u6,*) ' T :  translated functionals       '
      write(u6,*) ' FT:  fully translated functionals '
      write(u6,*) ' e.g. T:PBE for tPBE functional    '
      write(u6,*) ' **********************************'
      call abend()
    endif

    new%xc = get_base(new%otxc)

    if(.not.supported_functional(new%xc)) then
      call warningmessage(2,"functional type not supported")
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' Type-2 meta-GGA (functionals with '
      write(u6,*) ' orbital Laplacians as ingredients)'
      write(u6,*) ' are not supported in MC-PDFT.     '
      write(u6,*) ' **********************************'
      call abend()
    endif

    if(is_hybrid_xc(new%xc)) then
      call warningmessage(2,"Hybrid functionals not supported")
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' MC-PDFT does not translate hybrid '
      write(u6,*) ' functionals. If you want to run   '
      write(u6,*) ' hybrid MC-PDFT, use the LAMBda    '
      write(u6,*) ' keyword instead.                  '
      write(u6,*) '                                   '
      write(u6,*) ' EXAMPLE:                          '
      write(u6,*) '  tPBE0 = 75% tPBE + 25% MCSCF.    '
      write(u6,*) ' Usage:                            '
      write(u6,*) '  KSDFT=T:PBE                      '
      write(u6,*) '  LAMB =0.25                       '
      write(u6,*) ' **********************************'
      call abend()
    endif

    if (is_ft_meta(new%otxc)) then
      call warningmessage(2,"fully translated meta-GGA not supported")
      write(u6,*) ' ************* ERROR **************'
      write(u6,*) ' meta-GGA functional is currently  '
      write(u6,*) ' not compatible with full transla- '
      write(u6,*) ' tion as the full translation of   '
      write(u6,*) ' kinetic energy is not developed   '
      write(u6,*) ' **********************************'
      call abend()
    endif

  endfunction

  function valid_otxc(otxc)
    implicit none
    logical :: valid_otxc
    character(len=80),intent(in) :: otxc
    valid_otxc = otxc(1:2) == "T:" .or. otxc(1:3) == "FT:"
  endfunction

  function is_hybrid_xc(xc_base)
    implicit none
    logical :: is_hybrid_xc
    character(len=80),intent(in) :: xc_base
    real(kind=wp),external :: get_exfac
    is_hybrid_xc = abs(get_exfac(xc_base)) > 1.0d-8
  endfunction

  function is_hybrid(self)
    implicit none
    logical :: is_hybrid
    class(OTFNAL_t),intent(in) :: self
    is_hybrid = self%lambda > 0.0d0
  endfunction

  function func_type(xc_base)
    implicit none
    character(len=80),intent(in) :: xc_base
    integer(kind=iwp) :: func_type
    func_type=Get_Func_Type(xc_base)
  endfunction

  function get_base(otxc)
    implicit none
    character(len=80) :: get_base
    character(len=80),intent(in) :: otxc
    get_base = otxc(index(otxc,"T:")+2:)
  endfunction

  function is_ft_meta(otxc)
    implicit none
    character(len=80),intent(in) :: otxc
    logical :: is_ft_meta
    is_ft_meta = func_type(get_base(otxc))==meta_GGA_type1.and.get_ontop_type(otxc)==FullyTranslated
  endfunction

  function supported_functional(base_func)
    implicit none
    integer(kind=iwp) :: func_type
    character(len=80),intent(in) :: base_func
    logical :: supported_functional
    func_type=Get_Func_Type(base_func)
    if (func_type==meta_GGA_Type2.or.func_type==Other_Type) then
      supported_functional=.False.
    else
      supported_functional=.True.
    end if
  endfunction


  function get_ontop_type(otxc)
!   -1: invalid option
!    1: translated
!    2: fully-translated
!    3: converted (saved for future DCFT)
    implicit none
    integer(kind=iwp) :: get_ontop_type
    character(len=80),intent(in) :: otxc
    if (otxc(1:2)=="T:") then
      get_ontop_type=Translated
    else if (otxc(1:3)=="FT:") then
      get_ontop_type=FullyTranslated
    else
      get_ontop_type =Invalid_OTF
    end if
  endfunction

endmodule
