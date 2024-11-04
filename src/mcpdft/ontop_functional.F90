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

  public :: OTFNAL_t, Invalid_OTF, Translated, FullyTranslated, get_base

contains

  type(OTFNAL_t) function new(otxc,lambda)
    use mcpdft_output,only:lf
    implicit none

    real(kind=wp),intent(in) :: lambda
    character(len=80),intent(in) :: otxc

    new%lambda = lambda
    new%otxc = otxc
    call upcase(new%otxc)

    if(.not. valid_otxc(new%otxc)) then
      call warningmessage(2,"Wrong on-top functional for MC-PDFT")
      write(lf,*) ' ************* ERROR **************'
      write(lf,*) ' Current on-top functionals are:   '
      write(lf,*) ' T :  translated functionals       '
      write(lf,*) ' FT:  fully translated functionals '
      write(lf,*) ' e.g. T:PBE for tPBE functional    '
      write(lf,*) ' **********************************'
      call abend()
    endif

    new%xc = get_base(new%otxc)

    if(.not.supported_functional(new%xc)) then
      call warningmessage(2,"functional type not supported")
      write(lf,*) ' ************* ERROR **************'
      write(lf,*) ' Type-2 meta-GGA (functionals with '
      write(lf,*) ' orbital Laplacians as ingredients)'
      write(lf,*) ' are not supported in MC-PDFT.     '
      write(lf,*) ' **********************************'
      call abend()
    endif

    if(is_hybrid_xc(new%xc)) then
      call warningmessage(2,"Hybrid functionals not supported")
      write(lf,*) ' ************* ERROR **************'
      write(lf,*) ' MC-PDFT does not translate hybrid '
      write(lf,*) ' functionals. If you want to run   '
      write(lf,*) ' hybrid MC-PDFT, use the LAMBda    '
      write(lf,*) ' keyword instead.                  '
      write(lf,*) '                                   '
      write(lf,*) ' EXAMPLE:                          '
      write(lf,*) '  tPBE0 = 75% tPBE + 25% MCSCF.    '
      write(lf,*) ' Usage:                            '
      write(lf,*) '  KSDFT=T:PBE                      '
      write(lf,*) '  LAMB =0.25                       '
      write(lf,*) ' **********************************'
      call abend()
    endif

    if (is_ft_meta(new%otxc)) then
      call warningmessage(2,"fully translated meta-GGA not supported")
      write(lf,*) ' ************* ERROR **************'
      write(lf,*) ' meta-GGA functional is currently  '
      write(lf,*) ' not compatible with full transla- '
      write(lf,*) ' tion as the full translation of   '
      write(lf,*) ' kinetic energy is not developed   '
      write(lf,*) ' **********************************'
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

  function get_base(otxc)
    implicit none
    character(len=80) :: get_base
    character(len=80),intent(in) :: otxc
    get_base = otxc(index(otxc,"T:")+2:)
  endfunction

  function is_ft_meta(otxc)
    implicit none
    integer(kind=iwp) :: func_type, ontop_type
    character(len=80),intent(in) :: otxc
    logical :: is_ft_meta
    func_type=Get_Func_Type(get_base(otxc))
    ontop_type=get_ontop_type(otxc)
    is_ft_meta = func_type==meta_GGA_type1.and.ontop_type==FullyTranslated
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
