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
! Copyright (C) 1996-2006, Thorstein Thorsteinsson                     *
!               1996-2006, David L. Cooper                             *
!***********************************************************************

subroutine rules_cvb(chr)

implicit none
character(len=*), intent(in) :: chr

select case (chr)
  case ('MEM1')
    call chop1_cvb()
  case ('MEM2')
    call chop2_cvb()
  case ('MEM3')
    call chop3_cvb()
  case ('MEM4')
    call chop4_cvb()
  case ('MEM5')
    call chop5_cvb()
  case ('MEM6')
    call chop6_cvb()
  case ('MEM7')
    call chop7_cvb()
  case ('ORBFREE')
    call mkorbfree_cvb()
  case ('CIFREE')
    call mkcifree_cvb()
  case ('ICONFS')
    call mkiconfs_cvb()
  case ('GENDET')
    call mkciinfo_cvb()
    call mkvbinfo_cvb()
  case ('SYMELM')
    call mksymelm_cvb()
  case ('SYMINIT')
    call mksyminit_cvb()
  case ('CONSTRUC')
    call mkconstruc_cvb()
  case ('RDINT')
    ! do nothing
  case ('RDCAS')
    call mkrdcas_cvb()
  case ('SYMORBS')
    call mksymorbs_cvb()
  case ('SYMCVB')
    call mksymcvb_cvb()
  case ('GUESS')
    call mkguess_cvb()
  case ('ORBPERM')
    call mkorbperm_cvb()
  case ('TRNSPN')
    call mktrnspn_cvb()
  case ('STAT')
    call stat_cvb()
end select

return

end subroutine rules_cvb
