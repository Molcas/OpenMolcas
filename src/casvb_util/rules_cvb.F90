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

character*(*) chr

if (chr == 'MEM1') then
  call chop1_cvb()
else if (chr == 'MEM2') then
  call chop2_cvb()
else if (chr == 'MEM3') then
  call chop3_cvb()
else if (chr == 'MEM4') then
  call chop4_cvb()
else if (chr == 'MEM5') then
  call chop5_cvb()
else if (chr == 'MEM6') then
  call chop6_cvb()
else if (chr == 'MEM7') then
  call chop7_cvb()
else if (chr == 'ORBFREE') then
  call mkorbfree_cvb()
else if (chr == 'CIFREE') then
  call mkcifree_cvb()
else if (chr == 'ICONFS') then
  call mkiconfs_cvb()
else if (chr == 'GENDET') then
  call mkciinfo_cvb()
  call mkvbinfo_cvb()
else if (chr == 'SYMELM') then
  call mksymelm_cvb()
else if (chr == 'SYMINIT') then
  call mksyminit_cvb()
else if (chr == 'CONSTRUC') then
  call mkconstruc_cvb()
else if (chr == 'RDINT') then
  call mkrdint_cvb()
else if (chr == 'RDCAS') then
  call mkrdcas_cvb()
else if (chr == 'SYMORBS') then
  call mksymorbs_cvb()
else if (chr == 'SYMCVB') then
  call mksymcvb_cvb()
else if (chr == 'GUESS') then
  call mkguess_cvb()
else if (chr == 'ORBPERM') then
  call mkorbperm_cvb()
else if (chr == 'TRNSPN') then
  call mktrnspn_cvb()
else if (chr == 'STAT') then
  call stat_cvb()
end if

return

end subroutine rules_cvb
