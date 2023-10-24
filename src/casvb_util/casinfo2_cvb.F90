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

subroutine casinfo2_cvb()

use casvb_global, only: iorclos_c, iorcore_c, iorocc_c, istms2_c, istnel_c, istsy_c, mcore_c, mxirrep, mxstsy_ci, mxstt_ci, &
                        nstats_c, nstsym_c, weight_c
use Definitions, only: iwp

implicit none
integer(kind=iwp) :: i2s_c, isym_c, nel_c, neltot_c, norb_c

! Information from molcas interface file "JOBIPH":
call rdjobiph_cvb('JOBIPH')
call setjobiph_cvb(iorcore_c,iorclos_c,iorocc_c,mxirrep,nstsym_c,weight_c,istnel_c,istsy_c,istms2_c,nstats_c,mxstt_ci,mxstsy_ci, &
                   nel_c,norb_c,i2s_c,isym_c,mcore_c,neltot_c)

return

end subroutine casinfo2_cvb
