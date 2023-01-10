!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine fin_run_use()

use RunFile_data, only: i_run_CA_used, i_run_DA_used, i_run_DS_used, i_run_IA_used, i_run_IS_used, nTocCA, nTocDA, nTocDS, nTocIA, &
                        nTocIS
use Definitions, only: iwp

implicit none
logical(kind=iwp), external :: Reduce_Prt

if (Reduce_Prt()) return
!need_abend = 0
call check_use(nTocCA,i_run_CA_used,'cArray')
call check_use(nTocDA,i_run_DA_used,'dArray')
call check_use(nTocDS,i_run_DS_used,'dScalar')
call check_use(nTocIA,i_run_IA_used,'iArray')
call check_use(nTocIS,i_run_IS_used,'iScalar')
!if (need_abend == 1) call abend()

return

end subroutine fin_run_use
