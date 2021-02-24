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

subroutine gugadrt(ireturn)

use gugadrt_global, only: nci_dim
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp), intent(out) :: ireturn
real(kind=wp) :: sc, sc0, sc1
real(kind=r8), external :: seconds

sc0 = seconds()

call gugainit()

call gugadrt_mole_inf()
call gugadrt_paras_calculate()
call arrange_orbital_molcas()

call gugadrt_dbl_upwalk()       ! add by wyb 01.9.5
call gugadrt_ext_downwalk()     ! add by wyb 01.9.5
call gugadrt_active_drt()       ! add by wyb 01.9.5

call add_info('CI_DIM',[dble(nci_dim)],1,1)
call gugadrt_gugafinalize()
ireturn = 0

sc1 = seconds()
sc = sc1-sc0
write(6,910) sc

return

910 format(/,1x,'end of job, takes ',f10.2,1x,'seconds'/)

end subroutine gugadrt
