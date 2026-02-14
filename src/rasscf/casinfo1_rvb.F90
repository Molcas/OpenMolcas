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

subroutine casinfo1_rvb()

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: i2s_c, ireturn_rasscf, isym_c, nel_c, neltot_c, norb_c
logical(kind=iwp) :: iphex, oldex

! Information from molcas interface file 'JOBIPH' :
write(u6,'(a)') ' ------- Recover RASSCF-related information --------------------------------------'
call f_inquire('JOBIPH',iphex)
call f_inquire('JOBOLD',oldex)
if (iphex) then
  write(u6,'(/,a)') ' Using JOBIPH interface file.'
  !call systemf('cp -p JOBIPH JOBOLD')
  call Copy_JobIph('JOBIPH','JOBOLD')
else if (oldex) then
  write(u6,'(/,a)') ' Using JOBOLD interface file.'
  !call systemf('cp -p JOBOLD JOBIPH')
  call Copy_JobIph('JOBOLD','JOBIPH')
else
  write(u6,'(/,a)') ' Error: need either JOBOLD or JOBIPH file!'
  call abend_cvb()
end if

call rdjobiph_cvb('JOBIPH')
call setjobiph_cvb(nel_c,norb_c,i2s_c,isym_c,neltot_c)
call rasscf(ireturn_rasscf)
call clsfls_rasscf()
! rasscf will have overwritten jobiph ...
!call systemf('cp -p JOBOLD JOBIPH')
call Copy_JobIph('JOBOLD','JOBIPH')
write(u6,'(a)') ' ------- RASSCF-related information recovered ------------------------------------'

end subroutine casinfo1_rvb
