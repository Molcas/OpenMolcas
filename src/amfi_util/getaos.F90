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

subroutine getAOs(lhigh)
!bs get expansions of atomic orbitals in contracted functions

use AMFI_global, only: AOcoeffs, ncontrac, noccorb, occup
use Definitions, only: iwp

implicit none
integer(kind=iwp), intent(in) :: lhigh
integer(kind=iwp) :: icont, iorbital, lrun, Lu_33
logical(kind=iwp) :: EX
character(len=12) :: occread, occtext
character(len=18) :: textnorbmf, textnorbmf2
integer(kind=iwp), external :: IsFreeUnit

occtext = 'OCCUPATION: '
textnorbmf = 'Number of orbitals'
call f_inquire('AO-expansion',EX)
if (.not. EX) then
  !BS write(u6,*) 'get occupations from DATA-block'
  call getAOs2(lhigh)
else
  Lu_33 = IsFreeUnit(33)
  call molcas_open(Lu_33,'AO-expansion')
  !open(unit=Lu_33,file='AO-expansion',status='UNKNOWN')
  !BS write(u6,*) 'Orbitals for mean-field'
  do lrun=0,lhigh
    !BS write(u6,'(A3,I3)') 'L= ',lrun
    read(Lu_33,'(A18,I3)') textnorbmf2,noccorb(lrun)
    if (textnorbmf /= textnorbmf2) call SysAbendMsg('getAOs','wrong keyword for number of orbitals in getAOs',' ')
    !BS write(u6,*) 'number of orbitals ',noccorb(lrun)
    do iorbital=1,noccorb(lrun)
      read(Lu_33,'(A12,F5.3)') occread,occup(iorbital,lrun)
      !BS write(u6,'(A,F8.4)') occtext,occup(iorbital,lrun)
      if (occread /= occtext) call SysAbendMsg('getAOs','error reading AOs',' ')
      read(Lu_33,*) (AOcoeffs(icont,iorbital,lrun),icont=1,ncontrac(lrun))
      !BS write(u6,'(8F10.4)') (AOcoeffs(icont,iorbital,lrun),icont=1,ncontrac(lrun))
      !BS write(u6,*) ' '
      read(Lu_33,*)
    end do
  end do
  close(Lu_33)
end if

return

end subroutine getAOs
