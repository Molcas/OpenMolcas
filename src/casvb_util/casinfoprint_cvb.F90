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

subroutine casinfoprint_cvb()

use casvb_global, only: ipr, isym, isymv, ityp, mxirrep, nalf, nbet, nel, norb, nsym
use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: i, ii, incr, qisym(nsym)
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if ((ipr(1) >= 0) .and. (.not. up2date_cvb('CASPRINT'))) then
  write(u6,'(/,a,i4)') ' Number of active electrons :',nel
  write(u6,'(a,i4)') ' Number of active orbitals  :',norb
  write(u6,'(a,f4.1)') ' Total spin                 :',real(nalf-nbet,kind=wp)*Half
  if (nsym == 1) then
    write(u6,'(a,i4)') ' State symmetry             :',isym
  else
    incr = 0
    do i=1,mxirrep
      if (isymv(i) == 1) then
        incr = incr+1
        qisym(incr) = i
      end if
    end do
    write(u6,'(a,i4,7i3)') ' State symmetries           :',(qisym(ii),ii=1,nsym)
  end if
  write(u6,'(/,a,100i3)') ' Symmetries of active MOs   : ',(ityp(ii),ii=1,norb)
  call make_cvb('CASPRINT')
end if

return

end subroutine casinfoprint_cvb
