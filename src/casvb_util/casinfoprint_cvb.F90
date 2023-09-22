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

use Constants, only: Half
use Definitions, only: wp, iwp, u6

implicit none
#include "main_cvb.fh"
#include "print_cvb.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: i, ii, incr, iqisym
integer(kind=iwp), external :: mstacki_cvb
logical(kind=iwp), external :: up2date_cvb ! ... Make: up to date? ...

if ((ip(1) >= 0) .and. (.not. up2date_cvb('CASPRINT'))) then
  write(u6,'(/,a,i4)') ' Number of active electrons :',nel
  write(u6,'(a,i4)') ' Number of active orbitals  :',norb
  write(u6,'(a,f4.1)') ' Total spin                 :',real(nalf-nbet,kind=wp)*Half
  if (nsym == 1) then
    write(u6,'(a,i4)') ' State symmetry             :',isym
  else
    iqisym = mstacki_cvb(nsym)
    incr = 0
    do i=1,mxirrep
      if (isymv(i) == 1) then
        incr = incr+1
        iwork(incr+iqisym-1) = i
      end if
    end do
    write(u6,'(a,i4,7i3)') ' State symmetries           :',(iwork(ii+iqisym-1),ii=1,nsym)
    call mfreei_cvb(iqisym)
  end if
  write(u6,'(/,a,100i3)') ' Symmetries of active MOs   : ',(ityp(ii),ii=1,norb)
  call make_cvb('CASPRINT')
end if

return

end subroutine casinfoprint_cvb
