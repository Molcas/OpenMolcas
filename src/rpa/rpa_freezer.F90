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
! Copyright (C) 2013, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine RPA_Freezer()

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Figure out symmetry distribution of frozen occupied orbitals.

use Definitions, only: iwp

implicit none
#include "rpa_config.fh"
#include "rpa_data.fh"
#include "WrkSpc.fh"
integer(kind=iwp) :: iUHF, iSym, iSpin, ip_Fro, l_Fro
logical(kind=iwp) :: Freeze, Prnt
integer(kind=iwp), external :: RPA_iUHF

! set restricted(1)/unrestricted(2)
iUHF = RPA_iUHF()

! freeze orbitals (if requested)
iSpin = 1
Freeze = nFreeze(iSpin) > 0
do while (.not. Freeze .and. iSpin < iUHF)
  iSpin = iSpin+1
  Freeze = Freeze .or. nFreeze(iSpin) > 0
end do
if (Freeze) then
  Prnt = iPrint >= 0
  l_Fro = nSym
  call GetMem('OccFrz','Allo','Inte',ip_Fro,l_Fro)
  do iSpin=1,iUHF
    if (nFreeze(iSpin) > 0) then
      call RPA_Frz(nFreeze(iSpin),Prnt,nSym,Work(ip_OccEn(iSpin)),nFro(1,iSpin),nOcc(1,iSpin),iWork(ip_Fro))
      do iSym=1,nSym
        nFro(iSym,iSpin) = nFro(iSym,iSpin)+iWork(ip_Fro-1+iSym)
      end do
    end if
  end do
  call GetMem('OccFrz','Free','Inte',ip_Fro,l_Fro)
end if

! correct number of active occupied orbitals
do iSpin=1,iUHF
  do iSym=1,nSym
    nOcc(iSym,iSpin) = nOcc(iSym,iSpin)-nFro(iSym,iSpin)
  end do
end do

end subroutine RPA_Freezer
