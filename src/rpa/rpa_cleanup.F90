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

subroutine RPA_Cleanup(irc)

! Thomas Bondo Pedersen (CTCC,UiO), July 2013.
!
! Clean up after RPA run (deallocate etc.)

implicit none
integer irc
#include "rpa_config.fh"
#include "rpa_data.fh"

integer RPA_iUHF
external RPA_iUHF

integer i

irc = 0

! Set "Relax Method" on Runfile
call Put_cArray('Relax Method',RPAModel,8)

! Deallocate memory
do i=1,RPA_iUHF()
  if (l_CMO(i) > 0) then
    call GetMem('CMO(RPA)','Free','Real',ip_CMO(i),l_CMO(i))
  end if
  ip_CMO(i) = 0
  l_CMO(i) = 0
  if (l_EMO(i) > 0) then
    call GetMem('EMO(RPA)','Free','Real',ip_EMO(i),l_EMO(i))
  end if
  ip_EMO(i) = 0
  l_EMO(i) = 0
  if (l_OccEn(i) > 0) then
    call GetMem('OccEn','Free','Real',ip_OccEn(i),l_OccEn(i))
  end if
  ip_OccEn(i) = 0
  l_OccEn(i) = 0
  if (l_VirEn(i) > 0) then
    call GetMem('OccEn','Free','Real',ip_VirEn(i),l_VirEn(i))
  end if
  ip_VirEn(i) = 0
  l_VirEn(i) = 0
end do

end subroutine RPA_Cleanup
