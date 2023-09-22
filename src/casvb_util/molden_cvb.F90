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

subroutine molden_cvb()

use rctfld_module,only: lRF
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"

integer(kind=iwp) :: iDisk
real(kind=wp) :: Dummy(1)

call daname_cvb(JOBIPH,'JOBIPH')
iDisk = 0
call idafile(JOBIPH,2,iadr15,15,iDisk)
if (.not. lRF) then
  Dummy(1) = Zero
  call interf(0,Dummy,0,1)
end if

return

end subroutine molden_cvb
