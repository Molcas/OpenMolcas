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
! Copyright (C) 2017, Felix Plasser                                    *
!***********************************************************************

subroutine wfa(ireturn)

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(out) :: ireturn
#include "warnings.fh"

#ifdef _WFA_
integer(kind=iwp) :: iLU, ist, ien
character(len=180), external :: Get_Ln
character(len=180) :: Line
character(len=900) :: Inp
#endif

ireturn = _RC_NOT_AVAILABLE_

#ifdef _WFA_
! Read the input file and parse as a string to libwfa
ist = 1
call SpoolInp(iLU)
do while (Line(1:4) /= 'END ')
  Line = Get_ln(iLU)

  ien = ist+len(trim(Line))
  Inp(ist:ien) = trim(Line)
  ist = ien+1

  call Normal(Line)
end do

write(u6,*) 'Starting wavefunction analysis ...'
call wfa_driver(ireturn,Inp)
#else
write(u6,*) 'WFA module not installed!'
#endif

end subroutine wfa
