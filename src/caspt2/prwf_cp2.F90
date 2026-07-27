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
! Copyright (C) 1994, Per Ake Malmqvist                                *
!***********************************************************************
!--------------------------------------------*
! 1994  PER-AAKE MALMQUIST                   *
! DEPARTMENT OF THEORETICAL CHEMISTRY        *
! UNIVERSITY OF LUND                         *
! SWEDEN                                     *
!--------------------------------------------*

subroutine PRWF_CP2(ISYCI,NCO,CI,THR)

use sguga_states, only: CIS
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ISYCI, NCO
real(kind=wp), intent(in) :: CI(NCO), THR
integer(kind=iwp), parameter :: istate=1

write(u6,'(A)') repeat('-',80)
write(u6,'(a,es9.2)') ' CI COEFFICIENTS LARGER THAN ',THR
call PRWF1_CP2(CIS(istate)%NOCSF,CIS(istate)%IOCSF,CIS(istate)%NOW,CIS(istate)%IOW,ISYCI,CI,nCO,THR,CIS(istate)%nMidV)

end subroutine PRWF_CP2
