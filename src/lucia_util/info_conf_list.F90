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
! Copyright (C) 2001, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine INFO_CONF_LIST(NCONF_PER_OPEN,MAXOP,NEL,LENGTH_LIST,NCONF_TOT,IB_REO,IB_OCC)
! Info on configuration list form NCONF_PER_OPEN
!
! IB_REO : Offset for configuration  with given number of
!          open orbitals in of configuration reordering
! IB_OCC : Offset for configuration  with given number of
!          open orbitals in list of configuration occupations
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MAXOP, NCONF_PER_OPEN(MAXOP+1), NEL
integer(kind=iwp), intent(out) :: LENGTH_LIST, NCONF_TOT, IB_REO(MAXOP+1), IB_OCC(MAXOP+1)
integer(kind=iwp) :: JB_OCC, JB_REO, NOCOB, NOPEN

JB_REO = 1
JB_OCC = 1
!write(u6,*) ' MAXOP, NEL = ',MAXOP,NEL
do NOPEN=0,MAXOP
  IB_REO(NOPEN+1) = JB_REO
  IB_OCC(NOPEN+1) = JB_OCC
  if (mod(NEL-NOPEN,2) == 0) then
    NOCOB = NOPEN+(NEL-NOPEN)/2
    JB_OCC = JB_OCC+NOCOB*NCONF_PER_OPEN(NOPEN+1)
    JB_REO = JB_REO+NCONF_PER_OPEN(NOPEN+1)
  end if
end do

LENGTH_LIST = JB_OCC-1
NCONF_TOT = JB_REO-1

#ifdef _DEBUGPRINT_
write(u6,*) ' NCONF_PER_OPEN list'
call IWRTMA(NCONF_PER_OPEN,1,MAXOP+1,1,MAXOP+1)
write(u6,*) ' Length of configuration list :',LENGTH_LIST
write(u6,*) ' Total number of configurations : ',NCONF_TOT
#endif

end subroutine INFO_CONF_LIST
