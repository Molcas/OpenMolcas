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

subroutine WRT_CONF_LIST(ICONF,NCONF_FOR_OPEN,MAXOP,NELEC)
! Write list of configurations, given in packed form
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: ICONF(*), MAXOP, NCONF_FOR_OPEN(MAXOP+1), NELEC
integer(kind=iwp) :: IB, IOPEN, JCONF, NCONF_OP, NOCC_ORB

IB = 1
do IOPEN=0,MAXOP
  NCONF_OP = NCONF_FOR_OPEN(IOPEN+1)
  if (NCONF_OP /= 0) then
    write(u6,*) ' Number of configurations with ',IOPEN,' open orbitals is ',NCONF_OP

    NOCC_ORB = IOPEN+(NELEC-IOPEN)/2
    do JCONF=1,NCONF_OP
      call IWRTMA(ICONF(IB),1,NOCC_ORB,1,NOCC_ORB)
      IB = IB+NOCC_ORB
    end do
  end if
end do

end subroutine WRT_CONF_LIST
