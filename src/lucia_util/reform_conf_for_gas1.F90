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

subroutine REFORM_CONF_FOR_GAS1(ICONF_GAS,ICONF,IBORB,IBEL,MXPORB,NEL)
! Reform between local and global numbering of
! configuration for given GAS space
!
! Global => Local
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IBEL, NEL, ICONF(IBEL-1+NEL), IBORB, MXPORB
integer(kind=iwp), intent(inout) :: ICONF_GAS(MXPORB)
integer(kind=iwp) :: NTEST

ICONF_GAS(1:NEL) = ICONF(IBEL:IBEL+NEL-1)-IBORB+1

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Global => Local reform of conf'
  write(u6,*) ' ICONF_GAS :'
  call IWRTMA(ICONF_GAS,1,NEL,1,NEL)
  write(u6,*) ' Accessed part of ICONF'
  call IWRTMA(ICONF,1,IBEL-1+NEL,1,IBEL-1+NEL)
end if

end subroutine REFORM_CONF_FOR_GAS1
