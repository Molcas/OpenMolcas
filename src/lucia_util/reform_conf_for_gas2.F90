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
subroutine REFORM_CONF_FOR_GAS2(ICONF_GAS,ICONF,IBORB,IBEL,MXPORB,NEL)
! Reform between local and global numbering of
! configuration for given GAS space
!
! Local => GLobal
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: MXPORB, ICONF_GAS(MXPORB), IBORB, IBEL, NEL
integer(kind=iwp), intent(inout) :: ICONF(IBEL-1+NEL)

ICONF(IBEL:IBEL+NEL-1) = ICONF_GAS(1:NEL)+IBORB-1

#ifdef _DEBUGPRINT_
write(u6,*) ' Local => Global reform of conf'
write(u6,*) ' ICONF_GAS :'
call IWRTMA(ICONF_GAS,1,NEL,1,NEL)
write(u6,*) ' Accessed part of ICONF'
call IWRTMA(ICONF,1,IBEL-1+NEL,1,IBEL-1+NEL)
#endif

end subroutine REFORM_CONF_FOR_GAS2
