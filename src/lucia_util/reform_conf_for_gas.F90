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

subroutine REFORM_CONF_FOR_GAS(ICONF_GAS,ICONF,IBORB,IBEL,MXPORB,NEL,IWAY)
! Reform between local and global numbering of
! configuration for given GAS space
!
! IWAY = 1 : Global => Local
! IWAY = 2 : Local => GLobal
!
! Jeppe Olsen, November 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IBORB, IBEL, MXPORB, NEL, IWAY
integer(kind=iwp), intent(inout) :: ICONF_GAS(MXPORB), ICONF(IBEL-1+NEL)
integer(kind=iwp) :: IEL, NTEST

if (IWAY == 1) then
  do IEL=1,NEL
    ICONF_GAS(IEL) = ICONF(IBEL-1+IEL)-IBORB+1
  end do
else if (IWAY == 2) then
  do IEL=1,NEL
    ICONF(IBEL-1+IEL) = ICONF_GAS(IEL)+IBORB-1
  end do
else
  write(u6,*) ' Problem in REFORM_CONF ..., IWAY = ',IWAY
  !stop ' Problem in REFORM_CONF ..., IWAY ='
  call SYSABENDMSG('lucia_util/reform_conv','Internal error','')
end if

NTEST = 0
if (NTEST >= 100) then
  if (IWAY == 1) then
    write(u6,*) ' Global => Local reform of conf'
  else
    write(u6,*) ' Local => Global reform of conf'
  end if
  write(u6,*) ' ICONF_GAS :'
  call IWRTMA(ICONF_GAS,1,NEL,1,NEL)
  write(u6,*) ' Accessed part of ICONF'
  call IWRTMA(ICONF,1,IBEL-1+NEL,1,IBEL-1+NEL)
end if

end subroutine REFORM_CONF_FOR_GAS
