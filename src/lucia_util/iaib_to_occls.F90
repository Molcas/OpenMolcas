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

subroutine IAIB_TO_OCCLS(IAGRP,IATP,IBGRP,IBTP,IOC)
! Find address of the occupation class corresponding to given types
! of alpha and beta strings
!
! Jeppe Olsen, December 2001

use lucia_data, only: IBSPGPFTP, IOCLS, MXPNGAS, NELFSPGP, NGAS, NMXOCCLS
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: IAGRP, IATP, IBGRP, IBTP
integer(kind=iwp), intent(out) :: IOC
integer(kind=iwp) :: IABOCC(MXPNGAS), IATP_ABS, IBTP_ABS, INUM

IATP_ABS = IATP+IBSPGPFTP(IAGRP)-1
IBTP_ABS = IBTP+IBSPGPFTP(IBGRP)-1
!write(u6,*) ' IATP, IBTP, IAGRP, IBGRP = ',ATP,IBTP,IAGRP,IBGRP
!write(u6,*) ' IATP_ABS, IBTP_ABS ',IATP_ABS,IBTP_ABS

IABOCC(1:NGAS) = NELFSPGP(1:NGAS,IATP_ABS)+NELFSPGP(1:NGAS,IBTP_ABS)
! And the address of this occupation class
call CMP_IVEC_ILIST(IABOCC,IOCLS,NGAS,NMXOCCLS,INUM)

IOC = INUM

if (INUM == 0) then
  write(u6,*) ' Combination of alpha and beta string not found as occ-class'
  write(u6,*) ' Occ of alpha, Occ of beta, Occ of alpha+beta'
  call IWRTMA(NELFSPGP(1,IATP_ABS),1,NGAS,1,NGAS)
  call IWRTMA(NELFSPGP(1,IBTP_ABS),1,NGAS,1,NGAS)
  call IWRTMA(IABOCC,1,NGAS,1,NGAS)
  !stop ' Combination of alpha and beta string not found as occ-class'
  call SYSABENDMSG('lucia_util/iaib_to_occls','Internal error','')
end if

end subroutine IAIB_TO_OCCLS
