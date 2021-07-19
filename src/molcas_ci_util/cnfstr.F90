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
! Copyright (C) 1989, Jeppe Olsen                                      *
!***********************************************************************

subroutine CNFSTR(ICONF,ITYP,IASTR,IBSTR,NORB,NAEL,NBEL,IDET,IPRODT,ISCR,SGN,IPREXH)
! An orbital configuration ICONF is given
! Obtain the corresponding alpha strings, IASTR
!        the corresponding beta  strings, IBSTR
!        the corresponding sign  array  , SGN
!
! Jeppe Olsen, Summer of '89

use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICONF(*), ITYP, NORB, NAEL, NBEL, IDET, IPRODT(*) !IFG
integer(kind=iwp), intent(inout) :: IPREXH
integer(kind=iwp), intent(_OUT_) :: IASTR(*), IBSTR(*), ISCR(*) !IFG
real(kind=wp), intent(out) :: SGN(IDET)
integer(kind=iwp) :: ICLOS, IOCC, IOPEN, IP, ISGN, JDET, JTYP, KLDETS, KLFREE, NEL, NTEST
!PAM06 NOTE: NAEL and NBEL can legally be=0.
!PAM06 ... IASTR(NAEL,*), IBSTR(NBEL,*) ...
#include "spinfo.fh"
#include "ciinfo.fh"

NEL = NAEL+NBEL
IOPEN = ITYP-1+MINOP
ICLOS = (NAEL+NBEL-IOPEN)/2
IOCC = IOPEN+ICLOS

! Spin orbital occupations of determinants of configuration

KLFREE = 1
KLDETS = KLFREE
KLFREE = KLDETS+IDET*(NAEL+NBEL)
! Pointer for determinants of prototype ITYP
IP = 1
do JTYP=1,ITYP-1
  IP = IP+NDTFTP(JTYP)*(JTYP-1+MINOP)
end do
call CNDET(ICONF,IPRODT(IP),IDET,NAEL+NBEL,NORB,IOPEN,ICLOS,ISCR(KLDETS),IPREXH)

! Separate determinants into strings and determine sign change

do JDET=1,IDET
  call DETSTR2(ISCR(KLDETS+(JDET-1)*NEL),IASTR(1+NAEL*(JDET-1)),IBSTR(1+NBEL*(JDET-1)),NEL,NAEL,NBEL,ISGN,ISCR(KLFREE),IPREXH)
  !PAM06 ... ,IASTR(1,JDET),IBSTR(1,JDET), ...
  SGN(JDET) = real(ISGN,kind=wp)
end do

NTEST = 0
if (NTEST >= 1) then
  write(u6,*) ' Output from CNFSTR'
  write(u6,*) ' =================='
  write(u6,*) ' Input configuration'
  call IWRTMA(ICONF,1,IOCC,1,IOCC)
  write(u6,*) ' Corresponding alpha and beta strings'
  call IWRTMA(IASTR,NAEL,IDET,NAEL,IDET)
  call IWRTMA(IBSTR,NBEL,IDET,NBEL,IDET)
  write(u6,*) ' SIGN ARRAY'
  call WRTMAT(SGN,1,IDET,1,IDET)
end if

return

end subroutine CNFSTR
