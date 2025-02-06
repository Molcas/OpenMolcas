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

function ILEX_FOR_CONF_NEW(ICONF,NOCC_ORB,NORB,NEL,IARCW,IDOREO,IREO_new,nconf_op,ib_occls)
! A configuration ICONF of NOCC_ORB orbitals are given
! ICONF(I) = IORB implies IORB is singly occupied
! ICONF(I) = -IORB implies that IORB is doubly occupied
!
! Find lexical address
!
! IF IDOREO /= 0, IREO is used to reorder lexical number
! Jeppe Olsen, November 2001

use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: ILEX_FOR_CONF_NEW
integer(kind=iwp) :: NOCC_ORB, ICONF(NOCC_ORB), NORb, NEL, IARCW(NORB,NEL,2), IDOREO, IREO_new(*), nconf_op, ib_occls
integer(kind=iwp) :: i, IEL, ILEX, IOCC, length, n_ave, n_end, n_fin, n_ini, NTEST

IEL = 0
ILEX = 1

do IOCC=1,NOCC_ORB
  if (ICONF(IOCC) > 0) then
    IEL = IEL+1
    ILEX = ILEX+IARCW(ICONF(IOCC),IEL,1)
  else if (ICONF(IOCC) < 0) then
    IEL = IEL+2
    ILEX = ILEX+IARCW(-ICONF(IOCC),IEL,2)
  end if
end do

if (IDOREO /= 0) then

  !length = nconf_per_open(2*nocc_orb-nel+1,1)
  length = nconf_op
  ntest = 0
  if (ntest >= 10) then
    write(u6,*) 'in ilex_for_conf_new, check if all is right'
    write(u6,*) 'iconf:',(iconf(i),i=1,nocc_orb)
    write(u6,*) 'ilex and offset:',ilex,ib_occls
    write(u6,*) 'nconf_op =',nconf_op
    write(u6,*) 'ireo_new:'
    call iwrtma(ireo_new,1,length,1,length)
  end if
  !call iwrtma(nconf_per_open,1,4,1,4)
  ! look for the position where the value is ilex
  if (ireo_new(1) == ilex+ib_occls-1) then
    n_fin = 1
  elseif (ireo_new(length) == ilex+ib_occls-1) then
    n_fin = length
  else
    n_ini = 1
    n_end = length
666 continue
    n_ave = (n_ini+n_end)/2
    if (ireo_new(n_ave) == ilex+ib_occls-1) goto 100
    if (ireo_new(n_ave) < ilex+ib_occls-1) n_ini = n_ave
    if (ireo_new(n_ave) > ilex+ib_occls-1) n_end = n_ave
    goto 666
100 continue
    n_fin = n_ave
  end if
  ilex_for_conf_new = n_fin

else
  ILEX_FOR_CONF_NEW = ILEX
end if

NTEST = 0
if (NTEST >= 100) then
  write(u6,*) ' Configuration'
  call IWRTMA(ICONF,1,NOCC_ORB,1,NOCC_ORB)
  write(u6,*) ' new Lexical number = ',ILEX
  if (IDOREO /= 0) write(u6,*) ' new Reordered number = ',ILEX_FOR_CONF_NEW
end if

end function ILEX_FOR_CONF_NEW
