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
! Copyright (C) 1993, Jeppe Olsen                                      *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NEWTYP_MCLR(INCLS,INTP,IACOP,ITPOP,NOP,OUTCLS,OUTTP)
! an input string of given group and occupation type are given
! apply an string of elementary operators to this group and
! obtain group and type of new string
!
! Jeppe Olsen, October 1993
!
! ------
! Input
! ------
!
! INGRP : Group of input string
! INTP  : Type of input string (occupation type, # e in ras1,ras3)
! IACOP(I) = 1 : operator I is an annihilation operator
!          = 2 : operator I is a  creation   operator
! ITPOP(I) : orbitals space of operator I
! NOP : Number of operators
!
! Output
! ------
! OUTCLS : group of resulting string
! OUTTP  : Type of resulting string

use Definitions, only: iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
integer(kind=iwp), intent(in) :: INCLS, INTP, NOP, IACOP(NOP), ITPOP(NOP)
integer(kind=iwp), intent(out) :: OUTCLS, OUTTP
integer(kind=iwp) :: IDELTA, IEL1, IEL3, IOP

! Number of electrons in RAS1,RAS3
call GTSTTP_2(INCLS,IEL1,IEL3,INTP)

IDELTA = 0
do IOP=1,NOP
  ! Change in number of orbitals
  if (IACOP(IOP) == 1) then
    IDELTA = IDELTA-1
  else
    IDELTA = IDELTA+1
  end if
  ! Change in RAS1, RAS3
  if (ITPOP(IOP) == 1) then
    if (IACOP(IOP) == 1) then
      IEL1 = IEL1-1
    else
      IEL1 = IEL1+1
    end if
  else if (ITPOP(IOP) == 3) then
    if (IACOP(IOP) == 1) then
      IEL3 = IEL3-1
    else
      IEL3 = IEL3+1
    end if
  end if
end do
! Out class
OUTCLS = INCLS-IDELTA
!write(u6,*) ' OUTCLS,IEL1,IEL3 ',OUTCLS,IEL1,IEL3
! out type
call GTSTTP_1(OUTCLS,IEL1,IEL3,OUTTP)

#ifdef _DEBUGPRINT_
write(u6,*) ' NEWTYP, OUTCLS, OUTTP ',OUTCLS,OUTTP
#endif

return

end subroutine NEWTYP_MCLR
