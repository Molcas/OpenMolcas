!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

!#define _DEBUGPRINT_
subroutine NUMST4_MCLR(NEL,NORB1,NEL1MN,NEL1MX,NORB2,NORB3,NEL3MN,NEL3MX,NSTTP)
! Number of strings per type for class of strings with
!
! Between NEL1MN AND NEL1MX electrons in the first NORB1 orbitals
! Between NEL3MN AND NEL3MX electrons in the last  NORB3 orbitals

integer NSTTP(*)

NSTRIN = 0

#ifdef _DEBUGPRINT_
write(6,*) ' NUMST4'
write(6,*) ' NEL NEL1MN NEL1MX NEL3MN NEL3MX'
write(6,*) NEL,NEL1MN,NEL1MX,NEL3MN,NEL3MX
#endif

ITPMAX = 0
do IEL1=NEL1MN,min(NEL1MX,NORB1,NEL)
  NSTIN1 = IBION(NORB1,IEL1)
  IEL3MN = max(NEL3MN,NEL-(IEL1+NORB2))
  IEL3MX = min(NEL3MX,NEL-IEL1)
  do IEL3=IEL3MN,IEL3MX
    IEL2 = NEL-IEL1-IEL3
    NSTINT = NSTIN1*IBION(NORB2,IEL2)*IBION(NORB3,IEL3)
    NSTRIN = NSTRIN+NSTINT
    ITP = (NEL1MX-IEL1)*(NEL3MX-NEL3MN+1)+IEL3-NEL3MN+1
    NSTTP(ITP) = NSTINT
    ITPMAX = max(ITPMAX,ITP)
  end do
end do

#ifdef _DEBUGPRINT_
write(6,'(/A,I6)') '  Number of strings generated ... ',NSTRIN
write(6,*)
write(6,*) ' Largest type number ',ITPMAX
write(6,*) ' Number of strings per type'
call IWRTMA(NSTTP,1,ITPMAX,1,ITPMAX)
#endif

return

end subroutine NUMST4_MCLR
