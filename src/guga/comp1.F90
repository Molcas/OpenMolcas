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
! Copyright (C) 1986, Per E. M. Siegbahn                               *
!***********************************************************************

subroutine COMP1(LJ,ITYP,L,IT2,II,IID,JJ,JJD,JTYP,ITAI)

use guga_global, only: COUP, IADD10, ICASE, ICOUP, ICOUP1, IOUT, IX, JNDX, LN, Lu_10, NBUF, NMAT
use guga_util_global, only: COP, ICOP1, nCOP
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LJ, ITYP, L, IT2, II, IID, JJ, JJD, JTYP, ITAI(*)
integer(kind=iwp) :: IC1, IC2, ICT, IN_, IN2, IND, ITAIL, JND1, JND2, JOJ, KK1, KTYP
real(kind=wp) :: FAC
integer(kind=iwp), external :: ICUNP

FAC = One
ITAIL = IX(IT2+LJ)
do IN_=1,ITAIL
  IC1 = ICOUP(1)+IN_
  JND1 = JNDX(II+IC1)
  if (JND1 == 0) cycle
  IN2 = ITAI(IN_)
  if (IN2 == 0) cycle
  IC2 = ICOUP1(1)+IN2
  if (ITYP == 1) then
    KK1 = (JND1-1)*LN+L
    JOJ = ICUNP(ICASE,KK1)
    if (JOJ > 1) JOJ = JOJ-1
    FAC = JOJ
    if (JOJ == 0) cycle
  end if
  IC1 = JND1-IID
  JND2 = JNDX(JJ+IC2)
  if (JND2 == 0) cycle
  IC2 = JND2-JJD
  IOUT = IOUT+1
  KTYP = JTYP
  if (JTYP > 3) then
    ICT = IC1
    IC1 = IC2
    IC2 = ICT
    KTYP = JTYP-3
  end if
  !IND = KTYP+2**6*IC2
  IND = ior(KTYP,ishft(IC2,6))
  !ICOP1(IOUT) = IND+2**19*IC1
  ICOP1(IOUT) = ior(IND,ishft(IC1,19))
  COP(IOUT) = FAC*COUP(1)
  if (IOUT < NBUF) cycle
  ICOP1(nCOP+1) = NBUF
  !write(u6,*) 'WRITING BUFFER IN COMP1'
  !write(u6,*) '======================='
  !write(u6,'(A,I6)') 'ICOP1(nCOP+1) ',ICOP1(nCOP+1)
  !write(u6,'(A,I6)') 'IADD10     ',IADD10
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
end do

return

end subroutine COMP1
