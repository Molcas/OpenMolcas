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
!               2021, Ignacio Fdez. Galvan                             *
!***********************************************************************
! 2021: Remove GOTOs

subroutine COMP(I,LJ,ITYP,L,IT1,IT2)

use guga_global, only: COUP, IADD10, ICASE, ICOUP, ICOUP1, IOUT, IRC, IV0, IWAY, IX, J2, JNDX, JRC, LN, Lu_10, NBUF, NMAT
use guga_util_global, only: COP, ICOP1, nCOP
use Constants, only: One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: I, LJ, ITYP, L, IT1, IT2
integer(kind=iwp) :: IC1, IC2, II1, IN_, IND, ISTOP, ITAIL, IVL, IVL0, IVV, JJ, JJD, JND1, JND2, JOJ, KM
real(kind=wp) :: FAC
logical(kind=iwp) :: first
integer(kind=iwp), external :: ICUNP

if (IT1 /= IT2) then
  write(u6,*) 'Comp: IT1 /= IT2'
  write(u6,*) 'IT1,IT2=',IT1,IT2
  write(u6,*) 'ITYP,L=',ITYP,L
  call Abend()
end if
FAC = One
KM = I
first = .true.
do
  if (first) then
    KM = KM-1
    if (KM /= 0) IWAY(KM) = 1
    first = .false.
  end if
  if (KM /= 0) then
    call PATH(KM,ISTOP,IT1,IT2)
    if (ISTOP == 0) then
      first = .true.
      cycle
    end if
    KM = KM+1
    if (KM == I) exit
  else
    IVL = J2(1)
    ITAIL = IX(IT1+LJ)
    IVV = IV0-IVL
    if (IVV == 0) then
      JJ = 0
      JJD = 0
    else
      JJ = IRC(IVV)
      JJD = JRC(IVV)
    end if
    do IN_=1,ITAIL
      IC1 = ICOUP(1)+IN_
      IC2 = ICOUP1(1)+IN_
      JND1 = JNDX(JJ+IC1)
      if (JND1 == 0) cycle
      if (ITYP == 1) then
        II1 = (JND1-1)*LN+L
        JOJ = ICUNP(ICASE,II1)
        if (JOJ > 1) JOJ = JOJ-1
        FAC = JOJ
        if (JOJ == 0) cycle
      end if
      IC1 = JND1-JJD
      JND2 = JNDX(JJ+IC2)
      if (JND2 == 0) cycle
      IC2 = JND2-JJD
      IOUT = IOUT+1
      IVL0 = IV0-IVL
      !IND = IVL0+2**6*IC2
      !ICOP1(IOUT) = IND+2**19*IC1
      IND = ior(IVL0,ishft(IC2,6))
      ICOP1(IOUT) = ior(IND,ishft(IC1,19))

      COP(IOUT) = FAC*COUP(I)
      if (IOUT < NBUF) cycle
      ICOP1(nCOP+1) = NBUF
      call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
      call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
      NMAT = NMAT+NBUF
      IOUT = 0
    end do
    if (I == 1) exit
    KM = 1
  end if
end do

return

end subroutine COMP
