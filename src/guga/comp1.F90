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

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LJ, ITYP, L, IT2, II, IID, JJ, JJD, JTYP, ITAI(*)
#include "SysDef.fh"
#include "real_guga.fh"
#include "integ.fh"
#include "files_addr.fh"
#include "d.fh"
integer(kind=iwp) :: IC1, IC2, ICT, IN_, IN2, IND, ITAIL, JND1, JND2, JOJ, KK1, KTYP
real(kind=wp) :: FAC
integer(kind=iwp), external :: ICUNP
! statement function
integer(kind=iwp) :: JO, I
JO(I) = ICUNP(ICASE,I)

FAC = D1
ITAIL = IX(IT2+LJ)
do IN_=1,ITAIL
  IC1 = ICOUP(1)+IN_
  JND1 = JNDX(II+IC1)
  if (JND1 == 0) GO TO 90
  IN2 = ITAI(IN_)
  if (IN2 == 0) GO TO 90
  IC2 = ICOUP1(1)+IN2
  if (ITYP /= 1) GO TO 91
  KK1 = (JND1-1)*LN+L
  JOJ = JO(KK1)
  if (JOJ > 1) JOJ = JOJ-1
  FAC = JOJ
  if (JOJ == 0) GO TO 90
91 IC1 = JND1-IID
  JND2 = JNDX(JJ+IC2)
  if (JND2 == 0) GO TO 90
  IC2 = JND2-JJD
  IOUT = IOUT+1
  KTYP = JTYP
  if (JTYP <= 3) GO TO 92
  ICT = IC1
  IC1 = IC2
  IC2 = ICT
  KTYP = JTYP-3
92 continue
  !IND = KTYP+2**6*IC2
  IND = ior(KTYP,ishft(IC2,6))
  !ICOP1(IOUT) = IND+2**19*IC1
  ICOP1(IOUT) = ior(IND,ishft(IC1,19))
  COP(IOUT) = FAC*COUP(1)
  if (IOUT < NBUF) GO TO 90
  ICOP1(nCOP+1) = NBUF
  !write(u6,*) 'WRITING BUFFER IN COMP1'
  !write(u6,*) '======================='
  !write(u6,'(A,I6)') 'ICOP1(nCOP+1) ',ICOP1(nCOP+1)
  !write(u6,'(A,I6)') 'IADD10     ',IADD10
  call dDAFILE(Lu_10,1,COP,NCOP,IADD10)
  call iDAFILE(Lu_10,1,iCOP1,NCOP+1,IADD10)
  NMAT = NMAT+NBUF
  IOUT = 0
90 continue
end do

return

end subroutine COMP1
