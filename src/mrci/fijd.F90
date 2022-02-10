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

subroutine FIJD(INTSYM,INDX,C,DMO,JREFX,AREF)

implicit real*8(A-H,O-Z)
#include "SysDef.fh"
#include "mrci.fh"
dimension INTSYM(*), INDX(*), C(*), DMO(*), JREFX(*), AREF(*)
!Statement function
JSYM(L) = JSUNP(INTSYM,L)

ICHK = 0
IK = 0
ENPINV = 1.0d00/ENP
IADD10 = IAD10(8)
100 call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
LEN = ICOP1(nCOP+1)
if (LEN == 0) GO TO 100
if (LEN < 0) return
do IN=1,LEN
  IND = ICOP1(IN)
  if (ICHK /= 0) then
    ICHK = 0
    INDI = IND
    !NI = mod(INDI,2**10)
    !NK = mod(INDI/2**10,2**10)
    NI = ibits(INDI,0,10)
    NK = ibits(INDI,10,10)
    IK = IROW(NK)+NI
    GO TO 10
  end if
  if (IND == 0) then
    ICHK = 1
    GO TO 10
  end if
  !IVL = mod(IND,2**6)
  !IC2 = mod(IND/2**6,2**13)
  !IC1 = mod(IND/2**19,2**13)
  IVL = ibits(IND,0,6)
  IC2 = ibits(IND,6,13)
  IC1 = ibits(IND,19,13)
  if (IVL /= IVVER) GO TO 13
  DMO(IK) = DMO(IK)+COP(IN)*C(IC1)*C(IC2)*ENPINV
  if (ICPF == 0) GO TO 10
  IRC1 = JREFX(IC1)
  if (IRC1 == 0) GO TO 10
  IRC2 = JREFX(IC2)
  if (IRC2 == 0) GO TO 10
  DMO(IK) = DMO(IK)+COP(IN)*AREF(IRC1)*AREF(IRC2)*(1.0d00-ENPINV)
  GO TO 10
13 INDA = IRC(IVL)+IC1
  INDB = IRC(IVL)+IC2
  NA = INDX(INDA)
  NB = INDX(INDB)
  NS1 = JSYM(INDA)
  NS1L = MUL(NS1,LSYM)
  INUM = NVIR(NS1L)
  if (IVL >= 2) INUM = NVPAIR(NS1L)
  TERM = DDOT_(INUM,C(NA+1),1,C(NB+1),1)
  DMO(IK) = DMO(IK)+COP(IN)*TERM*ENPINV
10 continue
end do
GO TO 100

end subroutine FIJD
