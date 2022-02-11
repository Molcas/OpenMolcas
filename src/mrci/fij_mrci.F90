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

subroutine FIJ_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK,DBK)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: ICSPCK(*), INTSYM(*), INDX(*)
real(kind=wp) :: C(*), S(*), FC(*), A(*), B(*), FK(*), DBK(*)
#include "mrci.fh"
integer(kind=iwp) :: IADD25, IC1, IC2, ICHK, IIN, IK, ILEN, IND, INDA, INDB, INDI, INUM, IVL, NA, NB, NI, NK, NS1, NS1L
real(kind=wp) :: COPI
integer(kind=iwp), external :: JSUNP
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

ICHK = 0
IK = 0
IADD25 = 0
call dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
IADD10 = IAD10(8)
100 call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN == 0) GO TO 100
if (ILEN < 0) GO TO 200
do IIN=1,ILEN
  IND = ICOP1(IIN)
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
  COPI = COP(IIN)*FC(IK)
  if (IVL == IVVER) then
    S(IC1) = S(IC1)+COPI*C(IC2)
    S(IC2) = S(IC2)+COPI*C(IC1)
    GO TO 10
  end if
  INDA = IRC(IVL)+IC1
  INDB = IRC(IVL)+IC2
  NA = INDX(INDA)
  NB = INDX(INDB)
  NS1 = JSYM(INDA)
  NS1L = MUL(NS1,LSYM)
  INUM = NVIR(NS1L)
  if (IVL >= 2) INUM = NVPAIR(NS1L)
  call DAXPY_(INUM,COPI,C(NB+1),1,S(NA+1),1)
  call DAXPY_(INUM,COPI,C(NA+1),1,S(NB+1),1)
10 continue
end do
GO TO 100
200 if (ITER == 0) return
call AI_MRCI(INTSYM,INDX,C,S,FC,A,B,FK,DBK,0)
if ((ITER == 1) .and. (IREST == 0)) return
call AB_MRCI(ICSPCK,INTSYM,INDX,C,S,FC,A,B,FK)

return

end subroutine FIJ_MRCI
