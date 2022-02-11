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

subroutine FIJTD(INTSYM,INDX,C1,C2,TDMO)

use Definitions, only: wp, iwp, r8

implicit none
#include "mrci.fh"
integer(kind=iwp) :: INTSYM(*), INDX(*)
real(kind=wp) :: C1(*), C2(*), TDMO(NBAST,NBAST)
integer(kind=iwp) :: IC1, IC2, ICHK, IIN, ILEN, IND, INDA, INDB, INDI, INUM, IVL, NA, NB, ni, nk, NS1, NS1L
real(kind=wp) :: TERM
integer(kind=iwp), external :: JSUNP
real(kind=r8), external :: DDOT_
! Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

!------
! POW: Unnecessary but warning stopping initializations
ni = -1234567
nk = -1234567
!------
ICHK = 0
IADD10 = IAD10(8)
100 call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
ILEN = ICOP1(nCOP+1)
if (ILEN == 0) GO TO 100
if (ILEN < 0) return
do IIN=1,ILEN
  IND = ICOP1(IIN)
  if (ICHK /= 0) then
    ICHK = 0
    INDI = IND
    !NI = mod(INDI,2**10)
    !NK = mod(INDI/2**10,2**10)
    NI = ibits(INDI,0,10)
    NK = ibits(INDI,10,10)
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
  TDMO(NI,NK) = TDMO(NI,NK)+COP(IIN)*C1(IC1)*C2(IC2)
  if (NI /= NK) TDMO(NK,NI) = TDMO(NK,NI)+COP(IIN)*C2(IC1)*C1(IC2)
  GO TO 10
13 INDA = IRC(IVL)+IC1
  INDB = IRC(IVL)+IC2
  NA = INDX(INDA)
  NB = INDX(INDB)
  NS1 = JSYM(INDA)
  NS1L = MUL(NS1,LSYM)
  INUM = NVIR(NS1L)
  if (IVL >= 2) INUM = NVPAIR(NS1L)
  TERM = DDOT_(INUM,C1(NA+1),1,C2(NB+1),1)
  TDMO(NI,NK) = TDMO(NI,NK)+COP(IIN)*TERM
  if (NI == NK) goto 10
  TERM = DDOT_(INUM,C2(NA+1),1,C1(NB+1),1)
  TDMO(NK,NI) = TDMO(NK,NI)+COP(IIN)*TERM
10 continue
end do
GO TO 100

end subroutine FIJTD
