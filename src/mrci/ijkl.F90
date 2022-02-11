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

subroutine IJKL(INTSYM,INDX,C,S,FIJKL)

use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: INTSYM(*), INDX(*)
real(kind=wp) :: C(*), S(*), FIJKL(*)
#include "mrci.fh"
integer(kind=iwp) :: i, IADR, IC1, IC2, ICHK, IIN, ILEN, IND, INDA, INDB, INDI, INUM, IP, IVL, JP, KP, LENGTH, LP, NA, NB, NIJ, &
                     NIJKL, NKL, NS1, NS1L
real(kind=wp) :: COPI, fini
integer(kind=iwp), external :: JSUNP
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

!------
! POW: Unnecessary but warning stopping initialization
fini = huge(fini)
!------
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
call FZERO(FIJKL,NIJKL)
IADR = LASTAD(1)
201 call dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
call iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
LENGTH = INDSRT(NSRTMX+1)
IADR = INDSRT(NSRTMX+2)
!if (LENGTH > 0) call SCATTER(LENGTH,FIJKL,INDSRT,VALSRT)
do i=1,length
  FIJKL(INDSRT(i)) = VALSRT(i)
end do
if (IADR /= -1) GO TO 201
IADD10 = IAD10(5)
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
    !IP = mod(INDI,2**8)
    !JP = mod(INDI/2**8,2**8)
    !KP = mod(INDI/2**16,2**8)
    !LP = mod(INDI/2**24,2**8)
    IP = ibits(INDI,0,8)
    JP = ibits(INDI,8,8)
    KP = ibits(INDI,16,8)
    LP = ibits(INDI,24,8)
    NIJ = IROW(IP)+JP
    NKL = IROW(KP)+LP
    IND = NIJ*(NIJ-1)/2+NKL
    FINI = FIJKL(IND)
    goto 10
  end if
  if (IND == 0) then
    ICHK = 1
    goto 10
  end if
  !IVL = mod(IND,2**6)
  !IC2 = mod(IND/2**6,2**13)
  !IC1 = mod(IND/2**19,2**13)
  IVL = ibits(IND,0,6)
  IC2 = ibits(IND,6,13)
  IC1 = ibits(IND,19,13)
  COPI = COP(IIN)*FINI
  if (IVL == 0) then
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
200 continue

return

end subroutine IJKL
