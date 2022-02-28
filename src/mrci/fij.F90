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

subroutine FIJ(INTSYM,INDX,C,S,FC,A,B,FK,DBK)

use mrci_global, only: IRC, IREST, IROW, ITER, IVVER, LSYM, Lu_25, LUSYMB, NBTRI, NVIR, NVPAIR
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*), FC(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*), FK(*), DBK(*)
integer(kind=iwp) :: IADD10, IADD25, IC1, IC2, ICHK, IIN, IK, ILEN, IND, INDA, INDB, INDI, INUM, IVL, NA, NB, NI, NK, NS1, NS1L
real(kind=wp) :: COPI
integer(kind=iwp), external :: JSUNP

ICHK = 0
IK = 0
IADD25 = 0
call dDAFILE(Lu_25,2,FC,NBTRI,IADD25)
IADD10 = IAD10(8)
do
  call dDAFILE(LUSYMB,2,COP,nCOP,IADD10)
  call iDAFILE(LUSYMB,2,iCOP1,nCOP+1,IADD10)
  ILEN = ICOP1(nCOP+1)
  if (ILEN < 0) exit
  do IIN=1,ILEN
    IND = ICOP1(IIN)
    if (ICHK /= 0) then
      ICHK = 0
      INDI = IND
      NI = ibits(INDI,0,10)
      NK = ibits(INDI,10,10)
      IK = IROW(NK)+NI
    else if (IND == 0) then
      ICHK = 1
    else
      IVL = ibits(IND,0,6)
      IC2 = ibits(IND,6,13)
      IC1 = ibits(IND,19,13)
      COPI = COP(IIN)*FC(IK)
      if (IVL == IVVER) then
        S(IC1) = S(IC1)+COPI*C(IC2)
        S(IC2) = S(IC2)+COPI*C(IC1)
      else
        INDA = IRC(IVL)+IC1
        INDB = IRC(IVL)+IC2
        NA = INDX(INDA)
        NB = INDX(INDB)
        NS1 = JSUNP(INTSYM,INDA)
        NS1L = MUL(NS1,LSYM)
        INUM = NVIR(NS1L)
        if (IVL >= 2) INUM = NVPAIR(NS1L)
        S(NA+1:NA+INUM) = S(NA+1:NA+INUM)+COPI*C(NB+1:NB+INUM)
        S(NB+1:NB+INUM) = S(NB+1:NB+INUM)+COPI*C(NA+1:NA+INUM)
      end if
    end if
  end do
end do
if (ITER /= 0) then
  call AI_MRCI(INTSYM,INDX,C,S,FC,A,B,FK,DBK,0)
  if ((ITER /= 1) .or. (IREST /= 0)) call AB(INTSYM,INDX,C,S,FC,A,B,FK)
end if

return

end subroutine FIJ
