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

use mrci_global, only: INDSRT, IRC, IROW, LASTAD, Lu_70, LUSYMB, LN, LSYM, NSRTMX, NVIR, NVPAIR, VALSRT
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(in) :: C(*)
real(kind=wp), intent(inout) :: S(*)
real(kind=wp), intent(_OUT_) :: FIJKL(*)
integer(kind=iwp) :: i, IADD10, IADR, IC1, IC2, ICHK, IIN, ILEN, IND, INDA, INDB, INDI, INUM, IP, IVL, JP, KP, LENGTH, LP, NA, NB, &
                     NIJ, NIJKL, NKL, NS1, NS1L
real(kind=wp) :: COPI, FINI
integer(kind=iwp), external :: JSUNP

!------
! POW: Unnecessary but warning stopping initialization
FINI = huge(FINI)
!------
ICHK = 0
NIJ = IROW(LN+1)
NIJKL = NIJ*(NIJ+1)/2
FIJKL(1:NIJKL) = Zero
IADR = LASTAD(1)
do
  call dDAFILE(Lu_70,2,VALSRT,NSRTMX,IADR)
  call iDAFILE(Lu_70,2,INDSRT,NSRTMX+2,IADR)
  LENGTH = INDSRT(NSRTMX+1)
  IADR = INDSRT(NSRTMX+2)
  do i=1,length
    FIJKL(INDSRT(i)) = VALSRT(i)
  end do
  if (IADR == -1) exit
end do
IADD10 = IAD10(5)
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
      IP = ibits(INDI,0,8)
      JP = ibits(INDI,8,8)
      KP = ibits(INDI,16,8)
      LP = ibits(INDI,24,8)
      NIJ = IROW(IP)+JP
      NKL = IROW(KP)+LP
      IND = NIJ*(NIJ-1)/2+NKL
      FINI = FIJKL(IND)
    else if (IND == 0) then
      ICHK = 1
    else
      IVL = ibits(IND,0,6)
      IC2 = ibits(IND,6,13)
      IC1 = ibits(IND,19,13)
      COPI = COP(IIN)*FINI
      if (IVL == 0) then
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

return

end subroutine IJKL
