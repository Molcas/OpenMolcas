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

use Constants, only: One
use Definitions, only: wp, iwp, r8

implicit none
integer(kind=iwp) :: INTSYM(*), INDX(*), JREFX(*)
real(kind=wp) :: C(*), DMO(*), AREF(*)
#include "mrci.fh"
integer(kind=iwp) :: IC1, IC2, ICHK, IIN, IK, ILEN, IND, INDA, INDB, INDI, INUM, IRC1, IRC2, IVL, NA, NB, NI, NK, NS1, NS1L
real(kind=wp) :: ENPINV, TERM
integer(kind=iwp), external :: JSUNP
real(kind=r8), external :: DDOT_
!Statement function
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP(INTSYM,L)

ICHK = 0
IK = 0
ENPINV = One/ENP
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
      !NI = mod(INDI,2**10)
      !NK = mod(INDI/2**10,2**10)
      NI = ibits(INDI,0,10)
      NK = ibits(INDI,10,10)
      IK = IROW(NK)+NI
    else if (IND == 0) then
      ICHK = 1
    else
      !IVL = mod(IND,2**6)
      !IC2 = mod(IND/2**6,2**13)
      !IC1 = mod(IND/2**19,2**13)
      IVL = ibits(IND,0,6)
      IC2 = ibits(IND,6,13)
      IC1 = ibits(IND,19,13)
      if (IVL == IVVER) then
        DMO(IK) = DMO(IK)+COP(IIN)*C(IC1)*C(IC2)*ENPINV
        if (ICPF /= 0) then
          IRC1 = JREFX(IC1)
          if (IRC1 /= 0) then
            IRC2 = JREFX(IC2)
            if (IRC2 /= 0) DMO(IK) = DMO(IK)+COP(IIN)*AREF(IRC1)*AREF(IRC2)*(One-ENPINV)
          end if
        end if
      else
        INDA = IRC(IVL)+IC1
        INDB = IRC(IVL)+IC2
        NA = INDX(INDA)
        NB = INDX(INDB)
        NS1 = JSYM(INDA)
        NS1L = MUL(NS1,LSYM)
        INUM = NVIR(NS1L)
        if (IVL >= 2) INUM = NVPAIR(NS1L)
        TERM = DDOT_(INUM,C(NA+1),1,C(NB+1),1)
        DMO(IK) = DMO(IK)+COP(IIN)*TERM*ENPINV
      end if
    end if
  end do
end do

return

end subroutine FIJD
