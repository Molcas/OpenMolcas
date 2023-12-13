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

use mrci_global, only: IRC, IVVER, LSYM, LUSYMB, NBAST, NVIR, NVPAIR
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: INTSYM(*), INDX(*)
real(kind=wp), intent(in) :: C1(*), C2(*)
real(kind=wp), intent(inout) :: TDMO(NBAST,NBAST)
integer(kind=iwp) :: IADD10, IC1, IC2, ICHK, IIN, ILEN, IND, INDA, INDB, INDI, INUM, IVL, NA, NB, NI, NK, NS1, NS1L
real(kind=wp) :: TERM
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

!------
! POW: Unnecessary but warning stopping initializations
NI = -1234567
NK = -1234567
!------
ICHK = 0
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
    else if (IND == 0) then
      ICHK = 1
    else
      IVL = ibits(IND,0,6)
      IC2 = ibits(IND,6,13)
      IC1 = ibits(IND,19,13)
      if (IVL == IVVER) then
        TDMO(NI,NK) = TDMO(NI,NK)+COP(IIN)*C1(IC1)*C2(IC2)
        if (NI /= NK) TDMO(NK,NI) = TDMO(NK,NI)+COP(IIN)*C2(IC1)*C1(IC2)
      else
        INDA = IRC(IVL)+IC1
        INDB = IRC(IVL)+IC2
        NA = INDX(INDA)
        NB = INDX(INDB)
        NS1 = JSUNP(INTSYM,INDA)
        NS1L = MUL(NS1,LSYM)
        INUM = NVIR(NS1L)
        if (IVL >= 2) INUM = NVPAIR(NS1L)
        TERM = DDOT_(INUM,C1(NA+1),1,C2(NB+1),1)
        TDMO(NI,NK) = TDMO(NI,NK)+COP(IIN)*TERM
        if (NI /= NK) then
          TERM = DDOT_(INUM,C2(NA+1),1,C1(NB+1),1)
          TDMO(NK,NI) = TDMO(NK,NI)+COP(IIN)*TERM
        end if
      end if
    end if
  end do
end do

return

end subroutine FIJTD
