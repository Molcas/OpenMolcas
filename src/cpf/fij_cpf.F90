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
!               1986, Margareta R. A. Blomberg                         *
!***********************************************************************

subroutine FIJ_CPF(ICASE,JSY,INDX,C,S,FC,A,B,FK,DBK,ENP,EPP)

use cpf_global, only: IDENS, IRC, IREF0, IROW, ITER, IV0, LSYM, Lu_25, Lu_CIGuga, NNS, NORBT, NVIR
use guga_util_global, only: COP, IAD10, ICOP1, nCOP
use Symmetry_Info, only: Mul
use Definitions, only: wp, iwp

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: ICASE(*), JSY(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*), FC(*), FK(*), EPP(*)
real(kind=wp), intent(_OUT_) :: A(*), B(*), DBK(*)
real(kind=wp), intent(in) :: ENP(*)
integer(kind=iwp) :: IADD10, IADD25, IC1, IC2, ICHK, IIN, IK, ILEN, IND, INDA, INDB, INDI, INUM, IVL, NA, NB, NI, NK, NOB2, NS1, &
                     NS1L
real(kind=wp) :: COPI, TERM
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

IK = 0 ! dummy initialize
NOB2 = IROW(NORBT+1)
!if (IDENS == 1) write(u6,876) (FC(I),I=1,NOB2)
ICHK = 0
if (IDENS /= 1) then
  NOB2 = IROW(NORBT+1)
  IADD25 = 0
  call dDAFILE(Lu_25,2,FC,NOB2,IADD25)
end if
if ((IDENS == 1) .or. (ITER /= 1)) then
  IADD10 = IAD10(8)
  do
    call dDAFILE(Lu_CIGuga,2,COP,nCOP,IADD10)
    call iDAFILE(Lu_CIGuga,2,iCOP1,nCOP+1,IADD10)
    ILEN = ICOP1(nCOP+1)
    if (ILEN == 0) cycle
    if (ILEN < 0) exit
    do IIN=1,ILEN
      IND = ICOP1(IIN)
      if (ICHK == 0) then
        if (IND /= 0) then
          IVL = ibits(IND,0,6)
          IC2 = ibits(IND,6,13)
          IC1 = ibits(IND,19,13)
          COPI = COP(IIN)*FC(IK)
          if (IVL == IV0) then
            if (IC1 == IREF0) then
              if (IDENS /= 1) then
                COPI = COPI/sqrt(ENP(IC2))
                S(IC2) = S(IC2)+COPI
                if (ITER /= 1) EPP(IC2) = EPP(IC2)+COPI*C(IC2)
              else
                FC(IK) = FC(IK)+COP(IIN)*C(IC1)*C(IC2)/ENP(IC2)
              end if
            else if (IC2 == IREF0) then
              if (IDENS /= 1) then
                COPI = COPI/sqrt(ENP(IC1))
                S(IC1) = S(IC1)+COPI
                if (ITER /= 1) EPP(IC1) = EPP(IC1)+COPI*C(IC1)
              else
                FC(IK) = FC(IK)+COP(IIN)*C(IC1)*C(IC2)/ENP(IC1)
              end if
            else if (IDENS /= 1) then
              S(IC1) = S(IC1)+COPI*C(IC2)
              S(IC2) = S(IC2)+COPI*C(IC1)
            else
              FC(IK) = FC(IK)+COP(IIN)*C(IC1)*C(IC2)/(sqrt(ENP(IC1))*sqrt(ENP(IC2)))
            end if
          else
            INDA = IRC(IVL)+IC1
            INDB = IRC(IVL)+IC2
            NA = INDX(INDA)
            NB = INDX(INDB)
            NS1 = JSUNP(JSY,INDA)
            NS1L = MUL(NS1,LSYM)
            INUM = NVIR(NS1L)
            if (IVL >= 2) INUM = NNS(NS1L)
            if (IDENS /= 1) then
              S(NA+1:NA+INUM) = S(NA+1:NA+INUM)+COPI*C(NB+1:NB+INUM)
              S(NB+1:NB+INUM) = S(NB+1:NB+INUM)+COPI*C(NA+1:NA+INUM)
            else
              TERM = DDOT_(INUM,C(NA+1),1,C(NB+1),1)
              FC(IK) = FC(IK)+COP(IIN)*TERM/(sqrt(ENP(INDA))*sqrt(ENP(INDB)))
            end if
          end if
        else
          ICHK = 1
        end if
      else
        ICHK = 0
        INDI = IND
        NI = ibits(INDI,0,10)
        NK = ibits(INDI,10,10)
        IK = IROW(NK)+NI
      end if
    end do
  end do
end if
!if (IDENS == 1) write(u6,876) (FC(I),I=1,NOB2)
call AI_CPF(JSY,INDX,C,S,FC,C,A,B,FK,DBK,ENP,EPP,0)
if (ITER /= 1) call AB_CPF(ICASE,JSY,INDX,C,S,FC,A,B,FK,ENP)

return

!876 format(1X,'FIJ',5F12.6)

end subroutine FIJ_CPF
