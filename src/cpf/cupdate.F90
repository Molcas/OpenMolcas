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

subroutine CUPDATE(JSY,INDX,C,S,AP,BST,ENP)

use cpf_global, only: IAD25S, IADDP, IPRINT, IRC, IREF0, ITPUL, LSYM, Lu_25, Lu_30, Lu_CI, NCONF, NNS, NVIR, WLEV
use guga_util_global, only: nCOP
use Symmetry_Info, only: Mul
use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6

#include "intent.fh"

implicit none
integer(kind=iwp), intent(in) :: JSY(*), INDX(*)
real(kind=wp), intent(inout) :: C(*), S(*), BST(*)
real(kind=wp), intent(in) :: AP(*), ENP(*)
integer(kind=iwp) :: I, IAD, III, IIN, INUM, IP, IST, JJJ, NS1, NSIL
real(kind=wp) :: A, APW, EMP, W
integer(kind=iwp), external :: JSUNP
real(kind=wp), external :: DDOT_

W = WLEV

! VALENCE
IP = IRC(1)
do I=1,IP
  C(I) = AP(I)*C(I)-S(I)
end do

! SINGLES
IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  C(IST:IST+INUM-1) = C(IST:IST+INUM-1)*AP(IIN+I)-S(IST:IST+INUM-1)
end do

! DOUBLES
IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  C(IST:IST+INUM-1) = C(IST:IST+INUM-1)*AP(IIN+I)-S(IST:IST+INUM-1)
end do

! WRITE GRADIENT ONTO DISK, UNIT=30
IAD = IADDP(ITPUL)
call dDAFILE(Lu_30,1,C,NCONF,IAD)

! Reuse array S for HCOUT that was written in IJIJ.
IAD = IAD25S
do III=1,NCONF,nCOP
  JJJ = min(nCOP,NCONF+1-III)
  call dDAFILE(Lu_25,2,S(III),JJJ,IAD)
end do

! VALENCE
IP = IRC(1)
do I=1,IP
  APW = W-AP(I)
  C(I) = C(I)/(S(I)+APW)
  EMP = sqrt(ENP(I))
  C(I) = C(I)*EMP
  if (I == IREF0) C(I) = Zero
end do

! SINGLES
IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  APW = W-AP(IIN+I)
  C(IST:IST+INUM-1) = C(IST:IST+INUM-1)/(S(IST:IST+INUM-1)+APW)
  EMP = sqrt(ENP(IIN+I))
  C(IST:IST+INUM-1) = EMP*C(IST:IST+INUM-1)
end do

! DOUBLES
IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  APW = W-AP(IIN+I)
  C(IST:IST+INUM-1) = C(IST:IST+INUM-1)/(S(IST:IST+INUM-1)+APW)
  EMP = sqrt(ENP(IIN+I))
  C(IST:IST+INUM-1) = EMP*C(IST:IST+INUM-1)
end do

IAD = IADDP(ITPUL+1)
call dDAFILE(Lu_CI,1,C,NCONF,IAD)
IADDP(ITPUL+2) = IAD
if (IPRINT >= 15) write(u6,999) (C(I),I=1,NCONF)
A = DDOT_(NCONF,C,1,C,1)
if (A > Two) then
  write(u6,*) 'CUPDATE Error: A>2.0 (See code.)'
  call Abend()
end if
if (ITPUL == 1) BST(1) = A

return

999 format(6X,'C(UPD)',5F10.6)

end subroutine CUPDATE
