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

subroutine CUPDATE(JSY,INDX,C,S,AP,BST,T2,ENP)

use Constants, only: Zero, Two
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: JSY(*), INDX(*)
real(kind=wp) :: C(*), S(*), AP(*), BST(*), T2(*), ENP(*)
#include "cpfmcpf.fh"
#include "files_cpf.fh"
integer(kind=iwp) :: I, IAD, III, IIN, INUM, IP, IST, JJJ, NS1, NSIL
real(kind=wp) :: A, APW, EMP, W
integer(kind=iwp), external :: JSUNP_CPF
real(kind=r8), external :: DDOT_

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
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  call VSMSB(C(IST),1,AP(IIN+I),S(IST),1,C(IST),1,INUM)
end do

! DOUBLES
IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  call VSMSB(C(IST),1,AP(IIN+I),S(IST),1,C(IST),1,INUM)
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
  T2(I) = S(I)+APW
  C(I) = C(I)/T2(I)
  EMP = sqrt(ENP(I))
  C(I) = C(I)*EMP
  if (I == IREF0) C(I) = Zero
end do

! SINGLES
IP = IRC(2)-IRC(1)
IIN = IRC(1)
do I=1,IP
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NVIR(NSIL)
  IST = INDX(IIN+I)+1
  APW = W-AP(IIN+I)
  call VSADD(S(IST),1,APW,T2,1,INUM)
  call VDIV(T2,1,C(IST),1,C(IST),1,INUM)
  EMP = sqrt(ENP(IIN+I))
  call VSMUL(C(IST),1,EMP,C(IST),1,INUM)
end do

! DOUBLES
IP = IRC(4)-IRC(2)
IIN = IRC(2)
do I=1,IP
  NS1 = JSUNP_CPF(JSY,IIN+I)
  NSIL = MUL(NS1,LSYM)
  INUM = NNS(NSIL)
  IST = INDX(IIN+I)+1
  APW = W-AP(IIN+I)
  call VSADD(S(IST),1,APW,T2,1,INUM)
  call VDIV(T2,1,C(IST),1,C(IST),1,INUM)
  EMP = sqrt(ENP(IIN+I))
  call VSMUL(C(IST),1,EMP,C(IST),1,INUM)
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
