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

subroutine DIIS_CPF(DPI,DPJ,BST,MIT,BIJ,ITP,CN)

implicit real*8(A-H,O-Z)
dimension DPI(*), DPJ(*), BST(MIT,MIT), BIJ(ITP,ITP), CN(*), WHS(50)
#include "SysDef.fh"
#include "cpfmcpf.fh"
#include "files_cpf.fh"

if (ITPUL == 1) GO TO 26

ITM = ITPUL-1
do I=1,ITM
  do J=1,ITM
    BIJ(J,I) = BST(J,I)
  end do
end do
do I=1,ITPUL
  BIJ(ITP,I) = -1.0d0
  BIJ(I,ITP) = -1.0d0
end do
BIJ(ITP,ITP) = 0.0d0

do I=1,ITM
  IAD = IADDP(I+1)
  call dDAFILE(Lu_CI,2,DPJ,NCONF,IAD)
  T = DDOT_(NCONF,DPI,1,DPJ,1)
  BIJ(I,ITPUL) = T
  BIJ(ITPUL,I) = T
  BST(I,ITPUL) = T
  BST(ITPUL,I) = T
  if (I == 1) then
    T = DDOT_(NCONF,DPJ,1,DPJ,1)
    BIJ(1,1) = T
    BST(1,1) = T
  end if
end do
BIJ(ITPUL,ITPUL) = DDOT_(NCONF,DPI,1,DPI,1)
BST(ITPUL,ITPUL) = BIJ(ITPUL,ITPUL)

if (IPRINT < 10) GO TO 26
do I=1,ITP
  write(6,16) (BIJ(J,I),J=1,ITP)
  call XFLUSH(6)
16 format(6X,'BIJ ',6F12.6)
end do

26 if (IDIIS == 1) GO TO 25
do I=1,ITPUL
  IAD = IADDP(I)
  call dDAFILE(Lu_CI,2,DPJ,NCONF,IAD)
  do J=1,NCONF
    DPI(J) = DPI(J)+DPJ(J)
  end do
end do
if (IPRINT >= 15) write(6,14) (DPI(I),I=1,NCONF)
14 format(6X,'C(DIIS)',5F10.6)
return

25 call DECOMP(ITP,BIJ,BIJ)
do I=1,ITPUL
  WHS(I) = 0.0d00
end do
WHS(ITP) = -1.0d00
call SOLVE(ITP,BIJ,WHS,CN)

! UPDATE P AND DELTA P

call NEXT(DPI,DPJ,CN)

ITPUL = 0

return

end subroutine DIIS_CPF
