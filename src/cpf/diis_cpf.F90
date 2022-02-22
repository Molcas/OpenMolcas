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

use cpf_global, only: IADDP, IDIIS, IPRINT, ITPUL, Lu_CI, NCONF
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6, r8

implicit none
integer(kind=iwp) :: MIT, ITP
real(kind=wp) :: DPI(*), DPJ(*), BST(MIT,MIT), BIJ(ITP,ITP), CN(*)
integer(kind=iwp) :: I, IAD, ITM, J
real(kind=wp) :: T, WHS(50)
real(kind=r8), external :: DDOT_

if (ITPUL /= 1) then

  ITM = ITPUL-1
  do I=1,ITM
    do J=1,ITM
      BIJ(J,I) = BST(J,I)
    end do
  end do
  do I=1,ITPUL
    BIJ(ITP,I) = -One
    BIJ(I,ITP) = -One
  end do
  BIJ(ITP,ITP) = Zero

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

  if (IPRINT >= 10) then
    do I=1,ITP
      write(u6,16) (BIJ(J,I),J=1,ITP)
    end do
  end if

end if

if (IDIIS /= 1) then
  do I=1,ITPUL
    IAD = IADDP(I)
    call dDAFILE(Lu_CI,2,DPJ,NCONF,IAD)
    do J=1,NCONF
      DPI(J) = DPI(J)+DPJ(J)
    end do
  end do
  if (IPRINT >= 15) write(u6,14) (DPI(I),I=1,NCONF)
else
  call DECOMP(ITP,BIJ,BIJ)
  do I=1,ITPUL
    WHS(I) = Zero
  end do
  WHS(ITP) = -One
  call SOLVE(ITP,BIJ,WHS,CN)

  ! UPDATE P AND DELTA P

  call NEXT(DPI,DPJ,CN)

  ITPUL = 0
end if

return

14 format(6X,'C(DIIS)',5F10.6)
16 format(6X,'BIJ ',6F12.6)

end subroutine DIIS_CPF
