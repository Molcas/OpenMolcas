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

subroutine MPSQ2(C,S,W,MUL,INDX,JSY,NDIAG,INUM,IRC3,LSYM,NVIRT,SQ2)

use Definitions, only: wp, iwp

implicit none
real(kind=wp) :: C(*), S(*), W(*), SQ2
integer(kind=iwp) :: MUL(8,8), INDX(*), JSY(*), NDIAG(*), INUM, IRC3, LSYM, NVIRT
integer(kind=iwp) :: I, II1, MA, NA, NS1, NS1L
integer(kind=iwp), external :: JSUNP_CPF
! Statement function
!PAM97      INTEGER UNPACK
!PAM97      EXTERNAL UNPACK
!RL   JSYM(L)=IAND(ISHFT(JSY((L+19)/20),-3*((L+19)/20*20-L)),7)+1
!PAM96      JSYM(L)=UNPACK(JSY((L+9)/10),3*MOD(L-1,10)+1,3)+1
integer(kind=iwp) :: JSYM, L
JSYM(L) = JSUNP_CPF(JSY,L)

do I=1,INUM
  II1 = IRC3+I
  NS1 = JSYM(II1)
  NS1L = MUL(NS1,LSYM)
  if (NS1L /= 1) GO TO 10
  NA = INDX(II1)
  do MA=1,NVIRT
    C(NA+NDIAG(MA)) = SQ2*C(NA+NDIAG(MA))
    S(NA+NDIAG(MA)) = S(NA+NDIAG(MA))/SQ2
    W(NA+NDIAG(MA)) = W(NA+NDIAG(MA))/SQ2
  end do
10 continue
end do

return

end subroutine MPSQ2
