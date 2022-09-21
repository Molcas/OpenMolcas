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

subroutine decomp2ind(W,IDM,no,NF)

use Constants, only: Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: IDM, no, NF
real(kind=wp) :: W(IDM,*)
integer(kind=iwp) :: IJ, I, IJO, II, K, L, KL, LK, J, JIO
real(kind=wp) :: DD

!mp if (LLtrace) then
!mp   write(u6,*) 'Entering DECOMP2IND'
!mp   call xflush(u6)
!mp endif

! symmetrizes the upper index

!write(u6,*) 'decomp2ind:',IDM,no,NF
do I=1,no
  II = I*(I+1)/2
  do K=2,NF
    do L=1,K-1
      KL = (K-1)*NF+L
      LK = (L-1)*NF+K
      DD = Half*(W(KL,II)+W(LK,II))
      W(KL,II) = DD
      W(LK,II) = DD
    end do
  end do
end do
if (NO > 2) then
  do I=NO,2,-1
    IJ = I*(I-1)/2+1
    IJO = (I-1)*no+1
    call DCOPY_(IDM*I,W(1,IJ),1,W(1,IJO),1)
  end do
else if (NO == 2) then
  call DCOPY_(IDM,W(1,3),1,W(1,4),1)
  call DCOPY_(IDM,W(1,2),1,W(1,3),1)
end if
do I=2,no
  do J=1,I-1
    IJO = (I-1)*no+J
    JIO = (J-1)*no+I
    call map2_21_t3(W(1,IJO),W(1,JIO),NF,NF)
  end do
end do
!stop
!mp if (LLtrace) then
!mp   write(u6,*) 'Leaving DECOMP2IND'
!mp   call xflush(u6)
!mp end if

end subroutine decomp2ind
