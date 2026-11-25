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

subroutine SA_PREC2(rdia,S,CI,ENE)
! S = <i|H0-E|^-1|j>, inverted, over the SA states i, j, for the state of
! energy ENE (H0: ExpHinvv_SA, as in DMinvCI_SA)

use CI_Lib_Info, only: CI_Lib_Space
use input_mclr, only: nCSF, nRoots
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), parameter :: iState = 1
real(kind=wp), intent(in) :: rdia(nCSF(CI_Lib_Space(iState)%STSym)), CI(nCSF(CI_Lib_Space(iState)%STSym),nRoots), ENE
real(kind=wp), intent(out) :: S(nroots,nroots)
integer(kind=iwp) :: i, j, n
real(kind=wp), allocatable :: Tmp(:)
real(kind=wp), external :: DDot_

n = nCSF(CI_Lib_Space(iState)%STSym)
call mma_allocate(Tmp,n,Label='Tmp')
do j=1,nroots
  call ExpHinvv_SA(n,rdia,CI(:,j),Tmp,ENE)
  do i=1,nroots
    S(i,j) = DDot_(n,CI(:,i),1,Tmp,1)
  end do
end do
call mma_deallocate(Tmp)
call MatInvert(S,nroots)

end subroutine SA_PREC2
