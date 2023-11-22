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

subroutine two2mean13(carteSO,occup,AOcoeffs,onecart,ncontmf,norbsum,noccorb)
!bs gives the two first contributions
!bs < i M | j M >  with Malpha  and Mbeta
!bs the other orbit parts cancel

use AMFI_global, only: MxcontL
use Constants, only: Zero
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: ncontmf, norbsum, noccorb
real(kind=wp), intent(in) :: carteSO(ncontmf,ncontmf,norbsum,norbsum), occup(*), AOcoeffs(MxcontL,*)
real(kind=wp), intent(inout) :: onecart(MxcontL,MxcontL)
integer(kind=iwp) :: icartleft, icartright, Mrun
real(kind=wp) :: coeff

do icartleft=1,norbsum
  do icartright=1,norbsum
    coeff = Zero
    do Mrun=1,noccorb
      coeff = coeff+occup(Mrun)*AOcoeffs(icartleft,Mrun)*AOcoeffs(icartright,Mrun)
    end do
    onecart(1:ncontmf,1:ncontmf) = onecart(1:ncontmf,1:ncontmf)+coeff*carteSO(1:ncontmf,1:ncontmf,icartleft,icartright)
  end do
end do
!write(u6,*) 'effective integrals'
!do jrun=1,ncontmf
!  write(u6,'(4ES21.14)') (onecart(irun,jrun),irun=1,ncontmf)
!end do

return

end subroutine two2mean13
