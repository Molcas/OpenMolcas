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

subroutine dkh_prop(n,s,t,v,w,X,pXp,clight,xord,dkparam)
! Apply the arbitrary order DKH transformation to property integral

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One
use Definitions, only: wp, iwp

implicit none
! Input :
!   X     matrix of property operator
!   pXp   matrix representation of <pxXpx>+<pyXpy>+<pzXpz>
!   w     aka pVp
! Output :
!   X     store the transformed property integral
integer(kind=iwp), intent(in) :: n, xord, dkparam
real(kind=wp), intent(in) :: s(n,n), t(n,n), v(n,n), w(n,n), clight
real(kind=wp), intent(inout) :: X(n,n), pXp(n,n)
integer(kind=iwp) :: vord
real(kind=wp), allocatable :: Tr(:,:), Bk(:,:), EL(:,:), ES(:,:), OL(:,:), OS(:,:), Ep(:), E0(:), KC(:,:), Ws(:,:,:), Co(:), &
                              SCo(:), MW(:,:,:), ZW(:,:,:,:)

! Transform to free-particle FW picture

call mma_allocate(Tr,n,n,label='Tr')
call mma_allocate(Bk,n,n,label='Back')
call mma_allocate(EL,n,n,label='mEL')
call mma_allocate(ES,n,n,label='mES')
call mma_allocate(OL,n,n,label='mOL')
call mma_allocate(OS,n,n,label='mOS')
call mma_allocate(Ep,n,label='Ep')
call mma_allocate(E0,n,label='E0')
call mma_allocate(KC,n,3,label='KC')
call XDR_fpFW(n,s,t,v,w,Tr,Bk,EL,ES,OL,OS,Ep,E0,KC(:,1),KC(:,2),KC(:,3),clight)

! Calculate the DKH unitary transformation ( in terms of a series of W matrices 1..xord )

vord = xord*2
call mma_allocate(Ws,n,n,xord*2,label='Wsav')
call mma_allocate(Co,max(4,vord),label='Cof')
call dkh_cofu(vord,dkparam,Co)

call mma_allocate(Sco,max(4,vord),label='Cof2')
call mma_allocate(MW,n,n,6,label='Mat')
call mma_allocate(ZW,n,n,vord,10,label='Mat2')
call dkh_ham(n,vord,xord,vord,EL,ES,OL,OS,Ep,E0,Co,SCo,MW(:,:,1),MW(:,:,2),MW(:,:,3),MW(:,:,4),MW(:,:,5),MW(:,:,6),ZW(:,:,:,1), &
             ZW(:,:,:,2),ZW(:,:,:,3),ZW(:,:,:,4),ZW(:,:,:,5),ZW(:,:,:,6),ZW(:,:,:,7),ZW(:,:,:,8),ZW(:,:,:,9),ZW(:,:,:,10),Ws)

! Apply W[1..xord] determined transformation to property operator X

! convert X to fpFW picture
call XDR_fpFWprop(n,Tr,X,pXp,KC(:,1),KC(:,2),KC(:,3),EL,ES,OL,OS,MW(:,:,1))
call dkh_xpx(n,vord,xord,vord,EL,ES,OL,OS,Ep,E0,Co,SCo,MW(:,:,1),MW(:,:,2),MW(:,:,3),MW(:,:,4),MW(:,:,5),MW(:,:,6),ZW(:,:,:,1), &
             ZW(:,:,:,2),ZW(:,:,:,3),ZW(:,:,:,4),ZW(:,:,:,5),ZW(:,:,:,6),ZW(:,:,:,7),ZW(:,:,:,8),ZW(:,:,:,9),ZW(:,:,:,10),Ws)
call mma_deallocate(SCo)
call mma_deallocate(MW)
call mma_deallocate(ZW)

! Back transform to original non-orthogonal basis picture

call dmxma(n,'C','N',Bk,EL,ES,One)
call dmxma(n,'N','N',ES,Bk,X,One)

! Free temp memories

call mma_deallocate(Co)
call mma_deallocate(Ws)
call mma_deallocate(Tr)
call mma_deallocate(Bk)
call mma_deallocate(EL)
call mma_deallocate(ES)
call mma_deallocate(OL)
call mma_deallocate(OS)
call mma_deallocate(Ep)
call mma_deallocate(E0)
call mma_deallocate(KC)

return

end subroutine dkh_prop
