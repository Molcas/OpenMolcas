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

subroutine dkh_ts1e(n,s,t,v,w,ul,us,clight,dkord,xord,dkparam)
! Evaluate the arbitrary order DKH Hamiltonian ( and transform matrices ul & us )

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
! v : store the transformed relativistic one-electron Hamiltonian
integer(kind=iwp), intent(in) :: n, dkord, xord, dkparam
real(kind=wp), intent(in) :: s(n,n), t(n,n), w(n,n), clight
real(kind=wp), intent(inout) :: v(n,n)
real(kind=wp), intent(out) :: ul(n,n), us(n,n)
integer(kind=iwp) :: word, vord, n2
real(kind=wp), allocatable :: Tr(:,:), Bk(:,:), El(:,:), ES(:,:), OL(:,:), OS(:,:), Ep(:), E0(:), KC(:,:), Ws(:,:,:), Co(:), &
                              MW(:,:,:), ZW(:,:,:,:), SCo(:), XL(:,:), XS(:,:)

! Transform Hamiltonian matrix to the free-particle Foldy-Wouthuysen picture

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

! Call DKH transformation routines

! number of W matrices needed
word = max(dkord/2,xord)
! order of DKH needed, since high Xorder may need more transformation than DKHorder (for Hamiltonian)
vord = max(dkord,word*2)
call mma_allocate(Ws,n,n,xord*2,label='Wsav')
call mma_allocate(Co,max(4,vord),label='Cof')
! calculate expansion coefficient of general unitary transformation ( in terms of anti-Hermitian W )
call dkh_cofu(vord,dkparam,Co)

if (dkparam == 2) then
  ! special routine for EXP parameterization ( with fewer matrix multiplication than general routine )
  call mma_allocate(MW,n,n,5,label='NWork')
  call mma_allocate(ZW,n,n,vord,3,label='NNWork')
  MW(:,:,:) = Zero
  ZW(:,:,:,:) = Zero
  call AODKHEXP(n,vord,xord,dkord,Ep,E0,EL,ES,OL,MW(:,:,1),MW(:,:,2),MW(:,:,3),MW(:,:,4),MW(:,:,5),ZW(:,:,:,1),ZW(:,:,:,2), &
                ZW(:,:,:,3),Ws)
  call mma_deallocate(MW)
  call mma_deallocate(ZW)
else
  ! general parameterization routine
  call mma_allocate(SCo,vord,label='Cof2')
  call mma_allocate(MW,n,n,6,label='Mat')
  call mma_allocate(ZW,n,n,vord,10,label='Mat2')
  call dkh_ham(n,dkord,xord,vord,EL,ES,OL,OS,Ep,E0,Co,Sco,MW(:,:,1),MW(:,:,2),MW(:,:,3),MW(:,:,4),MW(:,:,5),MW(:,:,6),ZW(:,:,:,1), &
               ZW(:,:,:,2),ZW(:,:,:,3),ZW(:,:,:,4),ZW(:,:,:,5),ZW(:,:,:,6),ZW(:,:,:,7),ZW(:,:,:,8),ZW(:,:,:,9),ZW(:,:,:,10),Ws)
  call mma_deallocate(SCo)
  call mma_deallocate(MW)
  call mma_deallocate(ZW)
end if

! Calculate the transform matrices

if (xord > 0) then
  call mma_allocate(XL,n,n,label='fpUL')
  call mma_allocate(XS,n,n,label='fpUS')
  n2 = n+n
  call mma_allocate(MW,n2,n2,3,label='TmpZ')
  ! obtain transform matrices ( XL and XS )in fpFW picture
  call dkh_geneu(n,n2,xord,Co,Ws,XL,XS,MW(:,:,1),MW(:,:,2),MW(:,:,3))
  call mma_deallocate(MW)

  call mma_allocate(MW,n,n,4,label='TmpM')
  ! convert to original basis picture
  call XDR_mkutls(n,XL,XS,Tr,Bk,KC(:,1),KC(:,2),KC(:,3),ul,us,MW(:,:,1),MW(:,:,2),MW(:,:,3),MW(:,:,4))
  call mma_deallocate(MW)
  call mma_deallocate(XL)
  call mma_deallocate(XS)
end if

! Back transform Hamiltonian matrix to original non-orthogonal basis picture

call dmxma(n,'C','N',Bk,EL,ES,One)
call dmxma(n,'N','N',ES,Bk,v,One)

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

end subroutine dkh_ts1e
