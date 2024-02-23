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

subroutine project_exch(N1,N2,S1,S2,M1,M2,E1,E2,HEXCH,Jpar,Jc)
! this function determines the local pseudospins and rotates the hamiltonian
! to the local pseudospin basis
! E1, E2 : spin-orbit energies on each site
! S1, S2 : spin matrices on each site
! M1, M2 : magnetic moment matrices on each site
! HEXCH, HEXCH2, HEXCH3  : exchange hamiltonian

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N1, N2
complex(kind=wp), intent(in) :: S1(3,N1,N1), S2(3,N2,N2), M1(3,N1,N1), M2(3,N2,N2), HEXCH(N1,N1,N2,N2)
real(kind=wp), intent(in) :: E1(N1), E2(N2)
complex(kind=wp), intent(out) :: Jpar(N1-1,-(N1-1):N1-1,N2-1,-(N2-1):N2-1)
real(kind=wp), intent(out) :: Jc(3,3)
integer(kind=iwp) :: iprint, is1, is2, ns1, ns2
real(kind=wp) :: Ethr, gtens(2,3), maxes(2,3,3)
complex(kind=wp), allocatable :: HEXCH2(:,:,:,:), HEXCH3(:,:,:,:), MR1(:,:,:), MR2(:,:,:), SR1(:,:,:), SR2(:,:,:), TMP(:,:), &
                                 Z1(:,:), Z2(:,:)

!write(u6,'(A)') 'J parameters in the initial ab intio basis:'
!call JKQPar(N1,N2,HEXCH,Jpar)
!call tensor2cart(1,1,Jpar(1,-1:1,1,-1:1),Jc)

! determine the pseudospin on each site (Z1 and Z2):
! threshold for determination of the local pseudospin main anisotropy axis
Ethr = 0.2_wp
ns1 = 0
ns2 = 0
do is1=1,N1
  if (E1(is1) < Ethr) ns1 = ns1+1
end do
do is1=1,N2
  if (E2(is1) < Ethr) ns2 = ns2+1
end do
write(u6,'(A,i3)') 'size of local pseudospin, site 1  =',ns1
write(u6,'(A,i3)') 'size of local pseudospin, site 2  =',ns2
gtens(:,:) = Zero
maxes(:,:,:) = Zero
call atens(M1(:,1:ns1,1:ns1),ns1,gtens(1,:),maxes(1,:,:),2)
call atens(M2(:,1:ns2,1:ns2),ns2,gtens(2,:),maxes(2,:,:),2)
! rotate the magnetic moment to the coordinate system of main magnetic axes on Ln
call mma_allocate(SR1,3,N1,N1,label='SR1')
call mma_allocate(MR1,3,N1,N1,label='MR1')
call mma_allocate(SR2,3,N2,N2,label='SR2')
call mma_allocate(MR2,3,N2,N2,label='MR2')
call rotmom2(S1,N1,maxes(1,:,:),SR1)
call rotmom2(M1,N1,maxes(1,:,:),MR1)
call rotmom2(S2,N2,maxes(2,:,:),SR2)
call rotmom2(M2,N2,maxes(2,:,:),MR2)
call mma_deallocate(SR1)
call mma_deallocate(SR2)
iprint = 1
call mma_allocate(Z1,N1,N1,label='Z1')
call mma_allocate(Z2,N2,N2,label='Z2')
call pseudospin(MR1,N1,Z1,3,1,iprint)
call pseudospin(MR2,N2,Z2,3,1,iprint)
call mma_deallocate(MR1)
call mma_deallocate(MR2)
! rewrite the exchange matrix in the basis of local pseudospins:
call mma_allocate(HEXCH2,N1,N1,N2,N2,label='HEXCH2')
call mma_allocate(HEXCH3,N1,N1,N2,N2,label='HEXCH3')
call mma_allocate(TMP,max(N1,N2),max(N1,N2),label='TMP')
do is1=1,N2
  do is2=1,N2
    call ZGEMM_('C','N',N1,N1,N1,cOne,Z1,N1,HEXCH(:,:,is1,is2),N1,cZero,TMP(1:N1,1:N1),N1)
    call ZGEMM_('N','N',N1,N1,N1,cOne,TMP(1:N1,1:N1),N1,Z1,N1,cZero,HEXCH2(:,:,is1,is2),N1)
  end do
end do
do is1=1,N1
  do is2=1,N1
    call ZGEMM_('C','N',N2,N2,N2,cOne,Z2,N2,HEXCH2(is1,is2,:,:),N2,cZero,TMP(1:N2,1:N2),N2)
    call ZGEMM_('N','N',N2,N2,N2,cOne,TMP(1:N2,1:N2),N2,Z2,N2,cZero,HEXCH3(is1,is2,:,:),N2)
  end do
end do
call mma_deallocate(Z1)
call mma_deallocate(Z2)
call mma_deallocate(HEXCH2)
call mma_deallocate(TMP)
call JKQPar(N1,N2,HEXCH3,Jpar)
call mma_deallocate(HEXCH3)
call tensor2cart(Jpar(1,-1:1,1,-1:1),Jc)

return

end subroutine project_exch
