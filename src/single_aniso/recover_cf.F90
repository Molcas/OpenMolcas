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

subroutine recover_CF(N,HAM,Akq,B,C,Bstev)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: n
complex(kind=wp), intent(in) :: HAM(n,n), Akq(n-1,-(n-1):n-1)
real(kind=wp), intent(in) :: B(n,0:n), C(n,0:n), Bstev(n,-n:n)
integer(kind=iwp) :: i, info, j, k, q
real(kind=wp) :: c0, tdiff
complex(kind=wp) :: redME
real(kind=wp), allocatable :: w1(:), w2(:)
complex(kind=wp), allocatable :: Cm(:,:), Cp(:,:), HCF(:,:), O(:,:), W(:,:), Z(:,:)
real(kind=wp), external :: dznrm2_

do k=2,n-1
  do q=-k,k
    write(u6,'(A,i2,A,i3,A,2ES20.10)') 'Akq(',k,',',q,') = ',Akq(k,q)
  end do
end do
!=======================================================================
write(u6,'(A,ES20.10)') 'recover from Akq parameters'
call mma_allocate(HCF,n,n,label='HCF')
HCF(:,:) = cZero
call mma_allocate(Cm,n,n,label='Cm')
call mma_allocate(Cp,n,n,label='Cp')
do k=1,n-1
  do q=0,k
    ! generate the operator matrix K=ik, Q=iq, dimension=na
    call ITO(n,k,q,C0,Cp,Cm)
    if (q == 0) then
      HCF(:,:) = HCF(:,:)+Akq(k,q)*Cp(:,:)
    else
      HCF(:,:) = HCF(:,:)+Akq(k,q)*Cp(:,:)+Akq(k,-q)*Cm(:,:)
    end if
  end do !q
end do !k
call mma_deallocate(Cm)
call mma_deallocate(Cp)
tdiff = dznrm2_(n*n,HAM-HCF,1)
write(u6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
do i=1,n
  do j=1,n
    write(u6,'(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))') 'HAM(',i,',',j,')=',HAM(i,j),'      ','HCF(',i,',',j,')=',HCF(i,j), &
                                                           ' diff=',HAM(i,j)-HCF(i,j)
  end do
end do
call mma_allocate(w1,n,label='w1')
call mma_allocate(w2,n,label='w2')
call mma_allocate(Z,n,n,label='Z')
call diag_c2(HAM,n,info,w1,Z)
call diag_c2(HCF,n,info,w2,Z)
do i=1,n
  write(u6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))') 'W1(',i,')=',w1(i)-w1(1),'      ','W2(',i,')=',w2(i)-w2(1),' diff=',w1(i)-w2(i)
end do

!=======================================================================
write(u6,'(A,ES20.10)') 'recover from B and C parameters'
HCF(:,:) = cZero
call mma_allocate(O,n,n,label='O')
call mma_allocate(W,n,n,label='W')
do k=1,n-1
  do q=0,k
    call Liviu_ESO(n,k,q,O,W,redME)
    if (q == 0) then
      HCF(:,:) = HCF(:,:)+B(k,0)*O(:,:)
    else
      HCF(:,:) = HCF(:,:)+B(k,q)*O(:,:)+C(k,q)*W(:,:)
    end if
  end do
end do
tdiff = dznrm2_(n*n,(HAM-HCF),1)
write(u6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
do i=1,n
  do j=1,n
    write(u6,'(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))') 'HAM(',i,',',j,')=',HAM(i,j),'      ','HCF(',i,',',j,')=',HCF(i,j), &
                                                           ' diff=',HAM(i,j)-HCF(i,j)
  end do
end do
call diag_c2(HAM,n,info,w1,Z)
call diag_c2(HCF,n,info,w2,Z)
do i=1,n
  write(u6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))') 'W1(',i,')=',w1(i)-w1(1),'      ','W2(',i,')=',w2(i)-w2(1),' diff=',w1(i)-w2(i)
end do

!=======================================================================
write(6,'(A,ES20.10)') 'recover from Bstev'
HCF(:,:) = cZero
do k=1,n-1
  do q=0,k
    call ESO(n,k,q,O,W,redME)
    if (q == 0) then
      HCF(:,:) = HCF(:,:)+Bstev(k,0)*O(:,:)
    else
      HCF(:,:) = HCF(:,:)+Bstev(k,q)*O(:,:)+Bstev(k,-q)*W(:,:)
    end if
  end do
end do
tdiff = dznrm2_(n*n,(HAM-HCF),1)
write(u6,'(A,ES20.10)') 'total difference between HAM-HCF=',tdiff
do i=1,n
  do j=1,n
    write(u6,'(2(A,i2,A,i2,A,2ES20.10,A),2(2ES20.10,5x))') 'HAM(',i,',',j,')=',HAM(i,j),'      ','HCF(',i,',',j,')=',HCF(i,j), &
                                                           ' diff=',HAM(i,j)-HCF(i,j)
  end do
end do
call diag_c2(HAM,n,info,w1,Z)
call diag_c2(HCF,n,info,w2,Z)
do i=1,n
  write(u6,'(2(A,i2,A,ES20.10,A),2(2ES20.10,5x))') 'W1(',i,')=',w1(i)-w1(1),'      ','W2(',i,')=',w2(i)-w2(1),' diff=',w1(i)-w2(i)
end do
!=======================================================================

call mma_deallocate(w1)
call mma_deallocate(w2)
call mma_deallocate(Z)
call mma_deallocate(O)
call mma_deallocate(W)
call mma_deallocate(HCF)

return

end subroutine recover_CF
