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

use Constants, only: cZero
use Definitions, only: u6

implicit none
integer, intent(in) :: n
complex(kind=8), intent(in) :: HAM(n,n)
complex(kind=8), intent(in) :: Akq((n-1),-(n-1):(n-1))
real(kind=8), intent(in) :: B(n,0:n), C(n,0:n), Bstev(n,-n:n)
integer :: k, q, i, j, info
real(kind=8) :: tdiff
complex(kind=8) :: Cp(n,n), Cm(n,n), redME
complex(kind=8) :: O(n,n), W(n,n)
complex(kind=8) :: HCF(n,n), Z(n,n)
real(kind=8) :: w1(n), w2(n), c0, dznrm2_
external :: dznrm2_

do k=2,n-1
  do q=-k,k
    write(u6,'(A,i2,A,i3,A,2ES20.10)') 'Akq(',k,',',q,') = ',Akq(k,q)
  end do
end do
!=======================================================================
write(u6,'(A,ES20.10)') 'recover from Akq parameters'
HCF(:,:) = cZero
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
tdiff = dznrm2_(n*n,HAM-HCF,1)
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
write(u6,'(A,ES20.10)') 'recover from B and C parameters'
HCF(:,:) = cZero
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
tdiff = dznrm2_(n*n,(HAM(1:n,1:n)-HCF(1:n,1:n)),1)
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
tdiff = dznrm2_(n*n,(HAM(1:n,1:n)-HCF(1:n,1:n)),1)
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

return

end subroutine recover_CF
