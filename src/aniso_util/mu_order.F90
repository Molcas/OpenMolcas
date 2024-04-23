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

subroutine mu_order(d,MS,MM,gtens,order,HCF2,AMM,AMS,Z,iprint)
! This Subroutine receives the moment matrix dipso(3,d,d) and Returns the matrix re-builted using only the 1-st order operators.
!   MS: initial spin moment
!   MM: initial magnetic moment
!  AMS: transformed spin moment
!  AMM: transformed magnetic moment

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: One, Half, cZero, Onei
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: d, order, iprint
complex(kind=wp), intent(in) :: MS(3,d,d), MM(3,d,d)
real(kind=wp), intent(out) :: gtens(3)
complex(kind=wp), intent(out) :: HCF2(d,3,d,d), AMM(3,d,d), AMS(3,d,d), Z(d,d)
integer(kind=iwp) :: i, i1, i2, j, l, m, n
real(kind=wp) :: m_fact, maxes(3,3)
complex(kind=wp) :: O1, O2, SP_DIPO(3), SP_DIPW(3), SP_MOW
complex(kind=wp), allocatable :: AMS_TMP(:,:), B(:,:,:), BNMC(:,:,:), BNMS(:,:,:), DIP_O(:,:), DIP_W(:,:), TMP(:,:,:)
complex(kind=wp), external :: trace
!------------------------------------------------------------

Z(:,:) = cZero
! get the local pseudospin:
call pseudospin(MM,d,Z,3,1,iprint)
! re-write MM and MS to the new pseudospin basis:
AMS(:,:,:) = cZero
AMM(:,:,:) = cZero
do i=1,d
  do j=1,d
    do i1=1,d
      do i2=1,d
        AMM(:,i,j) = AMM(:,i,j)+MM(:,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
        AMS(:,i,j) = AMS(:,i,j)+MS(:,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
      end do
    end do
  end do
end do
if (iprint > 2) then
  call prMom('MU_ORDER:   AMM(l,i,j):',AMM,d)
  call prMom('MU_ORDER:   AMS(l,i,j):',AMS,d)
end if
! project the moment on ITO:
! obtain the b3m and c3m coefficients:
call mma_allocate(B,[1,3],[1,d],[-d,d],label='B')
call mma_allocate(BNMC,[1,3],[1,d],[0,d],label='BNMC')
call mma_allocate(BNMS,[1,3],[1,d],[0,d],label='BNMS')
call mma_allocate(DIP_O,d,d,label='DIP_O')
call mma_allocate(DIP_W,d,d,label='DIP_W')
call mma_allocate(AMS_TMP,d,d,label='AMS_TMP')
B(:,:,:) = cZero
do N=1,d-1
  do M=0,N
    DIP_O(:,:) = cZero
    DIP_W(:,:) = cZero
    call Stewens_matrixel(N,M,d,DIP_O,DIP_W,iprint)
    if (iprint > 5) then
      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_O, n = ',N,', m =',m
      write(u6,*)
      do i=1,d
        write(u6,'(20(2F10.6,2x))') (DIP_O(i,j),j=1,d)
      end do
      write(u6,*)
      write(u6,'( 5x,a,i2,a,i3)') 'DIP_W, n = ',N,', m =',m
      write(u6,*)
      do i=1,d
        write(u6,'(20(2F10.6,2x))') (DIP_W(i,j),j=1,d)
      end do
    end if

    SP_DIPO(:) = czero
    SP_DIPW(:) = czero
    SP_MOW = trace(d,DIP_O,DIP_W)
    do l=1,3
      AMS_TMP(:,:) = AMS(l,:,:)
      SP_DIPO(l) = trace(d,AMS_TMP,DIP_O)
      SP_DIPW(l) = trace(d,AMS_TMP,DIP_W)

      B(l,n,-m) = SP_DIPO(l)/SP_MOW
      B(l,n,m) = SP_DIPW(l)/SP_MOW
    end do ! l
  end do ! m
end do ! n

do n=1,d-1
  BNMC(:,n,0) = B(:,n,0)
  do m=1,N
    m_fact = (-One)**M
    BNMC(:,n,m) = B(:,n,m)+m_fact*B(:,n,-m)
    BNMS(:,n,m) = -Onei*(B(:,n,m)-m_fact*B(:,n,-m))
  end do
end do !n

if (iprint > 2) then
  write(u6,'(2A)') repeat('-',47),'|'
  write(u6,'(A)') '  n  |  m  |   |       B       |       C       |'
  do N=1,d-1
    write(u6,'(A)') '-----|-----|---|---------------|---------------|'
    do M=0,N
      if (M /= 0) write(u6,'(A)') '     |-----|---|---------------|---------------|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','X','|',real(BNMC(1,N,M)),'|',real(BNMS(1,N,M)),'|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Y','|',real(BNMC(2,N,M)),'|',real(BNMS(2,N,M)),'|'
      write(u6,'(2(2x,I1,2x,A),1x,A,1x,A,2(F14.9,1x,A))') N,'|',M,'|','Z','|',real(BNMC(3,N,M)),'|',real(BNMS(3,N,M)),'|'
    end do
  end do
  write(u6,'(2A)') repeat('-',47),'|'
end if

HCF2(:,:,:,:) = cZero
do N=1,d-1
  do M=0,N
    DIP_O(:,:) = cZero
    DIP_W(:,:) = cZero

    call Stewens_matrixel(N,M,d,DIP_O,DIP_W,iprint)

    do l=1,3
      do i=1,d
        do j=1,d
          if (M == 0) then
            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*DIP_O(i,j)
          else
            m_fact = (-One)**M
            O1 = Half*(m_fact*DIP_W(i,j)+DIP_O(i,j))
            O2 = -Half*Onei*(m_fact*DIP_W(i,j)-DIP_O(i,j))

            HCF2(N,l,i,j) = HCF2(N,l,i,j)+BNMC(l,N,M)*O1+BNMS(l,N,M)*O2
          end if
        end do
      end do
    end do
  end do
end do ! n

call mma_deallocate(B)
call mma_deallocate(BNMC)
call mma_deallocate(BNMS)
call mma_deallocate(DIP_O)
call mma_deallocate(DIP_W)
call mma_deallocate(AMS_TMP)

if (iprint > 2) then
  do N=1,d-1,2
    write(u6,*)
    write(u6,'( 5x,a,I2,a)') 'HCF2(',N,',l,i,j)'
    do l=1,3
      write(u6,*)
      write(u6,'(a,i3)') 'PROJECTION =',l
      do i=1,d
        write(u6,'(20(2F12.8,2x))') (HCF2(N,l,i,j),j=1,d)
      end do
    end do
    write(u6,*)
  end do
end if

call mma_allocate(TMP,3,d,d,label='TMP')
TMP(:,:,:) = HCF2(order,:,:,:)
call ATENS(TMP,d,gtens,maxes,1)
call mma_deallocate(TMP)

return

end subroutine mu_order
