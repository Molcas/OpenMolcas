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

subroutine moments(N,MS,MM,iprint)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, gElectron
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: N, iprint
complex(kind=wp), intent(in) :: MS(3,N,N), MM(3,N,N)
integer(kind=iwp) :: I, i1, i2, iDir, J, L
complex(kind=wp) :: Mf(3,3), Sf(3,3), Lf(3,3)
complex(kind=wp), allocatable :: AML(:,:,:), AMM(:,:,:), AMS(:,:,:), Z(:,:)
real(kind=wp), parameter :: g_e = -gElectron
!-----------------------------------------------------------------------

if (N < 1) return
call mma_allocate(Z,N,N,'Z')
call mma_allocate(AMS,3,N,N,'AMS')
call mma_allocate(AML,3,N,N,'AML')
call mma_allocate(AMM,3,N,N,'AMM')
!-----------------------------------------------------------------------
do iDir=1,3
  Z(:,:) = cZero
  AMM(:,:,:) = cZero
  AML(:,:,:) = cZero
  AMS(:,:,:) = cZero

  call pseudospin(MM,N,Z,iDir,1,1)

  do l=1,3
    do i=1,N
      do j=1,N
        do i1=1,N
          do i2=1,N
            AMM(l,i,j) = AMM(l,i,j)+MM(l,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
            AMS(l,i,j) = AMS(l,i,j)+MS(l,i1,i2)*conjg(Z(i1,i))*Z(i2,j)
          end do
        end do
        AML(l,i,j) = -AMM(l,i,j)-g_e*AMS(l,i,j)
      end do
    end do
  end do
  Mf(iDir,:) = AMM(:,1,1)
  Sf(iDir,:) = AMS(:,1,1)
  Lf(iDir,:) = AML(:,1,1)

  if (iprint > 3) then
    write(u6,*)
    write(u6,'(2x,a)') 'MOMENTS:  AMM(l,:,:)'
    do l=1,3
      write(u6,*)
      write(u6,'(a,i3)') 'PROJECTION2 =',l
      do i=1,N
        write(u6,'(20(2F12.8,2x))') (AMM(l,i,j),j=1,N)
      end do
    end do
    write(u6,*)
    write(u6,'(2x,a)') 'MOMENTS:  AMS(l,:,:)'
    do l=1,3
      write(u6,*)
      write(u6,'(a,i3)') 'PROJECTION2 =',l
      do i=1,N
        write(u6,'(20(2F12.8,2x))') (AMS(l,i,j),j=1,N)
      end do
    end do
    write(u6,*)
    write(u6,'(2x,a)') 'MOMENTS:  AML(l,:,:)'
    do l=1,3
      write(u6,*)
      write(u6,'(a,i3)') 'PROJECTION2 =',l
      do i=1,N
        write(u6,'(20(2F12.8,2x))') (AML(l,i,j),j=1,N)
      end do
    end do
  end if
end do !iDir

write(u6,'(A)') '--------------|--- H || Xm --|--- H || Ym --|--- H || Zm --|'
write(u6,'(8A)') repeat('-',14),'|',(repeat('-',14),'|',j=1,3)
write(u6,'(A,3(F13.9,1x,A))') ' <1| mu_X |1> |',(real(Mf(iDir,1)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') ' <1| mu_Y |1> |',(real(Mf(iDir,2)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') ' <1| mu_Z |1> |',(real(Mf(iDir,3)),'|',iDir=1,3)
write(u6,'(8a)') repeat('-',14),'|',(repeat('-',14),'|',j=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| L_X |1> |',(real(Lf(iDir,1)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| L_Y |1> |',(real(Lf(iDir,2)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| L_Z |1> |',(real(Lf(iDir,3)),'|',iDir=1,3)
write(u6,'(8a)') repeat('-',14),'|',(repeat('-',14),'|',j=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| S_X |1> |',(real(Sf(iDir,1)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| S_Y |1> |',(real(Sf(iDir,2)),'|',iDir=1,3)
write(u6,'(A,3(F13.9,1x,A))') '  <1| S_Z |1> |',(real(Sf(iDir,3)),'|',iDir=1,3)
write(u6,'(2a)') repeat('-',59),'|'

!-----------------------------------------------------------------------
call mma_deallocate(Z)
call mma_deallocate(AMS)
call mma_deallocate(AML)
call mma_deallocate(AMM)

return

end subroutine moments
