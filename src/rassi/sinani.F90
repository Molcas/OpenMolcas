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

subroutine SINANI(KDGN,IFUNCT,NSS,DIPSOm,SPNSFS,DIPSOm_SA)

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: cZero, cOne
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: KDGN, IFUNCT, NSS
complex(kind=wp), intent(in) :: DIPSOm(3,NSS,NSS)
complex(kind=wp), intent(in) :: SPNSFS(3,NSS,NSS)
real(kind=wp), intent(in) :: DIPSOm_SA
integer(kind=iwp) :: i, Ico1, Iso1, j, jCo1, Jso2, l
real(kind=wp) :: gtens(3), maxes(3,3)
character :: angm
real(kind=wp), allocatable :: UMATI(:,:), UMATR(:,:)
complex(kind=wp), allocatable :: DIPSOmSA(:,:,:), MATL(:,:), SPNSO(:,:,:), SPNSOSA(:,:,:), TEMP(:,:), Z(:,:)

if (.false.) then
  write(u6,'(/)')
  write(u6,'(A)') repeat('#',120)
  write(u6,'(25X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT IN THE BASIS OF SPIN ORBIT STATES'
  write(u6,'(A)') repeat('#',120)

  do l=1,3
    if (l == 1) angm = 'X'
    if (l == 2) angm = 'Y'
    if (l == 3) angm = 'Z'
    write(u6,'(/)')
    write(u6,'(4X,A12,A2)') 'PROJECTION: ',angm
    write(u6,'(/)')
    do Iso1=1,NSS
      write(u6,'(20(2X,2F10.6))') (DIPSOm(l,Iso1,Jso2),Jso2=1,NSS)
    end do
  end do
  write(u6,'(/)')

end if !if (IPGLOB >= 4)

call mma_allocate(DIPSOmSA,3,KDGN,KDGN,Label='DIPSOmSA')
do l=1,3
  do Ico1=1,KDGN
    do Jco1=1,KDGN
      DIPSOmSA(l,Ico1,Jco1) = cZero
      !S_SOM(L,I,J) = cZero
    end do
  end do
end do

do Iso1=1,KDGN
  do Jso2=1,KDGN
    Ico1 = Iso1+IFUNCT
    Jco1 = Jso2+IFUNCT
    do l=1,3
      !write(u6,*) 'DIPSOm',DIPSOm(l,Ico1,Jco1)
      DIPSOmSA(l,Iso1,Jso2) = -DIPSOm(l,Ico1,Jco1)
      !S_SOM(l,i,j) = S_SO(l,ic1,ic2)
    end do
  end do
end do

if (.false.) then
  write(u6,*)
  write(u6,'(10X,A)') 'MATRIX ELEMENTS OF THE MAGNETIC MOMENT IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
  do l=1,3
    write(u6,'(/)')
    write(u6,'(5X,A6,I3)') 'AXIS= ',l
    write(u6,*)
    do Ico1=1,KDGN
      write(u6,'(16(2F12.8,2x))') (DIPSOmSA(l,Ico1,Jco1),Jco1=1,KDGN)
    end do
  end do

end if

call ATENS_RASSI(DIPSOmSA,KDGN,gtens,maxes,3)
call mma_deallocate(DIPSOmSA)

if (.false.) then

  call mma_allocate(UMATR,NSS,NSS,Label='UMATR')
  call mma_allocate(UMATI,NSS,NSS,Label='UMATI')
  call mma_allocate(SPNSO,3,NSS,NSS,Label='SPNSO')
  call mma_allocate(Z,NSS,NSS,Label='Z')
  call get_dArray('UMATR_SINGLE',UMATR,NSS**2)
  call get_dArray('UMATI_SINGLE',UMATI,NSS**2)
  write(u6,'(/)')
  write(u6,'(5x,a)') 'umatr and umati'
  write(u6,'(/)')
  do i=1,NSS
    write(u6,'(5x,10(2f14.10,2x))') (UMATR(i,j),UMATI(i,j),j=1,NSS)
  end do

  do I=1,NSS
    do J=1,NSS
      do L=1,3
        SPNSO(L,I,J) = cZero
      end do
      Z(I,J) = cZero
    end do
  end do

  do i=1,NSS
    do j=1,NSS
      Z(i,j) = Z(i,j)+cmplx(UMATR(i,j),UMATI(i,j),kind=wp)
    end do
  end do
  call mma_deallocate(UMATR)
  call mma_deallocate(UMATI)

  call mma_allocate(TEMP,NSS,NSS,Label='TEMP')
  call mma_allocate(MATL,NSS,NSS,Label='MATL')
  do l=1,3
    do i=1,NSS
      do j=1,NSS
        MATL(i,j) = SPNSFS(L,i,j)
      end do
    end do

    call ZGEMM('C','N',NSS,NSS,NSS,cOne,Z,NSS,MATL,NSS,cZero,TEMP,NSS)
    call ZGEMM('N','N',NSS,NSS,NSS,cOne,TEMP,NSS,Z,NSS,cZero,MATL,NSS)

    do i=1,NSS
      do j=1,NSS
        SPNSO(L,i,j) = MATL(i,j)
      end do
    end do
  end do !l
  call mma_deallocate(TEMP)
  call mma_deallocate(MATL)

  write(u6,'(/)')
  write(u6,'(A)') repeat('#',120)
  write(u6,'(30X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT IN THE BASIS OF SPIN ORBIT STATES'
  write(u6,'(A)') repeat('#',120)
  write(u6,'(/)')
  do l=1,3
    if (l == 1) angm = 'X'
    if (l == 2) angm = 'Y'
    if (l == 3) angm = 'Z'
    write(u6,'(/)')
    write(u6,'(4X,A,A)') 'PROJECTION: ',angm
    write(u6,'(/)')
    do Iso1=1,NSS
      write(u6,'(20(2F10.6,2X))') (SPNSO(l,Iso1,Jso2),Jso2=1,NSS)
    end do
  end do

  call mma_deallocate(Z)
  call mma_allocate(SPNSOSA,3,KDGN,KDGN,Label='SPNSOSA')

  do l=1,3
    do Ico1=1,KDGN
      do Jco1=1,KDGN
        SPNSOSA(l,Ico1,Jco1) = cZero
      end do
    end do
  end do

  do Iso1=1,KDGN
    do Jso2=1,KDGN
      Ico1 = Iso1+IFUNCT
      Jco1 = Jso2+IFUNCT
      do l=1,3
        SPNSOSA(l,Iso1,Jso2) = SPNSO(l,Ico1,Jco1)
      end do
    end do
  end do

  write(u6,*)
  write(u6,'(10X,A)') 'MATRIX ELEMENTS OF THE SPIN MOMENT IN THE BASIS OF SPIN-ORBIT FUNCTIONS'
  do l=1,3
    write(u6,'(/)')
    write(u6,'(5X,A6,I3)') 'AXIS= ',l
    write(u6,*)
    do Ico1=1,KDGN
      write(u6,'(16(2F12.8,2x))') (SPNSOSA(l,Ico1,Jco1),Jco1=1,KDGN)
    end do
  end do
  call mma_deallocate(SPNSOSA)
  call mma_deallocate(SPNSO)

end if

!do l=1,3
!  do Iso1=1,KDGN
!    do Jso2=1,KDGN
!      Ico1 = Iso1+IFUNCT
!      Jco1 = Jso2+IFUNCT
!      write(u6,*) 'DIPSOm',DIPSOm(l,Ico1,Jco1)
!      IPSOm_SA(l,Iso1,Jso2) = -DIPSOm(l,Ico1,Jco1)
!      S_SOM( l,i,j) = S_SO(l,ic1,ic2)
!    end do
!  end do
!end do

! Avoid unused argument warnings
#include "macros.fh"
unused_var(DIPSOm_SA)

end subroutine SINANI
