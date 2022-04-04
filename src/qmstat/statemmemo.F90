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

! MO-basis route.
subroutine StateMMEmo(nAObas,nMObas,nState,nTyp,MME,iCent,Cha,Dip,Qua)

use qmstat_global, only: AvRed, BigT, MxMltp
use Index_Functions, only: nTri3_Elem, nTri_Elem
use Data_Structures, only: Alloc1DArray_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Half
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAObas, nMObas, nState, nTyp, iCent(nTri_Elem(nAObas))
type(Alloc1DArray_Type), intent(in) :: MME(nTri3_Elem(MxMltp))
real(kind=wp), intent(inout) :: Cha(nTri_Elem(nState),*), Dip(nTri_Elem(nState),3,*), Qua(nTri_Elem(nState),6,*)
integer(kind=iwp) :: i, iB1, iB2, iS1, iS2, iTyp, j, kaunta, kaunter, nSizeA, nSizeM
real(kind=wp) :: PerAake
real(kind=wp), allocatable :: AOG(:), AOG_S(:,:), MOG(:), MOG_s(:,:), O(:), TEMP(:,:)
! The reason why 8 and 7 are interchanged is that
! QMSTAT uses the ordering xx,xy,yy,xz,yz,zz while
! Seward uses the ordering xx,xy,xz,yy,yz,zz.
integer(kind=iwp), parameter :: xTyp(10) = [1,2,3,4,5,6,8,7,9,10]

kaunter = 0
nSizeA = nTri_Elem(nAObas)
nSizeM = nTri_Elem(nMObas)
call mma_allocate(MOG,nSizeM,label='Transition')
call mma_allocate(MOG_s,nMObas,nMObas,label='SqMO')
call mma_allocate(TEMP,nAObas,nMObas,label='TEMP')
call mma_allocate(AOG_s,nAObas,nAObas,label='SqAO')
call mma_allocate(AOG,nSizeA,label='TransitionA')
call mma_allocate(O,nTyp,label='OnTheWay')

! Loop over state pairs.

do iS1=1,nState
  do iS2=1,iS1
    kaunter = kaunter+1

    ! Collect the proper piece of the TDM in MO-basis.

    MOG(:) = BigT(1:nSizeM,kaunter)

    ! Additional transformation step from MO to AO.

    call Square(MOG,MOG_s,1,nMObas,nMObas)
    do i=1,nMObas
      do j=1,nMObas
        if (i == j) cycle
        MOG_s(j,i) = Half*MOG_s(j,i)
      end do
    end do
    call Dgemm_('N','N',nAObas,nMObas,nMObas,One,AvRed,nAObas,MOG_s,nMObas,Zero,TEMP,nAObas)
    call Dgemm_('N','T',nAObas,nAObas,nMObas,One,TEMP,nAObas,AvRed,nAObas,Zero,AOG_s,nAObas)
    do i=1,nAObas
      do j=1,nAObas
        if (i == j) cycle
        AOG_s(j,i) = Two*AOG_s(j,i)
      end do
    end do
    call SqToTri_Q(AOG_s,AOG,nAObas)

    ! Loop over AO-basis pairs.

    kaunta = 0
    do iB1=1,nAObas
      do iB2=1,iB1
        kaunta = kaunta+1
        PerAake = AOG(kaunta)
        do iTyp=1,nTyp
          O(iTyp) = MME(xTyp(iTyp))%A(kaunta)*PerAake
        end do
        Cha(kaunter,iCent(kaunta)) = Cha(kaunter,iCent(kaunta))+O(1)
        Dip(kaunter,:,iCent(kaunta)) = Dip(kaunter,:,iCent(kaunta))+O(2:4)
        Qua(kaunter,:,iCent(kaunta)) = Qua(kaunter,:,iCent(kaunta))+O(5:10)
      end do
    end do
  end do
end do
call mma_deallocate(MOG)
call mma_deallocate(MOG_s)
call mma_deallocate(TEMP)
call mma_deallocate(AOG_s)
call mma_deallocate(AOG)
call mma_deallocate(O)

return

end subroutine StateMMEmo
