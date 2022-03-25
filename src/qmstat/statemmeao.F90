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

! AO-basis route.
subroutine StateMMEao(nAObas,nState,nTyp,iBigT,iMME,iCent,Cha,Dip,Qua)

use qmstat_global, only: MxMltp
use Index_Functions, only: nTri3_Elem, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp) :: nAObas, nState, nTyp, iBigT, iMME(nTri3_Elem(MxMltp)), iCent(nTri_Elem(nAObas))
real(kind=wp) :: Cha(nTri_Elem(nState),*), Dip(nTri_Elem(nState),3,*), Qua(nTri_Elem(nState),6,*)
#include "WrkSpc.fh"
integer(kind=iwp) :: iB1, iB2, iS1, iS2, iTyp, kaunta, kaunter, nSize
real(kind=wp) :: PerAake
real(kind=wp), allocatable :: AOG(:), O(:)

kaunter = 0
nSize = nTri_Elem(nAObas)
call mma_allocate(AOG,nSize,label='Transition')
call mma_allocate(O,nTyp,label='OnTheWay')
! Loop over state pairs.
do iS1=1,nState
  do iS2=1,iS1
    kaunter = kaunter+1
    ! Collect this piece of the TDM in AO-basis.
    call dCopy_(nSize,Work(iBigT+nSize*(kaunter-1)),1,AOG,1)
    kaunta = 0
    ! Loop over AO-basis pairs and transform them as well as
    ! distribute their multipoles. Observe that the array iCent
    ! keeps track on where a certain AO-basis pair belongs.
    do iB1=1,nAObas
      do iB2=1,iB1
        kaunta = kaunta+1
        PerAake = AOG(kaunta)
        do iTyp=1,nTyp
          O(iTyp) = Work(iMME(iTyp)-1+kaunta)*PerAake
        end do
        Cha(kaunter,iCent(kaunta)) = Cha(kaunter,iCent(kaunta))+O(1)
        Dip(kaunter,1,iCent(kaunta)) = Dip(kaunter,1,iCent(kaunta))+O(2)
        Dip(kaunter,2,iCent(kaunta)) = Dip(kaunter,2,iCent(kaunta))+O(3)
        Dip(kaunter,3,iCent(kaunta)) = Dip(kaunter,3,iCent(kaunta))+O(4)
        Qua(kaunter,1,iCent(kaunta)) = Qua(kaunter,1,iCent(kaunta))+O(5)
        Qua(kaunter,2,iCent(kaunta)) = Qua(kaunter,2,iCent(kaunta))+O(6)
        ! The reason why 8 and 7 are interchanged is that
        ! QMSTAT uses the ordering xx,xy,yy,xz,yz,zz while
        ! Seward uses the ordering xx,xy,xz,yy,yz,zz.
        Qua(kaunter,3,iCent(kaunta)) = Qua(kaunter,3,iCent(kaunta))+O(8)
        Qua(kaunter,4,iCent(kaunta)) = Qua(kaunter,4,iCent(kaunta))+O(7)
        Qua(kaunter,5,iCent(kaunta)) = Qua(kaunter,5,iCent(kaunta))+O(9)
        Qua(kaunter,6,iCent(kaunta)) = Qua(kaunter,6,iCent(kaunta))+O(10)
      end do
    end do
  end do
end do
call mma_deallocate(AOG)
call mma_deallocate(O)

return

end subroutine StateMMEao
