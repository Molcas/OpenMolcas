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

subroutine RasRasTrans(nB,nStatePrim,Eig2,iPrint)

use qmstat_global, only: BigT, nState, RassiM
use Index_Functions, only: iTri, nTri_Elem
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nB, nStatePrim, iPrint
real(kind=wp), intent(in) :: Eig2(nStatePrim,nStatePrim)
integer(kind=iwp) :: i, iB, iBas, iDisk, iiS, indx, indypop, iS, j, jB, jBas, jjS, jS, kaunt, kaunter, LuIn, MEMMAX, nSize, &
                     nSizeBig, nSizeBigPrim, nTriS, nTriSP
character(len=30) :: OutLine
integer(kind=iwp), allocatable :: iTocBig(:)
real(kind=wp), allocatable :: AOG(:), BigV(:,:), Int1(:), Int2(:), Snt1(:,:), Snt2(:,:), Snt3(:,:)
#include "warnings.h"

!Guten Tag.

write(u6,*) '     ----- Transform from non-orthogonal RASSCF states to orthogonal RASSI states.'

! Set zeros and decide if transformation is at all possible.

LuIn = 66
kaunt = 0
iDisk = 0
call DaName(LuIn,RassiM)
nTriSP = nTri_Elem(nStatePrim)
call mma_allocate(iTocBig,nTriSP,label='iTocBig')
call iDaFile(LuIn,2,iTocBig,nTriSP,iDisk)
nSize = nTri_Elem(nB)
nTriS = nTri_Elem(nState)
nSizeBig = nSize*nTriS
nSizeBigPrim = nSize*nTriSP
call mma_maxDBLE(MEMMAX)

if (MEMMAX <= nSizeBig) then
  ! This means that we do not have memory enough for TDM in contracted
  ! form. Then there is no use to proceed at all.

  write(u6,*)
  write(u6,*) 'The transition density matrix is too big to put in memory!'
  write(u6,*) 'Either,'
  write(u6,*) '       (1) increase MOLCAS_MEM,'
  write(u6,*) '       (2) contract number of states further.'
  call Quit(_RC_GENERAL_ERROR_)

else if (MEMMAX >= (nSizeBig+nSizeBigPrim+nTriSP+nTriS+nStatePrim**2+nState*nStatePrim+nState**2)) then
  ! Here we go if there is enough memory for an in core transformation.

  call mma_allocate(BigT,nSize,nTriS,label='ALLES')
  call mma_allocate(BigV,nSize,nTriSP,label='ALLESin')
  call mma_allocate(Int1,nTriSP,label='Int1')
  call mma_allocate(Int2,nTriS,label='Int2')
  call mma_allocate(Snt1,nStatePrim,nStatePrim,label='Square1')
  call mma_allocate(Snt2,nState,nStatePrim,label='Square2')
  call mma_allocate(Snt3,nState,nState,label='Square3')
  BigT(:,:) = Zero
  kaunt = 0
  do i=1,nStatePrim
    do j=1,i
      kaunt = kaunt+1
      iDisk = iTocBig(kaunt)
      call dDaFile(LuIn,2,BigV(:,kaunt),nSize,iDisk)
    end do
  end do

  ! A lot of printing of TDM if requested.

  if (iPrint >= 25) then
    kaunt = 0
    do i=1,nStatePrim
      do j=1,i
        kaunt = kaunt+1
        write(OutLine,'(A,I3,I3)') 'TDM, Piece ',i,j
        call TriPrt(OutLine,' ',BigV(:,kaunt),nB)
      end do
    end do
  end if

  ! Proceed with transformation.

  kaunt = 0
  do iBas=1,nB
    do jBas=1,iBas
      kaunt = kaunt+1
      Int1(:) = BigV(kaunt,:)
      call Square(Int1,Snt1,1,nStatePrim,nStatePrim)
      call Dgemm_('T','N',nState,nStatePrim,nStatePrim,One,Eig2,nStatePrim,Snt1,nStatePrim,Zero,Snt2,nState)
      call Dgemm_('N','N',nState,nState,nStatePrim,One,Snt2,nState,Eig2,nStatePrim,Zero,Snt3,nState)
      call SqToTri_Q(Snt3,Int2,nState)
      BigT(kaunt,:) = Int2
    end do
  end do
  call mma_deallocate(BigV)
  call mma_deallocate(Int1)
  call mma_deallocate(Int2)
  call mma_deallocate(Snt1)
  call mma_deallocate(Snt2)
  call mma_deallocate(Snt3)

else
  ! Here we go if both TDM's can not be put in memory. Might be a bit
  ! slow due to its nested nature with repeated IO.

  call mma_allocate(BigT,nSize,nTriS,label='ALLES')
  call mma_allocate(AOG,nSize,label='AOGamma')
  BigT(:,:) = Zero
  do iiS=1,nStatePrim
    do jjS=1,nStatePrim
      indypop = iTri(iiS,jjS)
      iDisk = iTocBig(indypop)
      call dDaFile(LuIn,2,AOG,nSize,iDisk)
      kaunter = 0
      do iB=1,nB
        do jB=1,iB
          kaunter = kaunter+1
          do iS=1,nState
            do jS=1,iS
              indx = iTri(iS,jS)
              BigT(kaunter,indx) = BigT(kaunter,indx)+Eig2(iiS,iS)*Eig2(jjS,jS)*AOG(kaunter)
            end do
          end do
        end do
      end do
    end do
  end do
  call mma_deallocate(AOG)
end if

! Deallocations and finish up.

call mma_deallocate(iTocBig)
call DaClos(LuIn)

return

end subroutine RasRasTrans
