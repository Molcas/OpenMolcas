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

! This subroutine should be in a module, to avoid explicit interfaces
#ifndef _IN_MODULE_
#error "This file must be compiled inside a module"
#endif

subroutine Build_Mp2Dens(TriDens,nTriDens,MP2X_e,CMO,mSym,nOrbAll,Diagonalize)

use Index_Functions, only: iTri, nTri_Elem
use Data_Structures, only: V2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nTriDens, mSym, nOrbAll(8)
real(kind=wp), intent(inout) :: TriDens(nTriDens)
type(V2), intent(in) :: MP2X_e(8)
real(kind=wp), intent(in) :: CMO(*)
logical(kind=iwp), intent(in) :: Diagonalize
#include "corbinf.fh"
integer(kind=iwp) :: i, iLen, indx, ipSymLin(8), ipSymRec(8), ipSymTri(8), iSym, iUHF, lRecTot, LuMP2, nOrbAllMax, nOrbAllTot
character(len=30) :: Note
real(kind=wp), allocatable :: AORecBlock(:), AOTriBlock(:), EigenValBlock(:), EigenValTot(:), EigenVecBlock(:), EigenVecTot(:), &
                              Energies(:), MOTriBlock(:), TmpRecBlock(:)
integer(kind=iwp), allocatable :: IndT(:,:)
integer(kind=iwp), external :: IsFreeUnit

!                                                                      *
!***********************************************************************
!                                                                      *
nOrbAllTot = nOrbAll(1)
nOrbAllMax = nOrbAll(1)
lRecTot = nOrbAll(1)*nOrbAll(1)
do iSym=2,mSym
  nOrbAllTot = nOrbAllTot+nOrbAll(iSym)
  nOrbAllMax = max(nOrbAllMax,nOrbAll(iSym))
  lRecTot = lRecTot+nOrbAll(iSym)*nOrbAll(iSym)
end do

! A blockmatrix of the size Orb(iSym) X Orb(iSym) is
! allocated and set to zero
call mma_allocate(AORecBlock,nOrbAllMax**2,Label='AORecBlock')
call mma_allocate(TmpRecBlock,nOrbAllMax**2,Label='TmpRecBlock')
call mma_allocate(AOTriBlock,nTri_Elem(nOrbAllMax),Label='AOTriBlock')

if (Diagonalize) then
  call mma_allocate(MOTriBlock,nTri_Elem(nOrbAllMax),Label='MOTriBlock')
  call mma_allocate(EigenVecBlock,nOrbAllMax**2,Label='EigenVecBlock')
  call mma_allocate(EigenValBlock,nOrbAllMax,Label='EigenValBlock')
  call mma_allocate(EigenVecTot,lRecTot,Label='EigenVecTot')
  call mma_allocate(EigenValTot,nOrbAllTot,Label='EigenValTot')
  call mma_allocate(Energies,nOrbAllTot,Label='Energies')
  call mma_allocate(IndT,7,mSym,Label='IndT')
  Energies(:) = Zero
end if

AORecBlock(:) = Zero
TmpRecBlock(:) = Zero
AOTriBlock(:) = Zero

! Setup a pointer to a symmetryblock in rect or tri representation.
ipSymRec(1) = 0
ipSymTri(1) = 0
ipSymLin(1) = 0
do iSym=2,8
  ipSymTri(iSym) = ipSymTri(iSym-1)+nTri_Elem(nOrbAll(iSym-1))
  ipSymRec(iSym) = ipSymRec(iSym-1)+(nOrbAll(iSym-1))**2
  ipSymLin(iSym) = ipSymLin(iSym-1)+(nOrbAll(iSym-1))
end do

do iSym=1,mSym

  if (nOrbAll(iSym) /= 0) then

    ! Set initial values for the eigenvectors to the CMO-matrix.
    ! If we transform the CMO-matrix in the same way as the MO-matrix
    ! needs to be treated to be diagonal we get a transformation from
    ! AO to canonical basis. Also sets the energies to zero since they
    ! have no physical relevance in this basis.
    if (Diagonalize) then
      EigenVecBlock(1:nOrbAll(iSym)**2) = CMO(ipSymRec(iSym)+1:ipSymRec(iSym)+nOrbAll(iSym)**2)
    end if
    ! Transform the symmetryblock to AO-density
    call DGEMM_('N','N',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),One,CMO(ipSymRec(iSym)+1),nOrbAll(iSym),MP2X_e(iSym)%A, &
                nOrbAll(iSym),Zero,TmpRecBlock,nOrbAll(iSym))
    call DGEMM_('N','T',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),One,TmpRecBlock,nOrbAll(iSym),CMO(ipSymRec(iSym)+1), &
                nOrbAll(iSym),Zero,AORecBlock,nOrbAll(iSym))
    !call RecPrt('AODens:','(20F8.5)',AORecBlock,nOrb(iSym),nOrb(iSym))
    call Fold_Mat(1,nOrbAll(iSym),AORecBlock,AOTriBlock)
    iLen = nTri_Elem(nOrbAll(iSym))
    TriDens(ipSymTri(iSym)+1:ipSymTri(iSym)+iLen) = AOTriBlock(1:iLen)

    if (Diagonalize) then
      ! Make a normal folded matrix

      indx = 0
      do i=1,nOrbAll(iSym)
        MOTriBlock(indx+1:indx+i) = MP2X_e(iSym)%A(1:i,i)
        indx = indx+i
      end do

      call NIDiag(MOTriBlock,EigenVecBlock,nOrbAll(iSym),nOrbAll(iSym))

      do i=1,nOrbAll(iSym)
        EigenValBlock(i) = MOTriBlock(iTri(i,i))
      end do

      call SortEig(EigenValBlock,EigenVecBlock,nOrbAll(iSym),nOrbAll(iSym),-1,.false.)

#     ifdef _DEBUGPRINT_
      write(u6,*) 'The sorted eigenvalues are'
      do i=1,nOrbAll(iSym)
        write(u6,*) EigenValBlock(i)
      end do
      write(u6,*) 'Eigenvectors sorted'
      do i=1,nOrbAll(iSym)**2
        write(u6,*) EigenVecBlock(i)
      end do
#     endif

      iLen = nOrbAll(iSym)**2
      EigenVecTot(ipSymRec(iSym)+1:ipSymRec(iSym)+iLen) = EigenVecBlock(1:iLen)
      iLen = nOrbAll(iSym)
      EigenValTot(ipSymLin(iSym)+1:ipSymLin(iSym)+iLen) = EigenValBlock(1:iLen)

    end if

  end if
end do

! Put the MP2 natural orbitals we just produced on disk to the file MP2ORB.

if (Diagonalize) then
  LuMP2 = 50
  LuMP2 = IsFreeUnit(LuMP2)
  ! Build the TypeIndex array
  IndT(1,1:mSym) = nFro(1:mSym)
  IndT(2,1:mSym) = nOcc(1:mSym)
  IndT(3,1:mSym) = 0
  IndT(4,1:mSym) = 0
  IndT(5,1:mSym) = 0
  IndT(6,1:mSym) = nOrb(1:mSym)-nFro(1:mSym)-nOcc(1:mSym)-nDel(1:mSym)
  IndT(7,1:mSym) = nDel(1:mSym)
  Note = '*  Natural MP2 orbitals'
  call WrVec('MP2ORB',LuMP2,'COEI',mSym,nOrbAll,nOrbAll,EigenVecTot,EigenValTot,Energies,IndT,Note)
  ! Create a molden-file
  iUHF = 0
  call Molden_Interface(iUHF,'MP2ORB','MD_MP2')

end if

call mma_deallocate(AORecBlock)
call mma_deallocate(TmpRecBlock)
call mma_deallocate(AOTriBlock)

if (Diagonalize) then
  call mma_deallocate(MOTriBlock)
  call mma_deallocate(EigenVecBlock)
  call mma_deallocate(EigenValBlock)
  call mma_deallocate(EigenVecTot)
  call mma_deallocate(EigenValTot)
  call mma_deallocate(Energies)
  call mma_deallocate(IndT)
end if

return

end subroutine Build_Mp2Dens
