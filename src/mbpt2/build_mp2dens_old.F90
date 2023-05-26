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

subroutine Build_Mp2Dens_Old(TriDens,Density,CMO,mSym,nOrbAll,Diagonalize)

#include "intent.fh"

use Data_Structures, only: DSBA_Type
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(_OUT_) :: TriDens(*)
type(DSBA_Type), intent(in) :: Density
integer(kind=iwp), intent(in) :: mSym, nOrbAll(8)
real(kind=wp), intent(in) :: CMO(*)
logical(kind=iwp), intent(in) :: Diagonalize
integer(kind=iwp) :: i, idx, ipSymLin(8), ipSymRec(8), ipSymTri(8), iSym, iUHF, j, lRecTot, LuMP2, nOrbAllMax, nOrbAllTot
character(len=30) :: Note
integer(kind=iwp), allocatable :: IndT(:,:)
real(kind=wp), allocatable :: AORecBlock(:), AOTriBlock(:), EigenValBlock(:), EigenValTot(:), EigenVecBlock(:), EigenVecTot(:), &
                              Energies(:), MOTriBlock(:), TmpRecBlock(:)
integer(kind=iwp), external :: IsFreeUnit
#include "corbinf.fh"

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
call mma_allocate(AORecBlock,nOrbAllMax**2,label='AORecBlock')
call mma_allocate(TmpRecBlock,nOrbAllMax**2,label='TmpRecBlock')
call mma_allocate(AOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2,label='AOTriBlock')
if (Diagonalize) then
  call mma_allocate(MOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2,label='MOTriBlock')
  call mma_allocate(EigenVecBlock,nOrbAllMax**2,label='EigenVecBlock')
  call mma_allocate(EigenValBlock,nOrbAllMax,label='EigenValBlock')
  call mma_allocate(EigenVecTot,lRecTot,label='EigenVectors')
  call mma_allocate(EigenValTot,nOrbAllTot,label='EigenValues')
  call mma_allocate(Energies,nOrbAllTot,label='Energies')
  call mma_allocate(IndT,7,mSym,label='IndT')
  Energies(:) = Zero
end if
AORecBlock(:) = Zero
TmpRecBlock(:) = Zero
AOTriBlock(:) = Zero

! Setup a pointer to a symmetryblock in rect or tri representation.
ipSymRec(1) = 0
ipSymTri(1) = 1
ipSymLin(1) = 0
do iSym=2,8
  ipSymTri(iSym) = ipSymTri(iSym-1)+(nOrbAll(iSym-1))*(nOrbAll(iSym-1)+1)/2
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
      do i=1,nOrbAll(iSym)**2
        EigenVecBlock(i) = CMO(ipSymRec(iSym)+i)
      end do
    end if
    ! Transform the symmetryblock to AO-density
    call DGEMM_('N','N',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),One,CMO(ipSymRec(iSym)+1),nOrbAll(iSym),Density%SB(iSym)%A1, &
                nOrbAll(iSym),Zero,TmpRecBlock,nOrbAll(iSym))
    call DGEMM_('N','T',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),One,TmpRecBlock,nOrbAll(iSym),CMO(ipSymRec(iSym)+1), &
                nOrbAll(iSym),Zero,AORecBlock,nOrbAll(iSym))
    !call RecPrt('AODens:','(20F8.5)',AORecBlock,nOrb(iSym),nOrb(iSym))
    call Fold_Mat(1,nOrbAll(iSym),AORecBlock,AOTriBlock)
    call dcopy_(nOrbAll(iSym)*(nOrbAll(iSym)+1)/2,AOTriBlock,1,TriDens(ipSymTri(iSym)),1)

    if (Diagonalize) then
      ! Make a normal folded matrix

      idx = 1
      do i=1,nOrbAll(iSym)
        do j=1,i
          MOTriBlock(idx) = Density%SB(iSym)%A2(j,i)
          idx = idx+1
        end do
      end do

      call NIDiag(MOTriBlock,EigenVecBlock,nOrbAll(iSym),nOrbAll(iSym))

      do i=1,nOrbAll(iSym)
        EigenValBlock(i) = MOTriBlock(i*(i+1)/2)
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

      call dcopy_(nOrbAll(iSym)**2,EigenVecBlock,1,EigenVecTot(ipSymRec(iSym)+1),1)
      call dcopy_(nOrbAll(iSym),EigenValBlock,1,EigenValTot(ipSymLin(iSym)+1),1)

    end if

  end if
end do

! Put the MP2 natural orbitals we just produced on disk to the file MP2ORB.

if (Diagonalize) then
  LuMP2 = 50
  LuMP2 = IsFreeUnit(LuMP2)
  ! Build the TypeIndex array
  do iSym=1,mSym
    IndT(1,iSym) = nFro(iSym)
    IndT(2,iSym) = nOcc(iSym)
    IndT(3,iSym) = 0
    IndT(4,iSym) = 0
    IndT(5,iSym) = 0
    IndT(6,iSym) = nOrb(iSym)-nFro(iSym)-nOcc(iSym)-nDel(iSym)
    IndT(7,iSym) = nDel(iSym)
  end do
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

end subroutine Build_Mp2Dens_Old
