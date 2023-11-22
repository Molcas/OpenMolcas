!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2006, Thomas Bondo Pedersen                            *
!***********************************************************************

subroutine Cho_Stat_ParentDiag()
!
! Thomas Bondo Pedersen, February 2006.
!
! Purpose: print statistics about the diagonals from which the
!          Cholesky vectors are computed. I.e., whether they are
!          one-center or two-center diagonals. Does not work with
!          symmetry!

use Cholesky, only: InfVec, LuPri, nBasT, nnBstRT, nSym, NumCho, NumChT
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), parameter :: numAt = 6
integer(kind=iwp) :: i, i1, i2, iA, iAt, iAt0, iAt1, iAt2, iAtA, iAtB, iAtom, iB, iBatch, iClass(4), iRS1, iVec, nAt, nAtom, &
                     nBatch, nChk, nPseudo, nTot1, nTot2
real(kind=wp) :: R, Ratio(numAt), Rave, RClass(3), Rmax, Rmin
logical(kind=iwp) :: Debug
character(len=LenIn8) :: AtomLabel(MxBas)
integer(kind=iwp), allocatable :: iBF2Atom(:), mapRS2F(:,:), nBas_per_Atom(:), nBas_Start(:), nPC1(:)
real(kind=wp), allocatable :: Coord(:,:), RC2(:)
character(len=*), parameter :: SecNam = 'Cho_Stat_ParentDiag'

! Set debug flag.
! ---------------

#ifdef _DEBUGPRINT_
Debug = .true.
#else
Debug = .false.
#endif

! Symmetry not allowed.
! ---------------------

if (nSym /= 1) return

! Check if anything to do.
! ------------------------

if (NumCho(1) /= NumChT) call SysAbendMsg(SecNam,'NumChT /= NumCho(1)',' ')
if (NumChT < 1) return

! Get number of atoms.
! --------------------

!call Get_nAtoms_All(nAtom)
call Get_iScalar('Bfn Atoms',nAtom)
call Get_iScalar('Pseudo atoms',nPseudo)

! Get atomic labels and basis function labels.
! --------------------------------------------

call Get_cArray('Unique Basis Names',AtomLabel,LenIn8*nBasT)

! Allocate and get index arrays for indexation of basis functions on
! each atom.
! ------------------------------------------------------------------

call mma_allocate(nBas_per_Atom,nAtom,Label='nBas_per_Atom')
call mma_allocate(nBas_Start,nAtom,Label='nBas_Start')
call BasFun_Atom(nBas_per_Atom,nBas_Start,AtomLabel,nBasT,nAtom,Debug)

! Allocate and get nuclear coordinates.
! -------------------------------------

call mma_allocate(Coord,3,nAtom,Label='Coord')
!call Get_dArray('Unique Coordinates',Coord,3*nAtom)
call Get_dArray('Bfn Coordinates',Coord,3*nAtom)

! Allocate and set mapping from basis function to atom.
! -----------------------------------------------------

call mma_Allocate(iBF2Atom,nBasT,Label='iBF2Atom')
do iAtom=1,nAtom
  i1 = nBas_Start(iAtom)
  i2 = i1+nBas_per_Atom(iAtom)-1
  do i=i1,i2
    iBF2Atom(i) = iAtom
  end do
end do

if (Debug) then
  write(Lupri,*)
  write(Lupri,*) SecNam,': mapping from basis function to atom:'
  do i=1,nBasT
    write(Lupri,*) 'Basis function ',i,' is centered on atom ',iBF2Atom(i),' labeled ',AtomLabel(i)(1:LenIn)
  end do
end if

! Allocate and set mapping from 1st reduced set to full storage.
! --------------------------------------------------------------

call mma_allocate(mapRS2F,2,nnBstRT(1),Label='mapRS2F')
call Cho_RStoF(mapRS2F,2,nnBstRT(1),1)

! Allocate counters.
! ------------------

call mma_allocate(nPC1,nAtom,Label='nPC1')
call mma_allocate(RC2,NumChT,Label='RC2')

! Compute statistics.
! -------------------

Rmin = 1.0e15_wp
Rmax = -1.0e15_wp
Rave = Zero
nPC1(:) = 0
RC2(:) = Zero
do iVec=1,NumChT
  iRS1 = InfVec(iVec,1,1)
  iA = mapRS2F(1,iRS1)
  iB = mapRS2F(2,IRS1)
  iAtA = iBF2Atom(iA)
  iAtB = iBF2Atom(iB)
  if (iAtA == iAtB) then
    nPC1(iAtA) = nPC1(iAtA)+1
  else
    R = sqrt((Coord(1,iAtA)-Coord(1,iAtB))**2+(Coord(2,iAtA)-Coord(2,iAtB))**2+(Coord(3,iAtA)-Coord(3,iAtB))**2)
    RC2(iVec) = R
    Rmin = min(Rmin,R)
    Rmax = max(Rmax,R)
    Rave = Rave+R
  end if
end do
nTot1 = sum(nPC1(:))
nTot2 = NumChT-nTot1
if (nTot2 > 0) then
  Rave = Rave/real(nTot2,kind=wp)
else if (nTot2 < 0) then
  call SysAbendMsg(SecNam,'Setup error!','nTot2 < 0')
end if

! Print overall statisctics.
! --------------------------

write(Lupri,'(//,2X,A,/,2X,A)') 'Parent Diagonals','----------------'
write(Lupri,'(/,A,I9,A,F7.2,A)') 'Number of vectors from 1-center diagonals:',nTot1,' (', &
                                 1.0e2_wp*real(nTot1,kind=wp)/real(NumChT,kind=wp),'%)'
write(Lupri,'(A,I9,A,F7.2,A)') 'Number of vectors from 2-center diagonals:',nTot2,' (', &
                               1.0e2_wp*real(nTot2,kind=wp)/real(NumChT,kind=wp),'%)'

! Print statistics for 1-center vectors.
! --------------------------------------

write(Lupri,'(/,1X,A)') 'Vectors from 1-center diagonals:'
nBatch = (nAtom-nPseudo-1)/numAt+1 ! exclude pseudo-atoms
do iBatch=1,nBatch
  if (iBatch == nBatch) then
    ! exclude pseudo-atoms
    nAt = nAtom-nPseudo-numAt*(nBatch-1)
  else
    nAt = numAt
  end if
  iAt0 = numAt*(iBatch-1)
  iAt1 = iAt0+1
  iAt2 = iAt0+nAt
  do i=1,nAt
    iAt = iAt0+i
    if (nBas_per_Atom(iAt) < 1) then
      if (nPC1(iAt) > 0) call SysAbendMsg(SecNam,'No basis functions, but >0 vectors !?!?',' ')
      Ratio(i) = Zero
    else
      Ratio(i) = real(nPC1(iAt),kind=wp)/real(nBas_per_Atom(iAt),kind=wp)
    end if
  end do
  write(Lupri,'(/,A,6(6X,A))') 'Label              ',(AtomLabel(nBas_Start(i))(1:LenIn),i=iAt1,iAt2)
  write(Lupri,'(A,6(1X,I9))') 'Center no.         ',(i,i=iAt1,iAt2)
  write(Lupri,'(A,6(1X,I9))') 'Vectors (M)        ',(nPC1(i),i=iAt1,iAt2)
  write(Lupri,'(A,6(1X,I9))') 'Basis functions (N)',(nBas_per_Atom(i),i=iAt1,iAt2)
  write(Lupri,'(A,6(1X,F9.2))') 'Ratio (M/N)        ',(Ratio(i),i=1,nAt)
  if (iBatch /= nBatch) write(Lupri,*)
end do

! Print statistics for 2-center vectors.
! --------------------------------------

if (nTot2 > 0) then
  write(Lupri,'(/,1X,A)') 'Vectors from 2-center diagonals:'
  RClass(1) = Rave-(Rave-Rmin)*Half
  RClass(2) = Rave
  RClass(3) = Rave+(Rmax-Rave)*Half
  iClass(:) = 0
  nChk = 0
  do iVec=1,NumChT
    if (RC2(iVec) > Zero) then
      nChk = nChk+1
      if (RC2(iVec) <= RClass(1)) then
        iClass(1) = iClass(1)+1
      else if (RC2(iVec) <= RClass(2)) then
        iClass(2) = iClass(2)+1
      else if (RC2(iVec) <= RClass(3)) then
        iClass(3) = iClass(3)+1
      else
        iClass(4) = iClass(4)+1
      end if
    end if
  end do
  if (nChk /= nTot2) then
    ! FAQ
    ! --- Replace the call to SysAbendMsg with just a warning.
    !
    !call SysAbendMsg(SecNam,'nChk /= nTot2',' ')
    !
    ! --- TODO/FIX  figure out if with ghost atoms is just a mistmatch
    !               or there is really a bug

    write(u6,*) SecNam//': Warning! (nChk /= nTot2); could be due to the presence of ghost atoms.'
  end if
  write(Lupri,'(/,A,3ES15.5)') 'Min, average, and max center distance: ',Rmin,Rave,Rmax
  write(Lupri,'(A,ES12.2,A,I9)') '#vectors with center distance                R <= ',RClass(1),': ',iClass(1)
  write(Lupri,'(A,ES12.2,A,ES12.2,A,I9)') '#vectors with center distance ',RClass(1),' < R <= ',RClass(2),': ',iClass(2)
  write(Lupri,'(A,ES12.2,A,ES12.2,A,I9)') '#vectors with center distance ',RClass(2),' < R <= ',RClass(3),': ',iClass(3)
  write(Lupri,'(A,ES12.2,A,12X,A,I9)') '#vectors with center distance ',RClass(3),' < R    ',': ',iClass(4)
end if

! Deallocations.
! --------------

call mma_deallocate(RC2)
call mma_deallocate(nPC1)
call mma_deallocate(mapRS2F)
call mma_deallocate(iBF2Atom)
call mma_deallocate(Coord)
call mma_deallocate(nBas_Start)
call mma_deallocate(nBas_per_Atom)

end subroutine Cho_Stat_ParentDiag
