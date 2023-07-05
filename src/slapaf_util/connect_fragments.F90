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

subroutine Connect_Fragments(nAtoms,iTabBonds,nBondMax,nBonds,Coor,iTabAtoms,nMax,iANr)

use slapaf_parameters, only: rFuzz

implicit real*8(a-h,o-z)
#include "stdalloc.fh"
real*8 Coor(3,nAtoms)
integer iTabAtoms(2,0:nMax,nAtoms), iANr(nAtoms), iTabBonds(3,nBondMax)
integer, allocatable, dimension(:) :: iStack
integer, allocatable :: nSet(:)
real*8, allocatable :: SetDist(:)
#include "bondtypes.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(nSet,nAtoms,Label='nSet')
dR_Thr = rFuzz
HH_Thr = 5.0d0

call mma_allocate(iStack,nAtoms,label='iStack')

do
  Not_Defined = -1
  nSet(:) = Not_Defined

  nStack = 1
  iAtom = 1
  iStack(nStack) = iAtom
  iSet = 1
  nSet(iAtom) = iSet

# ifdef _DEBUGPRINT_
  write(6,*) '>>>> Connect Fragments <<<<<'
  call RecPrt('Coor',' ',Coor,3,nAtoms)
# endif
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  outer: do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! Pick up an atom from the top of the stack.

    iAtom = iStack(nStack)
#   ifdef _DEBUGPRINT_
    write(6,*)
    write(6,*) 'iAtom=',iAtom
    write(6,*) 'iSet=',iSet
    write(6,*) 'nStack=',nStack
    write(6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
#   endif

    ! Loop over all atoms to which this atom has a bond

    ! Decrement the stack size one step

    nStack = nStack-1
    do iNeighbor=1,iTabAtoms(1,0,iAtom)
      jAtom = iTabAtoms(1,iNeighbor,iAtom)
      iBond = iTabAtoms(2,iNeighbor,iAtom)
#     ifdef _DEBUGPRINT_
      write(6,*) 'jAtom=',jAtom
      write(6,*) 'iBond=',iBond
      write(6,*) 'iTabBonds(3,iBond)=',iTabBonds(3,iBond)
      write(6,*) 'nSet(jAtom)=',nSet(jAtom)
#     endif
      ! Skip if vdW bond.
      if (iTabBonds(3,iBond) == vdW_Bond) cycle
      if (nSet(jAtom) == Not_Defined) then

        ! If not defined add to the stack
        nStack = nStack+1
        iStack(nStack) = jAtom
        ! Assign which set the atom belongs to.
        nSet(jAtom) = iSet
      else if (nSet(jAtom) /= iSet) then

        ! Error if atom belongs to the wrong set.
        call WarningMessage(2,' Error in Connect_Fragments')
        write(6,*) 'Connect_Fragments:'
        write(6,*) 'Inconsistent set indices!'
        write(6,*) 'iSet,jSet=',iSet,nSet(jAtom)
        call Abend()
      end if
    end do
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! If there are more atoms in the stack process the next one.

    if (nStack /= 0) cycle outer
    !                                                                  *
    !*******************************************************************
    !                                                                  *
    ! If all atoms now belong to a set we are done. Otherwise
    ! we have a new set.

    ! Loop over all atoms. If an atom is found which is not assigned
    ! to a set add it to the stack.

    do iAtom=1,nAtoms
      if (nSet(iAtom) == Not_Defined) then
        ! Add atom to the stack
        nStack = 1
        iStack(nStack) = iAtom
        ! Set the new set index and assign it to the atom.
        iSet = iSet+1
        nSet(iAtom) = iSet
        exit
      end if
      if (iAtom == nAtoms) exit outer
    end do
  end do outer
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Now all atoms belong to a set
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! If all atoms belong to the same set then we are done.

  if (iSet == 1) exit
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Now we have to connect the disjoint sets.
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Find the shortest distance between every two sets

  call mma_allocate(SetDist,iSet*iSet,Label='SetDist')
  SetDist(:) = 1.0d6
  do kAtom=1,nAtoms
    kSet = nSet(kAtom)
    do lAtom=kAtom+1,nAtoms
      lSet = nSet(lAtom)
      if (lSet == kSet) cycle

      ! For l < k:
      ! In (k,l), the minimum distance between the sets k and l
      ! In (l,k), the minimum distance excluding H-H

      iOff = (min(kSet,lSet)-1)*iSet+max(kSet,lSet)
      iOffHH = (max(kSet,lSet)-1)*iSet+min(kSet,lSet)
      RTest = (Coor(1,kAtom)-Coor(1,lAtom))**2+(Coor(2,kAtom)-Coor(2,lAtom))**2+(Coor(3,kAtom)-Coor(3,lAtom))**2
      if ((iAnr(kAtom) /= 1) .or. (iANr(lAtom) /= 1)) SetDist(iOff) = min(SetDist(iOff),RTest)
      SetDist(iOffHH) = min(SetDist(iOffHH),RTest)
    end do
  end do

  ! Find the shortest distance between any two sets.

  rShort = 1.0d5
  do kSet=1,iSet
    do lSet=kSet+1,iSet

      ! H-H distances are only considered if the minimum distance
      ! (excluding H-H) is greater than HH_Thr
      iOff = (kSet-1)*iSet+lSet
      if (SetDist(iOff) > HH_Thr**2) iOff = (lSet-1)*iSet+kSet-1
      if (SetDist(iOff) < rShort) then
        rShort = SetDist(iOff)
      end if
    end do
  end do

  ! Add bonds between sets if they are shorter than rShort+dR_Thr

  do kAtom=1,nAtoms
    kSet = nSet(kAtom)
    inner: do lAtom=kAtom+1,nAtoms
      lSet = nSet(lAtom)
      if (lSet == kSet) cycle inner

      ! Again, add H-H bonds only if necessary
      iOff = (min(kSet,lSet)-1)*iSet+max(kSet,lSet)
      if ((SetDist(iOff) <= HH_Thr**2) .and. ((iANr(kAtom) == 1) .and. (iANr(lAtom) == 1))) cycle inner

      RTest = (Coor(1,kAtom)-Coor(1,lAtom))**2+(Coor(2,kAtom)-Coor(2,lAtom))**2+(Coor(3,kAtom)-Coor(3,lAtom))**2
      RTest = sqrt(RTest)
      if (RTest <= sqrt(rShort)+dR_Thr) then
        do iBond=1,nBonds

          ! Look through the bond list and find if it is there.
          if (((iTabBonds(1,iBond) == kAtom) .and. (iTabBonds(2,iBond) == lAtom)) .or. &
              ((iTabBonds(1,iBond) == lAtom) .and. (iTabBonds(2,iBond) == kAtom))) then
            iTabBonds(3,iBond) = Fragments_Bond

            ! If the bond is already in the list that is it.
            cycle inner
          end if
        end do

        ! Add the bond to the bond list
        if (nBonds+1 > nBondMax) then
          call WarningMessage(2,' Error in Connect_Fragments')
          write(6,*) 'Connect_Fragments: nBonds+1 > nBondMax'
          call Abend()
        end if

        nBonds = nBonds+1
        iTabBonds(1,nBonds) = lAtom
        iTabBonds(2,nBonds) = kAtom
        iTabBonds(3,nBonds) = Fragments_Bond

        ! Update atoms list

        nNeighbor_k = iTabAtoms(1,0,kAtom)+1
        if (nNeighbor_k > nMax) then
          call WarningMessage(2,' Error in Connect_Fragments')
          write(6,*) 'Connect_Fragments: nNeighbor_k > nMax'
          write(6,*) 'kAtom=',kAtom
          write(6,*) 'nNeighbor_k=',nNeighbor_k
          write(6,*) 'nMax=',nMax
          call Abend()
        end if
        iTabAtoms(1,0,kAtom) = nNeighbor_k
        iTabAtoms(1,nNeighbor_k,kAtom) = lAtom
        iTabAtoms(2,nNeighbor_k,kAtom) = nBonds

        nNeighbor_l = iTabAtoms(1,0,lAtom)+1
        if (nNeighbor_l > nMax) then
          call WarningMessage(2,' Error in Connect_Fragments')
          write(6,*) 'Connect_Fragments: nNeighbor_l > nMax'
          write(6,*) 'lAtom=',lAtom
          write(6,*) 'nNeighbor_l=',nNeighbor_l
          write(6,*) 'nMax=',nMax
          call Abend()
        end if
        iTabAtoms(1,0,lAtom) = nNeighbor_l
        iTabAtoms(1,nNeighbor_l,lAtom) = kAtom
        iTabAtoms(2,nNeighbor_l,lAtom) = nBonds

      end if
    end do inner
  end do
  call mma_deallocate(SetDist)

  ! Try again!
end do

call mma_deallocate(iStack)
call mma_deallocate(nSet)
!                                                                      *
!***********************************************************************
!                                                                      *

end subroutine Connect_Fragments
