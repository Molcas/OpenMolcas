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

subroutine Magic_Bonds(Coor,nAtoms,iTabBonds,nBondMax,nBonds,iTabAtoms,nMax)

use Slapaf_Info, only: Covalent_Bond, Magic_Bond, vdW_Bond
use Constants, only: deg2rad
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nBondMax, nMax
real(kind=wp), intent(in) :: Coor(3,nAtoms)
integer(kind=iwp), intent(inout) :: iTabBonds(3,nBondMax), nBonds, iTabAtoms(2,0:nMax,nAtoms)
integer(kind=iwp) :: iAtom, iBond, iBondType, iCase, jAtom, jBond, jBondType, kAtom, kBond, kNeighbor, nNeighbor_i, nNeighbor_j, &
                     nNeighbor_k
real(kind=wp) :: CosFi, CosFi_limit, Fi_limit, rij, rkj, xij, xkj, yij, ykj, zij, zkj

!                                                                      *
!***********************************************************************
!                                                                      *
!define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *

Fi_limit = 165.0_wp*deg2rad
CosFi_limit = cos(Fi_limit)
do iBond=1,nBonds
  iBondType = iTabBonds(3,iBond)
  if (iBondType == vdW_Bond) cycle
  if (iBondType > Magic_Bond) cycle

  do iCase=1,2

    if (iCase == 1) then
      iAtom = iTabBonds(1,iBond)
      jAtom = iTabBonds(2,iBond)
    else
      iAtom = iTabBonds(2,iBond)
      jAtom = iTabBonds(1,iBond)
    end if

    xij = Coor(1,iAtom)-Coor(1,jAtom)
    yij = Coor(2,iAtom)-Coor(2,jAtom)
    zij = Coor(3,iAtom)-Coor(3,jAtom)
    rij = sqrt(xij**2+yij**2+zij**2)

    nNeighbor_j = iTabAtoms(1,0,jAtom)
    neighbor: do kNeighbor=1,nNeighbor_j
      kAtom = iTabAtoms(1,kNeighbor,jAtom)
      jBond = iTabAtoms(2,kNeighbor,jAtom)
      if (jBond >= iBond) cycle neighbor
      jBondType = iTabBonds(3,jBond)
      if (jBondType == vdW_Bond) cycle neighbor
      if (jBondType > Magic_Bond) cycle neighbor

      ! If this is close to a linear system we will form
      ! a magic bond between iAtom and kAtom

      xkj = Coor(1,kAtom)-Coor(1,jAtom)
      ykj = Coor(2,kAtom)-Coor(2,jAtom)
      zkj = Coor(3,kAtom)-Coor(3,jAtom)
      rkj = sqrt(xkj**2+ykj**2+zkj**2)

      CosFi = (xij*xkj+yij*ykj+zij*zkj)/(rij*rkj)
#     ifdef _DEBUGPRINT_
      write(u6,*) 'iAtom,jAtom,kAtom=',iAtom,jAtom,kAtom
      write(u6,*) 'CosFi,CosFi_limit=',CosFi,CosFi_limit
#     endif

      ! Observe that this limit should be coordinated with
      ! the limit for excluding torsions in torsion_list!

      if (CosFi <= CosFi_Limit) then

#       ifdef _DEBUGPRINT_
        write(u6,*) 'Forming a "magic" bond'
        write(u6,*) 'iAtom,kAtom=',iAtom,kAtom
#       endif
        ! Double check that this bond is not included already.
        ! If that is the case just update the bond type if it
        ! is not classified already as a covalent bond.

        do kBond=1,nBonds
          if (((iTabBonds(1,kBond) == iAtom) .and. (iTabBonds(2,kBond) == kAtom)) .or. &
              ((iTabBonds(1,kBond) == kAtom) .and. (iTabBonds(2,kBond) == iAtom)) .and. (iTabBonds(3,kBond) /= Covalent_Bond)) then
            iTabBonds(3,kBond) = jAtom+Magic_Bond
            cycle neighbor
          end if
        end do

        if (nBonds+1 > nBondMax) then
          call WarningMessage(2,'Error in Magic_Bonds')
          write(u6,*) 'Magic_Bonds: nBonds > nBondMax'
          write(u6,*) 'iTabBonds:'
          do kBond=1,nBonds
            write(u6,*)
            write(u6,*) 'kBond=',kBond
            write(u6,*)
            write(u6,*) 'Atoms=',iTabBonds(1,kBond),iTabBonds(2,kBond)
            write(u6,*) 'Bondtype:',iTabBonds(3,kBond)
          end do
          call Abend()
        end if
        nBonds = nBonds+1
        iTabBonds(1,nBonds) = max(iAtom,kAtom)
        iTabBonds(2,nBonds) = min(iAtom,kAtom)
        iTabBonds(3,nBonds) = jAtom+Magic_Bond

        nNeighbor_i = iTabAtoms(1,0,iAtom)+1
        if (nNeighbor_i > nMax) then
          call WarningMessage(2,'Error in Magic_Bonds')
          write(u6,*) 'Magic_Bonds: nNeighbor_i > nMax'
          write(u6,*) 'iAtom=',iAtom
          write(u6,*) 'nNeighbor_i=',nNeighbor_i
          write(u6,*) 'nMax=',nMax
          call Abend()
        end if
        iTabAtoms(1,0,iAtom) = nNeighbor_i
        iTabAtoms(1,nNeighbor_i,iAtom) = kAtom
        iTabAtoms(2,nNeighbor_i,iAtom) = nBonds

        nNeighbor_k = iTabAtoms(1,0,kAtom)+1
        if (nNeighbor_k > nMax) then
          call WarningMessage(2,'Error in Magic_Bonds')
          write(u6,*) 'Magic_Bonds: nNeighbor_k > nMax'
          write(u6,*) 'kAtom=',kAtom
          write(u6,*) 'nNeighbor_k=',nNeighbor_k
          write(u6,*) 'nMax=',nMax
          call Abend()
        end if
        iTabAtoms(1,0,kAtom) = nNeighbor_k
        iTabAtoms(1,nNeighbor_k,kAtom) = iAtom
        iTabAtoms(2,nNeighbor_k,kAtom) = nBonds
      end if

    end do neighbor  ! kNeighbor

  end do     ! iCase

end do       ! iBond

return

end subroutine Magic_Bonds
