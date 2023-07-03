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

subroutine Find_Bonds(Coor,nAtoms,iTab,nMax,nx,ny,nz,iBox,iANr,iTabBonds,nBonds,nBondMax,iTabAtoms,ThrB)

implicit real*8(a-h,o-z)
real*8 Coor(3,nAtoms)
integer iTab(0:nMax,nx,ny,nz), iBox(3,nAtoms), iANr(nAtoms), iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Initiate

#ifdef _DEBUGPRINT_
call RecPrt('Coor',' ',Coor,3,nAtoms)
write(6,*) 'Find_Bonds: ThrB=',ThrB
write(6,*) 'Initialize iTabAtoms'
#endif

call ICopy(nBondMax*3,[0],0,iTabBonds,1)
call ICopy(2*(nMax+1)*nAtoms,[0],0,iTabAtoms,1)
nBonds = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Find all covalent bonds (type 0).

do iAtom=1,nAtoms

  iRow = iTabRow(iANr(iAtom))
  if (iRow == 0) Go To 99

# ifdef _DEBUGPRINT_
  write(6,*)
  write(6,*) 'iAtom, iAnr=',iAtom,iANr(iAtom)
  write(6,*)
# endif
  ix = iBox(1,iAtom)
  iy = iBox(2,iAtom)
  iz = iBox(3,iAtom)

  ! Loop over all adjacent boxes.

  do jx=ix-1,ix+1
    do jy=iy-1,iy+1
      do jz=iz-1,iz+1
        call Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,jx,jy,jz,iAtom,iRow,iANr,iTabBonds,nBonds,nBondMax,iTabAtoms,nMax,ThrB,1.0D+99)

      end do
    end do
  end do
99 continue
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'After Covalent Bonds'
write(6,*)
write(6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(6,*)
  write(6,*) 'iAtom=',iAtom
  write(6,*)
  write(6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
  write(6,*) ' nCoBond :',nCoBond(iAtom,nAtoms,nMax,iTabBonds,nBondMax,iTabAtoms)

end do
write(6,*)
write(6,*)
write(6,*) 'Bonds:'
do iBond=1,nBonds
  write(6,*)
  write(6,*) 'iBond=',iBond
  write(6,*)
  write(6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Find all vdW bonds (type 1).

! Loop over all adjacent boxes. Do vdW bonds.

ThrB_vdW = 1.0D-4*ThrB
do iAtom=1,nAtoms

  iRow = iTabRow(iANr(iAtom))
  if (iRow == 0) Go To 98

  ix = iBox(1,iAtom)
  iy = iBox(2,iAtom)
  iz = iBox(3,iAtom)

  ! Loop over all adjacent boxes.

  do jx=ix-1,ix+1
    do jy=iy-1,iy+1
      do jz=iz-1,iz+1
        call Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,jx,jy,jz,iAtom,iRow,iANr,iTabBonds,nBonds,nBondMax,iTabAtoms,nMax,ThrB,ThrB_vdW)

      end do
    end do
  end do
98 continue
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'After vdW Bonds'
write(6,*)
write(6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(6,*)
  write(6,*) 'iAtom=',iAtom
  write(6,*)
  write(6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(6,*)
write(6,*)
write(6,*)
write(6,*) 'Bonds:'
do iBond=1,nBonds
  write(6,*)
  write(6,*) 'iBond=',iBond
  write(6,*)
  write(6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Test connectivity in case of fragmentation of the system.

! These bonds are of type 2.

call Connect_Fragments(nAtoms,iTabBonds,nBondMax,nBonds,Coor,iTabAtoms,nMax,iANr)
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'After Connecting Fragments'
write(6,*)
write(6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(6,*)
  write(6,*) 'iAtom=',iAtom
  write(6,*)
  write(6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(6,*)
write(6,*)
write(6,*) 'Bonds:'
do iBond=1,nBonds
  write(6,*)
  write(6,*) 'iBond=',iBond
  write(6,*)
  write(6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
! Finally, for partially linear system we need "extra long bond",
! this to make sure that we include torsions over linear subsystems.
! This "magic bonds" are only included when forming torsions.

call Magic_Bonds(Coor,nAtoms,iTabBonds,nBondMax,nBonds,iTabAtoms,nMax)

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(6,*)
write(6,*) 'After Magic Bonds'
write(6,*)
write(6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(6,*)
  write(6,*) 'iAtom=',iAtom
  write(6,*)
  write(6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(6,*)
write(6,*)
write(6,*) 'Bonds:'
do iBond=1,nBonds
  write(6,*)
  write(6,*) 'iBond=',iBond
  write(6,*)
  write(6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif

return

end subroutine Find_Bonds
