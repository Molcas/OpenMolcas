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

#ifdef _DEBUGPRINT_
use Slapaf_Info, only: BondType
#endif
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nMax, nx, ny, nz, iTab(0:nMax,nx,ny,nz), iBox(3,nAtoms), iANr(nAtoms), nBondMax
real(kind=wp), intent(in) :: Coor(3,nAtoms), ThrB
integer(kind=iwp), intent(out) :: iTabBonds(3,nBondMax), nBonds, iTabAtoms(2,0:nMax,nAtoms)
integer(kind=iwp) :: iAtom, iRow, ix, iy, iz, jx, jy, jz
real(kind=wp) :: ThrB_vdW
integer(kind=iwp), external :: iTabRow

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
write(u6,*) 'Find_Bonds: ThrB=',ThrB
write(u6,*) 'Initialize iTabAtoms'
#endif

iTabBonds(:,:) = 0
iTabAtoms(:,:,:) = 0
nBonds = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Find all covalent bonds (type 0).

do iAtom=1,nAtoms

  iRow = iTabRow(iANr(iAtom))
  if (iRow == 0) cycle

# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'iAtom, iAnr=',iAtom,iANr(iAtom)
  write(u6,*)
# endif
  ix = iBox(1,iAtom)
  iy = iBox(2,iAtom)
  iz = iBox(3,iAtom)

  ! Loop over all adjacent boxes.

  do jx=ix-1,ix+1
    do jy=iy-1,iy+1
      do jz=iz-1,iz+1
        call Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,jx,jy,jz,iAtom,iRow,iANr,iTabBonds,nBonds,nBondMax,iTabAtoms,nMax,ThrB,1.0e99_wp)

      end do
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'After Covalent Bonds'
write(u6,*)
write(u6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(u6,*)
  write(u6,*) 'iAtom=',iAtom
  write(u6,*)
  write(u6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(u6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(u6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(u6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
  write(u6,*) ' nCoBond :',nCoBond(iAtom,nAtoms,nMax,iTabBonds,nBondMax,iTabAtoms)

end do
write(u6,*)
write(u6,*)
write(u6,*) 'Bonds:'
do iBond=1,nBonds
  write(u6,*)
  write(u6,*) 'iBond=',iBond
  write(u6,*)
  write(u6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(u6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Find all vdW bonds (type 1).

! Loop over all adjacent boxes. Do vdW bonds.

ThrB_vdW = 1.0e-4_wp*ThrB
do iAtom=1,nAtoms

  iRow = iTabRow(iANr(iAtom))
  if (iRow == 0) cycle

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
end do
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'After vdW Bonds'
write(u6,*)
write(u6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(u6,*)
  write(u6,*) 'iAtom=',iAtom
  write(u6,*)
  write(u6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(u6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(u6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(u6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(u6,*)
write(u6,*)
write(u6,*)
write(u6,*) 'Bonds:'
do iBond=1,nBonds
  write(u6,*)
  write(u6,*) 'iBond=',iBond
  write(u6,*)
  write(u6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(u6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
! Test connectivity in case of fragmentation of the system.

! These bonds are of type 2.

call Connect_Fragments(nAtoms,iTabBonds,nBondMax,nBonds,Coor,iTabAtoms,nMax,iANr)
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'After Connecting Fragments'
write(u6,*)
write(u6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(u6,*)
  write(u6,*) 'iAtom=',iAtom
  write(u6,*)
  write(u6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(u6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(u6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(u6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(u6,*)
write(u6,*)
write(u6,*) 'Bonds:'
do iBond=1,nBonds
  write(u6,*)
  write(u6,*) 'iBond=',iBond
  write(u6,*)
  write(u6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(u6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
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
write(u6,*)
write(u6,*) 'After Magic Bonds'
write(u6,*)
write(u6,*) 'iTabAtoms:'
do iAtom=1,nAtoms
  write(u6,*)
  write(u6,*) 'iAtom=',iAtom
  write(u6,*)
  write(u6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
  nn = iTabAtoms(1,0,iAtom)
  write(u6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
  write(u6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
  write(u6,*) ' Bondtype:',(iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
end do
write(u6,*)
write(u6,*)
write(u6,*) 'Bonds:'
do iBond=1,nBonds
  write(u6,*)
  write(u6,*) 'iBond=',iBond
  write(u6,*)
  write(u6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
  write(u6,*) 'Bondtype:',BondType(min(3,iTabBonds(3,iBond)))
end do
#endif

return

end subroutine Find_Bonds
