************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Find_Bonds(Coor,nAtoms,iTab,nMax,nx,ny,nz,iBox,iANr,
     &                      iTabBonds,nBonds,
     &                      nBondMax,iTabAtoms,ThrB)
      Implicit Real*8 (a-h,o-z)
      Real*8 Coor(3,nAtoms)
      Integer iTab(0:nMax,nx,ny,nz), iBox(3,nAtoms), iANr(nAtoms),
     &        iTabBonds(3,nBondMax), iTabAtoms(2,0:nMax,nAtoms)
#include "bondtypes.fh"
*                                                                      *
************************************************************************
*                                                                      *
!#define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
*     Initiate
*
#ifdef _DEBUGPRINT_
      Call RecPrt('Coor',' ',Coor,3,nAtoms)
      Write (6,*) 'Find_Bonds: ThrB=',ThrB
      Write (6,*) 'Initialize iTabAtoms'
#endif
*
      Call ICopy(nBondMax*3,[0],0,iTabBonds,1)
      Call ICopy(2*(nMax+1)*nAtoms,[0],0,iTabAtoms,1)
      nBonds = 0
*                                                                      *
************************************************************************
*                                                                      *
*     Find all covalent bonds (type 0).
*
      Do iAtom = 1, nAtoms
*
         iRow = iTabRow(iANr(iAtom))
         If (iRow.eq.0) Go To 99
*
#ifdef _DEBUGPRINT_
         Write (6,*)
         Write (6,*) 'iAtom, iAnr=',iAtom,iANr(iAtom)
         Write (6,*)
#endif
         ix=iBox(1,iAtom)
         iy=iBox(2,iAtom)
         iz=iBox(3,iAtom)
*
*------- Loop over all adjacent boxes.
*
         Do jx = ix-1, ix+1
            Do jy = iy-1, iy+1
               Do jz = iz-1, iz+1
                  Call Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,jx,jy,jz,
     &                             iAtom,iRow,iANr,
     &                             iTabBonds,nBonds,nBondMax,iTabAtoms,
     &                             nMax,ThrB,1.0D+99)

               End Do
            End Do
         End Do
 99      Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'After Covalent Bonds'
      Write (6,*)
      Write (6,*) 'iTabAtoms:'
      Do iAtom = 1, nAtoms
         Write (6,*)
         Write (6,*) 'iAtom=',iAtom
         Write (6,*)
         Write (6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
         nn=iTabAtoms(1,0,iAtom)
         Write (6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
         Write (6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
         Write (6,*) ' Bondtype:',(
     &               iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
         Write (6,*) ' nCoBond :', nCoBond(iAtom,nAtoms,nMax,iTabBonds,
     &                                     nBondMax,iTabAtoms)

      End Do
      Write (6,*)
      Write (6,*)
      Write (6,*) 'Bonds:'
      Do iBond = 1, nBonds
         Write (6,*)
         Write (6,*) 'iBond=',iBond
         Write (6,*)
         Write (6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
         Write (6,*) 'Bondtype:',BondType(Min(3,iTabBonds(3,iBond)))
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Find all vdW bonds (type 1).
*
*     Loop over all adjacent boxes. Do vdW bonds.
*
      ThrB_vdW=1.0D-4*ThrB
      Do iAtom = 1, nAtoms
*
         iRow = iTabRow(iANr(iAtom))
         If (iRow.eq.0) Go To 98
*
         ix=iBox(1,iAtom)
         iy=iBox(2,iAtom)
         iz=iBox(3,iAtom)
*
*------- Loop over all adjacent boxes.
*
         Do jx = ix-1, ix+1
            Do jy = iy-1, iy+1
               Do jz = iz-1, iz+1
                  Call Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,jx,jy,jz,
     &                             iAtom,iRow,iANr,
     &                             iTabBonds,nBonds,nBondMax,iTabAtoms,
     &                             nMax,ThrB,ThrB_vdW)

               End Do
            End Do
         End Do
 98      Continue
      End Do
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'After vdW Bonds'
      Write (6,*)
      Write (6,*) 'iTabAtoms:'
      Do iAtom = 1, nAtoms
         Write (6,*)
         Write (6,*) 'iAtom=',iAtom
         Write (6,*)
         Write (6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
         nn=iTabAtoms(1,0,iAtom)
         Write (6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
         Write (6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
         Write (6,*) ' Bondtype:',(
     &               iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
      End Do
      Write (6,*)
      Write (6,*)
      Write (6,*)
      Write (6,*) 'Bonds:'
      Do iBond = 1, nBonds
         Write (6,*)
         Write (6,*) 'iBond=',iBond
         Write (6,*)
         Write (6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
         Write (6,*) 'Bondtype:',BondType(Min(3,iTabBonds(3,iBond)))
      End Do
#endif
*                                                                      *
************************************************************************
*                                                                      *
*
*     Test connectivity in case of fragmentation of the system.
*
*     These bonds are of type 2.
*
      Call Connect_Fragments(nAtoms,iTabBonds,nBondMax,
     &                       nBonds,Coor,iTabAtoms,nMax,iANr)
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'After Connecting Fragments'
      Write (6,*)
      Write (6,*) 'iTabAtoms:'
      Do iAtom = 1, nAtoms
         Write (6,*)
         Write (6,*) 'iAtom=',iAtom
         Write (6,*)
         Write (6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
         nn=iTabAtoms(1,0,iAtom)
         Write (6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
         Write (6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
         Write (6,*) ' Bondtype:',(
     &               iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
      End Do
      Write (6,*)
      Write (6,*)
      Write (6,*) 'Bonds:'
      Do iBond = 1, nBonds
         Write (6,*)
         Write (6,*) 'iBond=',iBond
         Write (6,*)
         Write (6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
         Write (6,*) 'Bondtype:',BondType(Min(3,iTabBonds(3,iBond)))
      End Do
#endif
*
*                                                                      *
************************************************************************
*                                                                      *
*     Finally, for partially linear system we need "extra long bond",
*     this to make sure that we include torsions over linear subsystems.
*     This "magic bonds" are only included when forming torsions.
*
      Call Magic_Bonds(Coor,nAtoms,iTabBonds,nBondMax,nBonds,
     &                 iTabAtoms,nMax)
*
*                                                                      *
************************************************************************
*                                                                      *
#ifdef _DEBUGPRINT_
      Write (6,*)
      Write (6,*) 'After Magic Bonds'
      Write (6,*)
      Write (6,*) 'iTabAtoms:'
      Do iAtom = 1, nAtoms
         Write (6,*)
         Write (6,*) 'iAtom=',iAtom
         Write (6,*)
         Write (6,*) 'nNeighbor=',iTabAtoms(1,0,iAtom)
         nn=iTabAtoms(1,0,iAtom)
         Write (6,*) ' Neigbors:',(iTabAtoms(1,i,iAtom),i=1,nn)
         Write (6,*) ' Bond    :',(iTabAtoms(2,i,iAtom),i=1,nn)
         Write (6,*) ' Bondtype:',(
     &               iTabBonds(3,iTabAtoms(2,i,iAtom)),i=1,nn)
      End Do
      Write (6,*)
      Write (6,*)
      Write (6,*) 'Bonds:'
      Do iBond = 1, nBonds
         Write (6,*)
         Write (6,*) 'iBond=',iBond
         Write (6,*)
         Write (6,*) 'Atoms=',iTabBonds(1,iBond),iTabBonds(2,iBond)
         Write (6,*) 'Bondtype:',BondType(Min(3,iTabBonds(3,iBond)))
      End Do
#endif
*
      Return
      End
