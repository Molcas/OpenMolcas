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
      Subroutine Sort_to_Box(Coor,nAtoms,iTab,nMax,nx,ny,nz,iBox,iANr,
     &                       xmin,ymin,zmin,Box_Size)
      Implicit Real*8 (a-h,o-z)
      Real*8 Coor(3,nAtoms)
      Integer iTab(0:nMax,nx,ny,nz), iBox(3,nAtoms), iANr(nAtoms)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUG_
*                                                                      *
************************************************************************
*                                                                      *
*
      Call ICopy((nMax+1)*nx*ny*nz,[0],0,iTab,1)
*
      Do iAtom = 1, nAtoms
*
         If (iTabRow(iAnr(iAtom)).eq.0) Go To 99
*
         x = Coor(1,iAtom)-xmin
         ix = Int(x/Box_Size)+1
         y = Coor(2,iAtom)-ymin
         iy = Int(y/Box_Size)+1
         z = Coor(3,iAtom)-zmin
         iz = Int(z/Box_Size)+1
         iBox(1,iAtom)=ix
         iBox(2,iAtom)=iy
         iBox(3,iAtom)=iz
*
         Nr = iTab(0,ix,iy,iz) + 1
         If (Nr.gt.nMax) Then
            Call WarningMessage(2,'Sort_to_Box: Nr.gt.nMax')
            Call Abend()
         End If
#ifdef _DEBUG_
         Write (6,*) 'Sort_to_Box: ix,iy,iz,Nr,iAtom=',
     &                             ix,iy,iz,Nr,iAtom
#endif
*
         iTab(0, ix,iy,iz)=Nr
         iTab(Nr,ix,iy,iz)=iAtom
*
 99      Continue
      End Do
*
      Return
      End
