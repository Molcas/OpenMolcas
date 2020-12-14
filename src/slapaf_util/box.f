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
      Subroutine Box(Coor,nAtoms,iANr,iOptC,TabB,TabA,nBonds,nMax)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
      Real*8 Coor(3,nAtoms)
      Integer iANr(nAtoms)
      Integer, Allocatable:: TabB(:,:), TabA(:,:,:), iBox(:,:),
     &                       Tab(:,:,:,:)
*                                                                      *
************************************************************************
*                                                                      *
*define _DEBUGPRINT_
*                                                                      *
************************************************************************
*                                                                      *
      If (nAtoms.lt.2) Then
         Write (6,*) 'Too few atoms to relax: nAtoms=',nAtoms
         Call WarningMessage(2,'nAtoms.lt.2')
         Call Abend()
      End If
*
      ThrB=0.40D0
#ifdef _DEBUGPRINT_
      Call RecPrt('Box: Coor',' ',Coor,3,nAtoms)
      Write (6,*) 'Box: ThrB=',ThrB
#endif
*
      xmin= 1.0D+10
      ymin= 1.0D+10
      zmin= 1.0D+10
      xmax=-1.0D+10
      ymax=-1.0D+10
      zmax=-1.0D+10
*
*---- Establish boundaries
*
      Do iAtom = 1, nAtoms
         xmin=Min(xmin,Coor(1,iAtom))
         xmax=Max(xmax,Coor(1,iAtom))
         ymin=Min(ymin,Coor(2,iAtom))
         ymax=Max(ymax,Coor(2,iAtom))
         zmin=Min(zmin,Coor(3,iAtom))
         zmax=Max(zmax,Coor(3,iAtom))
      End Do
      xmin=xmin-1.0D-2
      xmax=xmax+1.0D-2
      ymin=ymin-1.0D-2
      ymax=ymax+1.0D-2
      zmin=zmin-1.0D-2
      zmax=zmax+1.0D-2
*
      Box_Size = 8.0D0  !   a.u.
      nx = Max(1,INT((xmax-xmin)/Box_Size)+1)
      adjust = (DBLE(nx)*Box_size - (xmax-xmin))/Two
      xmin=xmin-adjust
      xmax=xmax+adjust
      ny = MAX(1,INT((ymax-ymin)/Box_Size)+1)
      adjust = (DBLE(ny)*Box_size - (ymax-ymin))/Two
      ymin=ymin-adjust
      ymax=ymax+adjust
      nz = MAX(1,INT((zmax-zmin)/Box_Size)+1)
      adjust = (DBLE(nz)*Box_size - (zmax-zmin))/Two
      zmin=zmin-adjust
      zmax=zmax+adjust
#ifdef _DEBUGPRINT_
      Write (6,*) 'nx,ny,nz=',nx,ny,nz
#endif
*
      nMax=100
cnf      nMax=40
c AOM Fixed this size to account for double VdW counting
      nBondMax=nAtoms*(nAtoms+1)
c AOM
      Call mma_allocate(TabB,3,nBondMax,Label='TabB')
      Call mma_allocate(TabA,[1,2],[0,nMax],[1,nAtoms],Label='TabA')
      Call mma_allocate(Tab,[0,nMax],[1,nx],[1,ny],[1,nz],Label='Tab')
      Call mma_allocate(iBox,3,nAtoms,Label='iBox')
*
      Call Sort_to_Box(Coor,nAtoms,Tab,nMax,nx,ny,nz,
     &                 iBox,iANr,xmin,ymin,zmin,Box_Size)
*
      Call Find_Bonds(Coor,nAtoms,Tab,nMax,nx,ny,nz,iBox,iANr,
     &                iOptC,TabB,nBonds,nBondMax,TabA,ThrB)
*
      Call mma_deallocate(iBox)
      Call mma_deallocate(Tab)
*
      Return
      End
