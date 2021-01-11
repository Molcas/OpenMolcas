************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2013, Ignacio Fdez. Galvan                             *
************************************************************************
*  Align
*
*> @brief
*>   Align two structures.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Align a molecular structure with the reference, using the stored weights
*> (masses by default). This is sometimes needed to ensure that an optimal
*> structure is found when there are constraints expressed in weighted space
*> (e.g. `sphere` or `transverse`).
*>
*> @param[in,out] Coord Cartesian coordinates to align
*> @param[in]     Ref   Cartesian coordinates of the reference structure
*> @param[in]     nAtom Number of symmetry-unique atoms
************************************************************************
      Subroutine Align(Coord,Ref,nAtom)
      use Symmetry_Info, only: nIrrep, iOper
      use Slapaf_Info, only: Weights
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "sbs.fh"
#include "stdalloc.fh"
#include "weighting.fh"
      Real*8 Coord(3*nAtom), Ref(3*nAtom)
      Logical Invar
      Real*8, Allocatable:: Coor_All(:,:), Ref_All(:,:)
      Integer, Allocatable:: iStab(:)

*---- Do nothing if the energy is not rot. and trans. invariant
      Invar=(iAnd(iSBS,2**7).eq.0).and.(iAnd(iSBS,2**8).eq.0)
      If (.Not.Invar) Return

      Call mma_allocate(Coor_All,3,nAtom*8,Label='Coor_All')
      Call Expand_Coor(Coord,nAtom,Coor_All,mAtom)
      Call mma_allocate(Ref_All ,3,nAtom*8,Label='Ref_All')
      Call Expand_Coor(Ref,  nAtom,Ref_All,mAtom)

c     Call RecPrt('Coord before align',' ',Coor_All,3,mAtom)

      Call Superpose_w(Coor_All,Ref_All,Weights,mAtom,RMS,RMSMax)

*---- Get the stabilizers for each atom (to keep the symmetry)
*     (code copied from init_slapaf)
      Call mma_allocate(iStab,nAtom,Label='iStab')
      Do iAt=1,nAtom
        iAdr=(iAt-1)*3+1
        iChxyz=0
        Do i=0,2
          If (Ref(iAdr+i).ne.Zero) Then
            Do iIrrep=0,nIrrep-1
              If (iAnd(2**i,iOper(iIrrep)).ne.0)
     &           iChxyz=iOr(iChxyz,2**i)
            End Do
          End If
        End Do
        nStb=0
        Do iIrrep=0,nIrrep-1
          If ((nStb.le.1).And.
     &       (iAnd(iChxyz,iOper(iIrrep)).eq.0)) Then
            iStab(iAt)=iOper(iIrrep)
            nStb=nStb+1
          End If
        End Do
      End Do
*
      Call Fix_Symmetry(Coor_All,nAtom,iStab)
      Call mma_deallocate(iStab)

      call dcopy_(3*nAtom,Coor_All,1,Coord,1)

c     Call RecPrt('Coord after align',' ',Coor_All,3,mAtom)

      Call mma_deallocate(Coor_All)
      Call mma_deallocate(Ref_All)

      End
