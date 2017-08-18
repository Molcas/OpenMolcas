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

      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "sbs.fh"
#include "WrkSpc.fh"
#include "info_slapaf.fh"
#include "weighting.fh"
      Real*8 Coord(3*nAtom), Ref(3*nAtom)
      Logical Invar

*---- Do nothing if the energy is not rot. and trans. invariant
      Invar=(iAnd(iSBS,2**7).eq.0).and.(iAnd(iSBS,2**8).eq.0)
      If (.Not.Invar) Return

      Call Allocate_Work(ipx,3*nAtom*8)
      Call Expand_Coor(Coord,nAtom,Work(ipx),mAtom,nSym,iOper)
      Call Allocate_Work(ipy,3*nAtom*8)
      Call Expand_Coor(Ref,  nAtom,Work(ipy),mAtom,nSym,iOper)

c     Call RecPrt('Coord before align',' ',Work(ipx),3,mAtom)

      Call Superpose_w(Work(ipx),Work(ipy),Work(ipWeights),mAtom,
     &                 RMS,RMSMax)

*---- Get the stabilizers for each atom (to keep the symmetry)
*     (code copied from init_slapaf)
      Call Allocate_iWork(ipStab,nAtom)
      Do iAt=1,nAtom
        iAdr=(iAt-1)*3+1
        iChxyz=0
        Do i=0,2
          If (Ref(iAdr+i).ne.Zero) Then
            Do iIrrep=0,nSym-1
              If (iAnd(2**i,iOper(iIrrep)).ne.0)
     &           iChxyz=iOr(iChxyz,2**i)
            End Do
          End If
        End Do
        nStb=0
        Do iIrrep=0,nSym-1
          If ((nStb.le.1).And.
     &       (iAnd(iChxyz,iOper(iIrrep)).eq.0)) Then
            iWork(ipStab+iAt-1)=iOper(iIrrep)
            nStb=nStb+1
          End If
        End Do
      End Do
      Call Fix_Symmetry(Work(ipx),nAtom,iWork(ipStab))
      Call Free_iWork(ipStab)

      call dcopy_(3*nAtom,Work(ipx),1,Coord,1)

c     Call RecPrt('Coord after align',' ',Work(ipx),3,mAtom)

      Call Free_Work(ipx)
      Call Free_Work(ipy)

      End
