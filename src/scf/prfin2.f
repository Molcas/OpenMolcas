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
* Copyright (C) 2017, Roland Lindh                                     *
************************************************************************
      SubRoutine PrFin2(Ovlp,nDT,OccNo,nEO,CMO,nCMO,note)
      Implicit Real*8 (a-h,o-z)
*
      Real*8 Ovlp(nDT),OccNo(nEO),CMO(nBB)
* PAM 2007: Changed dimension of CMO array from nCMO to NBB:
* The larger size is needed here, and the allocated size is nBB.
      Character*80 Note
*
#include "real.fh"
#include "mxdm.fh"
#include "infscf.fh"
#include "file.fh"
#include "stdalloc.fh"
#include "rctfld.fh"
#include "oneswi.fh"
*
*---- Define local variables
      Real*8, Dimension(:), Allocatable:: Scr2

#include "SysDef.fh"
*
*---- Write orbitals on the file (the case InVec=3 and nIter=0
*     is set up in RdInp)
      If (jVOut.ge.2) Then
*
*------- Allocate memory for expanded CMO's and occupation numbers
         Call mma_allocate(Scr2,nBB,Label='Scr2')
*
*------- Prepare CMO in symmetry blocks nBas x nBas
         iVec  = 0
         iCMO  = 0
*         iseed = 153759
         Do iSym = 1, nSym
            Do i = 1, nBas(iSym)*nOrb(iSym)
               Scr2(iVec + i) = CMO(iCMO + i)
            End Do
            iVec = iVec + nBas(iSym)*nOrb(iSym)
            Do i = 1, nBas(iSym)*(nBas(iSym) - nOrb(iSym))
*              Scr2(iVec + i) = Random(iseed) - Half
               Scr2(iVec + i) = Zero
            End Do
            iVec = iVec + nBas(iSym)*(nBas(iSym) - nOrb(iSym))
            iCMO = iCMO + nOrb(iSym)*nBas(iSym)
         End Do
         Do i = 1, nBB
            CMO(i)=Scr2(i)
         End Do
*
*------- Prepare occupation numbers
         iOr = 0
         iBs = 0
         Do iSym = 1, nSym
            Do j = 1, nOrb(iSym)
               Scr2(iBs + j) = OccNo(iOr + j)
            End Do
            Do j = nOrb(iSym) + 1, nBas(iSym)
               Scr2(iBs + j) = Zero
            End Do
            iOr = iOr + nOrb(iSym)
            iBs = iBs + nBas(iSym)
         End Do
         Do i = 1, nnB
            OccNo(i) = Scr2(i)
         End Do
*
         Do iSym = 1, nSym
            nOrb(iSym) = nBas(iSym)
         End Do
         Call SetUp()
*
*------- Orthogonalize vectors
         Call Ortho(CMO,nCMO,Ovlp,nBT)
*
*------- Write on the file
         If (KSDFT.eq.'SCF') Then
            If(iUHF.eq.0) Then
               Note='* SCF orbitals'
               If (kIVO.ne.0 ) Note='* SCF orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='* SCF orbitals + arbitrary occupations'
            Else
               Note='* UHF orbitals'
               If (kIVO.ne.0 ) Note='* UHF orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='* UHF orbitals + arbitrary occupations'
            End If
         Else
            If(iUHF.eq.0) Then
               Note='* RKS-DFT orbitals'
               If (kIVO.ne.0 ) Note='* RKS-DFT orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='* RKS-DFT orbitals + arbitrary occupations'
            Else
               Note='* UKS-DFT orbitals'
               If (kIVO.ne.0 ) Note='* UKS-DFT orbitals + IVO'
               If (iCoCo.ne.0)
     &            Note='* UKS-DFT orbitals + arbitrary occupations'
            End If
         End If
*
*------- Deallocate memory for expanded CMO's and occupation numbers
         Call mma_deallocate(Scr2)
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
