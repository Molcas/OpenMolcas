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
* Copyright (C) 2019, Ignacio Fdez. Galvan                             *
************************************************************************
*  Transpose_TDM
*
*> @brief
*>   Transpose a transition density matrix in place.
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Transpose a transition density matrix, stored in symmetry blocks,
*> replacing the original matrix. The matrices contain only the symmetry
*> blocks that match the total symmetry of the transition.
*>
*> @param[in,out] TDM       Transition density matrix
*> @param[in]     Symmetry  Symmetry of the transition
************************************************************************
      Subroutine Transpose_TDM(TDM,Symmetry)
      Implicit None
      Real*8, Intent(InOut) :: TDM(*)
      Integer, Intent(In) :: Symmetry
#include "rassi.fh"
#include "symmul.fh"
#include "stdalloc.fh"
      Integer :: iSym1,iSym2,nTot,i,j
      Integer :: iBlock(0:8)
      Real*8, Allocatable :: Tmp(:)
* Compute the location of all the stored symmetry blocks
      nTot=0
      iBlock(0)=0
      Do iSym1=1,nSym
        iSym2=Mul(Symmetry,iSym1)
        nTot=nTot+nBasF(iSym1)*nBasF(iSym2)
        iBlock(iSym1)=nTot
      End Do
* Make a copy so we can transpose in place
      Call mma_Allocate(Tmp,nTot,Label='Tmp')
      Call dCopy_(nTot,TDM,1,Tmp,1)
* Transpose symmetry block (a,b) onto symmetry block (b,a)
      Do iSym1=1,nSym
        iSym2=Mul(Symmetry,iSym1)
        Do i=1,nBasF(iSym2)
          Do j=1,nBasF(iSym1)
            TDM(iBlock(iSym2-1)+(j-1)*nBasF(iSym2)+i) =
     &      Tmp(iBlock(iSym1-1)+(i-1)*nBasF(iSym1)+j)
          End Do
        End Do
      End Do
      Call mma_deAllocate(Tmp)
      End Subroutine
