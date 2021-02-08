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
      Subroutine Force(nFix,GrdX,nAtom,nInter,BMx,Iter,Grad,Lbl,Degen)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "stdalloc.fh"
#include "Molcas.fh"
      Real*8 GrdX(3*nAtom), BMx(3*nAtom,3*nAtom),
     &       Grad(nInter,Iter), Degen(3*nAtom)
      Character Lbl(nInter)*8
      Real*8 Dummy(1)
      Real*8, Allocatable:: Frc(:)
*
#ifdef _DEBUGPRINT_
      Call RecPrt('In Force:BMx ',' ',BMx ,3*nAtom,nInter)
      Call RecPrt('In Force:Degen ',' ',Degen ,1,3*nAtom)
      Call RecPrt('In Force:GrdX',' ',GrdX,3,nAtom)
#endif
*
*-----Frozen internal coordinates are introduced as constraints to the
*     energy functional with the Lagrange multipliers method. The new
*     energy functional will have zero gradient with respect to the
*     frozen parameters.
*
      Call mma_allocate(Frc,3*nAtom,Label='Frc')
*
*     Compute the norm of the cartesian force vector.
*
*     |dE/dx|=Sqrt(dE/dx|u|dE/dx)
*
      Do i = 1, 3*nAtom
         Frc(i) = Degen(i)*GrdX(i)
      End Do
*                                                                      *
************************************************************************
*                                                                      *
*---- Solve the equation dq/dx dE/dq = dE/dx
*
*     B dE/dq = dE/dx
*
      M = 3*nAtom
      N = nInter
      NRHS=1
      Call Eq_Solver('N',M,N,NRHS,BMx,.TRUE.,Dummy,Frc,
     &               Grad(1,Iter))
#ifdef _DEBUGPRINT_
      Call RecPrt(' Internal Forces in au before FIXIC ',
     &            ' ',Grad(1,Iter),nInter,1)
#endif
*                                                                      *
************************************************************************
*                                                                      *
*     Remove gradient components in the constraints directions.
*
      If (nFix.ne.0) Call Fixic(nFix,Grad(1,Iter),nInter,BMx,nAtom*3,
     &                          GrdX,Lbl,Degen)
*
*-----Write cartesian symmetry distinct forces which will be relaxed.
*
#ifdef _DEBUGPRINT_
      BLOCK
      use Slapaf_Info: only: AtomLbl
      Call PrList('Cartesian forces which will be relaxed'
     &            //' hartree/bohr',
     &            AtomLbl,nAtom,GrdX,3,nAtom)
      END BLOCK
#endif
*
      Call mma_deallocate(Frc)
*
      Return
      End
