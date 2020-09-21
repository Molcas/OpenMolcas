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
      Subroutine DeDe_SCF(Dens,TwoHam,nDens,mDens)
      use k2_arrays
      use Basis_Info, only: nBas
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "setup.fh"
      Integer nDens, mDens
      Real*8, Target:: Dens(nDens), TwoHam(nDens)
      Logical Special_NoSym, DFT_Storage
*
#ifdef _DEBUG_
      Call qEnter('DeDe_SCF')
#endif
      nr_of_Densities=1  ! Hardwired option
*
      nIndij=nShlls*(nShlls+1)/2
      nField=2+nr_of_Densities
      Call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')
*
*     The array with desymmetrized densities contain two additional
*     fields.
*     ipD00 is a null matrix, which should simplify the logic.
*     ipDijS is an auxilliary memory if not the whole set of a
*      desymmetrized density could be used.
*
      nDeDe_tot = nDeDe + MaxDe*nIrrep + MxDij
      Call mma_allocate(DeDe,nDeDe_tot,Label='DeDe')
      ipDeDe = 1
      ipD00 = ipDeDe + nDeDe
      ipDijS= ipD00  + MaxDe*nIrrep
      DeDe(:)=Zero
*
      Special_NoSym=.True.
      DFT_Storage=.False.
      Call mk_DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,
     &             ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &             DeDe,nDeDe)
*                                                                      *
************************************************************************
*                                                                      *
*     In case of no symmetry do a temporary square copy of the
*     density matrix.
*
*---- Change the folded density to pure triangular form,
*     i.e. off-diagonal elements are divided by two.
*
      If (nIrrep.eq.1) Then
         Call DScal_(nDens,Half,Dens,1)
         ij=0
         Do i = 1, nBas(0)
            ij = ij + i
            Dens(ij)=Two*Dens(ij)
         End Do
         mDens=nbas(0)*nbas(0)
         Call mma_allocate(Dq,mDens,Label='Dq')
         Call Square(Dens,Dq,1,nbas(0),nbas(0))
         pDq => Dq(:)
*
         Call mma_allocate(Fq,mDens,Label='Fq')
         Fq(:)=Zero
         pFq => Fq(:)
      Else
         mDens=nDens
         pDq => Dens(:)
         pFq => Twoham(:)
      End If
#ifdef _DEBUG_
      Call qExit('DeDe_SCF')
#endif
*
      Return
      End
