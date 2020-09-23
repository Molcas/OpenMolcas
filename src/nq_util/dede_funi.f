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
      Subroutine DeDe_Funi(Dens,nDens,nr_of_Densities)
      use k2_arrays
      use Sizes_of_Seward, only: S
      use Symmetry_Info, only: nIrrep
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "setup.fh"
      Real*8 Dens(nDens,nr_of_Densities)
      Logical Special_NoSym, DFT_Storage
*
      nIndij=S%nShlls*(S%nShlls+1)/2
      nField=2+nr_of_Densities
*
*
      Call mma_allocate(ipOffD,nField,nIndij,label='ipOffD')
      Call mma_allocate(DeDe,nDeDe_DFT+MaxDe*nIrrep,Label='DeDe')
      ipDeDe= 1
      ipD00 = ipDeDe + nDeDe_DFT
      ipDijs = -1  ! Dummy value
      DeDe(:)=Zero
*
      Special_NoSym=.False.
      DFT_Storage=.True.
      Call mk_DeDe(Dens,nDens,nr_of_Densities,ipOffD,nIndij,ipDeDe,
     &             ipD00,MaxDe,mDeDe,mIndij,Special_NoSym,DFT_Storage,
     &             DeDe,nDeDe_DFT)
*                                                                      *
************************************************************************
*                                                                      *
      Return
      End
