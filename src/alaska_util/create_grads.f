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
* Copyright (C) 2015, Ignacio Fdez. Galvan                             *
************************************************************************
*  Create_Grads
*
*> @brief Create an empty gradients file
*> @author Ignacio Fdez. Galv&aacute;n
*>
*> @details
*> Creates an empty gradients file (GRADS), to contain all the computed
*> gradients and couplings for the different roots.
*>
*> Structure of the gradients file:
*>  - 0: \p TOC (5 integers, each a disk address)
*>  - \p TOC(1): \p nRoots
*>  - \p TOC(2): \p nGrad
*>  - \p TOC(3): \p i_grad (\p nRoots integers)
*>  - \p TOC(4): \p i_nac (max(1,\p nRoots * (\p nRoots - 1)/2) integers)
*>  - \p TOC(5): Next free address
*>  - \p i_grad(i): Gradient of state \p i (\p nGrad reals)
*>  - \p i_nac(m): Coupling between states \p i and \p j
*>                 (\p m = \p i * (\p i - 1)/2 + \p j, \p i > \p j) (\p nGrad reals)
*>
*> If \p i_grad(i) or \p i_nac(m) is 0, the vector is not present in the
*> file. If it is negative, the vector is not present and cannot be
*> computed (so it is pointless to ask for it again).
*>
*> @param[in] FN     Filename for the gradients file (prgm-translatable)
*> @param[in] nRoots Number of roots to include in the file
*> @param[in] nGrad  Length of each gradient vector
************************************************************************
      Subroutine Create_Grads(FN,nRoots,nGrad)
      Implicit None
#include "stdalloc.fh"
      Integer :: nRoots,nGrad,Lu,TOC(5),nCoup,iAd
      Character(LEN=*) :: FN
      Integer, Dimension(:), Allocatable:: i_grad, i_nac
*
      nCoup=Max(1,nRoots*(nRoots-1)/2)
      Call mma_Allocate(i_grad,nRoots)
      Call mma_Allocate(i_nac,nCoup)
      Call iCopy(Size(TOC),[0],0,TOC,1)
      Call iCopy(nRoots,[0],0,i_grad,1)
      Call iCopy(nCoup,[0],0,i_nac,1)
*
      Lu=20
      Call DaName(Lu,Trim(FN))
      iAd=0
      Call iDaFile(Lu,1,TOC,SIZE(TOC),iAd)
      TOC(1)=iAd
      Call iDaFile(Lu,1,[nRoots],1,iAd)
      TOC(2)=iAd
      Call iDaFile(Lu,1,[nGrad],1,iAd)
      TOC(3)=iAd
      Call iDaFile(Lu,1,i_grad,nRoots,iAd)
      TOC(4)=iAd
      Call iDaFile(Lu,1,i_nac,nCoup,iAd)
      TOC(5)=iAd
      iAd=0
      Call iDaFile(Lu,1,TOC,Size(TOC),iAd)
      Call DaClos(Lu)
*
      Call mma_Deallocate(i_grad)
      Call mma_Deallocate(i_nac)
*
      End Subroutine Create_Grads
