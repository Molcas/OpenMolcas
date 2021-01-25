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
      Subroutine SCF_Energy(FstItr,E1_,E2_,EV)
      Use SCF_Arrays
      Use Interfaces_SCF, Only: PMat_SCF
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
      Logical FstItr
      Real*8, Dimension (:), Allocatable :: XCf
*
*
      nD = iUHF + 1
*
*     Allocate memory for coefficients for minimized densities.
*
      nXCf = MxIter*nD
      Call mma_allocate(XCf,nXCf)
*                                                                      *
************************************************************************
*                                                                      *
*     Compute
*     (1) density matrix,
*     (2) two-electron part of Fock matrix, and
*     (3) variational energies.
*                                                                      *
************************************************************************
*                                                                      *
*     (1) 1-particle density
*
      Call DMat(Dens,TwoHam,nBT,nDens,CMO,nBO,OccNo,nnB,nD,Ovrlp,
     &          XCf,nXCf,Vxc)
*                                                                      *
************************************************************************
*                                                                      *
*     (2) 2-electron part of the Fock matrix, and the potential of the
*         external field (none linear or bi-linear).
*
*         Affects data in position iPsLst
*
      Call PMat_SCF(Dens,OneHam,TwoHam,nBT,nDens,nXCf,FstItr,XCf,nD,
     &              EDFT,MxIter,Vxc,Fock)
*
      Call mma_deallocate(XCf)
*                                                                      *
************************************************************************
*                                                                      *
*     (3) the energy.
*
      Call EneClc(E1_,E2_,EV,Dens,OneHam,TwoHam,nBT,nDens,nD,EDFT,
     &            MxIter)
*
      Return
      End
