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
      Subroutine CSFGRad(Grad,nGrad)
************************************************************************
*                                                                      *
* Object: to compute the CSF component of the non-adiabatic derivative *
*         coupling                                                     *
*                                                                      *
* This routine assumes C1 symmetry and needs the 'D1ao-' matrix stored *
* in the runfile                                                       *
*                                                                      *
************************************************************************
      Implicit None
#include "itmax.fh"
#include "info.fh"
#include "stdalloc.fh"
#include "real.fh"
#include "nac.fh"
      Integer nGrad
      Real*8 Grad(nGrad)
*
      Integer nD,nB,lOper(1)
      Real*8 CCoor(3)
      Real*8, Dimension(:), Allocatable :: aDAO
      Logical Found
      Character(Len=80) Label
      External OvrGrd, OvrMmG

      Call DCopy_(nGrad,[Zero],0,Grad,1)

      nB=nBas(0)
      Call Qpg_dArray('D1ao-',Found,nD)
      Call mma_allocate(aDAO,nD)
      Call Get_dArray('D1ao-',aDAO,nD)
*     Call TriPrt('DAO-','',aDAO,nB)

*IFG Compute the CSF contribution to the coupling vector.
*    Inner product of S[x] and D^A (antisymmetric component of transition density matrix)
*    This is the same as the product of S[x]^A and D
      isCSF=.True.
      Call DCopy_(3,[Zero],0,CCoor,1)
      lOper(1)=1
      Label='The CSF Contribution'
      Call OneEl_g(OvrGrd,OvrMmG,Grad,nGrad,.False.,CCoor,
     &           aDAO,nD,lOper,1,0,Label)
      isCSF=.False.

      Call mma_deallocate(aDAO)

      End
