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
      Subroutine G_Nrm(GrdX,nAtom,nInter,GNrm,Iter,
     &                 Grad,Degen,mIntEff)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
      Real*8 GrdX(3*nAtom),GNrm(Iter),Grad(nInter,Iter),Degen(3*nAtom)
*
*
*     Compute the norm of the cartesian force vector.
*
*     |dE/dx|=Sqrt(dE/dx|u|dE/dx)
*
      Fabs=0.0D0
      Do i = 1, 3*nAtom
         Fabs=Fabs+ Degen(i)*GrdX(i)**2
      End Do
      Fabs=Sqrt(Fabs)
*#define _DEBUG_
#ifdef _DEBUG_
         Write (6,42) Fabs
42       Format (/,' Norm of the force vector',F20.15)
#endif
      GNrm(iter)=Fabs
*
*-----Write out the internal force vector.
*
      mIntEff=0
      Do i = 1, nInter
         If (Abs(Grad(i,Iter)).gt.1.0d-6) mIntEff=mIntEff+1
      End Do
      If (mIntEff.eq.0) mIntEff=1
#ifdef _DEBUG_
      Write (6,*) ' mIntEff=',mIntEff
#endif
*
      Return
      End
