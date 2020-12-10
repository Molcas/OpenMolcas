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
      Subroutine G_Nrm(nAtom,nInter,GNrm,Iter,Grad,mIntEff)
      use Slapaf_Info, only: Gx, Degen
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
#include "Molcas.fh"
      Real*8 GNrm(Iter),Grad(nInter,Iter)
*
*
*     Compute the norm of the cartesian force vector.
*
*     |dE/dx|=Sqrt(dE/dx|u|dE/dx)
*
      Fabs=0.0D0
      Do i = 1, nAtom
         Do j = 1, 3
            Fabs=Fabs+ Degen(j,i)*Gx(j,i,Iter)**2
         End Do
      End Do
      Fabs=Sqrt(Fabs)
*#define _DEBUGPRINT_
#ifdef _DEBUGPRINT_
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
#ifdef _DEBUGPRINT_
      Write (6,*) ' mIntEff=',mIntEff
#endif
*
      Return
      End
