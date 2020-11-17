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
      Function d_cart(Ind,nStab,mxdc,nSym)
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Integer nStab(mxdc), Ind(1)
      Real*8 d_cart
*                                                                      *
************************************************************************
*                                                                      *
*---- Cartesian coordinate  (iAtom)
*
      D_Cart=Zero
*
      iAtom=Ind(1)
*
      nU_A=nStab(iAtom)
*
*-----Now evaluate the degeneracy of the cartesian
*
      iDeg=nSym/nU_A
      d_cart=DBLE(iDeg)
*
*     Write (*,*) ' d_cart=',d_cart
*
      Return
      End
