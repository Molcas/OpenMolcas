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
      REAL*8 FUNCTION GTIJKL_MCLR(I,J,K,L)
      use Arrays, only: Int2
*
* Obtain  integral (I J ! K L )
* where I,J,K and l refers to active orbitals in type ordering
*
      IMPLICIT None
      Integer, Intent(In):: I, J, K, L
#include "detdim.fh"
#include "orbinp_mclr.fh"
      Integer iAbs, jAbs, kAbs, lAbs, ij, kl
      Integer itri

      itri(i,j)=Max(i,j)*(Max(i,j)-1)/2+Min(i,j)

      IABS = IREOTS(I)
*
      JABS = IREOTS(J)
*
      KABS = IREOTS(K)
*
      LABS = IREOTS(L)
*
      IJ=itri(iABS,JABS)
      KL=itri(kABS,lABS)

      GTIJKL_MCLR = INT2(itri(IJ,KL))
*
      END FUNCTION GTIJKL_MCLR
