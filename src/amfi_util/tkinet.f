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
      REAL*8 function Tkinet(l,alpha1,alpha2)
cbs   calculates the matrix element of kinetic energy
cbs   for primitive normalized functions with the same angular momentum l
cbs   and exponents alpha1 and alpha2
cbs   works only, if r**l is assumed for an l-value
cbs   formular obtained from the symmetric expression (d/dr's to (')
cbs   the left and to the right.
cbs   Overlaps of the different powers are partially crossed out
cbs   with  the overlap of functions with angular momentum l
cbs   final formula:
cbs   Tkinet=0.5*alpha12 (2l+3) (alpha1*alpha2/alpha12*alpha12)**((2L+7)/4)
cbs   with alpha12=0.5*(alpha1+alpha2)
cbs   as alpha12 has the dimensions 1/length**2, this can not be that bad...
      Implicit REAL*8 (a-h,o-z)
Cbs   alpha12 is the effective exponent
      Alpha12=0.5d0*(alpha1+alpha2)
      alphpro=alpha1*alpha2
      ll3=l+l+3
      ll7=l+l+7
      Tkinet=0.5d0*alpha12*DBLE(ll3)*(alphpro/
     *(alpha12*alpha12))**(0.25D0*DBLE(ll7))
      return
      end
