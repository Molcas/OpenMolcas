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
* Copyright (C) 2000, Gunnar Karlstrom                                 *
*               2000, Roland Lindh                                     *
************************************************************************
      Real*8 Function Anal_Gitt(cordsi,latato)
************************************************************************
*                                                                      *
*     Object:                                                          *
*                                                                      *
*     Authors: G. Karlstroem                                           *
*              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
*                                                                      *
*              and                                                     *
*                                                                      *
*              R. Lindh                                                *
*              Dept. of Chem. Phys., Univ. of Lund, Sweden.            *
*                                                                      *
*              March 2000                                              *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "real.fh"
      Real*8 cordsi(3,latato)
*
      Anal_Gitt=Zero
*
*.....analyze the lattice
*
      gatom=Zero
      Do i=1,latato
         faktor=One
         x1=cordsi(1,i)+One
         y1=cordsi(2,i)
         z1=cordsi(3,i)
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do

         x1=cordsi(1,i)-One
         y1=cordsi(2,i)
         z1=cordsi(3,i)
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do

         x1=cordsi(1,i)
         y1=cordsi(2,i)+One
         z1=cordsi(3,i)
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do

         x1=cordsi(1,i)
         y1=cordsi(2,i)-One
         z1=cordsi(3,i)
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do

         x1=cordsi(1,i)
         y1=cordsi(2,i)
         z1=cordsi(3,i)+One
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do

         x1=cordsi(1,i)
         y1=cordsi(2,i)
         z1=cordsi(3,i)-One
         Do j=1,latato
            r2 = (x1-cordsi(1,j))**2
     &         + (y1-cordsi(2,j))**2
     &         + (z1-cordsi(3,j))**2
            If (r2.lt.0.01D0) faktor=faktor+One
         End Do
         gatom=gatom+One/faktor
      End Do
      Anal_Gitt=gatom
*
      Return
      End
