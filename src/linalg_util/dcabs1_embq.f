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
      real*8 function dcabs1_emb(z)
      complex*16 z,zz
      real*8 t(2)
      equivalence (zz,t(1))
      zz = z
      dcabs1_emb = abs(t(1)) + abs(t(2))
      return
      end
