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
      SubRoutine GT1DIA_MCLR(H1DIA)
      use Arrays, only: FIMO
      Implicit Real*8 (a-h,o-z)

#include "Input.fh"
#include "Pointers.fh"
#include "glbbas_mclr.fh"
      Real*8 H1DIA(*)
      i=1
      Do iS=1,nSym
       ii=ipCM(iS)+nOrb(iS)*(nIsh(iS)-1)+nIsh(iS)-1
       Do iAsh=1,nAsh(iS)
        ii=ii+nOrb(iS)+1
        H1DIA(i) = FIMO(ii)
        i=i+1
       End Do
      End Do
      Return
      End
