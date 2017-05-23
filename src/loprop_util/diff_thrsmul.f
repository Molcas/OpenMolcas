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
      Subroutine Diff_ThrsMul(ipMP,ThrsMul,ThrsMul_Clever,nAt,nij,lMax)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      dMMax=0.0d0
      kauntA=0
      Do iAtom=1,nAt
        Do jAtom=1,iAtom
          kaunt=0
          Do l=0,1
            kComp=(l+1)*(l+2)/2
            Do k=1,kComp
              dM=Work(ipMP+nij*kaunt+kauntA)
              If(abs(dM).gt.dMMax)dMMax=abs(dM)
              kaunt=kaunt+1
            Enddo
          Enddo
          kauntA=kauntA+1
        Enddo
      Enddo

      ThrsMul_Clever=dMMax*ThrsMul

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(lMax)
      End
