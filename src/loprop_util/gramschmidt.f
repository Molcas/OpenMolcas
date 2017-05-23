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
      Subroutine GramSchmidt(S,C,nDim,Type,Center,iRestrict)
*
      Implicit ReaL*8 (A-H,O-Z)
#include "real.fh"
      Real*8 S(nDim,nDim),C(nDim,nDim)
      Integer Type(nDim), Occ, Vir,
     &        Center(nDim)
      Parameter(Occ=1,Vir=0)
*
c     write (*,*) 'nDim,nDim', nDim,nDim
C     Write (6,*) Center
      Do iOrb=1,nDim
         If (iRestrict .eq. 1 .AND. Type(iOrb).eq.Vir) Go To 99
         F = Zero
         If (S(iOrb,iOrb) .gt. Zero) F = One/Sqrt(S(iOrb,iOrb))
C
         Do ibas=1,nDim
            C(Ibas,Iorb)=F*C(Ibas,Iorb)
C        write (6,*) C(Ibas,Iorb), IBas, Iorb
         End Do
C
         Do Jorb=1,nDim
            S(Iorb,Jorb)=F*S(Iorb,Jorb)
         End Do
         Do Jorb=1,nDim
            S(Jorb,Iorb)=F*S(Jorb,Iorb)
         End Do
*
         iStart = 1
         If (iRestrict .eq. 0) iStart = iOrb+1
         Do jOrb=iStart,nDim
            If (iRestrict .eq. 1 .AND. Type(jOrb).eq.Occ) Go To 98
C           If (Center(iOrb).eq.Center(jOrb)) Go To 98
            A=S(Iorb,Jorb)
            Do Ibas=1,nDim
               C(Ibas,Jorb)=C(Ibas,Jorb)-A*C(Ibas,Iorb)
            End Do
            Do Korb=1,nDim
               S(Jorb,Korb)=S(Jorb,Korb)-A*S(Iorb,Korb)
            End Do
            Do Korb=1,nDim
               S(Korb,Jorb)=S(Korb,Jorb)-A*S(Korb,Iorb)
            End Do
 98         Continue
         End Do
 99      Continue
      End Do
      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer_array(Center)
      End
