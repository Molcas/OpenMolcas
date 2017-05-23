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
      Subroutine OffAtom(C1,C2,C3,C4,C5)
      Implicit Real*8 (a-z)

#include<constants.fh>
#include "real.fh"
      Parameter (AuAng=1d10*CONST_BOHR_RADIUS_IN_SI_)
      Dimension C1(3),C2(3),C3(3),C4(3),C5(3)
      Dimension D(3),E(3),U(3),V(3)

      D(1)=(C2(1)+C3(1))/2-C1(1)
      D(2)=(C2(2)+C3(2))/2-C1(2)
      D(3)=(C2(3)+C3(3))/2-C1(3)
      LD=sqrt(D(1)**2+D(2)**2+D(3)**2)
      D(1)=D(1)/LD
      D(2)=D(2)/LD
      D(3)=D(3)/LD

      U(1)=C2(1)-C1(1)
      U(2)=C2(2)-C1(2)
      U(3)=C2(3)-C1(3)
      V(1)=C3(1)-C1(1)
      V(2)=C3(2)-C1(2)
      V(3)=C3(3)-C1(3)
      E(1)=U(2)*V(3)-V(2)*U(3)
      E(2)=U(3)*V(1)-V(3)*U(1)
      E(3)=U(1)*V(2)-V(1)*U(2)
      LE=sqrt(E(1)**2+E(2)**2+E(3)**2)
      E(1)=E(1)/LE
      E(2)=E(2)/LE
      E(3)=E(3)/LE

      C4(1)=C1(1)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(1)
     &           +(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(1)
      C4(2)=C1(2)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(2)
     &           +(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(2)
      C4(3)=C1(3)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(3)
     &           +(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(3)
      C5(1)=C1(1)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(1)
     &           -(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(1)
      C5(2)=C1(2)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(2)
     &           -(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(2)
      C5(3)=C1(3)+(0.2767d0/AuAng)*cos(36.72d0*2d0*Pi/360d0)*D(3)
     &           -(0.2767d0/AuAng)*sin(36.72d0*2d0*Pi/360d0)*E(3)

      Return
      End
