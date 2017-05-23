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
      Subroutine YouGetThis(nAt,EC,Pot_Expo,Pot_Point,Pot_Fac,Diffed
     &                     ,ipMP,lMax,lMaxF,nij,LuYou)
      Implicit Real*8 (a-h,o-z)

#include "WrkSpc.fh"

      Dimension EC(3,nij),Pot_Expo(nij*4),Pot_Point(nij),Pot_Fac(nij*4)

      Logical Diffed(nij*4)

*
*-- Number of centres and maximal angular momentum.
*
      Write(LuYou,101)nij
      Write(LuYou,102)lMax,lMaxF
      kauntA=0
*
*-- Loop over centres.
*
      Do i=1,nij
        kauntA=kauntA+1
*
*---- Coordinates of centre.
*
        Write(LuYou,103)(EC(k,i),k=1,3)
        kaunt=0
        Do l=0,lMax
          nS=l*(l+1)*(l+2)/6
          nT=(l+1)*(l+2)*(l+3)/6
          If(l.le.lMaxF) then
            If(Diffed(2*(kauntA-1)+l+1)) then
*
*---- Factor and exponent.
*
              Write(LuYou,104)2.0d0*Pot_Expo(2*(kauntA-1)+l+1)
              Write(LuYou,105)(Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
            Else
*
*---- Factor and dummy-exponent, if this multipole should not be
*     made diffuse.
*

              Write(LuYou,104)-7.91204d0
              Write(LuYou,105)(Pot_Fac(4*(kauntA-1)+kk),kk=nS+1,nT)
            Endif
          Else
*
*---- The pure multipole, which under no circumstance can be diffuse.
*
            Write(LuYou,104)-7.91204d0
            Write(LuYou,105)(Work(ipMP+nij*kk+kauntA-1),kk=nS,nT-1)
          Endif
        Enddo
*
*---- Throw in the point-part.
*
        Write(LuYou,104)Pot_Point(kauntA)
      Enddo

101   Format(I5)
102   Format(2I5)
103   Format(3(F20.14))
104   Format(F20.14)
105   Format(3(F20.14))

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nAt)
      End
