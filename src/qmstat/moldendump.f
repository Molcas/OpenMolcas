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
      Subroutine MoldenDump(iC,CooRef,nP,nA,nC)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "WrkSpc.fh"
#include<constants.fh>
#include "real.fh"

      Parameter (AuAng=1d10*CONST_BOHR_RADIUS_IN_SI_)

      Dimension CooRef(MxCen,3),Coo(MxCen,3)
      Dimension iC(3)

      Logical ValidOrNot

*
*-- Clarifying words.
*
      Write(6,*)
      Write(6,*)
      Write(6,*)'   * Coordinates given in form for Molden *'
      Write(6,*)
      Write(6,*)' Put everything within the lines in a separate file'
     &//' and view with Molden.'
      Write(6,*)' Observe that the identity of molecules that are not'
     &//' valid water molecules is unknown.'
      Write(6,*)
      Write(6,*)'------------------------------------------------------'
     &//'------------------------------'

*
*-- Print total number of particles.
*
      Write(6,*)'  Substitue this line with number of atoms.'
      Write(6,*)
      Do 101, iP=1,nP
        ind=nC*(iP-1)
        Do 102, jC=1,nC
          Coo(jC,1)=Work(iC(1)+ind+jC-1)
          Coo(jC,2)=Work(iC(2)+ind+jC-1)
          Coo(jC,3)=Work(iC(3)+ind+jC-1)
102     Continue
        Call IsItValid(Coo,CooRef,ValidOrNot)
        If(.not.ValidOrNot) then
          Do 103, jC=1,nC
            Write(6,92)'C  ',(AuAng*Coo(jC,kk),kk=1,3)
103       Continue
        Else
          Write(6,92)'O  ',(AuAng*Coo(1,kk),kk=1,3)
          Write(6,92)'H  ',(AuAng*Coo(2,kk),kk=1,3)
          Write(6,92)'H  ',(AuAng*Coo(3,kk),kk=1,3)
        Endif
101   Continue
      Write(6,*)
      Write(6,*)'------------------------------------------------------'
     &//'------------------------------'

*
*-- Formats
*
92    Format(A,3(F10.6))

      Return
c Avoid unused argument warnings
      If (.False.) Call Unused_integer(nA)
      End
