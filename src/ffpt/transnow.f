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
      Subroutine TransNow(iV,ipS)
      Implicit Real*8 (a-h,o-z)

#include "input.fh"
#include "WrkSpc.fh"


*
*-- Be informative.
*
      If(ComStk(2,1,1,0)) then
        Write(6,*)
        Write(6,*)'    The DIPO perturbation is translated.'
      Else
        Write(6,*)
        Write(6,*)'    No translation of the perturbation',
     &    ' (only implemented for DIPO.'
        Write(6,*)'OBSERVE! Your result can be origo dependent!'
      Endif

*
*-- Translate the perturbation. Why, oh why, is the sign of Trana
*   the opposite of what it should be according to the formula?
*   Answer: Our charge is positive in the MLTPL integrals, hence
*   the field strength should be of opposite sign, see ptdipo.f
*   for example. Hence the minus in the translation formula becomes
*   a plus. Oh yeah!
*
      If((.not.ComStk(2,1,1,1)).and.(.not.ComStk(2,1,1,2)).and.
     &     (.not.ComStk(2,1,1,3))) Then
        Write(6,*)
        Write(6,*)'A strange error has occured. ComStk modified?'
        Call Abend()
      Endif
      kaunter=0
c      Write(6,*)'Trancoo',(TranCoo(i),i=1,3)
      Do 1009, i=1,nBas(1)
        Do 1008, j=1,i
          Do 1007, iXYZ=1,3
            If(ComStk(2,1,1,iXYZ)) Then
              Trana=ComVal(2,1,1,iXYZ)*TranCoo(iXYZ)*Work(ipS+kaunter)
*-----------Sign of Trana? See source code comment above.
              Work(iV+kaunter)=Work(iV+kaunter)+Trana
            EndIf
 1007     Continue
          kaunter=kaunter+1
 1008   Continue
 1009 Continue
      Return
      End
