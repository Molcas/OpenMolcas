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
      Subroutine OutCoor(TEXT,Char,nDim,FI,N1,N2,lAngstroms)
************************************************************************
*                                                                      *
*     Object: To generate a cartesian output with atomic labels        *
*             N1 and N2 are the real limits of dummy FI                *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "angstr.fh"
      Character*(*) TEXT
      Character*(*) Char(nDim)
      Logical lAngstroms
      Real*8 FI(N1,N2)
*
      Call qEnter('OutCoor')
*
      Lu=6
      Write (Lu,*)
      Write (Lu,*) '******************************************'//
     &             '***************'
      Write (Lu,*) Text
      Write (Lu,*) '******************************************'//
     &             '***************'
      Write (Lu,*) ' ATOM              X               Y      '//
     &             '         Z     '
      Do 10 I = 1, NDIM
         If (lAngstroms) then
            Write (Lu,300) Char(I),(FI(J,I)*angstr,J=1,3)
         else
            Write (Lu,300) Char(I),(FI(J,I),J=1,3)
         EndIf
300      Format (2X,A,3X,3F16.6)
10    Continue
*
      Write (Lu,*)
      Call qExit('OutCoor')
      Return
      End
