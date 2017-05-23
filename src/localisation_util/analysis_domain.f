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
* Copyright (C) 2006, Thomas Bondo Pedersen                            *
************************************************************************
      SubRoutine Analysis_Domain(iDomain,QD,f,Coord,AtomLbl,
     &                           nBas_Start,nAtom,nBas,nOcc)
C
C     Thomas Bondo Pedersen, January 2006.
C
C     Purpose: analyze orbital domains.
C
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom,nOcc)
      Real*8  QD(nOcc), f(nOcc), Coord(3,nAtom)
      Character*4 AtomLbl(2,nBas)
      Integer nBas_Start(nAtom)
#include "WrkSpc.fh"

      If (nAtom.lt.1 .or. nOcc.lt.1) Return

      Call Cho_Head('Orbital domain analysis','=',80,6)
      Do i = 1,nOcc
         nAt = iDomain(0,i)
         Rmin = 1.0d15
         Rmax = -1.0d15
         Rave = 0.0d0
         nij = 0
         Do jAt = 1,nAt-1
            jAtom = iDomain(jAt,i)
            Do iAt = jAt+1,nAt
               iAtom = iDomain(iAt,i)
               R = sqrt((Coord(1,iAtom)-Coord(1,jAtom))**2
     &                 +(Coord(2,iAtom)-Coord(2,jAtom))**2
     &                 +(Coord(3,iAtom)-Coord(3,jAtom))**2)
               Rmin = min(Rmin,R)
               Rmax = max(Rmax,R)
               Rave = Rave + R
               nij = nij + 1
            End Do
         End Do
         If (nij .eq. 0) Then
            Rmax = 0.0d0
            Rmin = 0.0d0
         Else
            Rave = Rave/dble(nij)
         End If
         Write(6,'(/,A,I6,A,I6)')
     &   'Orbital domain',i,':  size:',nAt
         Write(6,'(A,1P,2(1X,D15.5))')
     &   '  Charge, completeness function:',QD(i),f(i)
         Write(6,'(A,1P,3(1X,D15.5))')
     &   '  Rmin, Rmax, Rave             :',Rmin,Rmax,Rave
         Do iAt = 1,nAt
            iAtom = iDomain(iAt,i)
            Write(6,'(A,I6,2X,A,1X,3(1X,F12.3))')
     &      '  Atom:',iAtom,AtomLbl(1,nBas_Start(iAtom)),
     &      (Coord(j,iAtom),j=1,3)
         End Do
      End Do

      End
