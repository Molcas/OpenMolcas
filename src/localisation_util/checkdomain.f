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
      SubRoutine CheckDomain(irc,iDomain,nAtom,nOcc)
C
C     Thomas Bondo Pedersen, January 2006.
C
C     Purpose: check domain definition.
C
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom,nOcc)

      irc = 0
      Do i = 1,nOcc
         If (iDomain(0,i).lt.1 .or. iDomain(0,i).gt.nAtom) Then
            Write(6,*) 'Dimension of domain ',i,': ',iDomain(0,i)
            irc = irc + 1
         Else
            Do iAt = 1,iDomain(0,i)
               iAtom = iDomain(iAt,i)
               If (iAtom.lt.1 .or. iAtom.gt.nAtom) Then
                  Write(6,*) 'Atom ',iAt,' of domain ',i,': ',iAtom
                  irc = irc + 1
               End If
            End Do
         End If
      End Do

      End
