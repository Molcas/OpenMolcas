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
      SubRoutine Domain_Histogram(iDomain,nAtom,nOcc,Title)
C
C     Thomas Bondo Pedersen, January 2006.
C
C     Purpose: print histogram of domain sizes.
C
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom,nOcc)
      Character*(*) Title
#include "WrkSpc.fh"

      If (nAtom.lt.1 .or. nOcc.lt.1) Return

      i_min = iDomain(0,1)
      i_max = iDomain(0,1)
      x_ave = dble(iDomain(0,1))
      Do i = 2,nOcc
         i_min = min(i_min,iDomain(0,i))
         i_max = max(i_max,iDomain(0,i))
         x_ave = x_ave + dble(iDomain(0,i))
      End Do
      x_ave = x_ave/dble(nOcc)

      l_iCount = i_max - i_min + 1
      Call GetMem('Dm_Histo','Allo','Inte',ip_iCount,l_iCount)

      Call Cho_Head(Title,"=",80,6)
      Write(6,'(/,A,3X,I10,/,A,3X,I10,/,A,F13.2)')
     & 'Minimum size:',i_min,
     & 'Maximum size:',i_max,
     & 'Average size:',x_ave
      Call Domain_Histo1(iDomain,nAtom,nOcc,iWork(ip_iCount),i_min,
     &                   i_max)
      Fac = 1.0d2/dble(nOcc)
      Pct = Fac*dble(iWork(ip_iCount))
      i = i_min
      Write(6,'(/,A,I10,A,I10,3X,F7.2,A)')
     & 'Number with size',i,':',iWork(ip_iCount),Pct,'%'
      Do iC = 2,l_iCount
         Pct = Fac*dble(iWork(ip_iCount-1+iC))
         i = i + 1
         Write(6,'(A,I10,A,I10,3X,F7.2,A)')
     &   'Number with size',i,':',iWork(ip_iCount-1+iC),Pct,'%'
      End Do

      Call GetMem('Dm_Histo','Free','Inte',ip_iCount,l_iCount)

      End
      SubRoutine Domain_Histo1(iDomain,nAtom,nOcc,iCount,i_min,i_max)
      Implicit Real*8 (a-h,o-z)
      Integer iDomain(0:nAtom,nOcc), iCount(*)

      nC = i_max - i_min + 1
      Call iCopy(nC,0,0,iCount,1)

      Do i = 1,nOcc
         iC = iDomain(0,i) - i_min + 1
         iCount(iC) = iCount(iC) + 1
      End Do

      End
