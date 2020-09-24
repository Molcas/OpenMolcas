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
      SubRoutine Cho_DiaSP()
C
C     Thomas Bondo Pedersen, March 2006.
C
C     Purpose: prescreening of diagonal.
C
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "choptr.fh"
#include "WrkSpc.fh"

      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j
      Tmax(i,j)=Work(ip_Tmax-1+nShell*(j-1)+i)

#if defined (_DEBUG_)
#endif

      If (Cho_PreScreen) Then ! prescreening with approx. diagonal

         l_Tmax = nShell**2
         Call GetMem('Cho_Tmax','Allo','Real',ip_Tmax,l_Tmax)

         Call Shell_MxSchwz(nShell,Work(ip_Tmax))
         Tmax_All = Tmax(1,1)
         Do i = 2,nShell
            Do j = 1,i
               Tmax_All = max(Tmax_All,Tmax(i,j))
            End Do
         End Do

         Tau = Thr_PreScreen
         nnShl = 0
         Do i = 1,nShell
            Do j = 1,i
               If (Tmax_All*Tmax(i,j) .gt. Tau) Then
                  nnShl = nnShl + 1
               End If
            End Do
         End Do
         l_iSP2F = nnShl
         Call GetMem('SP2F','Allo','Inte',ip_iSP2F,l_iSP2F)

         ij = 0
         Do i = 1,nShell
            Do j = 1,i
               If (Tmax_All*Tmax(i,j) .gt. Tau) Then
                  ij = ij + 1
                  iWork(ip_iSP2F-1+ij) = iTri(i,j)
               End If
            End Do
         End Do

         Call GetMem('Cho_Tmax','Free','Real',ip_Tmax,l_Tmax)

      Else ! no prescreening, include all shell pairs.

         nnShl = nnShl_Tot
         l_iSP2F = nnShl
         Call GetMem('SP2F','Allo','Inte',ip_iSP2F,l_iSP2F)

         Do ij = 1,nnShl
            iWork(ip_iSP2F-1+ij) = ij
         End Do

      End If

#if defined (_DEBUG_)
      If (.not.Cho_PreScreen) Tau = 0.0d0
      Write(LuPri,*) '>>> Exit from Cho_DiaSP:'
      Write(LuPri,*) '    Screening threshold               : ',Tau
      Write(LuPri,*) '    Total number of shell pairs       : ',
     &                nnShl_Tot
      Write(LuPri,*) '    Contributing number of shell pairs: ',
     &                nnShl
      If (nnShl_Tot .ne. 0) Then
         Write(LuPri,*) '    Screening-%: ',
     &                   1.0d2*DBLE(nnShl_Tot-nnShl)/DBLE(nnShl_Tot)
      End If
#endif

      End
