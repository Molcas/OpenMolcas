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
      use ChoArr, only: iSP2F
      Implicit Real*8 (a-h,o-z)
#include "cholesky.fh"
#include "stdalloc.fh"

      Real*8, Allocatable:: TMax(:,:)

      iTri(i,j)=max(i,j)*(max(i,j)-3)/2+i+j

      If (Cho_PreScreen) Then ! prescreening with approx. diagonal

         Call mma_allocate(Tmax,nShell,nShell,Label='nShell')

         Call Shell_MxSchwz(nShell,Tmax)
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
         Call mma_allocate(iSP2F,nnShl,Label='iSP2F')

         ij = 0
         Do i = 1,nShell
            Do j = 1,i
               If (Tmax_All*Tmax(i,j) .gt. Tau) Then
                  ij = ij + 1
                  iSP2F(ij) = iTri(i,j)
               End If
            End Do
         End Do

         Call mma_deallocate(TMax)

      Else ! no prescreening, include all shell pairs.

         nnShl = nnShl_Tot
         Call mma_allocate(iSP2F,nnShl,Label='iSP2F')

         Do ij = 1,nnShl
            iSP2F(ij) = ij
         End Do

      End If

#if defined (_DEBUGPRINT_)
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
