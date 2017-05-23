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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Real*8 Function LDF_InteractionRange(A)
C
C     Thomas Bondo Pedersen, October 2010.
C
C     Purpose: return the max interatomic distance (in bohr) for atom
C              pairs in which atom A takes part.
C
C     Note: Atom to atom pair map must have been set up
C           (call LDF_SetA2AP()); if not, -1.0d0 is returned.
C
      Implicit None
      Integer A
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "ldf_a2ap.fh"

#if defined (_DEBUG_)
      Character*20 SecNam
      Parameter (SecNam='LDF_InteractionRange')
#endif

      Real*8   LDF_AtomicDistance
      External LDF_AtomicDistance

      Real*8  Rmax
      Integer n, ip, iAB, AB, A_, B_

      Integer i, j
      Integer A2AP
      Integer AP_Atoms
      A2AP(i,j)=iWork(ip_A2AP-1+2*(j-1)+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      LDF_InteractionRange=-1.0d0
      If (l_A2AP.gt.0) Then
         Rmax=-1.0d0
         n=A2AP(1,A)
         ip=A2AP(2,A)-1
         Do iAB=1,n
            AB=iWork(ip+iAB)
            A_=AP_Atoms(1,AB)
            B_=AP_Atoms(2,AB)
#if defined (_DEBUG_)
            If (A_.ne.A .and. B_.ne.A) Then
               Call WarningMessage(2,
     &                          SecNam//': A2AP info not properly set!')
               Call LDF_Quit(1)
            End If
#endif
            Rmax=max(Rmax,LDF_AtomicDistance(A_,B_))
         End Do
         LDF_InteractionRange=Rmax
      End If

      End
