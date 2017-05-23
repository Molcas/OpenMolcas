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
      Subroutine LDF_UnsetAtomInfo(irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Unset atom info stored in ldf_atom_info.fh
C              I.e. free memory allocated for this information and
C              garble the rest.
C
C     Return codes:
C       irc=0: success
C       irc=1: not set, i.e. info already unset and nothing is done
C
      Implicit None
      Integer irc
#include "ldf_atom_info.fh"
#include "WrkSpc.fh"

      Character*6 Label

      Integer iAtom, ip, l

      Integer i, j
      Integer A_Shells, A_AuxShells
      A_Shells(i,j)=iWork(ip_A_Shells-1+2*(j-1)+i)
      A_AuxShells(i,j)=iWork(ip_A_AuxShells-1+2*(j-1)+i)

      ! init return code
      irc=0

      ! Check status
      If (LDF_AtomInfo_Status.eq.LDF_AtomInfo_Unset) Then
         Call WarningMessage(1,'LDF_UnsetAtomInfo: Info already unset!')
         irc=1
         Return
      End If

      ! Deallocate auxiliary shell lists
      Do iAtom=1,NumberOfAtoms
         l=A_AuxShells(1,iAtom)
         If (l.gt.0) Then
            Write(Label,'(A,I4.4)') 'AA',iAtom-1
            ip=A_AuxShells(2,iAtom)
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
      End Do

      ! Deallocate valence shell lists
      Do iAtom=1,NumberOfAtoms
         l=A_Shells(1,iAtom)
         If (l.gt.0) Then
            Write(Label,'(A,I4.4)') 'SA',iAtom-1
            ip=A_Shells(2,iAtom)
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
      End Do

      ! Deallocate A_AuxShells
      Call GetMem('A_AuxShells','Free','Inte',ip_A_AuxShells,
     &                                         l_A_AuxShells)
      ip_A_AuxShells=0
      l_A_AuxShells=0

      ! Deallocate A_Shells
      Call GetMem('A_Shells','Free','Inte',ip_A_Shells,l_A_Shells)
      ip_A_Shells=0
      l_A_Shells=0

      ! Deallocate A_Unique
      Call GetMem('A_Unique','Free','Inte',ip_A_Unique,l_A_Unique)
      ip_A_Unique=0
      l_A_Unique=0

      ! Deallocate atomic coordinates
      Call GetMem('LDF_Coord','Free','Real',ip_Coord,l_Coord)
      ip_Coord=0
      l_Coord=0

      ! Zero NumberOfAtoms
      NumberOfAtoms=0

      ! Mark info unset
      LDF_AtomInfo_Status=LDF_AtomInfo_Unset

      End
