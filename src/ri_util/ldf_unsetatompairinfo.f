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
      Subroutine LDF_UnsetAtomPairInfo(irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Unset atom pair info stored in ldf_atom_pair_info.fh
C              I.e. free memory allocated for this information and
C              garble the rest.
C
C     Return codes:
C       irc=0: success
C       irc=1: not set, i.e. info already unset and nothing is done
C
      Implicit None
      Integer irc
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*8 Label

      Integer iAtomPair
      Integer ip, l

      Integer  LDF_AtomPair_DiagDim
      External LDF_AtomPair_DiagDim

      Integer i, j
      Integer AP_2CFunctions, AP_1CLinDep
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)

      ! Init return code
      irc=0

      ! Check status
      If (LDF_AtomPairInfo_Status.eq.LDF_AtomPairInfo_Unset) Then
         Call WarningMessage(2,'LDF_UnsetAtomPairInfo: already unset!')
         irc=1
         Return
      End If

      ! Deallocate disk addresses for C
      Call GetMem('AP_DiskC','Free','Inte',ip_AP_DiskC,l_AP_DiskC)
      ip_AP_DiskC=0
      l_AP_DiskC=0

      ! Deallocate unique atom pairs
      Call GetMem('AP_Unique','Free','Inte',ip_AP_Unique,l_AP_Unique)
      ip_AP_Unique=0
      l_AP_Unique=0

      ! Deallocate AP_2CFunctions
      Do iAtomPair=1,NumberOfAtomPairs
         l=4*AP_2CFunctions(1,iAtomPair)
         If (l.gt.0) Then
            Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
            ip=AP_2CFunctions(2,iAtomPair)
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
      End Do
      Call GetMem('AP2CFN','Free','Inte',ip_AP_2CFunctions,
     &                                    l_AP_2CFunctions)
      ip_AP_2CFunctions=0
      l_AP_2CFunctions=0

      ! Deallocate AP_1CLinDep
      Do iAtomPair=1,NumberOfAtomPairs
         l=3*AP_1CLinDep(1,iAtomPair)
         If (l.gt.0) Then
            Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
            ip=AP_1CLinDep(2,iAtomPair)
            Call GetMem(Label,'Free','Inte',ip,l)
         End If
      End Do
      Call GetMem('AP1CLD','Free','Inte',ip_AP_1CLinDep,l_AP_1CLinDep)
      ip_AP_1CLinDep=0
      l_AP_1CLinDep=0

      ! Deallocate diagonal blocks
      Call LDF_DeallocateBlockMatrix('APD',ip_AP_Diag)
      ip_AP_Diag=0
      l_AP_Diag=0
      Call LDF_DeallocateBlockMatrix('APB',ip_AP_DiagBak)
      ip_AP_DiagBak=0
      l_AP_DiagBak=0

      ! Deallocate AP_Atoms
      Call GetMem('LDFAPA','Free','Inte',ip_AP_Atoms,l_AP_Atoms)
      ip_AP_Atoms=0
      l_AP_Atoms=0

      ! Zero number of atom pairs
      NumberOfAtomPairs=0

      ! Mark info unset
      LDF_AtomPairInfo_Status=LDF_AtomPairInfo_Unset

      End
