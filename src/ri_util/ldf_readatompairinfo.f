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
      Subroutine LDF_ReadAtomPairInfo(irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Read atom pair info from disk.
C
      Implicit None
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*8 Label

      Integer Lu
      Integer iAddr
      Integer iBuf(1)
      Integer ip, l
      Integer iAtomPair

      Integer i, j
      Integer AP_1CLinDep, AP_2CFunctions
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Init return code
      irc=0

      ! Open file
      Lu=7
      Call DAName_MF(Lu,'LDFAP')

      ! Init disk address
      iAddr=0

      ! Read info
      l=1
      Call iDAFile(Lu,2,iBuf,l,iAddr)
      NumberOfAtomPairs=iBuf(1)
      l_AP_Atoms=2*NumberOfAtomPairs
      Call GetMem('LDFAPA','Allo','Inte',ip_AP_Atoms,l_AP_Atoms)
      Call iDAFile(Lu,2,iWork(ip_AP_Atoms),l_AP_Atoms,iAddr)
      l_AP_Unique=NumberOfAtomPairs
      Call GetMem('AP_Unique','Allo','Inte',ip_AP_Unique,l_AP_Unique)
      Call iDAFile(Lu,2,iWork(ip_AP_Unique),l_AP_Unique,iAddr)
      l_AP_DiskC=NumberOfAtomPairs
      Call GetMem('AP_DiskC','Allo','Inte',ip_AP_DiskC,l_AP_DiskC)
      Call iDAFile(Lu,2,iWork(ip_AP_DiskC),l_AP_DiskC,iAddr)
      l_AP_1CLinDep=2*NumberOfAtomPairs
      Call GetMem('AP1CLD','Allo','Inte',ip_AP_1CLinDep,l_AP_1CLinDep)
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         Call iDAFile(Lu,2,iBuf,l,iAddr)
         iWork(ip_AP_1CLinDep+2*(iAtomPair-1))=iBuf(1)
         l=3*AP_1CLinDep(1,iAtomPair)
         If (l.lt.1) Then
            iWork(ip_AP_1CLinDep+2*(iAtomPair-1)+1)=0
         Else
            Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_AP_1CLinDep+2*(iAtomPair-1)+1)=ip
            Call iDAFile(Lu,2,iWork(ip),l,iAddr)
         End If
      End Do
      l_AP_2CFunctions=2*NumberOfAtomPairs
      Call GetMem('AP2CFN','Allo','Inte',ip_AP_2CFunctions,
     &                                    l_AP_2CFunctions)
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         Call iDAFile(Lu,2,iBuf,l,iAddr)
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=iBuf(1)
         l=4*AP_2CFunctions(1,iAtomPair)
         If (l.lt.1) Then
            iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=0
         Else
            Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip
            Call iDAFile(Lu,2,iWork(ip),l,iAddr)
         End If
      End Do
      Call LDF_AllocateBlockMatrix('APD',ip_AP_Diag)
      l_AP_Diag=NumberOfAtomPairs
      Call LDF_AllocateBlockMatrix('APB',ip_AP_DiagBak)
      l_AP_DiagBak=NumberOfAtomPairs
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         Call iDAFile(Lu,2,iBuf,l,iAddr)
         l=iBuf(1)
         If (l.lt.1) Then
            Call WarningMessage(1,
     &                 'LDF_ReadAtomPairInfo: zero diagonal dimension?')
            Write(6,'(A,I10)') 'Atom pair:',iAtomPair
            Call xFlush(6)
         Else
            ip=iWork(ip_AP_Diag-1+iAtomPair)
            Call dDAFile(Lu,2,Work(ip),l,iAddr)
            ip=iWork(ip_AP_DiagBak-1+iAtomPair)
            Call dDAFile(Lu,2,Work(ip),l,iAddr)
         End If
      End Do

      ! Close file
      Call DAClos(Lu)

      ! Mark atom pair info set
      LDF_AtomPairInfo_Status=LDF_AtomPairInfo_Set

      End
