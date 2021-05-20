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
      Subroutine LDF_CopyUniqueAtomPairs(irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Copy unique atom pair data
C
      Implicit None
      Integer irc
#include "ldf_atom_pair_info.fh"

      Integer iAtomPair

      ! Init return code
      irc=0

      ! Copy
      Do iAtomPair=1,NumberOfAtomPairs
         Call LDF_CopyUniqueAtomPair(iAtomPair)
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_CopyUniqueAtomPair(iAtomPair)
      Implicit None
      Integer iAtomPair
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim, LDF_DiskAddressOfC
      External LDF_AtomPair_DiagDim, LDF_DiskAddressOfC

      Character*8 Label

      Integer jAtomPair
      Integer ip, l

      Integer i,j
      Integer ip_D, AP_Unique, AP_1CLinDep, AP_2CFunctions
      ip_D(i)=iWork(ip_AP_Diag-1+i)
      AP_Unique(i)=iWork(ip_AP_Unique-1+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Get unique atom pair
      jAtomPair=AP_Unique(iAtomPair)

      ! Only copy if iAtomPair is not the unique one
      If (jAtomPair.ne.iAtomPair) Then
         ! Copy 1CLinDep
         iWork(ip_AP_1CLinDep+2*(iAtomPair-1))=AP_1CLinDep(1,jAtomPair)
         If (AP_1CLinDep(1,iAtomPair).gt.0) Then
            Write(Label,'(A,I5.5)') '1CL',iAtomPair-1
            l=3*AP_1CLinDep(1,iAtomPair)
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_AP_1CLinDep+2*(iAtomPair-1)+1)=ip
            Call iCopy(l,iWork(AP_1CLinDep(2,jAtomPair)),1,
     &                   iWork(AP_1CLinDep(2,iAtomPair)),1)
         End If
         ! Copy 2CFunctions
         iWork(ip_AP_2CFunctions+2*(iAtomPair-1))=AP_2CFunctions(1,
     &                                                        jAtomPair)
         If (AP_2CFunctions(1,iAtomPair).gt.0) Then
            Write(Label,'(A,I5.5)') '2CF',iAtomPair-1
            l=4*AP_2CFunctions(1,iAtomPair)
            Call GetMem(Label,'Allo','Inte',ip,l)
            iWork(ip_AP_2CFunctions+2*(iAtomPair-1)+1)=ip
            Call iCopy(l,iWork(AP_2CFunctions(2,jAtomPair)),1,
     &                   iWork(AP_2CFunctions(2,iAtomPair)),1)
         End If
         ! Copy updated diagonal block
         l=LDF_AtomPair_DiagDim(iAtomPair)
         Call dCopy_(l,Work(ip_D(jAtomPair)),1,Work(ip_D(iAtomPair)),1)
         ! Copy disk address of coefficients
         iWork(ip_AP_DiskC-1+iAtomPair)=LDF_DiskAddressOfC(jAtomPair)
      End If

      End
