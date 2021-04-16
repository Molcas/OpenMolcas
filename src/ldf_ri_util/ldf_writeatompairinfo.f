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
      Subroutine LDF_WriteAtomPairInfo(irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Write atom pair info to disk for later use.
C
      Implicit None
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim
      External LDF_AtomPair_DiagDim

      Integer Lu
      Integer iAddr
      Integer iBuf(1)
      Integer ip, l
      Integer iAtomPair

      Integer i, j
      Integer ip_D, ip_DB
      Integer AP_1CLinDep, AP_2CFunctions
      ip_D(i)=iWork(ip_AP_Diag-1+i)
      ip_DB(i)=iWork(ip_AP_DiagBak-1+i)
      AP_1CLinDep(i,j)=iWork(ip_AP_1CLinDep-1+2*(j-1)+i)
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      ! Init return code
      irc=0

      ! Open file
      Lu=7
      Call DAName_MF(Lu,'LDFAP')

      ! Init disk address
      iAddr=0

      ! Write info
      l=1
      iBuf(1)=NumberOfAtomPairs
      Call iDAFile(Lu,1,iBuf,l,iAddr)
      l=2*NumberOfAtomPairs
      Call iDAFile(Lu,1,iWork(ip_AP_Atoms),l,iAddr)
      l=NumberOfAtomPairs
      Call iDAFile(Lu,1,iWork(ip_AP_Unique),l,iAddr)
      l=NumberOfAtomPairs
      Call iDAFile(Lu,1,iWork(ip_AP_DiskC),l,iAddr)
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         iBuf(1)=AP_1CLinDep(1,iAtomPair)
         Call iDAFile(Lu,1,iBuf,l,iAddr)
         l=3*AP_1CLinDep(1,iAtomPair)
         If (l.gt.0) Then
            ip=AP_1CLinDep(2,iAtomPair)
            Call iDAFile(Lu,1,iWork(ip),l,iAddr)
         End If
      End Do
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         iBuf(1)=AP_2CFunctions(1,iAtomPair)
         Call iDAFile(Lu,1,iBuf,l,iAddr)
         l=4*AP_2CFunctions(1,iAtomPair)
         If (l.gt.0) Then
            ip=AP_2CFunctions(2,iAtomPair)
            Call iDAFile(Lu,1,iWork(ip),l,iAddr)
         End If
      End Do
      Do iAtomPair=1,NumberOfAtomPairs
         l=1
         iBuf(1)=LDF_AtomPair_DiagDim(iAtompair)
         Call iDAFile(Lu,1,iBuf,l,iAddr)
         l=iBuf(1)
         If (l.gt.0) Then
            ip=ip_D(iAtomPair)
            Call dDAFile(Lu,1,Work(ip),l,iAddr)
            ip=ip_DB(iAtomPair)
            Call dDAFile(Lu,1,Work(ip),l,iAddr)
         End If
      End Do

      ! Close file
      Call DAClos(Lu)

      End
