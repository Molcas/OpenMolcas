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
* Copyright (C) 2010,2012, Thomas Bondo Pedersen                       *
************************************************************************
      Subroutine LDF_SetAtomPairInfo(useUnique,Verbose,irc)
C
C     Thomas Bondo Pedersen, June 2010.
C
C     Purpose: Set atom pair info in ldf_atom_pair_info.fh
C
C     Return codes:
C       irc=0: success
C       irc=1: already set
C       irc=2: an error occurred
C
C     -----
C     NOTES
C     -----
C
C     Thomas Bondo Pedersen, February 2012:
C     useUnique introduced:
C     Unique mapping should be disabled, useUnique=.False.; a warning
C     message is printed if not.
C     This implies that unique atom pairs are disabled.
C     The problem is that linear dependencies and 2C
C     functions are not handled correctly. In particular, linear
C     dependency arrays point to the parent atom and through that the
C     absolute shell - which is obviously wrong for a duplicate pair.
C     TODO/FIXME: fix index arrays for duplicate pairs and turn unique
C     atoms and atom pairs back on. This should be done at the end of
C     the LDF where the information is set for the non-unique pairs.
C
      Implicit None
      Logical useUnique
      Logical Verbose
      Integer irc
#include "ldf_atom_pair_info.fh"
#include "WrkSpc.fh"

      Character*19 SecNam
      Parameter (SecNam='LDF_SetAtomPairInfo')

      Integer iAtomPair
      Integer ip0

      Logical FirstCall
      Save FirstCall
      Data FirstCall /.True./

      ! Init return code
      irc=0

      ! If already set, print warning message and return
      If (FirstCall) Then
         FirstCall=.False.
      Else
         If (LDF_AtomPairInfo_Status.eq.LDF_AtomPairInfo_Set) Then
            If (Verbose) Then
               Call WarningMessage(1,
     &                        SecNam//'LDF Atom Pair Info already set!')
            End If
            irc=1
            Return
         End If
      End If

      ! Find significant atom pairs.
      Call LDF_FindSignificantAtomPairs(irc)
      If (irc.ne.0) Then
         If (Verbose) Then
            Write(6,'(A,A,I8)')
     &      SecNam,': LDF_FindSignificantAtomPairs returned code',irc
         End If
         LDF_AtomPairInfo_Status=LDF_AtomPairInfo_Unset
         irc=2
         Return
      End If

      ! Allocate and get unique significant atom pairs
      l_AP_Unique=NumberOfAtomPairs
      Call GetMem('AP_Unique','Allo','Inte',ip_AP_Unique,l_AP_Unique)
      If (useUnique) Then
         ! tbp: see note above
         Call WarningMessage(1,SecNam//
     &': WARNING: setting unique atom pair list; this may cause errors')
         Call xFlush(6)
         Call LDF_GetAtomPairToUniqueAtomPairMap(iWork(ip_AP_Unique),
     &                                           l_AP_Unique)
      Else
         ! tbp: make unique list a trivial mapping => no duplicate atom
         ! pairs.
         ip0=ip_AP_Unique-1
         Do iAtomPair=1,NumberOfAtomPairs
            iWork(ip0+iAtomPair)=iAtomPair
         End Do
      End If

      ! Allocate and initialize disk addresses for C
      l_AP_DiskC=NUmberOfAtomPairs
      Call GetMem('AP_DiskC','Allo','Inte',ip_AP_DiskC,l_AP_DiskC)
      Do iAtomPair=0,NumberOfAtomPairs-1
         iWork(ip_AP_DiskC+iAtomPair)=-1
      End Do

      ! Mark info set
      LDF_AtomPairInfo_Status=LDF_AtomPairInfo_Set

      ! Print
      If (Verbose) Then
         Call LDF_PrintAtomPairInfo()
      End If

      End
