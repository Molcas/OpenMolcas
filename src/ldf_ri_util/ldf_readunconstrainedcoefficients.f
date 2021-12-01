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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_ReadUnconstrainedCoefficients(AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: read unconstrained coefficients (1C lin dep removed).
C              This routine is used only for debugging purposes, and
C              the file may not exist (in which case irc=-1 is
C              returned).
C
      Implicit None
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair

      Character*5 FileName
      Parameter (FileName='LDFUC')

      Logical FileExists
      Integer l
      Integer Lu_UC
      Integer iAddr
      Integer i
      Integer iOpt

      Integer j
      Integer Unique
      Logical isUnique
      Unique(j)=iWork(ip_AP_Unique-1+j)
      isUnique(j)=Unique(j).eq.j

      FileExists=.False.
      Call f_Inquire(FileName,FileExists)
      If (FileExists) Then
         l=LDF_AtomPair_DiagDim(AB)*LDF_nBasAux_Pair(AB)
         If (l_C.lt.l) Then
            irc=1
         Else
            Lu_UC=7
            Call DAName_mf_wa(Lu_UC,FileName)
            iAddr=0
            Do i=1,Unique(AB)-1
               If (isUnique(i)) Then
                  iAddr=iAddr
     &                 +LDF_AtomPair_DiagDim(i)*LDF_nBasAux_Pair(i)
               End If
            End Do
            iOpt=2
            Call dDAFile(Lu_UC,iOpt,C,l,iAddr)
            Call DAClos(Lu_UC)
            irc=0
         End If
      Else
         irc=-1
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_WriteUnconstrainedCoefficients(AB,l_C,C,irc)
C
C     Thomas Bondo Pedersen, March 2011.
C
C     Purpose: write unconstrained coefficients to disk (for debugging).
C
      Implicit None
      Integer AB
      Integer l_C
      Real*8  C(l_C)
      Integer irc
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim, LDF_nBasAux_Pair
      External LDF_AtomPair_DiagDim, LDF_nBasAux_Pair

      Character*5 FileName
      Parameter (FileName='LDFUC')

      Integer l
      Integer Lu_UC
      Integer iAddr
      Integer i
      Integer iOpt

      Integer j
      Logical isUnique
      isUnique(j)=iWork(ip_AP_Unique-1+j).eq.j

      If (.not.isUnique(AB)) Then
         irc=-1
      Else
         l=LDF_AtomPair_DiagDim(AB)*LDF_nBasAux_Pair(AB)
         If (l_C.lt.l) Then
            irc=1
         Else
            Lu_UC=7
            Call DAName_mf_wa(Lu_UC,FileName)
            iAddr=0
            Do i=1,AB-1
               If (isUnique(i)) Then
                  iAddr=iAddr
     &                 +LDF_AtomPair_DiagDim(i)*LDF_nBasAux_Pair(i)
               End If
            End Do
            iOpt=1
            Call dDAFile(Lu_UC,iOpt,C,l,iAddr)
            Call DAClos(Lu_UC)
            irc=0
         End If
      End If

      End
