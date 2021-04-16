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
      Subroutine LDF_UpdateDiagonal(iAtomPair,l_CBar,CBar,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Update integral diagonal for atom pair iAtomPair,
C                 (uv|uv)=(uv|uv)-sum(K) CBar[uv,K]^2
C              Check for (too) negative diagonals.
C
C     Return codes:
C        irc=0: all ok
C        irc>1: too negative diagonals encountered (irc will be the
C               number of too negative diagonals)
C
C
      Implicit None
      Integer iAtomPair
      Integer l_CBar
      Real*8  CBar(l_CBar)
      Integer irc
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBasAux_Pair, LDF_AtomPair_DiagDim
      External LDF_nBasAux_Pair, LDF_AtomPair_DiagDim

      Integer M
      Integer nuv
      Integer ip_D
      Integer K
      Integer uvK
      Integer uv

      ! Get dimensions
      M=LDF_nBasAux_Pair(iAtomPair)
      nuv=LDF_AtomPair_DiagDim(iAtomPair)

      ! Update diagonal
      ip_D=iWork(ip_AP_Diag-1+iAtomPair)-1
      Do K=0,M-1
         uvK=nuv*K
         Do uv=1,nuv
            Work(ip_D+uv)=Work(ip_D+uv)-CBar(uvK+uv)**2
         End Do
      End Do

      ! Check for too negative diagonals
      irc=0
      Do uv=1,nuv
         If (Work(ip_D+uv).lt.TooNegative) Then
            irc=irc+1
         End If
      End Do

      End
