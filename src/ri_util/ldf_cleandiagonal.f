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
      Subroutine LDF_CleanDiagonal(iAtomPair)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: zero negative diagonals for atom pair iAtomPair.
C              If any diagonal is too negative, execution is stopped!
C
      Implicit None
      Integer iAtomPair
#include "localdf.fh"
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_AtomPair_DiagDim
      External LDF_AtomPair_DiagDim

      Integer ip_D, nuv, uv
#if defined (_DEBUG_)
      Integer nClean
#endif

#if defined (_DEBUG_)
      Write(6,'(A,A,I9)')
     & 'LDF_CleanDiagonal: zeroing negative diagonals for ',
     & 'atom pair',iAtomPair
      Call xFlush(6)
      nClean=0
#endif

      nuv=LDF_AtomPair_DiagDim(iAtomPair)
      ip_D=iWork(ip_AP_Diag-1+iAtomPair)-1
      Do uv=1,nuv
         If (Work(ip_D+uv).lt.0.0d0) Then
            If (Work(ip_D+uv).lt.TooNegative) Then
               Call WarningMessage(2,
     &                      'LDF_CleanDiagonal: too negative diagonal!')
               Write(6,'(A,I9)') 'Atom Pair:',iAtomPair
               Write(6,'(A,I9,1X,1P,D15.6)')
     &         'Diagonal element (no. and value):',uv,Work(ip_D+uv)
               Write(6,'(A,1P,D15.6,A)')
     &         '(Too negative diagonals are those <',TooNegative,')'
               Call LDF_Quit(1)
            End If
            Work(ip_D+uv)=0.0d0
#if defined (_DEBUG_)
            nClean=nClean+1
#endif
         End If
      End Do

#if defined (_DEBUG_)
      Write(6,'(A,I9,A)')
     & '   - zeroed ',nClean,' diagonals.'
      Call xFlush(6)
#endif

      End
