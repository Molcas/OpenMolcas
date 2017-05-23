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
#if defined (_MOLCAS_MPP_)
      Subroutine LDF_P_AddAuxBasVector(ip_V)
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Add aux bas vector across nodes.
C
      Implicit None
      Integer ip_V
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"
#include "para_info.fh"

      Integer  LDF_nAtom, LDF_nBasAux_Atom
      External LDF_nAtom, LDF_nBasAux_Atom

      Integer iAtom, nAtom, iAtomPair
      Integer ipV0, ip, l

      Integer i, j
      Integer AP_2CFunctions
      AP_2CFunctions(i,j)=iWork(ip_AP_2CFunctions-1+2*(j-1)+i)

      If (nProcs.gt.1 .and. Is_Real_Par()) Then
         nAtom=LDF_nAtom()
         ipV0=ip_V-1
         Do iAtom=1,nAtom
            l=LDF_nBasAux_Atom(iAtom)
            If (l.gt.0) Then
               ip=iWork(ipV0+iAtom)
               Call GAdGOp(Work(ip),l,'+')
            End If
         End Do
         ipV0=ipV0+nAtom
         Do iAtomPair=1,NumberOfAtomPairs
            l=AP_2CFunctions(1,iAtomPair)
            If (l.gt.0) Then
               ip=iWork(ipV0+iAtomPair)
               Call GAdGOp(Work(ip),l,'+')
            End If
         End Do
      End If
      End
#elif defined (NAGFOR)
c Some compilers do not like empty files
      Subroutine empty_LDF_P_AddAuxBasVector()
      End
#endif
