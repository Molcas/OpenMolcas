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
      Subroutine LDF_WriteC(iAtomPair,l_C,C,LuC,iAddr)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Write LDF fitting coefficients for atom pair iAtomPair to
C     disk at address iAddr. Save disk address in atom pair info.
C
C     Note: iAddr is updated here (to disk address of next atom pair)!
C
      Implicit None
      Integer iAtomPair
      Integer l_C
      Real*8  C(l_C)
      Integer LuC
      Integer iAddr
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Character*10 SecNam
      Parameter (SecNam='LDF_WriteC')

      Integer  LDF_nBas_Atom, LDF_nBasAux_Pair
      External LDF_nBas_Atom, LDF_nBasAux_Pair

      Integer LenOfC, iOpt

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      ! Save disk address
      If (iAddr.lt.0) Then
         Call WarningMessage(2,SecNam//': Invalid disk address')
         Call LDF_Quit(1)
      Else
         iWork(ip_AP_DiskC-1+iAtomPair)=iAddr
      End If

      ! Compute dimension of C
      LenOfC=LDF_nBas_Atom(AP_Atoms(1,iAtomPair))
     &      *LDF_nBas_Atom(AP_Atoms(2,iAtomPair))
     &      *LDF_nBasAux_Pair(iAtomPair)

      ! Write C to disk
      If (LenOfC.gt.l_C) Then
         Call WarningMessage(2,SecNam//': LenOfC>l_C')
         Call LDF_Quit(1)
      Else
         iOpt=1 ! write option
         Call dDAFile(LuC,iOpt,C,LenOfC,iAddr)
      End If

      End
