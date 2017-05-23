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
      Subroutine LDF_Add2CenterFunctions(iAtomPair,ip_CBar,l_CBar,
     &                                   ip_Z,l_Z,Modified,irc)
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Purpose: Add two-center functions to the auxiliary basis of atom
C     pair iAtomPair (if needed), and update CBar and Z matrices
C     (reallocated and computed here). Note that the atoms constituting
C     atom pair iAtomPair need not be different. In this case, the
C     functions added are not two-center functions, of course, but taken
C     from the AO product space on the atom in question. For aCD/acCD
C     sets, there should be no need for this, but predefined sets may
C     need additional functions to achieve the required accuracy.
C
C     The addition of two-center functions is based on the updated
C     diagonal stored along with the atom pair info.
C
C     If the auxiliary basis was modified (i.e. if 2C functions were
C     added), Modified=.True. on exit.
C
      Implicit None
      Integer iAtomPair
      Integer ip_CBar, l_CBar
      Integer ip_Z, l_Z
      Logical Modified
      Integer irc
#include "localdf_print.fh"

      Integer n2CF_Added

      ! Init return code
      irc=0

      ! Set list of potential two-center functions from diagonals.
      ! Diagonals larger than the required threshold are included.
      n2CF_Added=0
      Call LDF_Set2CL(iAtomPair,n2CF_Added)
      If (iPrint.ge.Inf_AuxBas) Then
         Call Cho_Head('Auxiliary Basis Info after Initial 2C Addition',
     &                 '-',80,6)
         Call LDF_PrintAuxBasInfo(iAtomPair)
      End If

      ! Set Modified
      Modified=n2CF_Added.gt.0

      If (Modified) Then
         ! Cholesky decompose residual matrix to get minimal set of
         ! two-center functions
         Call LDF_ResidualCD(iAtomPair,ip_CBar,l_CBar,irc)
         If (irc.ne.0) Then
            Write(6,'(A,I8)')
     &      'LDF_Add2CenterFunctions: LDF_ResidualCD returned code',irc
            irc=1
            Return
         End If
         If (iPrint.ge.Inf_AuxBas) Then
            Call Cho_Head('Auxiliary Basis Info after Residual CD',
     &                    '-',80,6)
            Call LDF_PrintAuxBasInfo(iAtomPair)
         End If
         ! Deallocate CBar and Z
         Call GetMem('CBar','Free','Real',ip_CBar,l_CBar)
         ip_CBar=0
         l_CBar=0
         Call GetMem('ZVec','Free','Real',ip_Z,l_Z)
         ip_Z=0
         l_Z=0
         ! Compute new CBar and Z
         Call LDF_ComputeCBar(iAtomPair,ip_CBar,l_CBar,ip_Z,l_Z,irc)
         If (irc.ne.0) Then
            Write(6,'(A,I8)')
     &      'LDF_Add2CenterFunctions: LDF_ComputeCBar returned code',irc
            irc=1
            Return
         End If
      End If

      End
