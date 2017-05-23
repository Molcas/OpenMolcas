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
* Copyright (C) 2014, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_Fock_CoulombOnly_XIDI(Mode,tau,nD,FactC,
     *                                     ip_DBlocks,ip_FBlocks)
C
C     Thomas Bondo Pedersen, June 2014.
C
C     Compute correction from errors in the integral diagonal blocks,
C     adding them to the blocked Fock matrix.
C
      Implicit None
      Integer Mode
      Real*8  tau
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Logical  Rsv_Tsk
      External Rsv_Tsk
      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet
      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical IPI_set_here

      Integer TaskListID
      Integer AB
      Integer A, B
      Integer nAB
      Integer ip_Int, l_Int, n_Int
      Integer iD
      Integer ipD, ipF

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (.not.LDF_IntegralPrescreeningInfoIsSet()) Then
         Call LDF_SetIntegralPrescreeningInfo()
         IPI_set_here=.True.
      Else
         IPI_set_here=.False.
      End If

      Call Init_Tsk(TaskListID,NumberOfAtomPairs)
      Do While (Rsv_Tsk(TaskListID,AB))
         A=AP_Atoms(1,AB)
         B=AP_Atoms(2,AB)
         nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
         If (nAB.gt.0) Then
            n_Int=nAB**2
            l_Int=2*n_Int
            Call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
            Call LDF_ComputeValenceIntegrals(AB,AB,n_Int,Work(ip_Int))
            Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,AB,n_Int,
     &                                            Work(ip_Int+n_Int))
           Call dAXPY_(n_Int,-1.0d0,Work(ip_Int+n_Int),1,Work(ip_Int),1)
            Do iD=1,nD
               ipD=iWork(ip_DBlocks(iD)-1+AB)
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               Call dGeMV_('N',nAB,nAB,
     &                    FactC(iD),Work(ip_Int),max(nAB,1),
     &                    Work(ipD),1,1.0d0,Work(ipF),1)
            End Do
            Call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
         End If
      End Do
      Call Free_Tsk(TaskListID)

      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      End
