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
      Subroutine LDF_FVIFC(UsePartPermSym,Mode,tau,
     &                     nD,FactC,ip_DBlocks,ip_FBlocks)
C
C     Thomas Bondo Pedersen, January 2011.
C     - based on LDF_FCI by T.B. Pedersen.
C
C     Purpose: Compute Coulomb contribution to Fock matrix using
C              two-electron integrals computed from LDF coefficients
C              (debug code).
C
      Implicit None
      Logical UsePartPermSym
      Integer Mode
      real*8  tau
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Logical  LDF_IntegralPrescreeningInfoIsSet
      External LDF_IntegralPrescreeningInfoIsSet
      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Logical IPI_set_here

      Integer AB, CD
      Integer A, B, C, D
      Integer nAB, nCD
      Integer ip_Int, l_Int
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

      If (UsePartPermSym) Then ! use particle permutation symmetry
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
            Do CD=1,AB-1
               C=AP_Atoms(1,CD)
               D=AP_Atoms(2,CD)
               nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
               l_Int=nAB*nCD
               Call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
               Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,CD,
     &                                               l_Int,Work(ip_Int))
               Do iD=1,nD
                  ipD=iWork(ip_DBlocks(iD)-1+CD)
                  ipF=iWork(ip_FBlocks(iD)-1+AB)
                  Call dGeMV_('N',nAB,nCD,
     &                       FactC(iD),Work(ip_Int),max(nAB,1),
     &                       Work(ipD),1,1.0d0,Work(ipF),1)
               End Do
               Do iD=1,nD
                  ipD=iWork(ip_DBlocks(iD)-1+AB)
                  ipF=iWork(ip_FBlocks(iD)-1+CD)
                  Call dGeMV_('T',nAB,nCD,
     &                       FactC(iD),Work(ip_Int),max(nAB,1),
     &                       Work(ipD),1,1.0d0,Work(ipF),1)
               End Do
               Call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
            End Do
            l_Int=nAB**2
            Call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
            Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,AB,
     &                                            l_Int,Work(ip_Int))
            Do iD=1,nD
               ipD=iWork(ip_DBlocks(iD)-1+AB)
               ipF=iWork(ip_FBlocks(iD)-1+AB)
               Call dGeMV_('N',nAB,nAB,
     &                    FactC(iD),Work(ip_Int),max(nAB,1),
     &                    Work(ipD),1,1.0d0,Work(ipF),1)
            End Do
            Call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
         End Do
      Else ! do not use particle permutation symmetry
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            nAB=LDF_nBas_Atom(A)*LDF_nBas_Atom(B)
            Do CD=1,NumberOfAtomPairs
               C=AP_Atoms(1,CD)
               D=AP_Atoms(2,CD)
               nCD=LDF_nBas_Atom(C)*LDF_nBas_Atom(D)
               l_Int=nAB*nCD
               Call GetMem('FCIInt','Allo','Real',ip_Int,l_Int)
               Call LDF_ComputeValenceIntegralsFromC(Mode,tau,AB,CD,
     &                                               l_Int,Work(ip_Int))
               Do iD=1,nD
                  ipD=iWork(ip_DBlocks(iD)-1+CD)
                  ipF=iWork(ip_FBlocks(iD)-1+AB)
                  Call dGeMV_('N',nAB,nCD,
     &                       FactC(iD),Work(ip_Int),nAB,Work(ipD),1,
     &                       1.0d0,Work(ipF),1)
               End Do
               Call GetMem('FCIInt','Free','Real',ip_Int,l_Int)
            End Do
         End Do
      End If

      If (IPI_set_here) Then
         Call LDF_UnsetIntegralPrescreeningInfo()
      End If

      End
