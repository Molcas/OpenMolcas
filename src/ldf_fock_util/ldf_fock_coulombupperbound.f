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
      Subroutine LDF_Fock_CoulombUpperBoundNorm_Full(PrintNorm,
     &                                               PackedD,
     &                                               nD,FactC,ip_D,
     &                                               UBFNorm)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute norm of the upper bound to the Coulomb Fock
C              matrix error in LDF.
C
      Implicit None
      Logical PrintNorm
      Logical PackedD
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Real*8  UBFNorm(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer ip_DBlkP, l_DBlkP
      Integer iD

      If (nD.lt.1) Return
      If (NumberOfAtomPairs.lt.1) Return

      l_DBlkP=nD
      Call GetMem('CUBDBP','Allo','Inte',ip_DBlkP,l_DBlkP)
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
         Call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,
     &                         iWork(ip_DBlkP-1+iD))
         Call LDF_ScaleOffDiagonalMatrixBlocks(iWork(ip_DBlkP-1+iD),
     &                                         2.0d0)
      End Do
      Call LDF_Fock_CoulombUpperBoundNorm(PrintNorm,nD,FactC,
     &                                    iWork(ip_DBlkP),UBFNorm)
      Do iD=1,nD
         Call LDF_DeallocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
      End Do
      Call GetMem('CUBDBP','Free','Inte',ip_DBlkP,l_DBlkP)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombUpperBoundNorm(PrintNorm,nD,FactC,
     &                                          ip_DBlocks,UBFNorm)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute norm of the upper bound to the Coulomb Fock
C              matrix error in LDF.
C
      Implicit None
      Logical PrintNorm
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Real*8  UBFNorm(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer ip
      Integer ip_U, l_U
      Integer iD
      Integer AB
      Integer nAB
      Integer ipDel
      Integer uv

      Integer i, j
      Integer AP_Atoms
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      If (nD.lt.1) Return
      If (NumberOfAtomPairs.lt.1) Return

c-tbp Call LDF_GetQuadraticDiagonal(ip)
      ip=ip_AP_Diag
      l_U=nD
      Call GetMem('CUBNrmU','Allo','Real',ip_U,l_U)
      Call LDF_ComputeU(ip,nD,ip_DBlocks,Work(ip_U))
      Do iD=1,nD
         UBFNorm(iD)=0.0d0
         Do AB=1,NumberOfAtomPairs
            nAB=LDF_nBas_Atom(AP_Atoms(1,AB))
     &         *LDF_nBas_Atom(AP_Atoms(2,AB))
            ipDel=iWork(ip-1+AB)-1
            Do uv=1,nAB
               UBFNorm(iD)=UBFNorm(iD)+Work(ipDel+uv)
            End Do
         End Do
         UBFNorm(iD)=FactC(iD)*Work(ip_U-1+iD)*sqrt(UBFNorm(iD))
      End Do
      Call GetMem('CUBNrmU','Free','Real',ip_U,l_U)
C-tbp Call LDF_FreeQuadraticDiagonal(ip)

      If (PrintNorm) Then
         Do iD=1,nD
            Write(6,'(A,I10,A,1P,D20.10,1X,A,D20.10,A)')
     &      'Norm of upper bound Coulomb error for density',iD,':',
     &      UBFNorm(iD),
     &      '(BlockRMS=',
     &      sqrt(UBFNorm(iD)**2/dble(NumberOfAtomPairs)),')'
         End Do
         Call xFlush(6)
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombUpperBound_Full(PrintNorm,
     &                                           Add,PackedD,PackedF,
     &                                           nD,FactC,ip_D,ip_F)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute the upper bound correction to the LDF Fock
C              matrix (Coulomb only):
C
C              |Fexact(uv)-F(uv)| <= sqrt[(Delta(uv)|Delta(uv))]*U
C
C              U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|
C
      Implicit None
      Logical PrintNorm
      Logical Add
      Logical PackedD
      Logical PackedF
      Integer nD
      Real*8  FactC(nD)
      Integer ip_D(nD)
      Integer ip_F(nD)
#include "WrkSpc.fh"
#include "localdf_bas.fh"
#include "ldf_atom_pair_info.fh"

      Integer l
      Integer iD
      Integer ip_DBlkP, l_DBlkP
      Integer ip_FBlkP, l_FBlkP

      ! Return if nothing to do
      If (nD.lt.1) Return
      If (NumberOfAtomPairs.lt.1) Return

      ! Allocate and extract density matrix blocks
      l_DBlkP=nD
      Call GetMem('CUBFDBP','Allo','Inte',ip_DBlkP,l_DBlkP)
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
         Call LDF_Full2Blocked(Work(ip_D(iD)),PackedD,
     &                         iWork(ip_DBlkP-1+iD))
         Call LDF_ScaleOffdiagonalMatrixBlocks(iWork(ip_DBlkP-1+iD),
     &                                         2.0d0)
      End Do

      ! If not Add, initialize Fock matrices
      If (.not.Add) Then
         If (PackedF) Then
            l=nBas_Valence*(nBas_Valence+1)/2
         Else
            l=nBas_Valence**2
         End If
         Do iD=1,nD
            Call Cho_dZero(Work(ip_F(iD)),l)
         End Do
      End If

      ! Allocate and extract Fock matrix blocks
      l_FBlkP=nD
      Call GetMem('CUBFFBP','Allo','Inte',ip_FBlkP,l_FBlkP)
      Do iD=1,nD
         Call LDF_AllocateBlockMatrix('Fck',iWork(ip_FBlkP-1+iD))
         Call LDF_Full2Blocked(Work(ip_F(iD)),PackedF,
     &                         iWork(ip_FBlkP-1+iD))
      End Do

      ! Compute upper bound and add to blocked Fock matrices
      Call LDF_Fock_CoulombUpperBound(PrintNorm,nD,FactC,
     &                                iWork(ip_DBlkP),iWork(ip_FBlkP))

      ! Get full storage Fock matrices from blocked ones
      Do iD=1,nD
         Call LDF_Blocked2Full(iWork(ip_FBlkP-1+iD),PackedF,
     &                         Work(ip_F(iD)))
      End Do

      ! Deallocations
      Do iD=1,nD
         Call LDF_DeallocateBlockMatrix('Fck',iWork(ip_FBlkP-1+iD))
      End Do
      Call GetMem('CUBFFBP','Allo','Inte',ip_FBlkP,l_FBlkP)
      Do iD=1,nD
         Call LDF_DeallocateBlockMatrix('UBD',iWork(ip_DBlkP-1+iD))
      End Do
      Call GetMem('CUBFDBP','Free','Inte',ip_DBlkP,l_DBlkP)

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CoulombUpperBound(PrintNorm,nD,FactC,
     &                                      ip_DBlocks,ip_FBlocks)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: add the upper bound correction to the LDF Fock
C              matrix (Coulomb only):
C
C              |Fexact(uv)-F(uv)| <= sqrt[(Delta(uv)|Delta(uv))]*U
C
C              U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|
C
      Implicit None
      Logical PrintNorm
      Integer nD
      Real*8  FactC(nD)
      Integer ip_DBlocks(nD)
      Integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer ip_U, l_U, ip
      Integer ip_Norm, l_Norm
      Integer iD
      Integer AB

      Real*8 UBFNorm

      If (nD.lt.1) Return
      If (NumberOfAtomPairs.lt.1) Return

      l_U=nD
      Call GetMem('LDFCU','Allo','Real',ip_U,l_U)
C-tbp Call LDF_GetQuadraticDiagonal(ip)
      ip=ip_AP_Diag
      Call LDF_ComputeU(ip,nD,ip_DBlocks,Work(ip_U))
      Call LDF_Fock_CUB(ip,nD,FactC,Work(ip_U),ip_FBlocks)
C-tbp Call LDF_FreeQuadraticDiagonal(ip)
      Call GetMem('LDFCU','Free','Real',ip_U,l_U)

      If (PrintNorm) Then
         If (NumberOfAtomPairs.gt.0) Then
            l_Norm=NumberOfAtomPairs
            Call GetMem('UBFNrm','Allo','Real',ip_Norm,l_Norm)
            Do iD=1,nD
               Call LDF_BlockMatrixNorm(ip_FBlocks(iD),ip_Norm)
               UBFNorm=0.0d0
               Do AB=1,NumberOfAtomPairs
                  UBFNorm=UBFNorm+Work(ip_Norm-1+AB)**2
               End Do
               Write(6,'(A,A,I10,A,1P,D20.10,1X,A,D20.10,A)')
     &         'Norm of Fock matrix after adding Coulomb upper bound',
     &         ' for density',iD,':',
     &         sqrt(UBFNorm),
     &         '(BlockRMS=',sqrt(UBFNorm/dble(NumberOfAtomPairs)),')'
            End Do
            Call xFlush(6)
            Call GetMem('UBFNrm','Free','Real',ip_Norm,l_Norm)
         End If
      End If

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_Fock_CUB(ip_AP_QD,nD,FactC,U,ip_FBlocks)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute
C
C              sqrt[(Delta(uv)|Delta(uv))]*U
C
C              and add to Fock matrix.
C
C     Note: diagonal integrals for A=B must be stored quadratically.
C
      Implicit None
      Integer ip_AP_QD
      Integer nD
      Real*8  FactC(nD)
      Real*8  U(nD)
      Integer ip_FBlocks(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer iD
      Integer AB
      Integer A, B
      Integer nA, nB
      Integer uv
      Integer ipDel, ipFB

      Real*8  UU

      Integer i, j
      Integer ip_Delta
      Integer AP_Atoms
      ip_Delta(i)=iWork(ip_AP_QD-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      Do iD=1,nD
         UU=FactC(iD)*U(iD)
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            nA=LDF_nBas_Atom(A)
            nB=LDF_nBas_Atom(B)
            ipDel=ip_Delta(AB)-1
            ipFB=iWork(ip_FBlocks(iD)-1+AB)-1
            Do uv=1,nA*nB
               Work(ipFB+uv)=Work(ipFB+uv)+sqrt(Work(ipDel+uv))*UU
            End Do
         End Do
      End Do

      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Subroutine LDF_ComputeU(ip_AP_QD,nD,ip_DBlocks,U)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Purpose: compute
C
C              U = sum_uv sqrt[(Delta(uv)|Delta(uv))]*|D(uv)|
C
C     Note: diagonal integrals for A=B must be stored quadratically.
C
      Implicit None
      Integer ip_AP_QD
      Integer nD
      Integer ip_DBlocks(nD)
      Real*8  U(nD)
#include "WrkSpc.fh"
#include "ldf_atom_pair_info.fh"

      Integer  LDF_nBas_Atom
      External LDF_nBas_Atom

      Integer iD
      Integer AB
      Integer A, B
      Integer nA, nB
      Integer uv
      Integer ipDel, ipDB

      Integer i, j
      Integer ip_Delta
      Integer AP_Atoms
      ip_Delta(i)=iWork(ip_AP_QD-1+i)
      AP_Atoms(i,j)=iWork(ip_AP_Atoms-1+2*(j-1)+i)

      Do iD=1,nD
         U(iD)=0.0d0
         Do AB=1,NumberOfAtomPairs
            A=AP_Atoms(1,AB)
            B=AP_Atoms(2,AB)
            nA=LDF_nBas_Atom(A)
            nB=LDF_nBas_Atom(B)
            ipDel=ip_Delta(AB)-1
            ipDB=iWork(ip_DBlocks(iD)-1+AB)-1
            Do uv=1,nA*nB
               U(iD)=U(iD)+sqrt(Work(ipDel+uv))*abs(Work(ipDB+uv))
            End Do
         End Do
      End Do

      End
