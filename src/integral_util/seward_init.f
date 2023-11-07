***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 1990,2020, Roland Lindh                                *
!               1990, IBM                                              *
!***********************************************************************
      Subroutine Seward_Init()
!***********************************************************************
!                                                                      *
!     Object: to set data which is stored in common blocks             *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose      *
!             January 1990                                             *
!***********************************************************************
      use EFP_Module, only: lEFP, nEFP_fragments
      use k2_arrays, only: XMem
      use k2_structure, only: k2_processed
      use Basis_Info, only: Seward_Activated
      use RICD_Info, only: iRI_Type, Do_RI
      use rmat, only: RmatR, Epsabs, Epsrel, qCoul, Epsq, bParm, dipol,
     &                Dipol1, keyr, Quadpack, nagint, testint,
     &                RMat_On, lgamma
      use DCR_mod, only: DCR_Init
      use NAC, only: isNAC, isCSF
      Implicit None
      Logical, External:: Reduce_Prt
      Integer, External:: iPrintLevel
      Integer iPL
#include "twoswi.fh"
#include "print.fh"
      Character(LEN=180) Env
!                                                                      *
!***********************************************************************
!                                                                      *
!
!-----Info
!
      Seward_Activated=.False.
!
!-----LInfo
!
      Call GetEnvF('MOLCAS_NEW_DEFAULTS', Env)
      Call UpCase(Env)
      If (Env.eq.'YES') Then
         Do_RI=.True.
         iRI_Type=4
      End If
!
      iPL=iPrintLevel(-1)
      If (iPL.eq.2) Then
         iPL=5
      Else If (iPL.eq.3) Then
         iPL=6
      Else If (iPL.eq.4) Then
         iPL=7
      Else If (iPL.eq.5) Then
         iPL=49  ! 99 would be just too much
      End If
      nPrint(:)=iPL
      If ((Reduce_Prt().and.iPL.lt.6).or.iPL.eq.0) Then
         Show=.False.
      Else
         Show=.True.
      End If
!
      NDDO=.False.
!
      Seward_Activated=.True.
      XMem=.False.
      k2_processed=.False.
!
      Call Set_Binom()
      Call Set_CanInd()
!
!     Set some default value for RMAT type integration
!
!---- rmat.fh
!
      RmatR     = 10.0D0
      Epsabs    = 10.D-10
      Epsrel    = 1.D-14
      qCoul     = 0.0D0
      Epsq      = 1.D-8
      bParm     = 0.0D0
      dipol(1)  = 0.0D0
      dipol(2)  = 0.0D0
      dipol(3)  = 0.0D0
      Dipol1    = 0.0D0
      keyr      = 6
      Quadpack  = .True.
      nagint    = .False.
      testint   = .False.
      RMat_On   = .False.
      lgamma = 9
!
      Call DCR_Init()
!
      Call Set_Basis_Mode('Valence')
!
!     nac.fh
!
      isNAC = .False.
      isCSF = .False.
!
!     EFP stuff
!
      lEFP=.False.
      nEFP_fragments=0
!
      End Subroutine Seward_Init
