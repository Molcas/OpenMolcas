!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2005, Giovanni Ghigo                                   *
!***********************************************************************
!  Cho_GenC
!
!> @brief
!>   Routine for the generation of the ``A,B`` block of Coulomb integrals (in symmetries \p iSymA, \p iSymB)
!>   for occupied MO \p iI, \p iJ in symmetries \p iSymI, \p iSymJ.
!> @author Giovanni Ghigo
!>
!> @details
!> The routine generates the ``A,B`` block of integrals gathering 9
!> sub-blocks. These are combination of inactive, active, and
!> secondary ``A,B`` MO.
!>
!> @param[in] iSymI,iSymJ,iSymA,iSymB Symmetry block of the two-electrons integrals
!> @param[in] NumV                    Number of Cholesky vectors to transform in the current batch
!> @param[in] AddCou                  Array of the ``A,B`` integrals block
!> @param[in] LenCou                  Length of the ``A,B`` integrals block
!***********************************************************************
      Subroutine Cho_GenC(iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV,         &
     &                    AddCou,LenCou, LenEx)
!***********************************************************************
! Author :  Giovanni Ghigo                                             *
!           Lund University, Sweden & Torino University, Italy         *
!           July 2005                                                  *
!----------------------------------------------------------------------*
! This is the routine that really generates the A,B block of coulomb   *
! integrals (in symmetries iSymA,iSymB) for occupied MO iI,iJ in       *
! symmetries iSymI,iSymJ. The 3 x 3 A,B block is built gathering 9     *
! sub-blocks. These are combination of inactive, active, and secondary *
! A,B MO                                                               *
! OBS !!!!!  By now, it works only for iSymA .EQ. iSymB  !!!           *
!***********************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
      Integer LenCou, LenEx
      Real*8 AddCou(LenCou)
#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"

      Integer LenA(3), LenB(3)

      Type V1
        Real*8, Allocatable:: A(:)
      End Type V1
      Type (V1):: AddSB(3,3)

      Real*8, Allocatable:: AddSq(:)
!                                                                      *
!***********************************************************************
!                                                                      *
      Interface
        Subroutine MkCouSB11(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB11
!       Subroutine MkCouSB12(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!       Real*8, Allocatable:: A(:)
!       Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
!       End Subroutine MkCouSB12
!       Subroutine MkCouSB13(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!       Real*8, Allocatable:: A(:)
!       Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
!       End Subroutine MkCouSB13
        Subroutine MkCouSB21(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB21
        Subroutine MkCouSB22(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB22
!       Subroutine MkCouSB23(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!       Real*8, Allocatable:: A(:)
!       Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
!       End Subroutine MkCouSB23
        Subroutine MkCouSB31(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB31
        Subroutine MkCouSB32(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB32
        Subroutine MkCouSB33(A,iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        Real*8, Allocatable:: A(:)
        Integer iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV
        End Subroutine MkCouSB33
      End Interface
!                                                                      *
!***********************************************************************
!                                                                      *

!GG   ------------------------------------------------------------------
!      IfTest=.True.
!GG   ------------------------------------------------------------------

      LenA(1) = nIsh(iSymA)
      LenA(2) = nAsh(iSymA)
      LenA(3) = nSsh(iSymA)
      LenB(1) = nIsh(iSymB)
      LenB(2) = nAsh(iSymB)
      LenB(3) = nSsh(iSymB)

! --- GENERATION of SubBlocks
!GG   ------------------------------------------------------------------
      If (SubBlocks(1,1)) Call MkCouSB11(AddSB(1,1)%A,                  &
     &                  iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest.and.Allocated(AddSB(1,1)%A)) Then
      Write(6,*)'       SB_11 :',nIsh(iSymA),' x',nIsh(iSymB)
      Write(6,'(8F10.6)') AddSB(1,1)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(1,2)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
! excluded
!          Call MkCouSB12(AddSB(1,2)%A,
!     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB=nIsh(iSymA)*nAsh(iSymB)
          Call mma_Allocate(AddSB(1,2)%A,LenSB,Label='AddSB')
          AddSB(1,2)%A(:)=0.0D0
        Endif
      Endif
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(1,2)%A)) then
      Write(6,*)'       SB_12 :',nIsh(iSymA),' x',nAsh(iSymB)
      Write(6,'(8F10.6)')AddSB(1,2)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(1,3)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
! excluded
!          Call MkCouSB13(AddSB(1,3)%A,
!     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB=nIsh(iSymA)*nSsh(iSymB)
          Call mma_allocate(AddSB(1,3)%A,LenSB,Label='AddSB')
          AddSB(1,3)%A(:)=0.0D0
        Endif
      Endif
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(1,3)%A)) then
      Write(6,*)'       SB_13 :',nIsh(iSymA),' x',nSsh(iSymB)
      Write(6,'(8F10.6)')AddSB(1,3)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(2,1)) Call MkCouSB21(AddSB(2,1)%A,                  &
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(2,1)%A)) then
      Write(6,*)'       SB_21 :',nAsh(iSymA),' x',nIsh(iSymB)
      Write(6,'(8F10.6)')AddSB(2,1)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(2,2)) Call MkCouSB22(AddSB(2,2)%A,                  &
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(2,2)%A)) then
      Write(6,*)'       SB_22 :',nAsh(iSymA),' x',nAsh(iSymB)
      Write(6,'(8F10.6)')AddSB(2,2)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(2,3)) then
        If (iSymA.NE.iSymB) then
          Continue ! CGG
! excluded
!          Call MkCouSB23(AddSB(2,3)%A,
!     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
        else
          LenSB=nAsh(iSymA)*nSsh(iSymB)
          Call mma_allocate(AddSB(2,3)%A,LenSB,Label='AddSB')
          AddSB(2,3)%A(:)=0.0D0
        Endif
      Endif
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(2,3)%A)) then
      Write(6,*)'       SB_23 :',nAsh(iSymA),' x',nSsh(iSymB)
      Write(6,'(8F10.6)')AddSB(2,3)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(3,1)) Call MkCouSB31(AddSB(3,1)%A,                  &
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(3,1)%A)) then
      Write(6,*)'       SB_31 :',nSsh(iSymA),' x',nIsh(iSymB)
      Write(6,'(8F10.6)')AddSB(3,1)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If (SubBlocks(3,2)) Call MkCouSB32(AddSB(3,2)%A,                  &
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(3,2)%A)) then
      Write(6,*)'       SB_32 :',nSsh(iSymA),' x',nAsh(iSymB)
      Write(6,'(8F10.6)')AddSB(3,2)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If(SubBlocks(3,3)) Call MkCouSB33(AddSB(3,3)%A,                   &
     &     iSymI,iSymJ,iSymA,iSymB, iI,iJ, numV)
!GG   ------------------------------------------------------------------
      If(IfTest .and.Allocated(AddSB(3,3)%A)) then
      Write(6,*)'       SB_33 :',nSsh(iSymA),' x',nSsh(iSymB)
      Write(6,'(8F10.6)')AddSB(3,3)%A(:)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

!GG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)'     END GENERATION of SubBlocks'
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------
! --- END GENERATION of SubBlocks

! --- GATERING of SubBlocks
      Call mma_allocate(AddSq,LenEx,Label='AddSq')

      iAddCouSB = 1

      Do iSB_B = 1, 3
        Do iB = 1,LenB(iSB_B)
          Do iSB_A = 1, 3

            If (LenA(iSB_A)==0) Cycle

            ! SB(iSB_A,iSB_B)
            iAddSBi = 1 + LenA(iSB_A) * (iB-1)
            Call dCopy_(LenA(iSB_A),AddSB(iSB_B,iSB_A)%A(iAddSBi),1,    &
     &                  AddSq(iAddCouSB),1)
            iAddCouSB = iAddCouSB + LenA(iSB_A)

          EndDo ! iSB_B
        EndDo  ! iB
      EndDo ! iSB_B
      nOrbA = nOrb(iSymA)
!GG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Write(6,*)'        The Square Gatered matrix'
      Call PrintSquareMat(nOrbA,AddSq)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

      Call Local_Triang(nOrbA, AddSq)
      call daxpy_(LenCou,1.0d0,AddSq,1,AddCou,1)
      Call mma_deallocate(AddSq)
!GG   ------------------------------------------------------------------
      If(IfTest) then
      Write(6,*)
      Call XFlush(6)
      Write(6,*)'        The Triangular Integrals matrix'
      Call PrintDiagMat(nOrbA,AddCou)
      Call XFlush(6)
      EndIf
!GG   ------------------------------------------------------------------

! --- END GATERING of SubBlocks

      Do iSB_A = 1, 3
        Do iSB_B = 1, 3
          If (Allocated(AddSB(iSB_A,iSB_B)%A))                          &
     &        Call mma_deallocate(AddSB(iSB_A,iSB_B)%A)
        EndDo
      EndDo

!GG   ------------------------------------------------------------------
!      IfTest=.False.
!GG   ------------------------------------------------------------------
      Return
      End
