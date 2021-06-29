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
! Copyright (C) 2004,2005, Giovanni Ghigo                              *
!               2021, Roland Lindh                                     *
!***********************************************************************
!  Cho_TraA
!
!> @brief
!>   Routine for the transformation of the Cholesky vectors in MO-based TCVx for case \c Sym(i) &ne; \c Sym(j)
!> @author Giovanni Ghigo
!>
!> @details
!> In the inner batch the Cholesky Full Vectors are transformed and
!> stored in memory. Adresses (``1``) and length (``2``) are stored in the
!> matrix <tt>TCVX(iType, Sym(i), Sym(j)) of allocatable 2D arrays</tt>.
!>
!> - \c iType = ``1``: TCVA
!> - \c iType = ``2``: TCVB
!> - \c iType = ``3``: TCVC
!> - \c iType = ``4``: TCVD
!> - \c iType = ``5``: TCVE
!> - \c iType = ``6``: TCVF
!> - \c iType = ``7``: TCVG
!>
!> Types ``1``, ``2`` and ``4``--``7`` are generated only if \c DoTCVA = ``.True.``
!> TCVC is always generated.
!>
!> In the first half-transformation the vectors are contracted
!> only with the occupied (inactive and active) MO coefficients
!> for \c Sym(j). In the second half-transformation the vectors are
!> contracted with all MO coefficients.
!>
!> @note
!> The logical matrix \c TCVXist must be defined.
!>
!> @param[in]     iSymL       Symmetry of the Cholesky vector
!> @param[in]     iSym        Symmetry(``i``) of the Cholesky Full Vector
!> @param[in]     jSym        Symmetry(``j``) of the Cholesky Full Vector
!> @param[in]     NumV        Number of Cholesky vectors to transform in the current batch
!> @param[in]     CMO         MO coefficients
!> @param[in]     NCMO        Total number of MO coefficients
!> @param[in]     lUCHFV      Unit number of the Cholesky Full Vector to transform (``CHFV``)
!> @param[in,out] iStrtVec_AB Current initial disk pointer of the Cholesky Full Vector to transform (``CHFV``)
!> @param[in]     nFVec       Number of Cholesky vectors to transform in the inner batch procedure
!***********************************************************************
      Subroutine Cho_TraA(iSymL, iSym,jSym, NumV, CMO,NCMO,             &
     &                               lUCHFV, iStrtVec_AB, nFVec)
!***********************************************************************
!  This is the routine for the transformation from AO basis to MO      *
!  basis of the Cholesky Full Vectors when  iSym.NE.jSym.              *
!  The new Transformed Cholesky Full Vectors L are :                   *
!                                                                      *
!  IF DoFull.EQ.False. (CASPT2, MCLR,MBPT2)                            *
!   TCVA : L_ij, L_ji=T(L_ij) if DoTCVA=.True.                         *
!   TCVB : L_tj, L_ui         if DoTCVA=.True.                         *
!   TCVC : L_aj, L_bi                                                  *
!   TCVD : L_tu, L_ut=T(L_tu) if DoTCVA=.True.                         *
!   TCVE : L_au, L_bt         if DoTCVA=.True.                         *
!   TCVF : L_ab               if DoTCVA=.True.                         *
!   TCVG : L_jt, L_iu         if DoTCVA=.True.                         *
!  For generation of <pk|ql>  p,q: All MO, k,l: Occupied (i & t)       *
!  MO Indices  i,j: Inactive;   t,u: Active;   a,b: Secondary          *
!                                                                      *
!----------------------------------------------------------------------*
!  Author  :  Giovanni Ghigo                                           *
!             Lund University, Sweden & University di Torino, Italia   *
!  Written :  October 2004                                             *
!  Modified:  July 2005                                                *
!***********************************************************************
      use Cho_Tra
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)

#include "rasdim.fh"
#include "stdalloc.fh"
#include "SysDef.fh"
      Dimension CMO(NCMO)
      Logical TCVA,TCVB,TCVBt,TCVC,TCVCt,TCVD,TCVE,TCVEt,TCVF

      Real*8, Allocatable:: XAj(:), XAu(:), XAb(:), XBi(:), XBt(:)
      Real*8, Allocatable:: FAB(:,:)

! --- Memory to allocate & Nr. of Cholesky vectors transformable
!     A=Alpha(AO);  B=Beta(AO)
      TCVA = .False.
      TCVB = .False.
      TCVBt= .False.
      TCVC = .False.
      TCVCt= .False.
      TCVD = .False.
      TCVE = .False.
      TCVEt= .False.
      TCVF = .False.

      NFAB = 0
      Nij = 0  ! A
      Nji = 0  ! A"
      Ntj = 0  ! B
      Njt = 0  ! G
      Nui = 0  ! B"
      Niu = 0  ! G"
      Naj = 0  ! C
      Nbi = 0  ! C"
      Ntu = 0  ! D
      Nut = 0  ! D"
      Nau = 0  ! E
      Nbt = 0  ! E"
      Nab = 0  ! F

      Len_XAj = 0   ! A/A", B, C
      Len_XAu = 0   ! D/D", E
      Len_XAb = 0   ! F
      Len_XBi = 0   ! B",C"
      Len_XBt = 0   ! E"

      NFAB = nBas(iSym) * nBas(jSym)


!     Allocate memory for Transformed Cholesky Vectors - TCVx
! TCV-A :
      If (TCVXist(1,iSym,jSym)) Then
        TCVA = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Nij = nIsh(iSym) * nIsh(jSym)
        Nji = nIsh(jSym) * nIsh(iSym)
        Call mma_allocate(TCVX(1,iSym,jSym)%A,Nij,NumV,Label='TCVA')
        Call mma_allocate(TCVX(1,jSym,iSym)%A,Nji,NumV,Label='TCVA')
      EndIf

! TCV-B :
      If (TCVXist(2,iSym,jSym)) Then
        TCVB = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Ntj = nAsh(iSym) * nIsh(jSym)
        Njt = nIsh(jSym) * nAsh(iSym)
        Call mma_allocate(TCVX(2,iSym,jSym)%A,Ntj,NumV,Label='TCVB')
        Call mma_allocate(TCVX(7,jSym,iSym)%A,Njt,NumV,Label='TCVB')
      EndIf
      If (TCVXist(2,jSym,iSym)) Then
        TCVBt= .True.
        Len_XBi = nBas(jSym) * nIsh(iSym)
        Nui = nAsh(jSym) * nIsh(iSym)
        Niu = nIsh(iSym) * nAsh(jSym)
        Call mma_allocate(TCVX(2,jSym,iSym)%A,Nui,NumV,Label='TCVBt')
        Call mma_allocate(TCVX(7,iSym,jSym)%A,Niu,NumV,Label='TCVBt')
      EndIf

! TCV-C :
      If (TCVXist(3,iSym,jSym)) Then
        TCVC = .True.
        Len_XAj = nBas(iSym) * nIsh(jSym)
        Naj = nSsh(iSym) * nIsh(jSym)
        Call mma_allocate(TCVX(3,iSym,jSym)%A,Naj,NumV,Label='TCVC')
      EndIf
      If (TCVXist(3,jSym,iSym)) Then
        TCVCt= .True.
        Len_XBi = nBas(jSym) * nIsh(iSym)
        Nbi = nSsh(jSym) * nIsh(iSym)
        Call mma_allocate(TCVX(3,jSym,iSym)%A,Nbi,NumV,Label='TCVCt')
      EndIf

! TCV-D :
      If (TCVXist(4,iSym,jSym)) Then
        TCVD = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Ntu = nAsh(iSym) * nAsh(jSym)
        Nut = nAsh(jSym) * nAsh(iSym)
        Call mma_allocate(TCVX(4,iSym,jSym)%A,Ntu,NumV,Label='TCVD')
        Call mma_allocate(TCVX(4,jSym,iSym)%A,Nut,NumV,Label='TCVD')
      EndIf

! TCV-E :
      If (TCVXist(5,iSym,jSym)) Then
        TCVE = .True.
        Len_XAu = nBas(iSym) * nAsh(jSym)
        Nau = nSsh(iSym) * nAsh(jSym)
        Call mma_allocate(TCVX(5,iSym,jSym)%A,Nau,NumV,Label='TCVE')
      EndIf
      If (TCVXist(5,jSym,iSym)) Then
        TCVEt= .True.
        Len_XBt = nBas(jSym) * nAsh(iSym)
        Nbt = nSsh(jSym) * nAsh(iSym)
        Call mma_allocate(TCVX(5,jSym,iSym)%A,Nbt,NumV,Label='TCVEt')
      EndIf

! TCV-F :
      If (TCVXist(6,iSym,jSym)) Then
        TCVF = .True.
        Len_XAb = nBas(iSym) * nSsh(jSym)
        Nab = nSsh(iSym) * nSsh(jSym)
        Call mma_allocate(TCVX(6,jSym,iSym)%A,Nab,NumV,Label='TCVF')
      EndIf

      iStrt = 1
      Do i=1,iSym-1
         iStrt = iStrt + nBas(i) * nBas(i)
      EndDo

      jStrt = 1
      Do j=1,jSym-1
        jStrt = jStrt + nBas(j) * nBas(j)
      EndDo

! --- START LOOP iiVec   -----------------------------------------------
      DO iiVec = 1, NumV, nFVec
        NumFV=Max(nFVec,NumV-iiVec+1)
        iFBatch = (iiVec+nFVec-1)/nFVec

!       Allocate memory & Load Full Cholesky Vectors - CHFV

        iStrtVec_FAB = iStrtVec_AB + nFVec * (iFBatch-1)

        Call mma_allocate(FAB,NFAB,NumFV,Label='FAB')
        Call RdChoVec(FAB,NFAB,NumFV,iStrtVec_FAB,lUCHFV)

!  ---  Start Loop  jVec  ---
        Do jVec=iiVec,iiVec+NumFV-1   ! Loop  jVec
          iVec = jVec - iiVec + 1

!     --- 1st Half-Transformation  iBeta(AO) -> q(MO) only occupied
          jStrt0MO = jStrt + nFro(jSym) * nBas(jSym)

!         From CHFV A(Alpha,Beta) to XAj(Alpha,jMO)
          If ( TCVA .or. TCVB .or. TCVC ) then
            Call mma_allocate(XAj,Len_XAj,Label='XAj')
            Call ProdsA_2(FAB(:,iVec), nBas(iSym),nBas(jSym),           &
     &                  CMO(jStrt0MO),nIsh(jSym), XAj)
          EndIf
          jStrt0MO = jStrt0MO + nIsh(jSym) * nBas(jSym)

!         From CHFV A(Alpha,Beta) to XAu(Alpha,uMO)
          If ( TCVD .or. TCVE ) then
            Call mma_allocate(XAu,Len_XAu,Label='XAu')
            Call ProdsA_2(FAB(:,iVec), nBas(iSym),nBas(jSym),           &
     &                  CMO(jStrt0MO),nAsh(jSym), XAu)
          EndIf
          jStrt0MO = jStrt0MO + nAsh(jSym) * nBas(jSym)

!         From CHFV A(Alpha,Beta) to XAb(Alpha,bMO)
          If ( TCVF ) then
            Call mma_allocate(XAb,Len_XAb,Label='XAb')
            Call ProdsA_2(FAB(:,iVec), nBas(iSym),nBas(jSym),           &
     &                  CMO(jStrt0MO),nSsh(jSym), XAb)
          EndIf

          iStrt0MO = iStrt + nFro(iSym) * nBas(iSym)

!         From CHFV A(Alpha,Beta) to XBi(Beta,iMO)
          If ( TCVBt .or. TCVCt ) then
            Call mma_allocate(XBi,Len_XBi,Label='XBi')
            Call ProdsA_2t(FAB(:,iVec), nBas(iSym),nBas(jSym),          &
     &                  CMO(iStrt0MO),nIsh(iSym), XBi)
          EndIf
          iStrt0MO = iStrt0MO + nIsh(iSym) * nBas(iSym)

!         From CHFV A(Alpha,Beta) to XBt(Beta,tMO)
          If ( TCVEt ) then
            Call mma_allocate(XBt,Len_XBt,Label='XBt')
            Call ProdsA_2t(FAB(:,iVec), nBas(iSym),nBas(jSym),          &
     &                  CMO(iStrt0MO),nAsh(iSym), XBt)
          EndIf

!     --- 2nd Half-Transformation  iAlpha(AO) -> p(MO)
          iStrt0MO = iStrt + nFro(iSym) * nBas(iSym)

!         From XAj(Alpha,jMO) to ij(i,j)
          If ( TCVA ) then
            Call ProdsA_1(XAj, nBas(iSym),nIsh(jSym),                   &
     &                    CMO(iStrt0MO),nIsh(iSym),                     &
     &                    TCVX(1,iSym,jSym)%A(:,jVec))
            Call Trnsps(nIsh(iSym),nIsh(jSym),                          &
     &                  TCVX(1,iSym,jSym)%A(:,jVec),                    &
     &                  TCVX(1,jSym,iSym)%A(:,jVec))
          EndIf
          iStrt0MO = iStrt0MO + nIsh(iSym) * nBas(iSym)

!         From XAj(Alpha,jMO) to tj(t,j)
          If ( TCVB ) then
            Call ProdsA_1(XAj, nBas(iSym),nIsh(jSym),                   &
     &                  CMO(iStrt0MO),nAsh(iSym),                       &
     &                  TCVX(2,iSym,jSym)%A(:,jVec))
            Call Trnsps(nAsh(iSym),nIsh(jSym),                          &
     &                  TCVX(2,iSym,jSym)%A(:,jVec),                    &
     &                  TCVX(7,jSym,iSym)%A(:,jVec))
          EndIf

!         From XAu(Alpha,jMO) to tu(t,u)
          If ( TCVD ) then
            Call ProdsA_1(XAu, nBas(iSym),nAsh(jSym),                   &
     &                   CMO(iStrt0MO),nAsh(iSym),                      &
     &                   TCVX(4,iSym,jSym)%A(:,jVec))
            Call Trnsps(nAsh(iSym),nAsh(jSym),                          &
     &                  TCVX(4,iSym,jSym)%A(:,jVec),                    &
     &                  TCVX(4,jSym,iSym)%A(:,jVec))
          EndIf
          iStrt0MO = iStrt0MO + nAsh(iSym) * nBas(iSym)

!         From XAj(Alpha,jMO) to aj(a,j)
          If ( TCVC ) then
            Call ProdsA_1(XAj, nBas(iSym),nIsh(jSym),                   &
     &                   CMO(iStrt0MO),nSsh(iSym),                      &
     &                   TCVX(3,iSym,jSym)%A(:,jVec))
          EndIf

!         From XAu(Alpha,jMO) to au(a,u)
          If ( TCVE ) then
            Call ProdsA_1(XAu, nBas(iSym),nAsh(jSym),                   &
     &                    CMO(iStrt0MO),nSsh(iSym),                     &
     &                    TCVX(5,iSym,jSym)%A(:,jVec))
          EndIf

!         From XAb(Alpha,jMO) to ab(a,b)
          If ( TCVF ) then
            Call ProdsA_1(XAb, nBas(iSym),nSsh(jSym),                   &
     &                   CMO(iStrt0MO),nSsh(iSym),                      &
     &                   TCVX(6,jSym,iSym)%A(:,jVec))
          EndIf

          jStrt0MO = jStrt + ( nFro(jSym) + nIsh(jSym) ) * nBas(jSym)

!         From XBi(Beta,jMO) to ui(u,i)
          If ( TCVBt ) then
            Call ProdsA_1(XBi, nBas(jSym),nIsh(iSym),                   &
     &                    CMO(jStrt0MO),nAsh(jSym),                     &
     &                    TCVX(2,jSym,iSym)%A(:,jVec))
            Call Trnsps(nAsh(jSym),nIsh(iSym),                          &
     &                  TCVX(2,jSym,iSym)%A(:,jVec),                    &
     &                  TCVX(7,iSym,jSym)%A(:,jVec))
          EndIf
          jStrt0MO = jStrt0MO + nAsh(jSym) * nBas(jSym)

!         From XBi(Beta,jMO) to bi(b,i)
          If ( TCVCt ) then
            Call ProdsA_1(XBi, nBas(jSym),nIsh(iSym),                   &
     &                   CMO(jStrt0MO),nSsh(jSym),                      &
     &                   TCVX(3,jSym,iSym)%A(:,jVec))
          EndIf

!         From XBt(Beta,jMO) to bt(b,t)
          If ( TCVEt ) then
            Call ProdsA_1(XBt, nBas(jSym),nAsh(iSym),                   &
     &                   CMO(jStrt0MO),nSsh(jSym),                      &
     &                   TCVX(5,jSym,iSym)%A(:,jVec))
          EndIf

!     --- End of Transformations

          If (Allocated(XAj)) Call mma_deallocate(XAj)
          If (Allocated(XAu)) Call mma_deallocate(XAu)
          If (Allocated(XAb)) Call mma_deallocate(XAb)
          If (Allocated(XBi)) Call mma_deallocate(XBi)
          If (Allocated(XBt)) Call mma_deallocate(XBt)

        EndDo
!  ---  End Loop  jVec  ---

        Call mma_deallocate(FAB)
      ENDDO
! --- END LOOP iiVec   -------------------------------------------------

      Return
! Avoid unused argument warnings
      If (.False.) Call Unused_integer(iSymL)
      End

      Subroutine ProdsA_1(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8 AB(iA,iB), CMO(iA,nMO), Y(nMO,iB)
      Call DGEMM_('T','N',nMO,iB,iA,1.0d0,CMO,iA,AB,iA,0.0d0,Y,nMO)
      Return
      End

      Subroutine ProdsA_2(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8 AB(iA,iB), CMO(iB,nMO), Y(iA,nMO)
      Call DGEMM_('N','N',iA,nMO,iB,1.0d0,AB,iA,CMO,iB,0.0d0,Y,iA)
      Return
      End

      Subroutine ProdsA_2t(AB,iA,iB, CMO,nMO, Y)
      Implicit Real*8 (a-h,o-z)
      Implicit Integer (i-n)
      Real*8 AB(iA,iB), CMO(iA,nMO), Y(iB,nMO)
      Call DGEMM_('T','N',iB,nMO,iA,1.0d0,AB,iA,CMO,iA,0.0d0,Y,iB)
      Return
      End
