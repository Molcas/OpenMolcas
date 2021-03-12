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
* Copyright (C) 1996, Markus P. Fuelscher                              *
*               1998, Roland Lindh                                     *
************************************************************************
      Subroutine TraDrv(IPR,lSquare,
     &                  iSym,jSym,kSym,lSym,
     &                  iBas,jBas,kBas,lBas,
     &                  iOrb,jOrb,kOrb,lOrb,
     &                  iFro,jFro,kFro,lFro,
     &                  iIsh,jIsh,kIsh,lIsh,
     &                  iAsh,jAsh,kAsh,lAsh,
     &                  ij_Bas_pairs,kl_Bas_pairs,
     &                  ij_Orb_pairs,kl_Orb_pairs,
     &                  off_PUVX,off_sqMat,off_ltMat,mxSym,
     &                  CMO,PUVX,D1I,FI,D1A,FA,ExFac)
************************************************************************
*                                                                      *
*     control subsection for:                                          *
*     - transformation of ERIs from AO to MO basis                     *
*     - Fock matrix generation                                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     M.P. Fuelscher                                                   *
*     University of Lund, Sweden, 1996                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
*     Modified to only need unique symmetry blocks, R. Lindh, March '98*
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "stdalloc.fh"
#include "timers.fh"
*
      Integer off_PUVX(mxSym,mxSym,mxSym), off_sqMat(*), off_ltMat(*)
      Real*8 CMO(*), PUVX(*)
      Real*8 D1I(*), D1A(*), FI(*), FA(*)
      Integer case1, case2
      Logical Process_Twice, lSquare

      Real*8, Allocatable, Target:: Scrt1(:), PQVX(:), TURS(:),
     &                              InBuf(:)
      Real*8, Allocatable :: Buf2(:), Buf3(:)
      Real*8, Pointer :: Buf9(:)=>Null()
      Real*8, Pointer :: PQRS(:)=>Null(), PQRS_(:)=>Null()
*
      iTri(i,j)=i*(i-1)/2+j
*
*     generate offsets
      iiOff  = off_sqMat(iSym) + iFro*iBas + 1
      jjOff  = off_sqMat(jSym) + jFro*jBas + 1
      kkOff  = off_sqMat(kSym) + kFro*kBas + 1
      llOff  = off_sqMat(lSym) + lFro*lBas + 1
      iiiOff = iiOff + iIsh*iBas
      jjjOff = jjOff + jIsh*jBas
      kkkOff = kkOff + kIsh*kBas
      lllOff = llOff + lIsh*lBas
*
*     Generate logical flag for unique symmetry blocks
*
      Process_Twice = (iTri(iSym,jSym).ne.iTri(kSym,lSym)) .and.
     &                .Not.lSquare
*
*     find cases
      case1 = 4
      If ( iSym.eq.jSym ) case1 = case1-2
      If ( iSym.eq.kSym ) case1 = case1-1
      case2 = iAsh*jAsh*kAsh*lAsh
*
*     quit, if this integral block is not used
      If ( case1.eq.4 .and. case2.eq.0 ) Return
*
*     allocate memory
      nScrt1=0
      If (iSym.eq.jSym) nScrt1=iBas*jBas
      If (kSym.eq.lSym) nScrt1 = Max(nScrt1,kBas*lBas)
      nBuf2 = Max(jBas*iAsh,kBas*lAsh,iBas*jAsh)
      nBuf3 = Max(jOrb*iAsh,kAsh*lAsh,iOrb*jAsh)
*
*     PQVX and TURS are the halftransformed integrals
*
      nPQVX = ij_Bas_pairs*kl_Orb_pairs
      If (Process_Twice) Then
         nTURS = ij_Orb_pairs*kl_Bas_pairs
         nBuf2 = Max(nBuf2,ij_Orb_pairs,lBas*kAsh)
         nBuf3 = Max(nBuf3,kAsh*lOrb,kOrb*lAsh)
      Else
         nTURS = 0
      End If
*
      If (nScrt1.ne.0) Call mma_allocate(Scrt1,nScrt1,Label='Scrt1')
      Call mma_allocate(Buf2,nBuf2,Label='Buf2')
      Call mma_allocate(Buf3,nBuf3,Label='Buf3')
*
      Call mma_allocate(PQVX,nPQVX,Label='PQVX')
      PQVX(:)=Zero
      If (nTURS.gt.0) Then
         Call mma_allocate(TURS,nTURS,Label='TURS')
         TURS(:)=Zero
      End If
*
      Call mma_maxDBLE(nInBuf)
      nInBuf=Min(nInBuf,(ij_Bas_pairs*kl_Bas_pairs+1))
      nInBuf=Max(nInBuf,(kl_Bas_pairs+1))
      Call mma_allocate(InBuf,nInBuf,Label='InBuf')
*
      If ( IPR.ge.5 .and. case2.ne.0 ) then
         Write(6,'(1X,4I2,2X,4I4,2X,4I4,2X,4I4)')
     &                  iSym,jSym,kSym,lSym,
     &                  iBas,jBas,kBas,lBas,
     &                  iOrb,jOrb,kOrb,lOrb,
     &                  iAsh,jAsh,kAsh,lAsh
      End If
*
      nPairs  = 0
      ij_pair = 0
      nOff = 0
      Do i = 1,iBas
        jMax = jBas
        If ( iSym.eq.jSym ) jMax=i
        Do j = 1,jMax
          ij_pair = ij_pair+1
*
*         read a block of electron repulsion integrals
          If ( nPairs.eq.0 ) then
            iOpt = 1
            If ( ij_pair.ne.1 ) iOpt = 2
            Call RdOrd(iRc,iOpt,iSym,jSym,kSym,lSym,
     &                 InBuf,nInBuf,nPairs)
            nOff = 0
          End If
*
*         unwrap triangular matrix of electron repulsion integrals
          If ( kSym.eq.lSym ) then
            klBas=kBas*(kBas+1)/2
            Call Square(InBuf(nOff+1),Scrt1,1,kBas,kBas)
            PQRS(1:kBas**2) => Scrt1(1:kBas**2)
          Else
            klBas=kBas*lBas
            PQRS(1:klBas) => InBuf(nOff+1:nOff+klBas)
          End If
          PQRS_(1:klBas) => InBuf(nOff+1:nOff+klBas)
*
*         generate Fock matrices
          If ( case1.ne.4 ) then
            Call Timing(Piaget_1,Swatch,Swatch,Swatch)
            If ( case1.eq.2 ) then
              Call Ftwo(case1,ExFac,iSym,kSym,i,j,
     &                  off_sqMat,off_ltMat,
     &                  D1I,FI,D1A,FA,PQRS)
            Else
              Call Ftwo(case1,ExFac,iSym,jSym,i,j,
     &                  off_sqMat,off_ltMat,
     &                  D1I,FI,D1A,FA,PQRS)
            End If
            Call Timing(Piaget_2,Swatch,Swatch,Swatch)
            Piaget_2 = Piaget_2 - Piaget_1
            Piaget_3 = Piaget_3 + Piaget_2
          End If
*
*         first half transformation of electron repulsion integrals:
*         (ij!kl) --> (vx!ij)
          If ( case2.ne.0 ) then
            Call Timing(Candino_1,Swatch,Swatch,Swatch)
            Call Tra2A(ij_pair,ij_Bas_pairs,kl_Orb_pairs,
     &                 kSym,lSym,kBas,lBas,kAsh,lAsh,
     &                 CMO(kkkOff),CMO(lllOff),
     &                 PQRS,Buf2,Buf3,PQVX)
            If (Process_Twice) Then
                Call Tra2C(i,iSym,iBas,iAsh,j,jSym,jBas,jAsh,
     &                     kl_Bas_pairs,ij_Orb_pairs,
     &                     CMO(iiiOff),CMO(jjjOff),
     &                     PQRS_,Buf2,TURS)
            End If
            Call Timing(Candino_2,Swatch,Swatch,Swatch)
            Candino_2 = Candino_2 - Candino_1
            Candino_3 = Candino_3 + Candino_2
          End If
          PQRS =>Null()
          PQRS_=>Null()
*
          nPairs = nPairs-1
          nOff = nOff+kl_Bas_pairs
*
        End Do
      End Do
      Call Timing(Candino_1,Swatch,Swatch,Swatch)
      If (IPR.ge.99) Then
         Call RecPrt('PQVX',' ',PQVX,ij_Bas_Pairs,kl_Orb_Pairs)
         If (Process_Twice) Call RecPrt('TURS',' ',TURS,
     &                                  kl_Bas_Pairs,ij_Orb_Pairs)
      End If
*
*     Do second half transformation, skip if block is zero long
*
      If ( case2.ne.0 ) Then
         i1 = off_PUVX(iSym,jSym,kSym)
         i2 = off_PUVX(jSym,iSym,kSym)
*        Call FZero(PUVX(1+i1),iOrb*jAsh*kl_Orb_pairs)
*        Call FZero(PUVX(1+i2),jOrb*iAsh*kl_Orb_pairs)
         Do kl_pair = 1, kl_Orb_pairs
           iOff=(kl_pair-1)*ij_Bas_pairs
*
*          unwrap triangular matrix of electron repulsion integrals
           If ( iSym.eq.jSym ) then
             Call Square(PQVX(iOff+1),Scrt1,1,iBas,iBas)
             Buf9(1:iBas**2) => Scrt1(1:iBas**2)
           Else
             Buf9(1:iBas*jBas) => PQVX(iOff+1:iOff+iBas*jBas)
           End If
*
*          second half transformation of electron repulsion integrals:
*          (vx!ij) --> (pu!vx)
           Call Tra2B(iSym,jSym,iBas,jBas,iAsh,jAsh,iOrb,jOrb,
     &                kl_pair,kl_Orb_pairs,
     &                CMO(iiOff),CMO(jjOff),CMO(iiiOff),CMO(jjjOff),
     &                Buf9,Buf2,Buf3,Buf3,
     &                PUVX(1+i1),PUVX(1+i2) )
          Buf9=>Null()
*
         End Do
         If (IPR.ge.99) Then
             Call RecPrt('PUVX(i,j,k)',' ',PUVX(1+i1),iOrb,
     &                    jAsh*kl_Orb_pairs)
             Call RecPrt('PUVX(j,i,k)',' ',PUVX(1+i2),jOrb,
     &                    iAsh*kl_Orb_pairs)
         End If
*
         If (Process_Twice) Then
            i1 = off_PUVX(kSym,lSym,iSym)
            i2 = off_PUVX(lSym,kSym,iSym)
*           Call FZero(PUVX(1+i1),kOrb*iAsh*ij_Orb_pairs)
*           Call FZero(PUVX(1+i2),lOrb*kAsh*ij_Orb_pairs)
            Do ij_pair = 1, ij_Orb_pairs
               iOff=(ij_pair-1)*kl_Bas_pairs
*
*              unwrap triangular matrix of electron repulsion integrals
              If ( kSym.eq.lSym ) then
                 Call Square(TURS(iOff+1),Scrt1,1,kBas,kBas)
                 lBuf9 = ip_of_Work(Scrt1)
                 Buf9(1:kbas**2) => Scrt1(1:kbas**2)
              Else
                 Buf9(1:kBas*lBas) => TURS(iOff+1:iOff+kBas*lBas)
              End If
*
*             second half transformation of electron repulsion integrals
*             (vx!ij) --> (pu!vx)
              Call Tra2B(kSym,lSym,kBas,lBas,kAsh,lAsh,kOrb,lOrb,
     &                   ij_pair,ij_Orb_pairs,
     &                   CMO(kkOff),CMO(llOff),CMO(kkkOff),CMO(lllOff),
     &                   Buf9,Buf2,Buf3,Buf3,
     &                   PUVX(1+i1),PUVX(1+i2) )
              Buf9=>Null()
            End Do
            If (IPR.ge.99) Then
                Call RecPrt('PUVX(k,l,i)',' ',PUVX(1+i1),kOrb,
     &                       lAsh*ij_Orb_pairs)
                Call RecPrt('PUVX(l,k,i)',' ',PUVX(1+i2),lOrb,
     &                       kAsh*ij_Orb_pairs)
            End If
         End If
*
      End If
*
*     deallocate memory
      Call mma_deallocate(InBuf)
      If (nTURS.ne.0) Call mma_deallocate(TURS)
      Call mma_deallocate(PQVX)
      Call mma_deallocate(Buf3)
      Call mma_deallocate(Buf2)
      If (Allocated(Scrt1)) Call mma_deallocate(Scrt1)
      Call Timing(Candino_2,Swatch,Swatch,Swatch)
      Candino_2 = Candino_2 - Candino_1
      Candino_3 = Candino_3 + Candino_2
*
*
      Return
      End
