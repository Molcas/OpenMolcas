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
      Parameter(Zero=0.0d0)
#include "WrkSpc.fh"
#include "timers.fh"
*
      Integer off_PUVX(mxSym,mxSym,mxSym), off_sqMat(*), off_ltMat(*)
      Dimension CMO(*), PUVX(*)
      Dimension D1I(*), D1A(*), FI(*), FA(*)
      Integer case1, case2
      Logical Process_Twice, lSquare
*
      iTri(i,j)=i*(i-1)/2+j
*
      Call qEnter('TraDrv')
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
      If ( case1.eq.4 .and. case2.eq.0 ) then
        Call qExit('TraDrv')
        Return
      End If
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
      If (nScrt1.ne.0)
     &    Call GetMem('TraScr1','Allo','Real',ipScrt1,nScrt1)
      Call GetMem('TraScr2','Allo','Real',lBuf2,nBuf2)
      Call GetMem('TraScr3','Allo','Real',lBuf3,nBuf3)
*
      Call GetMem('TraScr4','Allo','Real',ipPQVX,nPQVX)
      call dcopy_(nPQVX,[Zero],0,Work(ipPQVX),1)
      If (nTURS.gt.0) Then
         Call GetMem('TraScr6','Allo','Real',ipTURS,nTURS)
         call dcopy_(nTURS,[Zero],0,Work(ipTURS),1)
      End If
*
      Call GetMem('TraScr5','Max ','Real',ipInBuf,nInBuf)
      nInBuf=Min(nInBuf,(ij_Bas_pairs*kl_Bas_pairs+1))
      nInBuf=Max(nInBuf,(kl_Bas_pairs+1))
      Call GetMem('TraScr5','Allo','Real',ipInBuf,nInBuf)
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
     &                 Work(ipInBuf),nInBuf,nPairs)
            nOff = 0
          End If
*
*         unwrap triangular matrix of electron repulsion integrals
          ipPQRS_ = ipInBuf+nOff
          If ( kSym.eq.lSym ) then
            Call Square(Work(ipInBuf+nOff),Work(ipScrt1),1,kBas,kBas)
            ipPQRS = ipScrt1
          Else
            ipPQRS = ipInBuf+nOff
          End If
*
*         generate Fock matrices
          If ( case1.ne.4 ) then
            Call Timing(Piaget_1,Swatch,Swatch,Swatch)
            If ( case1.eq.2 ) then
              Call Ftwo(case1,ExFac,
     &                  iSym,kSym,
     &                  i,j,
     &                  off_sqMat,off_ltMat,
     &                  D1I,FI,D1A,FA,Work(ipPQRS))
            Else
              Call Ftwo(case1,ExFac,
     &                  iSym,jSym,
     &                  i,j,
     &                  off_sqMat,off_ltMat,
     &                  D1I,FI,D1A,FA,Work(ipPQRS))
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
     &                 kSym,lSym,
     &                 kBas,lBas,
     &                 kAsh,lAsh,
     &                 CMO(kkkOff),CMO(lllOff),
     &                 Work(ipPQRS),Work(lBuf2),
     &                 Work(lBuf3),Work(ipPQVX))
            If (Process_Twice) Then
                Call Tra2C(i,iSym,iBas,iAsh,
     &                     j,jSym,jBas,jAsh,
     &                     kl_Bas_pairs,ij_Orb_pairs,
     &                     CMO(iiiOff),CMO(jjjOff),
     &                     Work(ipPQRS_),Work(lBuf2),Work(ipTURS))
            End If
            Call Timing(Candino_2,Swatch,Swatch,Swatch)
            Candino_2 = Candino_2 - Candino_1
            Candino_3 = Candino_3 + Candino_2
          End If
*
          nPairs = nPairs-1
          nOff = nOff+kl_Bas_pairs
*
        End Do
      End Do
      Call Timing(Candino_1,Swatch,Swatch,Swatch)
      If (IPR.ge.99) Then
         Call RecPrt('PQVX',' ',Work(ipPQVX),ij_Bas_Pairs,
     &                kl_Orb_Pairs)
         If (Process_Twice) Call RecPrt('TURS',' ',Work(ipTURS),
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
           iOff=ipPQVX+(kl_pair-1)*ij_Bas_pairs
*
*          unwrap triangular matrix of electron repulsion integrals
           If ( iSym.eq.jSym ) then
             Call Square(Work(iOff),Work(ipScrt1),1,iBas,iBas)
             lBuf9 = ipScrt1
           Else
             lBuf9 = iOff
           End If
*
*          second half transformation of electron repulsion integrals:
*          (vx!ij) --> (pu!vx)
           Call Tra2B(iSym,jSym,
     &                iBas,jBas,
     &                iAsh,jAsh,
     &                iOrb,jOrb,
     &                kl_pair,kl_Orb_pairs,
     &                CMO(iiOff),CMO(jjOff),CMO(iiiOff),CMO(jjjOff),
     &                Work(lBuf9),Work(lBuf2),
     &                Work(lBuf3),Work(lBuf3),
     &                PUVX(1+i1),PUVX(1+i2) )
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
               iOff=ipTURS+(ij_pair-1)*kl_Bas_pairs
*
*              unwrap triangular matrix of electron repulsion integrals
              If ( kSym.eq.lSym ) then
                 Call Square(Work(iOff),Work(ipScrt1),1,kBas,kBas)
                 lBuf9 = ipScrt1
              Else
                 lBuf9 = iOff
              End If
*
*             second half transformation of electron repulsion integrals
*             (vx!ij) --> (pu!vx)
              Call Tra2B(kSym,lSym,
     &                   kBas,lBas,
     &                   kAsh,lAsh,
     &                   kOrb,lOrb,
     &                   ij_pair,ij_Orb_pairs,
     &                   CMO(kkOff),CMO(llOff),CMO(kkkOff),CMO(lllOff),
     &                   Work(lBuf9),Work(lBuf2),
     &                   Work(lBuf3),Work(lBuf3),
     &                   PUVX(1+i1),PUVX(1+i2) )
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
      Call GetMem('TraScr5','Free','Real',ipInBuf,nInBuf)
      If (nTURS.ne.0) Call GetMem('TraScr6','Free','Real',ipTURS,nTURS)
      Call GetMem('TraScr4','Free','Real',ipPQVX,nPQVX)
      Call GetMem('TraScr3','Free','Real',lBuf3,nBuf3)
      Call GetMem('TraScr2','Free','Real',lBuf2,nBuf2)
      If (nScrt1.ne.0)
     &   Call GetMem('TraScr1','Free','Real',ipScrt1,nScrt1)
      Call Timing(Candino_2,Swatch,Swatch,Swatch)
      Candino_2 = Candino_2 - Candino_1
      Candino_3 = Candino_3 + Candino_2
*
      Call qExit('TraDrv')
*
      Return
      End
