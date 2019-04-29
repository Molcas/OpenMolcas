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
************************************************************************

      Subroutine Tra_Ctl2(CMO,PUVX,TUVX,D1I,FI,D1A,FA,IPR,lSquare,ExFac)

************************************************************************
*                                                                      *
*     main control section for:                                        *
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
************************************************************************
      Implicit Real*8 (A-H,O-Z)
      Parameter ( Zero=0.0d0 )

#include "rasdim.fh"
#include "general.fh"
#include "WrkSpc.fh"
*
      Dimension CMO(*), PUVX(*), TUVX(*)
      Dimension D1I(*), D1A(*), FI(*), FA(*)
      Integer state_symmetry, SymProd
      Integer off_PUVX(mxSym,mxSym,mxSym),
     &        off_sqMat(mxSym), off_ltMat(mxSym)
      Logical lSquare
*
      SymProd(i,j)=1+iEor(i-1,j-1)
*
      Call qEnter('Tra_Ctl2')
*
          IPR=0
      If ( IPR.gt.1 ) then
        Write(6,*)
        Write(6,*) ' Enter transformation section'
        Write(6,*) ' ============================'
        Write(6,*)
        Call GetMem('TRA_CTL','List','Real',junk,junk)
      End If
*
*     nasty, but necessary
      state_symmetry=lSym
*
*     generate offsets
*
      iStack = 0
      Do iSym = 1,nSym
         off_sqMat(iSym) = iStack
         iBas = nBas(iSym)
         iStack = iStack+ iBas*iBas
      End Do
*
      iStack = 0
      Do iSym = 1,nSym
         off_ltMat(iSym) = iStack
         iBas = nBas(iSym)
         iStack = iStack+ (iBas*iBas+iBas)/2
      End Do
*
      iStack = 0
      Do iSym = 1,nSym
        iOrb = nOrb(iSym)
        Do jSym = 1,nSym
          jAsh = nAsh(jSym)
          ijSym=SymProd(iSym,jSym)
          Do kSym = 1,nSym
            kAsh = nAsh(kSym)
            lSym=SymProd(ijSym,kSym)
            If (lSym.le.kSym) Then
              lAsh = nAsh(lSym)
              kl_Orb_pairs = kAsh*lAsh
              If ( kSym.eq.lSym ) kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
              off_PUVX(iSym,jSym,kSym)=iStack
              iStack = iStack + iOrb*jAsh*kl_Orb_pairs
            End If
          End Do
        End Do
      End Do
      nPUVX=iStack
*
*     Init Fock matrices
      Call dCopy_(nTot1,Zero,0,FI,1)
      Call dCopy_(nTot1,Zero,0,FA,1)
*
*     start transformation section
*
      If ( IPR.ge.5 ) then
        Write(6,*)
     &  ' Symmetry  Basis functions   total orbitals    active orbitals'
        Write(6,*)
     &  ' -------------------------------------------------------------'
      End If
      Do iSym = 1,nSym
        iBas = nBas(iSym)
        iOrb = nOrb(iSym)
        iFro = nFro(iSym)
        iIsh = nIsh(iSym)
        iAsh = nAsh(iSym)
        Do jSym = 1,iSym
          jBas = nBas(jSym)
          jOrb = nOrb(jSym)
          jFro = nFro(jSym)
          jIsh = nIsh(jSym)
          jAsh = nAsh(jSym)
          ijSym=SymProd(iSym,jSym)
          kSymMax = nSym
          If (.Not.lSquare) kSymMax=iSym
          Do kSym = 1,kSymMax
            kBas = nBas(kSym)
            kOrb = nOrb(kSym)
            kFro = nFro(kSym)
            kIsh = nIsh(kSym)
            kAsh = nAsh(kSym)
            lSym=SymProd(ijSym,kSym)
            If (lSym.le.kSym) Then
              lBas = nBas(lSym)
              lOrb = nOrb(lSym)
              lFro = nFro(lSym)
              lIsh = nIsh(lSym)
              lAsh = nAsh(lSym)
*
              If (iBas*jBas*kBas*lBas.ne.0 ) then
                ij_Bas_pairs = iBas*jBas
                ij_Orb_pairs = iAsh*jAsh
                If ( iSym.eq.jSym ) then
                  ij_Bas_pairs = (iBas*iBas+iBas)/2
                  ij_Orb_pairs = (iAsh*iAsh+iAsh)/2
                End If
                kl_Bas_pairs = kBas*lBas
                kl_Orb_pairs = kAsh*lAsh
                If ( kSym.eq.lSym ) then
                  kl_Bas_pairs = (kBas*kBas+kBas)/2
                  kl_Orb_pairs = (kAsh*kAsh+kAsh)/2
                End If
                Call TraDrv(IPR,lSquare,
     &                      iSym,jSym,kSym,lSym,
     &                      iBas,jBas,kBas,lBas,
     &                      iOrb,jOrb,kOrb,lOrb,
     &                      iFro,jFro,kFro,lFro,
     &                      iIsh,jIsh,kIsh,lIsh,
     &                      iAsh,jAsh,kAsh,lAsh,
     &                      ij_Bas_pairs,kl_Bas_pairs,
     &                      ij_Orb_pairs,kl_Orb_pairs,
     &                      off_PUVX,off_sqMat,off_ltMat,
     &                      mxSym,CMO,PUVX,
     &                      D1I,FI,D1A,FA,ExFac)
              End If
*
            End If
          End Do
        End Do
      End Do
      If ( IPR.ge.5 ) then
        Write(6,*)
     &  ' -------------------------------------------------------------'
      End If

C  Synchronize Fock matrices if running parallel:
      Call GADsum(FI,nTot1)
      Call GADsum(FA,nTot1)

*     print FI and FA
      If ( IPR.ge.10 ) then
        Write(6,*)
        Write(6,*) ' FI in AO-basis'
        Write(6,*) ' --------------'
        Write(6,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          if (iOrb.gt.0) Call TriPrt(' ',' ',FI(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
        Write(6,*)
        Write(6,*) ' FA in AO-basis'
        Write(6,*) ' --------------'
        Write(6,*)
        iOff = 1
        Do iSym = 1,nSym
          iOrb = nOrb(iSym)
          if (iOrb.gt.0) Call TriPrt(' ',' ',FA(iOff),iOrb)
          iOff = iOff + (iOrb*iOrb+iOrb)/2
        End Do
      End If

C  Synchronize PUVX if running parallel:
      Call GAdsum(PUVX,nPUVX)

*     select integrals TUVX
      Call Get_TUVX(PUVX,TUVX)

*     save integrals on disk
      iDisk=0
      Call DDaFile(LUINTM,1,PUVX,nPUVX,iDisk)

*     nasty, but necessary
      lSym=state_symmetry

      Call qExit('Tra_Ctl2')

      Return
      End
