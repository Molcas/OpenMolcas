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
* Copyright (C) 1992, Per-Olof Widmark                                 *
*               1992, Markus P. Fuelscher                              *
*               1992, Piotr Borowski                                   *
************************************************************************
      SubRoutine SetUp
************************************************************************
*                                                                      *
*     purpose: Set up needed parameters                                *
*                                                                      *
*     called from: ReadIn, TrGen, Freeze, PrFin                        *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     written by:                                                      *
*     P.O. Widmark, M.P. Fuelscher and P. Borowski                     *
*     University of Lund, Sweden, 1992                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (a-h,o-z)
#include "mxdm.fh"
#include "infscf.fh"
      Integer maxnOcc(MxSym),minnOcc(MxSym)
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*
*---- Set up global parameters used later
      nnOc   = 0
      nnFr   = 0
      nnB    = 0
      nnO    = 0
      nBT    = 0
      nOT    = 0
      nBO    = 0
      nBB    = 0
      nOO    = 0
      nOV    = 0
      nOFS   = 0
      nOFT   = 0
      MaxBas = 0
      MaxOrb = 0
      MaxOrF = 0
      MaxOrO = 0
      MaxBxO = 0
cvv      MaxOcc = 0
      MaxFro = 0
      MaxBOF = 0
      MaxBOO = 0
      Do iSym = 1, nSym
         if (iUHF.eq.0) then
             maxnOcc(iSym)=nOcc(iSym,1)
             minnOcc(iSym)=nOcc(iSym,1)
         else
             maxnOcc(iSym)=max(nOcc(iSym,1),nOcc(iSym,2))
             minnOcc(iSym)=min(nOcc(iSym,1),nOcc(iSym,2))
         endif
         If (nBas(iSym).gt.MxBas) Then
            Write (6,*) 'SetUp: nBas(iSym).gt.MxBas'
            Write (6,*) 'nBas(iSym),MxBas=',nBas(iSym),MxBas
            Call Abend()
         End If
         If (nOrb(iSym).gt.nBas(iSym)) Then
            Write (6,*) 'SetUp: nOrb(iSym).gt.nBas(iSym)'
            Write (6,*) 'nOrb(iSym),nBas(iSym)=',nOrb(iSym),nBas(iSym)
            Call Abend()
         End If
         If (maxnOcc(iSym).gt.nOrb(iSym)) Then
            Write (6,*) 'iSym=',iSym
            Write (6,*) 'SetUp: nOcc(iSym).gt.nOrb(iSym)'
            Write (6,*) 'nOcc(iSym),nOrb(iSym)=',
     &                          maxnOcc(iSym),nOrb(iSym)
            Call Abend()
         End If
         If (nFro(iSym).gt.minnOcc(iSym)) Then
            Write (6,*) 'SetUp: nFro(iSym).gt.nOcc(iSym)'
            Write (6,*) 'nFro(iSym),nOcc(iSym)=',
     &                          nFro(iSym),minnOcc(iSym)
            Call Abend()
         End If
         nnOc   = nnOc + nOcc(iSym,1)
          if(iUHF.eq.1) then
           nnOc = nnOc + nOcc(iSym,2)
          endif
         nnFr   = nnFr + nFro(iSym)
         nnB    = nnB  + nBas(iSym)
         nnO    = nnO  + nOrb(iSym)
         nBT    = nBT  + nBas(iSym)*(nBas(iSym) + 1)/2
         nOT    = nOT  + nOrb(iSym)*(nOrb(iSym) + 1)/2
         nBO    = nBO  + nBas(iSym)*nOrb(iSym)
         nBB    = nBB  + nBas(iSym)*nBas(iSym)
         nOO    = nOO  + nOrb(iSym)*nOrb(iSym)
         nOV    = nOV  + (maxnOcc(iSym)-nFro(iSym))*
     &                   (nOrb(iSym)-minnOcc(iSym))
         nOFS   = nOFS + (nOrb(iSym)-nFro(iSym))**2
         nOFT   = nOFT + (nOrb(iSym)-nFro(iSym))
     &                 * (nOrb(iSym)-nFro(iSym)+1)/2
         MaxBas = Max(MaxBas,nBas(iSym))
         MaxOrb = Max(MaxOrb,nOrb(iSym))
         MaxOrF = Max(MaxOrF,nOrb(iSym) - nFro(iSym))
         MaxOrO = Max(MaxOrO,nOrb(iSym) - minnOcc(iSym))
         MaxBxO = Max(MaxBxO,nBas(iSym)*nOrb(iSym))
cvv         MaxOcc = Max(MaxOcc,nOcc(iSym,1))
         MaxFro = Max(MaxFro,nFro(iSym))
         MaxBOF = Max(MaxBOF,nBas(iSym)*(nOrb(iSym)-nFro(iSym)))
         MaxBOO = Max(MaxBOO,nBas(iSym)*(nOrb(iSym)-minnOcc(iSym)))
      End Do
*
      If (nnB.gt.2*MxBas .and. .not.DSCF ) Then
         Write (6,*) 'SetUp: nnB.gt.2*MxBas .and. .not.DSCF'
         Write (6,*) 'nnB,MxBas=',nnB,MxBas
         Call Abend()
      Else If (nnB.gt.4*MxBas .and. DSCF ) Then
         Write (6,*) 'SetUp: nnB.gt.4*MxBas .and. DSCF'
         Write (6,*) 'nnB,MxBas=',nnB,MxBas
         Call Abend()
      End If
*
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
*
      Return
      End
