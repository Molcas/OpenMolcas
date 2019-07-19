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
* Copyright (C) Per-Olof Widmark                                       *
************************************************************************
************************************************************************
*                                                                      *
* This routine pads orbital energy vectors.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      SubRoutine PadEor(Eor1,Eor2,nSym,nBas,nOrb)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8  Eor1(*)
      Real*8  Eor2(*)
      Integer nSym
      Integer nBas(*)
      Integer nOrb(*)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer iFrom(8)
      Integer iTo(8)
      Integer iPtr
      Integer iSym
      Integer ndata
      Integer i
*----------------------------------------------------------------------*
* Transfer orbital energies.                                           *
*----------------------------------------------------------------------*
      iFrom(1) = nOrb(1)
      iTo(1)   = nOrb(1)
      Do iSym=1,nSym-1
         iFrom(iSym+1) = iFrom(iSym) + nOrb(iSym+1)
         iTo(iSym+1)   = iTo(iSym)   + nOrb(iSym+1)
     &                 + nBas(iSym)-nOrb(iSym)
      End Do
      Do iSym=nSym,1,-1
         ndata=nOrb(iSym)
         Do i=1,ndata
            Eor2(iTo(iSym)+1-i)=Eor1(iFrom(iSym)+1-i)
         End Do
         If(nBas(iSym).gt.nOrb(iSym)) Then
            ndata=nBas(iSym)-nOrb(iSym)
            iPtr=iTo(iSym)+1
            Call dCopy_(ndata,[0.0d0],0,Eor2(iPtr),1)
         End If
      End Do
*----------------------------------------------------------------------*
* Finish                                                               *
*----------------------------------------------------------------------*
      Return
      End
