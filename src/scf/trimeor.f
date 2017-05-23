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
* This routine trim orbital energy vectors.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      SubRoutine TrimEor(Eor1,Eor2,nSym,nBas,nOrb)
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
      Integer iSym
      Integer ndata
*----------------------------------------------------------------------*
* Transfer orbital energies.                                           *
*----------------------------------------------------------------------*
      iFrom(1) = 1
      iTo(1)   = 1
      Do iSym=1,nSym-1
         iFrom(iSym+1) = iFrom(iSym) + nBas(iSym)
         iTo(iSym+1)   = iTo(iSym)   + nOrb(iSym)
      End Do
      Do iSym=nSym,1,-1
         ndata=nOrb(iSym)
         call dcopy_(ndata,Eor1(iFrom(iSym)),1,Eor2(iTo(iSym)),1)
      End Do
*----------------------------------------------------------------------*
* Finish                                                               *
*----------------------------------------------------------------------*
      Return
      End
