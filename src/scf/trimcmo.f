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
* This routine trim CMO's from nBas x nBas to nBas x nOrb.             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      SubRoutine TrimCMO(CMO1,CMO2,nSym,nBas,nOrb)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8  CMO1(*)
      Real*8  CMO2(*)
      Integer nSym
      Integer nBas(*)
      Integer nOrb(*)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer iFrom(8)
      Integer iTo(8)
      Integer iSym
      Integer ndata, i
*----------------------------------------------------------------------*
* Transfer orbitals.                                                   *
*----------------------------------------------------------------------*
      iFrom(1) = 1
      iTo(1)   = 1
      Do iSym=1,nSym-1
         iFrom(iSym+1) = iFrom(iSym) + nBas(iSym)*nBas(iSym)
         iTo(iSym+1)   = iTo(iSym)   + nBas(iSym)*nOrb(iSym)
         If (iTo(iSym+1).gt.iFrom(iSym+1)) Then
            Write (6,*) 'Error in TrimCMO'
            Call Abend()
         End If
      End Do
      Do iSym=1,nSym
         ndata=nBas(iSym)*nOrb(iSym)
*
*        Note that CMO1 and CMO2 might overlap. Hence, we can not use
*        an ordinary call to DCopy!
*
         If (iFrom(iSym).ne.iTo(iSym)) Then
            Do i = 0, nData-1
               CMO2(iTo(iSym)+i) = CMO1(iFrom(iSym)+i)
            End Do
         End If
      End Do
*----------------------------------------------------------------------*
* Finish                                                               *
*----------------------------------------------------------------------*
      Return
      End
