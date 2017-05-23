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
      SubRoutine MulPop(CMO,mBB,nD,Ovrlp,mBT,OccNo,mmB)
************************************************************************
*                                                                      *
* This routine is a wrapper for Charge (which prints Mulliken popu-    *
* lation analyzes) since charge needs square CMO matrices.             *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund UNiversity, Sweden                                     *
*                                                                      *
************************************************************************
      Implicit None
#include "mxdm.fh"
#include "infscf.fh"
#include "stdalloc.fh"
      Integer mBB, nD, mBT, mmB
      Real*8 CMO(mBB,nD), Ovrlp(mBT), OccNo(mmB,nD)
************************************************************************
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer iSym
      Logical isOK
      Real*8, Dimension(:), Allocatable:: Aux1, Aux2
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      isOK=.true.
      Do iSym=1,nSym
         isOK=isOK .and. nBas(iSym).eq.nOrb(iSym)
      End Do
      If(isOK) Then
         If(iUHF.eq.0) Then
            Call Charge(nSym,nBas,Name,CMO(1,1),OccNo(1,1),
     &                  Ovrlp,2,.false.,.false.)
         Else
            Call Charge(nSym,nBas,Name,CMO(1,1),OccNo(1,1),
     &                  Ovrlp,0,.false.,.false.)
            Call Charge(nSym,nBas,Name,CMO(1,2),OccNo(1,2),
     &                  Ovrlp,1,.false.,.false.)
         End If
      Else
         Call mma_allocate(Aux1,nBB,Label='Aux1')
         Call mma_allocate(Aux2,nnB,Label='Aux2')
         If(iUHF.eq.0) Then
            Call PadCMO(CMO(1,1),Aux1,nSym,nBas,nOrb)
            Call PadEor(OccNo(1,1),Aux2,nSym,nBas,nOrb)
            Call Charge(nSym,nBas,Name,Aux1,Aux2,
     &                  Ovrlp,2,.false.,.false.)
         Else
            Call PadCMO(CMO(1,1),Aux1,nSym,nBas,nOrb)
            Call PadEor(OccNo(1,1),Aux2,nSym,nBas,nOrb)
            Call Charge(nSym,nBas,Name,Aux1,Aux2,
     &                  Ovrlp,0,.false.,.false.)
            Call PadCMO(CMO(1,2),Aux1,nSym,nBas,nOrb)
            Call PadEor(OccNo(1,2),Aux2,nSym,nBas,nOrb)
            Call Charge(nSym,nBas,Name,Aux1,Aux2,
     &                  Ovrlp,1,.false.,.false.)
         End If
         Call mma_deallocate(Aux1)
         Call mma_deallocate(Aux2)
      End If
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
