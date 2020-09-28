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
* Copyright (C) 1991, Manuela Merchan                                  *
************************************************************************
      Subroutine AutoCut
************************************************************************
*                                                                      *
*     Purpose:                                                         *
*     This routine counts the number of orbitals that have an          *
*     occupation number smaller than a given threshold and             *
*     marks them as deleted orbitals.                                  *
*                                                                      *
***** M. Merchan, University of Valencia, Spain, 1991 ******************
*
      Implicit Real*8 (A-H,O-Z)

#include "motra_global.fh"
#include "trafo_motra.fh"
*
*----------------------------------------------------------------------*
*     Start procedure                                                  *
*----------------------------------------------------------------------*
      ipBas=0
      Do iSym=1,nSym
        iDel=0
        Do iBas=1,nBas(iSym)
          If ( Occ(ipBas+iBas).le.abs(CutThrs(iSym)) ) iDel=iDel+1
        End Do
        ipBas=ipBas+nBas(iSym)
        If ( nDel(iSym).lt.iDel ) nDel(iSym)=iDel
        If ( (nDel(iSym)+nFro(iSym)).gt.nBas(iSym) ) Then
           Write (6,*) 'AutoCut:nDel(iSym)+nFro(iSym)).gt.nBas(iSym)'
           Write (6,*) 'iSym=',iSym
           Write (6,*) 'nDel(iSym)=',nDel(iSym)
           Write (6,*) 'nFro(iSym)=',nFro(iSym)
           Write (6,*) 'nBas(iSym)=',nBas(iSym)
           Call QTrace()
           Call Abend()
        End If
      End Do
*----------------------------------------------------------------------*
*     Normal termination                                               *
*----------------------------------------------------------------------*
      nOrbt=0
      nOrbtt=0
      Do iSym=1,nSym
        nOrb(iSym)=nBas(iSym)-nFro(iSym)-nDel(iSym)
        nOrbt=nOrbt+nOrb(iSym)
        nOrbtt=nOrbtt+nOrb(iSym)*(nOrb(iSym)+1)/2
      End Do
      Return
      End
