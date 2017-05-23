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
* This routine dumps orbital energies on the runfile, expands them if  *
* necessary.                                                           *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      SubRoutine DumpEor(Label,Eor,nSym,nBas,nOrb)
      Implicit None
#include "stdalloc.fh"
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Character*(*) Label
      Real*8  Eor(*)
      Integer nSym
      Integer nBas(*)
      Integer nOrb(*)
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer npDump
      Integer iFrom(8)
      Integer iTo(8)
      Integer iSym
      Integer ndata
      Real*8, Dimension(:), Allocatable:: Dump
*----------------------------------------------------------------------*
* Preliminaries                                                        *
*----------------------------------------------------------------------*
      npDump=0
      Do iSym=1,nSym
         npDump=npDump+nBas(iSym)
      End Do
      Call mma_allocate(Dump,npDump,Label='DumpOE')
*----------------------------------------------------------------------*
* Dump orbital energies                                                *
*----------------------------------------------------------------------*
      iFrom(1)=1
      iTo(1)=1
      Do iSym=1,nSym-1
         iFrom(iSym+1)=iFrom(iSym)+nOrb(iSym)
         iTo(iSym+1)=iTo(iSym)+nBas(iSym)
      End Do
      Do iSym=nSym,1,-1
         ndata=nOrb(iSym)
         call dcopy_(ndata,Eor(iFrom(iSym)),1,Dump(iTo(iSym)),1)
      End Do
      Call Put_dArray(Label,Dump,npDump)
*----------------------------------------------------------------------*
* Finish                                                               *
*----------------------------------------------------------------------*
      Call mma_deallocate(Dump)
      Return
      End
