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
* Copyright (C) 2006, Per-Olof Widmark                                 *
************************************************************************
      Subroutine Scram(CMO,nSym,nBas,nOrb,ScrFac)
************************************************************************
*                                                                      *
* This routine scrambles start orbitals in order to introduce symmetry *
* breaking where it is desirable.                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
* Written: September 2006                                              *
*                                                                      *
************************************************************************
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Real*8  CMO(*)
      Integer nSym
      Integer nBas(nSym)
      Integer nOrb(nSym)
      Real*8  ScrFac
*----------------------------------------------------------------------*
* External references                                                  *
*----------------------------------------------------------------------*
      Real*8 fRandom_Molcas
      External fRandom_Molcas
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer iSeed
      Save    iSeed
      Integer iSym
      Integer iOrb
      Integer jOrb
      Integer iBas
      Integer iOff
      Integer indx
      Integer jndx
      Real*8  p
      Real*8  q
      Real*8  u
      Real*8  v
      Data iSeed/13/
*----------------------------------------------------------------------*
* Do small rotations                                                   *
*----------------------------------------------------------------------*
      iOff=0
      Do iSym=1,nSym
*        Write(6,*) 'Scrambling irrep',iSym
         Do iOrb=1,nOrb(iSym)-1
            jOrb=iOrb+1
            q=ScrFac*(2.0d0*fRandom_Molcas(iSeed)-1.0d0)
            p=Sqrt(1.0d0-q*q)
*           Write(6,*) 'q=',q
            Do iBas=1,nBas(iSym)
               indx=iOff+(iOrb-1)*nBas(iSym)+iBas
               jndx=iOff+(jOrb-1)*nBas(iSym)+iBas
               u=p*CMO(indx)-q*CMO(jndx)
               v=q*CMO(indx)+p*CMO(jndx)
               CMO(indx)=u
               CMO(jndx)=v
            End Do
         End Do
         iOff=iOff+nBas(iSym)*nOrb(iSym)
      End Do
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
