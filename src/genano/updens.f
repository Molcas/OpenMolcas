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
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine UpDens
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      Do 100 i=1,nDsym
         pDsym(i)=0.0d0
100   Continue
      IndOcc=1
      IndCmo=1
      IndNme=1
      Do 200 iSym=1,nSym
*        Write(*,'(a,i5)') ' iSym:   ',iSym
         Do 210 iOrb=1,nBas(iSym)
*           Write(*,'(a,i5)') '   iOrb:   ',iOrb
*           Write(*,'(a,i5)') '   IndOcc: ',IndOcc
*           Write(*,'(a,i5)') '   IndCmo: ',IndCmo
            o=Occ(IndOcc)
            e=Eps(IndOcc)
            w=wc0-wc1*e
            w=Max(w,1.0d0)
*           Write(6,'(a,2i3,3f12.6)') 'iSym.iOrb,o,e,w',iSym,iOrb,o,e,w
            If(Abs(o).gt.thr) Then
               Call UpOrb(nBas(iSym),o,w,Cmo(IndCmo),
     &                    Name(IndNme))
            End If
            IndCmo=IndCmo+nBas(iSym)
            IndOcc=IndOcc+1
210      Continue
         IndNme=IndNme+nBas(iSym)
200   Continue
*     Write(*,*) '*** Partial density matrix in UpDens ***'
      iBlk=0
      Do 300 iLqn=0,MxLqn
         Do 310 iShell=-iLqn,iLqn
            iBlk=iBlk+1
            If(nPrim(iLqn).gt.0) Then
*              Write(*,'(a,2i5)') ' Block',iLqn,iShell
*              Call Triprt(' ','(6F12.6)',pDsym(iSymbk(iBlk)),
*    &                     nPrim(iLqn))
            End If
310      Continue
300   Continue
      Do 400 i=1,nDsym
         tDsym(i)=tDsym(i)+wSet(kSet)*pDsym(i)
400   Continue
*     Write(*,*) '*** Total density matrix in UpDens ***'
      iBlk=0
      Do 500 iLqn=0,MxLqn
         Do 510 iShell=-iLqn,iLqn
            iBlk=iBlk+1
            If(nPrim(iLqn).gt.0) Then
*              Write(*,'(a,2i5)') ' Block',iLqn,iShell
*              Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),
*    &                     nPrim(iLqn))
            End If
510      Continue
500   Continue
      Return
      End
