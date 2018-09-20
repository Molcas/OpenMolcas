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
      Subroutine CpOvlp(from,to)
      Implicit Real*8 (a-h,o-z)
      Dimension from(*),to(*)
#include "parm.fh"
#include "common.fh"
      logical oki,okj
#include "symlab.fh"
      Do 100 iBlk=1,(MxLqn+1)*(MxLqn+1)
         indx=0
         ijBas=0
         iBasX=0
         Do 110 iSym=1,nSym
         Do 110 iBas=1,nBas(iSym)
         iBasX=iBasX+1
         oki=Name(iBasX)(1:LENIN).eq.Center
         oki=oki .and. Name(iBasX)(LENIN1:LENIN8).eq.type(iBlk)
         If(oki) indx=indx+1
         jndx=0
         jBasX=iBasX-iBas
         Do 110 jBas=1,iBas
            jBasX=jBasX+1
            ijBas=ijBas+1
            okj=Name(jBasX)(1:LENIN).eq.Center
            okj=okj .and. Name(jBasX)(LENIN1:LENIN8).eq.type(iBlk)
            If(okj) jndx=jndx+1
            If(oki .and. okj) Then
               ind=jndx+indx*(indx-1)/2+iSymBk(iBlk)-1
               to(ind)=from(ijBas)
            End If
110      Continue
100   Continue
*     Write(*,*) '*** Overlap matrix in CpOvlp ***'
      iBlk=0
      Do 200 iLqn=0,MxLqn
         Do 210 iShell=-iLqn,iLqn
            iBlk=iBlk+1
            If(nPrim(iLqn).gt.0) Then
*              Write(*,'(a,2i5)') ' Block',iLqn,iShell
*              Call Triprt(' ','(6F12.6)',to(iSymBk(iBlk)),nPrim(iLqn))
            End If
210      Continue
200   Continue
      Return
      End
