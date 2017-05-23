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
      Subroutine SphAve
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
      iPtr=0
      Do 100 iLqn=0,MxLqn
         nDim=nPrim(iLqn)*(nPrim(iLqn)+1)/2
         Do 110 ij=1,nDim
            ind=iPtr+ij
            t=0.0d0
            Do 111 i=0,2*iLqn
               t=t+tDsym(ind+i*nDim)
111         Continue
            t=t/(2*iLqn+1)
            Do 112 i=0,2*iLqn
               tDsym(ind+i*nDim)=t
112         Continue
110      Continue
         iPtr=iPtr+(2*iLqn+1)*nDim
100   Continue
*      Write(*,*) '*** Density matrix in SphAve ***'
*      iBlk=0
*      Do 200 iLqn=0,MxLqn
*         Do 210 iShell=-iLqn,iLqn
*            iBlk=iBlk+1
*            If(nPrim(iLqn).gt.0) Then
*               Write(*,'(a,2i5)') ' Block',iLqn,iShell
*               Call Triprt(' ','(6F12.6)',tDsym(iSymbk(iBlk)),
*    &                      nPrim(iLqn))
*            End If
*210      Continue
*200   Continue
      Return
      End
