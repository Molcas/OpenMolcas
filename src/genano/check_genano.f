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
* This routine perform some consistency checks.                        *
*                                                                      *
*======================================================================*
*                                                                      *
* Author: Per-Olof Widmark                                             *
*         IBM Sweden                                                   *
*                                                                      *
************************************************************************
      Subroutine Check_genano
      Implicit Real*8 (a-h,o-z)
#include "parm.fh"
#include "common.fh"
#include "symlab.fh"
      Dimension icnt(0:MxLqn)
      Do 10 i=0,MxLqn
         icnt(i)=0
10    Continue
      ind=0
      Do 100 iSym=1,nSym
*        Write(*,'(a,i1,a)') ' <<< Symmetry ',iSym,' >>>'
         Do 110 iBas=1,nBas(iSym)
            ind=ind+1
*           Write(*,*) Name(1,ind),Name(2,ind)
            Do 111 l=0,MxLqn
               If(Name(ind)(1:LENIN).eq.Center) Then
                  If(Name(ind)(LENIN1:LENIN4).eq.type(l*(l+1)+1)) Then
                     icnt(l)=icnt(l)+1
                  End If
               End If
111         Continue
110      Continue
100   Continue
      Do 200 l=0,MxLqn
         If(icnt(l).ne.nPrim(l)) Then
            Write(6,*) 'Number of primitives do not match!'
            Write(6,'(1x,a,2i5)') type(l*(l+1)+1),nPrim(l),icnt(l)
            Call Quit_OnUserError()
         End If
200   Continue
      Do 300 i=1,nDsym
         pDsym(i)=0.0d0
300   Continue
      Return
      End
