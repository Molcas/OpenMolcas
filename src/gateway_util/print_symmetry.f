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
* Copyright (C) 2006, Roland Lindh                                     *
************************************************************************
      SubRoutine Print_Symmetry()
************************************************************************
*                                                                      *
*     Object: to write the output of seward            .               *
*                                                                      *
*                                                                      *
* Called from: Seward                                                  *
*                                                                      *
* Calling    : QEnter                                                  *
*              QExit                                                   *
*                                                                      *
*     Author: Roland Lindh, Dept Chem. Phys., Lund University, Sweden  *
*             September '06                                            *
************************************************************************
      use Symmetry_Info, only: nIrrep, iChTbl, iOper, lIrrep, lBsFnc,
     &                         SymLab
      Implicit Real*8 (A-H,O-Z)
#include "real.fh"
#include "print.fh"
#include "gateway.fh"
      Character SymOpr(0:7)*29, Format*80,  ChSymO(0:7)*5
      Data SymOpr/' Unit operation              ',
     &            ' Reflection in the yz-plane  ',
     &            ' Reflection in the xz-plane  ',
     &            ' Rotation around the z-axis  ',
     &            ' Reflection in the xy-plane  ',
     &            ' Rotation around the y-axis  ',
     &            ' Rotation around the x-axis  ',
     &            ' Inversion through the origin'/
      Data ChSymO/'  E  ','s(yz)','s(xz)','C2(z)',
     &            's(xy)','C2(y)','C2(x)','  i  '/
*                                                                      *
************************************************************************
*                                                                      *
      iRout=2
      iPrint=nPrint(iRout)
      If (iPrint.eq.0) Return
      Call QEnter('Print_Symmetry')
      LuWr=6
*                                                                      *
************************************************************************
*                                                                      *
      Write (LuWr,*)
      Call CollapseOutput(1,'   Symmetry information:')
      Write (LuWr,'(3X,A)') '   ---------------------'
      Write (LuWr,*)
*                                                                      *
************************************************************************
*                                                                      *
      If (nIrrep.ne.1) Then
         Write (LuWr,'(19X,A)') ' --- Group Generators ---'
         nOper=0
         If (nIrrep.eq.8) nOper=3
         If (nIrrep.eq.4) nOper=2
         If (nIrrep.eq.2) nOper=1
         Do i = 1, nOper
            j=i
            if (i.eq.3) j=4
            Write (LuWr,'(19X,A)') SymOpr(iOper(j))
         End Do
         Write (LuWr,*)
      ENd If
*                                                                      *
************************************************************************
*                                                                      *
      Write (LuWr,'(19X,A,A)') ' Character Table for ', SymLab
      Write (LuWr,*)
      Write (Format,'(A,I1,A)')   '(20X,A3,1X,',nIrrep,
     &                            '(1X,I5),2X,A)'
      Write (LuWr,'(27X,8(A5,1X))') (ChSymO(iOper(iIrrep)),iIrrep=
     &      0,nIrrep-1)
      Do iIrrep = 0, nIrrep-1
         LenlBs=Len(lBsFnc(iIrrep))
         Write (LuWr,Format)
     &         lIrrep(iIrrep),(iChTbl(iIrrep,jIrrep),
     &                         jIrrep=0,nIrrep-1),
     &         lBsFnc(iIrrep)(1:iCLast(lBsFnc(iIrrep),LenlBs))
      End Do
*                                                                      *
************************************************************************
*                                                                      *
      Call CollapseOutput(0,'  Symmetry information:')
      Write (LuWr,*)
      Call QExit('Print_Symmetry')
      Return
      End
