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
* Copyright (C) 1996, Anders Bernhardsson                              *
************************************************************************
       SubRoutine FCOOR(LUT,COOR)
********************************************************************
*                                                                  *
*      Transforms a symmetry adapted gradient to unsymmetric  form *
*                                                                  *
*       Written by Anders Bernhardsson                             *
*       960427                                                     *
*                                                                  *
********************************************************************
      use Basis_Info
      use Center_Info
      Implicit Real*8(a-h,o-z)
#include "itmax.fh"
#include "info.fh"
      Real*8 A(3),COOR(3,*), B(3)
      Character*(LENIN) Lab
      mdc=0
      iIrrep=0
*
      Write(LUT,'(A)') '*BEGIN COORDINATES'
      Write(LUT,'(A)') '*LABEL COORDINATES CHARGE '
      Do iCnttp=1,nCnttp
         Do iCnt=1,dbsc(iCnttp)%nCntr
            mdc=mdc+1
            call dcopy_(3,Coor(1,mdc),1,A,1)
            Do iCo=0,nIrrep/dc(mdc)%nStab-1
               kop=dc(mdc)%iCoSet(iCo,0)
               Call OA(kOp,A,B)
               ii=nint(dbsc(icnttp)%Charge)
               Lab=dc(mdc)%LblCnt(1:LENIN)
               call setLab(Lab,ico)
               write (LUT,'(1X,A,1X,3F20.10,1X,I3)')
     &                  Lab,B(1:3),ii
             End Do
         End Do
      End Do
      Write(LUT,'(A)') '*END COORDINATES'
      Return
      End
      Subroutine SetLab(label,j)
      Implicit Real*8 (a-h,o-z)
      Character*(*) Label
      Logical no
      no=.true.
      do i=1,LEN(label)
      If (no.and.label(i:i).eq.' ') then
      Write(Label(i:i),'(I1)') j
      no=.false.
      end if
      end do
      Return
      end
