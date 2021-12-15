************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      Subroutine Get_Grad(Grad,nGrad)
      Implicit None
      Integer nGrad
      Real*8 Grad(nGrad)
*     Local variables
      Integer :: mGrad=0
      Character(LEN=24), Parameter:: Label='GRAD'
      Logical :: Found=.False.

      Call qpg_dArray(Label,Found,mGrad)
      If(.not.Found .or. nGrad==0) Then
         Call SysAbendmsg('get_grad','Did not find:',Label)
      End If
      If (mGrad/=nGrad) Then
         Write (6,*) 'mGrad=',mGrad
         Write (6,*) 'nGrad=',nGrad
         Call SysAbendmsg('get_grad','mGrad/=nGrad:',Label)
      End If
      Call Get_dArray(Label,Grad,nGrad)

      Return
      End
