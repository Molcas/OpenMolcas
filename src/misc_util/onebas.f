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
      Subroutine OneBas(Label)
************************************************************************
*                                                                      *
*     Change nBas in OneDat.fh                                         *
*                                                                      *
************************************************************************
*
      Implicit Integer (A-Z)
      Integer IntBas(8)
      Character*(*) Label
*
#include "OneDat.fh"
*
      If (Label.eq.'CONT') Then
         Call Get_iArray('nBas',IntBas,nSym)
      Else If (Label.eq.'PRIM') Then
         Call Get_iArray('nBas_Prim',IntBas,nSym)
      Else
         Write (6,*) 'OneBas: Illegal Label value!'
         Write (6,*) 'Value: ',Label
         Call Abend()
      End If
      Call ICopy(nSym,IntBas,1,nBas,1)
      Return
      End
