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
* Copyright (C) 1991, Roland Lindh                                     *
************************************************************************
      Subroutine PXMem(nRys,MemPV,la,lb,lr)
      External NAMem, MltMem, EFMem, CntMem
#include "property_label.fh"
*
      If (PLabel.eq.'NAInt ') Then
         Call PVMem(nRys,MemPV,la,lb,lr,NAMem)
      Else If (PLabel.eq.'MltInt') Then
         Call PVMem(nRys,MemPV,la,lb,lr,MltMem)
      Else If (PLabel.eq.'EFInt ') Then
         Call PVMem(nRys,MemPV,la,lb,lr,EFMem)
      Else If (PLabel.eq.'CntInt') Then
         Call PVMem(nRys,MemPV,la,lb,lr,CntMem)
      Else
         Call WarningMessage(2,'PXMem: Illegal type!')
         Write(6,*) '       PLabel=',PLabel
         Call Abend()
      End If
*
      Return
      End
