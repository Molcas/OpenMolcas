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
* Copyright (C) 2017, Roland Lindh                                     *
************************************************************************
      Subroutine GET_EFP()
      use EFP_Module
      Integer CoorType
      Call Get_lScalar('EFP',EFP)
      If (EFP) Then
         Call Get_iScalar('nEFP_fragments',nEFP_fragments)
         Call Get_iScalar('nEFP_Coor',nEFP_Coor)
         Call Get_iScalar('Coor_Type',CoorType)
         Coor_Type=INT(CoorType,c_int)
*
         Allocate (FRAG_type(nEFP_fragments))
         Call Get_cArray('FRAG_Type',FRAG_Type,180*nEFP_fragments)
*
         Allocate (ABC(3,nEFP_fragments))
         Call Get_cArray('ABC',ABC,180*3*nEFP_fragments)
*
         Allocate (EFP_Coors(nEFP_Coor,nEFP_fragments))
         Call Get_dArray('EFP_COORS',EFP_COORS,nEFP_Coor*nEFP_fragments)
      End If
      Return
      End
