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
      Subroutine DMP_EFP()
      use EFP_Module
      Integer CoorType
      Call Put_lScalar('EFP',EFP)
      If (EFP) Then
         Call Put_iScalar('nEFP_fragments',nEFP_fragments)
         CoorType=Coor_Type
         Call Put_iScalar('Coor_Type',CoorType)
         Call Put_cArray('FRAG_Type',FRAG_Type,180*nEFP_fragments)
         Call Put_cArray('ABC',ABC,3*180*nEFP_fragments)
         Call Put_iScalar('nEFP_Coor',nEFP_Coor)
         Call Put_dArray('EFP_COORS',EFP_COORS,nEFP_Coor*nEFP_fragments)
      End If
      Return
      End
