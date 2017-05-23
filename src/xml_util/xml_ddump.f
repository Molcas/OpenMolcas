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
* This routine dumps doubles in xml format.                            *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Subroutine xml_dDump(Name,Appear,Units,Level,Data,nx,ny)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Character*(*) Name
      Character*(*) Appear
      Character*(*) Units
      Real*8        Data(*)
      Integer       nx
      Integer       ny
      Integer       Level
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      Call xml_ddumpc(Name,Len(Name),
     &                Appear,Len(Appear),
     &                Units,Len(Units),Level,
     &                Data,nx,ny)
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
