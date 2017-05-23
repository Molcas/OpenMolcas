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
* This routine opens an xml container.                                 *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Subroutine xml_open(Name,Appear,Units,Level,Value)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Character*(*) Name
      Character*(*) Appear
      Character*(*) Units
      Character*(*) Value
      Integer Level
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Character*16 myName
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*
      myName=Name
      Call Upcase(myName)
      If(myName.eq.'MODULE') Then
         Call poke_iScalar("xml opened",1)
      End If
      Call xml_openc(Name,Len(Name),
     &               Appear,Len(Appear),
     &               Units,Len(Units),Level,
     &               Value,Len(Value))
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*
      Return
      End
