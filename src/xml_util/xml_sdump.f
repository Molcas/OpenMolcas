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
* This routine dumps characters in xml format.                         *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
* Author:  Per-Olof Widmark                                            *
*          Lund University, Sweden                                     *
*                                                                      *
************************************************************************
      Subroutine xml_sDump(Name,opt)
      Implicit None
*----------------------------------------------------------------------*
* Dummy arguments                                                      *
*----------------------------------------------------------------------*
      Character*(*) Name
*----------------------------------------------------------------------*
* Local variables                                                      *
*----------------------------------------------------------------------*
      Integer opt
*----------------------------------------------------------------------*
*                                                                      *
*----------------------------------------------------------------------*

      if(opt.eq.0) then
      Call xml_cDumps(Name,Len(Name))
      else
      Call xml_cDumpc(Name,Len(Name))
      endif
*----------------------------------------------------------------------*
* Done                                                                 *
*----------------------------------------------------------------------*

      Return
      End
