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
************************************************************************
*                                                                      *
* This routine query the existence of array data on runfile.           *
*                                                                      *
************************************************************************
      Subroutine LookUp_label(i, Type, Label)
#include "pg_ca_info.fh"
*----------------------------------------------------------------------*
* Arguments                                                            *
*----------------------------------------------------------------------*
      Character*(*) Type
      Character* 16 Label
*----------------------------------------------------------------------*
* Define local variables                                               *
*----------------------------------------------------------------------*
      Character*16 RecLab(256)
*
      Integer      nTmp,iTmp
      Integer      i
*----------------------------------------------------------------------*
* Initialize local variables                                           *
*----------------------------------------------------------------------*
*----------------------------------------------------------------------*
* Read info from runfile.                                              *
*----------------------------------------------------------------------*
      Call ffRun(Type,nTmp,iTmp)
      Call cRdRun(Type,RecLab,16*256)
      Label=RecLab(i)
      Return
      End
