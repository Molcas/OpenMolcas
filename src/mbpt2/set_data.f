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
      Subroutine Set_Data
      Implicit Real*8 (a-h,o-z)
*
#include "mxdim.fh"
#include "cdtfaux.fh"
#include "files_mbpt2.fh"
*
CStart Molcas
*
*     Aux
*
      Thize=1.0d-6
      SIntTh=1.0d-14
      TfThre=1.0d-14
      nDisc =0
CEnd
*                                                                      *
************************************************************************
*                                                                      *
*     Set up some file names
*
      LUINTA=40
      FNINTA='ORDINT'
      LUHLF1=50
      FNHLF1='LUHLF1'
      LUHLF2=60
      FNHLF2='LUHLF2'
      LUHLF3=70
      FNHLF3='LUHLF3'
      LUINTM=80
      FNINTM='MOLINT'
*
      Return
      End
