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
      Subroutine init_dump
#include "fciqmc_global.fh"
#include "files_fciqmc.fh"
*----------------------------------------------------------------------*
*     Define  file names and unit numbers.                             *
*----------------------------------------------------------------------*
      FnInpOrb='INPORB'
      FnJobIph='JOBIPH'
      FnOneAO='ONEINT'
      FnTwoAO='ORDINT'
      FnOneMO='TRAONE'
      FnTwoMO='TRAINT'
      FnHalf='TEMP1'
      FnExt='EXTRACT'
      FnCom='COMFILE'
      Debug=0
      iPrint=0
      iOneOnly=0
      iCTonly=0
      iDoInt=0
      iVecTyp=2
      iAutoCut=0
      iRFpert=0
      LuInpOrb=10
      LuJobIph=15
      LuOneAO=20
      LuTwoAO=40
      LuOneMO=30
      LuTwoMO=50
      LuHalf=60
      LuExt=18
      LuCom=22
*----------------------------------------------------------------------*
      End
