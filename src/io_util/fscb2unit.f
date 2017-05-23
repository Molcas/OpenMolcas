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
* Copyright (C) Victor P. Vysotskiy                                    *
************************************************************************
       Subroutine FSCB2UNIT(cunit,LuP)

************************************************************
*
*   <DOC>
*     <Name>FSCB2UNIT</Name>
*     <Syntax>FSCB2UNIT(cunit)</Syntax>
*     <Arguments>
*       \Argument{cinit}{System (C)file descriptor}{integer}{inout}
*     </Arguments>
*     <Purpose>Translate system (C)file descriptor into internal Molcas's one </Purpose>
*     <Dependencies></Dependencies>
*     <Author>V. Vysotskiy</Author>
*     <Modified_by></Modified_by>
*     <Side_Effects></Side_Effects>
*     <Description>
*       Translate system (C)file descriptor into internal Molcas's one
*     </Description>
*    </DOC>
*
************************************************************
#include "fio.fh"
#include "pfio.fh"

       Integer cunit,i,Lu
*
       Lu=-1
       Do i=1,MxFile
C          print *,i,FSCB(i),cunit
          If(FSCB(i).eq.cunit) Then
          Lu=i
          End If
       End Do
#ifndef _I8_
       Lu=MPUnit(0,Lu)
#endif
       LuP=-1
       If(Lu.eq.-1) Call Abend()
       Do i=1,NProfFiles
C          print *, i,LuNameProf(i),LuName(Lu)
          If(LuNameProf(i).eq.LuName(Lu)) Then
             LuP=i
          End If
       End Do

       If(LuP.eq.-1) Call Abend()

       Return
       End
