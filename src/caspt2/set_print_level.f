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
      Subroutine Set_Print_Level
      use output, only:silent,terse,usual,verbose,debug,insane,iPrGlb
      Implicit None
      Integer, External :: iPrintLevel
      Logical, External :: Reduce_Prt
* Print levels request from molcas?
      IPRGLB=IPRINTLEVEL(-1)
* If inside an optimization loop, minimize output
* unless we *really* want a lot of output.
      If (Reduce_Prt()) THEN
       IPRGLB=MAX(IPRGLB-USUAL,SILENT)
      END IF
      END
