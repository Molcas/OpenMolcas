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
* Copyright (C) 2000-2016, Valera Veryazov                             *
************************************************************************
******************************************************************************
*                                                                            *
* Author:   Valera Veryazov 2000-2016                                        *
*           Theoretical Chemistry                                            *
*           Lund University                                                  *
*           Sweden                                                           *
*                                                                            *
******************************************************************************
      Subroutine systemf(c,rc)
      Implicit None
      Character*(*) C
      Character*1024 C2
      Integer LenC,StrnLn,i,rc

      LenC=StrnLn(C)
      if(LenC.gt.1024-1)then
        Write(6,*)' Error in systemf.f ! LenC :',LenC
        call abend()
      endif

      do 100 i=1,lenc
100   c2(i:i)=c(i:i)
      call systemc(c2,lenc,rc)
      Return
      End
      Subroutine systemf2(c,rc)
      Implicit None
      Character*(*) C
      Character*1024 C2
      Integer LenC,StrnLn,i,rc

      LenC=StrnLn(C)
      if(LenC.gt.1024-1)then
        Write(6,*)' Error in systemf.f ! LenC :',LenC
        call abend()
      endif

      do 100 i=1,lenc
100   c2(i:i)=c(i:i)
      call systemc2(c2,lenc,rc)
      Return
      End
