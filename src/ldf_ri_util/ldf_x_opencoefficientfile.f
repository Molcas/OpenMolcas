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
* Copyright (C) 2010, Thomas Bondo Pedersen                            *
************************************************************************
      Integer Function LDF_X_OpenCoefficientFile()
C
C     Thomas Bondo Pedersen, July 2010.
C
C     Opens LDF coefficient file and returns logical unit.
C
      Implicit None
      Integer Lu
      Lu=7
      Call DAName_MF_WA(Lu,'LDFC')
      LDF_X_OpenCoefficientFile=Lu
      End
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      Integer Function LDF_X_OpenC()
      Implicit None
      Integer LDF_X_OpenCoefficientFile
      LDF_X_OpenC=LDF_X_OpenCoefficientFile()
      End
