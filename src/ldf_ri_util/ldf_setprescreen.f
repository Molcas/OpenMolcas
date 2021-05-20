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
* Copyright (C) 2011, Thomas Bondo Pedersen                            *
************************************************************************
      Subroutine LDF_SetPrescreen(Value)
C
C     Thomas Bondo Pedersen, January 2011.
C
C     Set LDF prescreening threshold (for determining significant atom
C     pairs).
C
      Implicit None
      Real*8 Value
#include "localdf.fh"
      Thr_Prescreen=Value
      End
