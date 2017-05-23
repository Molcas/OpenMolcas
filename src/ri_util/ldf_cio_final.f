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
      Subroutine LDF_CIO_Final()
C
C     Thomas Bondo Pedersen, September 2010.
C
C     Purpose: Shut down LDF coefficient I/O.
C
      Implicit None
#include "ldf_cio.fh"

      ! Deallocate buffer
      If (l_LDFC_Buffer.gt.0) Then
         Call GetMem('CBuffer','Free','Real',
     &               ip_LDFC_Buffer,l_LDFC_Buffer)
         ip_LDFC_Buffer=0
         l_LDFC_Buffer=0
      End If
      If (l_LDFC_Blocks.gt.0) Then
         Call GetMem('LDFC_Blk','Free','Inte',
     &               ip_LDFC_Blocks,l_LDFC_Blocks)
         ip_LDFC_Blocks=0
         l_LDFC_Blocks=0
      End If

      ! No atom pairs in buffer
      LastAtomPair=0

      ! Close coefficient file
      If (Lu_LDFC.gt.0) Then
         Call DAClos(Lu_LDFC)
         Lu_LDFC=0
      End If

#if defined (_DEBUG_)
      Write(6,'(/,A)')
     & 'LDF_CIO_Final: coefficient I/O has been shut down!'
      Call xFlush(6)
#endif

      End
