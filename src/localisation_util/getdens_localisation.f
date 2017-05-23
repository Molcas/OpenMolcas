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
* Copyright (C) Thomas Bondo Pedersen                                  *
************************************************************************
      SubRoutine GetDens_Localisation(Dens,CMO,nBas,nOcc)
C
C     Author: T.B. Pedersen
C
C     Purpose: compute density from CMOs as Dens = CMO * (CMO)^T
C
      Implicit None
      Real*8  Dens(*), CMO(*)
      Integer nBas, nOcc

      Integer nTBs

      nTBs = max(nBas,1)
      Call DGEMM_('N','T',nBas,nBas,nOcc,1.0d0,CMO,nTBs,CMO,nTBs,
     &           0.0d0,Dens,nTBs)

      End
