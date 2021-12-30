!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2005, Thomas Bondo Pedersen                            *
!***********************************************************************
      SubRoutine RotateOrb_ER(R,CMO,nBasis,nOrb2Loc,Debug)
!
!     Thomas Bondo Pedersen, November 2005.
!
!     Purpose: rotate ER orbitals,
!              CMO -> CMO * U
!              U = R*[R^T*R]^(-1/2)
!
      Implicit Real*8 (a-h,o-z)
      Real*8  R(nOrb2Loc,nOrb2Loc), CMO(nBasis,nOrb2Loc)
      Logical Debug
#include "WrkSpc.fh"

      Character*12 SecNam
      Parameter (SecNam = 'RotateOrb_ER')

      If (nOrb2Loc.lt.1 .or. nBasis.lt.1) Return

!     Allocate transformation matrix U.
!     ---------------------------------

      lU = nOrb2Loc**2
      Call GetMem('Umat','Allo','Real',ipU,lU)

!     Compute U.
!     ----------

      Call GetU_ER(Work(ipU),R,nOrb2Loc)

!     Debug: check that U is unitary.
!     -------------------------------

      If (Debug) Then
         ThrU = 1.0d-10
         irc = -1
         Call Chk_Unitary(irc,Work(ipU),nOrb2Loc,ThrU)
         If (irc .ne. 0) Then
            Call SysAbendMsg(SecNam,'U matrix is not unitary!',' ')
         End If
      End If

!     Update C.
!     ---------

      lCMO = nBasis*nOrb2Loc
      Call GetMem('CMOscr','Allo','Real',ipCMO,lCMO)
      Call dCopy_(lCMO,CMO,1,Work(ipCMO),1)
      Call DGEMM_('N','N',nBasis,nOrb2Loc,nOrb2Loc,                     &
     &           1.0d0,Work(ipCMO),nBasis,Work(ipU),nOrb2Loc,           &
     &           0.0d0,CMO,nBasis)
      Call GetMem('CMOscr','Free','Real',ipCMO,lCMO)

!     De-allocate U.
!     --------------

      Call GetMem('Umat','Free','Real',ipU,lU)

      End
