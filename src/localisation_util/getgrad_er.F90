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
      SubRoutine GetGrad_ER(Functional,GradNorm,R,CMO,nBasis,nOrb2Loc,  &
     &                      Timing)
!
!     Thomas Bondo Pedersen, November 2005.
!
!     Purpose: compute ER Functional and its gradient norm.
!              The R matrix is computed as R(i,j) = (ij|jj).
!
!     Note: symmetry is NOT allowed (but is not tested!).
!
      use Data_structures, only: DSBA_Type
      use Data_structures, only: Allocate_DSBA, Deallocate_DSBA
      Implicit Real*8 (a-h,o-z)
      Real*8  R(nOrb2Loc,nOrb2Loc), CMO(nBasis,nOrb2Loc)
      Logical Timing

      Character(LEN=10), Parameter:: SecNam = 'GetGrad_ER'

      Character*80 Txt

      Integer, Parameter:: nSym = 1
      Integer nOcc(nSym)

      Type (DSBA_Type) CMOt

!     Initialization.
!     ---------------

      Functional = 0.0d0
      GradNorm = 0.0D0
      If (nOrb2Loc.lt.1 .or. nBasis.lt.1) Return

!     Transpose CMO (only the part to be localised).
!     ----------------------------------------------

      Call Allocate_DSBA(CMOt,[nOrb2Loc],[nBasis],nSym)
      Do i = 1,nOrb2Loc
         CMOt%SB(1)%A2(i,:) = CMO(:,i)
      End Do

!     Compute R.
!     ----------

      nOcc(1) = nOrb2Loc
      irc = -1
      Call Cho_Get_Rij(irc,CMOt,nOcc,R,Timing)
      If (irc .ne. 0) Then
         Write(Txt,'(A,I6)') 'Cho_Get_Rij returned',irc
         Call SysAbendMsg(SecNam,'Calculation of ER gradient failed:',  &
     &                    Txt)
      End If
      Call Deallocate_DSBA(CMOt)

!     Compute gradient norm and functional.
!     -------------------------------------

      Do i = 1,nOrb2Loc-1
         Functional = Functional + R(i,i)
         Do j = i+1,nOrb2Loc
            GradNorm = GradNorm + (R(i,j)-R(j,i))**2
         End Do
      End Do
      Functional = Functional + R(nOrb2Loc,nOrb2Loc)
      GradNorm = 4.0d0*sqrt(GradNorm)

      End
