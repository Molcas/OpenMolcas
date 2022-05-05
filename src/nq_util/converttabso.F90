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
! Copyright (C) 2021, Jie J. Bao                                       *
!***********************************************************************
! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Dec. 08, 2021, created this file.               *
! ****************************************************************
      Subroutine ConvertTabSO(TabSO2,TabSO,mAO,mGrid,nMOs)
      use nq_pdft, only: lft, lGGA

      INTEGER mAO,mGrid,nMOs,iGrid,nAOGrid
      Real*8 :: TabSO(mAO,mGrid,nMOs)
      Real*8 :: TabSO2(nMOs,mAO*mGrid)

      INTEGER :: iSt, iEnd, iAO, jAO, iOff

      nAOGrid=mAO*mGrid   ! TabSO : mAO*mGrid x nMOs
                          ! TabSO2: nMOs x mAO*nGrid

      ! loop over first and optionally second derivatives of the SOs
      ! this defines the length of nAO to 3 or 9.
      iSt = 1
      If (lft.and.lGGA) Then
         iEnd = 9
      Else
         iEnd = 3
      End If

      Do iGrid=1,mGrid


         Do jAO=iSt, iEnd

            iOff = (iGrid-1)*mAO + jAO

            iAO=jAO+1
            CALL DCopy_(nMOs,TabSO(iAO,iGrid,1),nAOGrid,                &
     &                       TabSO2(:,iOff),1)
         End Do
      End Do

      RETURN
      End Subroutine ConvertTabSO
