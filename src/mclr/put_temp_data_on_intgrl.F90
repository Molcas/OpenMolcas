!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************
      Subroutine put_temp_data_on_intgrl(LUINTMZ_, NSYMZ_, NORBZ_,      &
     &                                   NISHZ_, NASHZ_  )
      use Intgrl, only: NSYMZ,IAD2M,LUINTMZ,NORBZ
      Implicit None
      Integer ::       LUINTMZ_, NSYMZ_
      Integer ::       NORBZ_(8), NISHZ_(8), NASHZ_(8)
      Integer ::       nLength, iAddress, i
      iAddress=0
      IAD2M(1:3,1:36*36)=0
      nLength=3*36*36
      ! read the address list from the existing file
      Call iDaFile(LUINTMZ_,2,IAD2M,nLength,iAddress)
      NSYMZ=NSYMZ_
      LUINTMZ=LUINTMZ_
      Do i=1,NSYMZ_
        NORBZ(i) = NORBZ_(i)
      End Do
      Return
! Avoid unused argument warnings
      If (.False.) Then
        Call Unused_integer_array(NISHZ_)
        Call Unused_integer_array(NASHZ_)
      End If
      End Subroutine put_temp_data_on_intgrl
