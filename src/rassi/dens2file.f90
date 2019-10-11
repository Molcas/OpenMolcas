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
  subroutine dens2file(array1,array2,array3,adim,lu,adr,iEmpty,iOpt,iGo)
  implicit none

  integer, intent(in) :: adim, lu, iEmpty, iOpt, iGo
  integer, intent(inout) :: adr
  real*8 , intent(inout) :: array1(adim),array2(adim),array3(adim)

    If (IAND(iGo,1).ne.0) Then
       If (IAND(iEmpty,1).ne.0) Then
          call ddafile(lu,iOpt,array1,adim,adr)
       Else
          If (iOpt.eq.2) array1(:)=0.0D0
       End If
    Else
       If (IAND(iEmpty,1).ne.0) Then
          call ddafile(lu,0,array1,adim,adr)
       End If
    End If
    If (IAND(iGo,2).ne.0) Then
       If (IAND(iEmpty,2).ne.0) Then
          call ddafile(lu,iOpt,array2,adim,adr)
       Else
          If (iOpt.eq.2) array2(:)=0.0D0
       End If
    Else
       If (IAND(iEmpty,2).ne.0) Then
          call ddafile(lu,0,array2,adim,adr)
       End If
    End If
    If (IAND(iGo,4).ne.0) Then
       If (IAND(iEmpty,4).ne.0) Then
          call ddafile(lu,iOpt,array3,adim,adr)
       Else
          If (iOpt.eq.2) array3(:)=0.0D0
       End If
    End If
!
  end subroutine dens2file
