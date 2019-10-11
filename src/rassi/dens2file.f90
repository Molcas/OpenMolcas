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
  subroutine dens2file(array1,array2,array3,adim,lu,adr,iEmpty,iOpt,iGo, &
                       iState,jState)
  use rassi_aux, Only : AO_Mode, Job_Index, nasht_save, CMO1, CMO2
  implicit none

  integer, intent(in) :: adim, lu, iEmpty, iOpt, iGo, iState, jState
  integer, intent(inout) :: adr
  real*8 , intent(inout) :: array1(adim),array2(adim),array3(adim)
  Integer JOB1, JOB2, bdim, IRC

    bdim=adim
    If (.NOT.AO_Mode .AND. iOpt.eq.2) bdim=nasht_save**2

    If (IAND(iGo,1).ne.0) Then
       If (IAND(iEmpty,1).ne.0) Then
          call ddafile(lu,iOpt,array1,bdim,adr)
       Else
          If (iOpt.eq.2) array1(:)=0.0D0
       End If
    Else
       If (IAND(iEmpty,1).ne.0) Then
          call ddafile(lu,0,array1,bdim,adr)
       End If
    End If
    If (IAND(iGo,2).ne.0) Then
       If (IAND(iEmpty,2).ne.0) Then
          call ddafile(lu,iOpt,array2,bdim,adr)
       Else
          If (iOpt.eq.2) array2(:)=0.0D0
       End If
    Else
       If (IAND(iEmpty,2).ne.0) Then
          call ddafile(lu,0,array2,bdim,adr)
       End If
    End If
    If (IAND(iGo,4).ne.0) Then
       If (IAND(iEmpty,4).ne.0) Then
          call ddafile(lu,iOpt,array3,bdim,adr)
       Else
          If (iOpt.eq.2) array3(:)=0.0D0
       End If
    End If
!
    If (.Not.AO_Mode) Then
!
!      Expand the TDMs to AO basis.
!
       Write (6,*) 'Not implemented yet.'
       JOB1=JOB_Index(iState)
       JOB2=JOB_Index(jState)
       CALL RDCMO_RASSI(JOB1,CMO1)
       CALL RDCMO_RASSI(JOB2,CMO2)
!
       If (IAND(iGo,1).ne.0.AND.iAND(iEmpty,1).ne.0) Then
          SCR(:)=array1(1:NASHT**2)
          CALL MKTDZZ(CMO1,CMO2,SCR,array1,iRC)
       End If
!
       If (IAND(iGo,2).ne.0.AND.iAND(iEmpty,2).ne.0) Then
          SCR(:)=array2(1:NASHT**2)
          CALL MKTDZZ(CMO1,CMO2,SCR,array2,iRC)
       End If
!
       If (IAND(iGo,4).ne.0.AND.iAND(iEmpty,4).ne.0) Then
          SCR(:)=array3(1:NASHT**2)
          CALL MKTDZZ(CMO1,CMO2,SCR,array3,iRC)
       End If
!
    End If
!
  end subroutine dens2file
