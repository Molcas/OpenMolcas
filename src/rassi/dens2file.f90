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
  use rassi_aux, Only : AO_Mode, Job_Index, nasht_save, CMO1, CMO2,      &
                        DMAB, mTRA
  implicit none
  External DDot_
  Real*8 DDot_
#include "stdalloc.fh"
  integer, intent(in) :: adim, lu, iEmpty, iOpt, iGo, iState, jState
  integer, intent(inout) :: adr
  real*8 , intent(inout) :: array1(adim),array2(adim),array3(adim)
  Integer JOB1, JOB2, bdim, IRC, J1, J2
  Real*8, Allocatable:: TRA1(:), TRA2(:)

    bdim=adim
    If (.NOT.AO_Mode .AND. iOpt.eq.2) bdim=nasht_save**2+1
    If (bdim.gt.adim) Then
       Write (6,*) 'Dens2file: bdim.gt.adim'
       Call Abend()
    End If

!   Write (*,*) 'iState,jState=',iState,jState
    If (IAND(iGo,1).ne.0) Then
       If (IAND(iEmpty,1).ne.0) Then
!         If (iOpt.eq.1) Write (6,*) 'array1=',DDot_(bdim-1,Array1,1,Array1,1), Array1(bdim)
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
!         If (iOpt.eq.1) Write (6,*) 'array2=',DDot_(bdim-1,Array2,1,Array2,1), Array2(bdim)
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
!         If (iOpt.eq.1) Write (6,*) 'array3=',DDot_(bdim-1,Array3,1,Array3,1), Array3(bdim)
          call ddafile(lu,iOpt,array3,bdim,adr)
       Else
          If (iOpt.eq.2) array3(:)=0.0D0
       End If
    End If
!
    If (.Not.AO_Mode.AND.iOpt.eq.2) Then
!
!      Expand the TDMs to AO basis.
!
       JOB1=JOB_Index(iState)
       JOB2=JOB_Index(jState)
!      Write (6,*) 'iState,Job1=',iState,Job1
!      Write (6,*) 'jState,Job2=',jState,Job2
       J1=Max(JOB1,JOB2)
       J2=Min(JOB1,JOB2)
       CALL RDCMO_RASSI(J1,CMO1)
       CALL RDCMO_RASSI(J2,CMO2)
       Call mma_Allocate(TRA1,mTra,Label='TRA1')
       Call mma_Allocate(TRA2,mTra,Label='TRA2')
       CALL FINDT(CMO1,CMO2,TRA1,TRA2)
       Call mma_deallocate(TRA2)
       Call mma_deallocate(TRA1)

!
       If (IAND(iGo,1).ne.0.AND.iAND(iEmpty,1).ne.0) Then
!         Write (*,*) 'array1=',DDot_(nasht_save**2,Array1,1,Array1,1),Array1(bdim)
          CALL MKTDAB(array1(bdim),array1,DMAB,iRC)
          CALL MKTDZZ(CMO1,CMO2,DMAB,array1,iRC)
       End If
!
       If (IAND(iGo,2).ne.0.AND.iAND(iEmpty,2).ne.0) Then
!         Write (*,*) 'array2=',DDot_(nasht_save**2,Array2,1,Array2,1),Array2(bdim)
          CALL MKTDAB(array2(bdim),array2,DMAB,iRC)
          CALL MKTDZZ(CMO1,CMO2,DMAB,array2,iRC)
       End If
!
       If (IAND(iGo,4).ne.0.AND.iAND(iEmpty,4).ne.0) Then
!         Write (*,*) 'array3=',DDot_(nasht_save**2,Array3,1,Array3,1),Array3(bdim)
          CALL MKTDAB(array3(bdim),array3,DMAB,iRC)
          CALL MKTDZZ(CMO1,CMO2,DMAB,array3,iRC)
       End If
!
    End If
!
  end subroutine dens2file
