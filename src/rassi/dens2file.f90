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
! Copyright (C) 2019, Roland Lindh                                     *
!***********************************************************************
  subroutine dens2file(array1,array2,array3,adim,lu,adr,iEmpty,iOpt,iGo, &
                       iState,jState)
  use rassi_aux, Only : AO_Mode, Job_Index, nasht_save, CMO1, CMO2,      &
                        DMAB, mTRA, Job1_Old, Job2_Old
  implicit none
  External DDot_
  Real*8 DDot_
#include "stdalloc.fh"
  integer, intent(in) :: adim, lu, iEmpty, iOpt, iGo, iState, jState
  integer, intent(inout) :: adr
  real*8 , intent(inout) :: array1(adim),array2(adim),array3(adim)
  Integer JOB1, JOB2, bdim, IRC, J1, J2
  Real*8, Allocatable:: TRA1(:), TRA2(:)
  Logical Redo_binat

!**********************************************************************
!   This is the generalized disk interface for TDMs.
!   Some of the parameters and variables are explained here.
!
!   AO_Mode: true if TDMs are in the AO basis, otherwise the TDMs are
!            stored in the basis of the active orbitals only (no sym).
!   iEmtpy: the three lowest bits are set if the TDMAB, TSDMAB, and
!           WDMAB, are stored on disk, respectively. That is, for
!           example, if iEmpty=5 the code will write only the first
!           and the last matrix. On read of the second matrix the
!           routine will generate a zero matrix.
!   iOpt: values are 1, or 2, for write and read, respectively.
!   iGo: the three lowest bits are set to tell which matrices are
!        requested. For example, iGo=2, means that only spin-densities
!        are requested.
!
!**********************************************************************
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
       J1=Max(JOB1,JOB2)
       J2=Min(JOB1,JOB2)
       Redo_Binat=.False.
       If (J1.ne.JOB1_Old.or.J2.ne.JOB2_Old) Then
          CALL RDCMO_RASSI(J1,CMO1)
          CALL RDCMO_RASSI(J2,CMO2)
          JOB1_OLD=J1
          JOB2_OLD=J2
          Redo_Binat=.True.
       End If
!
       If (Redo_Binat) Then
          Call mma_Allocate(TRA1,mTra,Label='TRA1')
          Call mma_Allocate(TRA2,mTra,Label='TRA2')
          CALL FINDT(CMO1,CMO2,TRA1,TRA2)
          Call mma_deallocate(TRA2)
          Call mma_deallocate(TRA1)
       End If

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
