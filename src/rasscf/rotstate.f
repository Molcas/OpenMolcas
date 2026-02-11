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
! Copyright (C) 2020, Jie J. Bao                                       *
!***********************************************************************
      Subroutine RotState()
      use stdalloc, only: mma_allocate, mma_deallocate
      use rasscf_global, only: ICMSP, ITER, IXMSP, LROOTS, IADR15,      &
     &                         Ener
      use PrintLevel, only: DEBUG,USUAL
      use output_ras, only: LF,IPRLOC
      use general_data, only: JOBIPH,NCONF
      use Molcas, only: MxRoot
      use RASDim, only: MxIter
      Implicit None


! ****************************************************************
! history:                                                       *
! Jie J. Bao, on Mar. 13, 2020, created this file.               *
! ****************************************************************


      Integer NHrot                ! storing info in H0_Rotate.txt
      Integer NRState            ! storing info in Do_Rotate.txt
      Integer rcidisk
      INTEGER JRoot,IPRLEV
      CHARACTER(Len=18)::MatInfo
      Real*8, Allocatable:: CIVEC(:,:), CIScr(:,:), HScr(:), State(:),  &
     &                      HRot(:,:)
      Integer i, iad15

      IPRLEV=IPRLOC(3)

      IF(IPRLEV.ge.USUAL) THEN
      write(LF,*)
      write(LF,*) repeat('=',71)
      write(LF,*)
      write(LF,'(11X,A)')'Do_Rotate.txt is found in scratch directory.'
      IF(IXMSP.eq.1) THEN
       write(LF,'(11X,A)')                                              &
     & 'Following properties are for XMS intermediate states.'
      ELSE IF(ICMSP.eq.1) THEN
       write(LF,'(11X,A)')                                              &
     & 'Following properties are for CMS intermediate states.'
      ELSE
       write(LF,'(11X,A)')                                              &
     & 'Following properties are for intermediate states'
       write(LF,'(11X,A)')                                              &
     & ' obtained from the user-supplied rotation matrix'
      ENDIF
      ENDIF

      NRState=lRoots**2
      NHRot=NRState

      CALL mma_allocate(CIVec,nConf,lRoots,Label='CIVec')
      CALL mma_allocate(CIScr,nConf,lRoots,Label='CIScr')
      CALL mma_allocate(HScr,NHRot,Label='HScr')
      CALL mma_allocate(State,NRState,Label='State')
      CALL mma_allocate(HRot,lRoots,lRoots,Label='HRot')


!JB   read rotation matrix in Do_Rotate.txt
      CALL ReadMat2('ROT_VEC',MatInfo,State,lRoots,lRoots,              &
     &              7,18,'T')
      iF(IPRLEV.GE.DEBUG) Then
        write(LF,*)'rotation matrix'
        CALL RecPrt(' ',' ',State,lRoots,lRoots)
      eND iF
      HRot(:,:)=0.0D0
      NHRot=lRoots**2
      Do I=1,lRoots
        HRot(I,I)=ENER(I,ITER)
      End Do
      Call DGEMM_('t','n',lRoots,lRoots,lRoots,1.0D0,State,             &
     &     lRoots,HRot,lRoots,0.0D0,HScr,lRoots)
      Call DGEMM_('n','n',lRoots,lRoots,lRoots,1.0D0,HScr,              &
     &     lRoots,State,lRoots,0.0D0,HRot,lRoots)
      CALL PrintMat2('ROT_HAM',MatInfo,HRot,lRoots,lRoots,              &
     &              7,18,'T')
      if(IPRLEV.GE.DEBUG) Then
       write(LF,'(6X,A)') 'Rotated Hamiltonian matrix '
       Call RecPrt('HRot',' ',hRot,lRoots,lRoots)
      End if

!JB   read CI vector from jobiph
      rcidisk=IADR15(4)
      Do jRoot = 1,lRoots
        Call DDafile(JOBIPH,2,CIScr(:,jRoot),nConf,rcidisk)
      End Do
      Call DGEMM_('n','n',NConf,lRoots,lRoots,1.0D0,CIScr,              &
     &     nConf,State,lRoots,0.0D0,CIVec,nConf)

!     updating final energies as those for rotated states
      rcidisk=IADR15(4)
      Do I=1,lRoots
        Call DDafile(JOBIPH,1,CIVec(:,I),nConf,rcidisk)
        ENER(I,ITER)=HRot(I,I)
      End Do
      IAD15 = IADR15(6)
      CALL DDAFILE(JOBIPH,1,ENER,mxRoot*mxIter,IAD15)


      IF(IPRLEV.GE.DEBUG) Then
      write(LF,'(2A)')'Printing the coeff of the first CSF',            &
     &' for each state'
      Do I=1,lRoots
        write(LF,*) CIVec(1,I)
      End Do
      End If

      Call mma_deallocate(HScr)
      Call mma_deallocate(CIScr)
      Call mma_deallocate(State)
      Call mma_deallocate(CIVec)
      Call mma_deallocate(HRot)

      IF(IPRLEV.ge.USUAL) THEN
      write(LF,*)
      write(LF,*) repeat('=',71)
      END IF

      End Subroutine RotState
