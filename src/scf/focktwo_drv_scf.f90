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
      Subroutine FockTwo_Drv_scf(nSym,nBas,nAux,Keep,DLT,DSQ,FLT,nFLT,ExFac,nBSQT,nBMX,nD,nOcc,lOcc,iDummy_run)
      use OFembed, only: Do_OFemb,OFE_first,FMaux
      use OFembed, only: Rep_EN
      use ChoSCF, only: ALGO
      use Cholesky, only: timings
      use Constants, only: Zero
      use stdalloc, only: mma_allocate, mma_deallocate
      Implicit None
      Integer nSym, nBas(8), nAux(8), Keep(8), nFLT, nBSQT, nBMX, nD, lOcc, iDummy_run
      Real*8 FLT(nFLT,nD)
      Real*8 DLT(nFLT,nD)
      Real*8 DSQ(nBSQT,nD)
      Integer nOcc(lOcc,nD)
      Real*8 ExFac

      Logical DoCholesky,GenInt
      Character(LEN=50) CFmt
      Integer iRC, lBuf
      Real*8 TotCPU, TotCPU1, TotCPU2
      Real*8 TotWall, TotWall1, TotWall2

      Real*8, Allocatable :: FSQ(:,:)
      Real*8, Allocatable :: W1(:), W2(:)
      Real*8, Allocatable :: tFLT(:,:)
!
      GenInt=.false.
      DoCholesky=.false.
      if(ALGO.eq.0) GenInt=.true. !use GenInt to regenerate integrals
      Call DecideOnCholesky(DoCholesky)

!      write(6,*)'*************************'
!      write(6,*)'ONLY COULOMB CONTRIBUTION'
!      write(6,*)'*************************'
!      exFac=0.d0
!      write(6,*)'ExFac= ',ExFac
!
      If (Do_OFemb) Then ! Coul. potential from subsys B
         If (OFE_first) Call mma_allocate(FMaux,nFlt,Label='FMaux')

         Call Coul_DMB(OFE_first,nD,Rep_EN,FMaux,DLT(:,1),DLT(:,nD),nFlt)
         OFE_first=.false.
      End If
!
      Call mma_allocate(FSQ,nBSQT,nD,Label='FSQ')
      FSQ(:,:)=Zero

      if ((.not.DoCholesky).or.(GenInt)) then
         Call mma_Allocate(W2,NBMX*NBMX,Label='W2')
      end if
!
! nFlt is the total dimension of the LT fock matrix
      Call mma_allocate(tFLT,nFLT,nD,Label='tFLT')
      tFLT(:,:)=Zero
!
!
      Call mma_maxDBLE(LBUF)
!
! Standard building of the Fock matrix from Two-el integrals
! or
! Building of the Fock matrix regenerating the integrals on the fly
!
      Call CWTIME(TotCPU1,TotWALL1)

      IF (.not.DoCholesky .or. (DoCholesky.and.GenInt) ) THEN

         IF (DoCholesky.and.GenInt) Then
            ! save some space for GenInt
            LBUF = MAX(LBUF-LBUF/10,0)
            ! Make sure that the ri/ch vectors are in reordered mode
            Call Cho_X_ReOVec(irc)
         End If
         Call mma_allocate(W1,LBUF,Label='W1')
!
       If (LBUF.LT.NBMX**2) Then
         WRITE(6,*)'FockTwo_Drv_SCF Error: Too little memory remains for the call to FOCKTWO_SCF.'
         WRITE(6,*)' Largest allocatable array size LBUF=',LBUF
         WRITE(6,*)' Max nr of bf in any symmetry,  NBMX=',NBMX
         WRITE(6,*)' Required minimum size       NBMX**2=',NBMX**2
         WRITE(6,*)'    (All in Real*8-size words)'
         Call  ABEND()
       End If
!
       Call FOCKTWO_scf(nSym,nBas,nAux,Keep,DLT,DSQ,tFLT,nFlt,FSQ,W1,Size(W1),W2,Size(W2),ExFac,nD,nBSQT)

      ENDIF
!
      Call CWTIME(TotCPU2,TotWALL2)
      TOTCPU  = TotCPU2 - TotCPU1
      TOTWALL = TotWALL2 - TotWALL1

!---- Timings information for conventional or Cholesky with ALGO=0
      IF ( .not.DoCholesky .or. (DoCholesky.and.GenInt) ) THEN
      if(timings)then
      CFmt='(2x,A)'
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      if(DoCholesky)then
      Write(6,CFmt)'---    Cholesky SCF - Integral regeneration   ---'
      else
      Write(6,CFmt)'-----------     Conventional SCF     ------------'
      endif
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,CFmt)'Fock matrix construction        CPU       WALL   '
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,'(2x,A26,2f10.2)')'TOTAL                                     ',TOTCPU,TOTWALL
      Write(6,CFmt)'- - - - - - - - - - - - - - - - - - - - - - - - -'
      Write(6,*)
      endif
      ENDIF
!
! Building of the Fock matrix directly from Cholesky vectors
!

      IF (DoCholesky .and. .not.GenInt.and.iDummy_run.eq.1) THEN
       Write(6,*) '*** Warning: missing feature in Cholesky code'
       Write(6,*) 'Use the results with extra care!'
      endif

      IF (DoCholesky .and. .not.GenInt.and.iDummy_run.eq.0) THEN
!
         CALL CHOscf_drv(nBSQT,nD,nSym,nBas,DSQ(:,1),DLT(:,1),DSQ(:,nD),DLT(:,nD),tFLT(:,1),tFLT(:,nD),nFLT,ExFac,     &
                         FSQ,nOcc(:,1),nOcc(:,nD))

      ENDIF
!
      FLT(:,:)=FLT(:,:)+tFLT(:,:)
!
      Call mma_deallocate(tFLT)
!
      If (Do_OFemb) Then ! add FM from subsystem B
        FLT(:,1)=FLT(:,1)+FMaux(:)
        If (nD==2) FLT(:,2)=FLT(:,2)+FMaux(:)
      EndIf
!
      IF ((.not.DoCholesky).or.(GenInt)) THEN
          Call mma_deallocate(W1)
          Call mma_deallocate(W2)
      END IF

      Call mma_deallocate(FSQ)
!
      Return
      End Subroutine FockTwo_Drv_scf
