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
      Subroutine GetGDMat(GDMat)
      use lucia_data, only: DStmp, Dtmp
      use stdalloc, only: mma_allocate, mma_deallocate
      use Lucia_Interface, only: Lucia_Util
      use rasscf_global, only: lRoots, nAc, iAdr15
      use general_data, only: JOBIPH,NCONF
      Implicit None

!     Output
      Real*8,DIMENSION(lRoots*(lRoots+1)/2,NAC,NAC)::GDMat

!     Auxiliary qunatities
      INTEGER CIDisk1,CIDisk2
      INTEGER NIJ2, jRoot, kRoot, iOrb, jOrb
      Real*8, Allocatable:: SDtmp(:), TmpD(:)
      Real*8, Allocatable:: VecL(:), VecR(:)


      Call mma_allocate(VecL,NConf,Label='VecL')
      Call mma_allocate(VecR,NConf,Label='VecR')
      Call mma_allocate(TmpD,NAC**2,Label='TmpD')
      Call mma_allocate(SDtmp,NAC**2,Label='SDtmp')
      SDtmp(:)=DStmp(:)
      DStmp(:)=0.0D0
      TmpD(:)=Dtmp(:)
      Dtmp(:)=0.0D0
      CIDisk1=IADR15(4)
      Do jRoot=1,lRoots
       Call DDafile(JOBIPH,2,VecL,nConf,CIDisk1)
       CIDisk2=IADR15(4)
       Do kRoot=1,jRoot
        Call DDafile(JOBIPH,2,VecR,nConf,CIDisk2)
!        write(6,*) 'VecL and VecR for states',jRoot,kRoot
!        write(6,*)(VecL(I),I=0,NConf-1)
!        write(6,*)(VecR(I),I=0,NConf-1)
        Call Lucia_Util('Densi',                                        &
     &                  CI_Vector=VecL(:),                              &
     &                  RVEC=VECR(:))
!        write(6,*)'GDMat for states',jRoot,kRoot
         dO IOrb=1,NAC
          do JOrb=1,NAC
          NIJ2=jRoot*(jRoot-1)/2+kRoot
          GDMat(NIJ2,JOrb,IOrb)=Dtmp(JOrb+(IOrb-1)*NAC)
          end do
!          write(6,'(10(F8.4,2X))')(GDMat(NIJ2,IOrb,JOrb),JOrb=1,NAC)
         eND dO
       End Do
      End DO
      DStmp(:)=SDtmp(:)
      DTmp(:)=TmpD(:)
      Call mma_deallocate(SDtmp)
      Call mma_deallocate(TmpD)
      Call mma_deallocate(VecL)
      Call mma_deallocate(VecR)

      END Subroutine GetGDMat
