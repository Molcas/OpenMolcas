************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2021, Rulin Feng                                       *
************************************************************************
*       ****************************************************
*                  Do SVD for SO-TDM in AO basis
*       ****************************************************
*        This routine is made to do single value decompositon
*        to spin-orbit coupled transition density matrices in AO
*        basis.
*        Input: TDMZZ, transition densitry matrix
*               TSDMZZ,transition spin density matrix(x,y,z)
*        Remember that TDMZZ and TSDMZZ contain six TDM's each
*        See sonatorbm_full.f
*        TDMZZ(3,:) and TSDMZZ(1-3,:), for the real part of TDM
*        The transition density matrix TDMZZ does not depend on
*        spin matrices, thus TDMZZ(1-3,:) are the same, so are the
*        imaginary part TDMZZ(4-6,:).
*        TDMZZ(6,:) and TSDMZZ(4-6,:), for the imaginary part of TDM
*
*                                                       -RF 8/18,2021

      SUBROUTINE DO_AOTDMNTO(TDMZZ,TSDMZZ,ANTSIN,ISTATE,JSTATE,nb,nb2)
      use OneDat, only: sNoNuc, sNoOri, sOpSiz
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: IfArgu

      IMPLICIT None
#include "rassi.fh"
      Integer ISTATE,JSTATE,nb,nb2
      REAL*8 TDMZZ(6,nb2)
      REAL*8 TSDMZZ(6,nb2)
      REAL*8 ANTSIN(6,nb2)
      REAL*8, ALLOCATABLE:: TDMZZL(:,:), TSDMZZL(:,:)
      COMPLEX*16, ALLOCATABLE:: YMAT(:,:)
      COMPLEX*16, ALLOCATABLE:: TDMZZLC(:),TDMZZC(:),BUFF(:),DIPsC(:)
      COMPLEX*16, ALLOCATABLE:: SVDU(:),SVDVH(:),RESI(:)
      REAL*8, ALLOCATABLE:: SVDS(:)
      COMPLEX*16, ALLOCATABLE:: BUFF1(:),BUFF2(:),SumofYdiag(:)
      COMPLEX*16  Transition_Dipole
      Integer i,j,info,lwork,di,icmp,iopt,irc,isylab
      REAL*8 NumofEc, Sumofeigen, eigen_print_limit,Zero,Two,pi
      REAL*8 SumofTDMZZLC
      REAL*8 Dummy(1)
      Integer, DIMENSION(1)::SIZ
      COMPLEX*16, ALLOCATABLE:: SIZC(:)
      Integer LU, isfreeunit, iDummy(7,8)
c start Phase factor stuff
c trace of transition dipole real and imaginary (x,y,and z)
      REAL*8 ttdr(3),ttdi(3)
      REAL*8 phi,sd
c end
      CHARACTER(LEN=8) LABEL
      CHARACTER(LEN=7) STATENAME,STATENAMETMP
      CHARACTER(LEN=128) FNAME
      CHARACTER(LEN=72) NOTE
      Real*8, Allocatable:: Dips(:), Dip(:)
      Real*8, Allocatable:: TMPR(:), TMPI(:)
      Real*8, Allocatable:: SZZs(:), SZZ(:)
      Real*8, Allocatable:: EIG(:), EIGM(:)
      Real*8, Allocatable:: TMP(:)
      Real*8, Allocatable:: BFF(:)
      Real*8, Allocatable:: RESIX(:), RESIR(:)
      Real*8, Allocatable:: SM(:), SMI(:)
      Real*8, Allocatable:: SVDUR(:), SVDUI(:)
      Integer, Allocatable:: PIV(:)
      Real*8, Allocatable:: SVDVHR(:), SVDVHI(:)
      Real*8, Allocatable:: SVDVR(:), SVDVI(:)

      Zero=0.0D0
      Two=2.0D0
      pi=ACOS(-1.0D0)

c ANTISYMMETRIC matrix needs a little fixing
      do i=0, nb-1
        do j=1, nb
          if(i.LT.j-1) then
            ANTSIN(3,i*nb+j)=-ANTSIN(3,i*nb+j)
          else if(i.EQ.j-1) then
            ANTSIN(3,i*nb+j)=zero
          endif
        enddo
      enddo
c The imaginary part may need a negative sign
      TDMZZ(4,:) = -TDMZZ(4,:)
      TDMZZ(5,:) = -TDMZZ(5,:)
      TDMZZ(6,:) = -TDMZZ(6,:)
c     'HERMSING' ITYPE=1
c     'ANTISING' ITYPE=2
c     'HERMTRIP' ITYPE=3
c     'ANTITRIP' ITYPE=4

C Thus we obtained the AO based transition density matrix TDMZZ
c and the transition spin density matrix TSDMZZ
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C we can do testing here we calculate the oscillator strength
c by dot with dipole moment AO matrix and trace
c do the test before the Lowdin orthogonalization
c      If (TestPrint) then
      Call MMA_ALLOCATE(BUFF,nb2,LABEL="LBUFF")
      Call MMA_ALLOCATE(TDMZZC,nb2,LABEL="TDMZZC")
      do di=1, 3
        Call mma_allocate(DIPs,nb2,Label='Dips')
        LABEL(1:8)='MLTPL  1'
        IRC = -1
        ICMP = di
        ISYLAB = 1
        IOPT = ibset(0,sOpSiz) ! Only read the size of the array
        CALL IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
        !no nuclear contrib, no origin of operator
        IOPT = ibset(ibset(0,sNoOri),sNoNuc)
        Call mma_allocate(DIP,SIZ(1),Label='DIP')
        CALL RDONE(IRC,IOPT,LABEL,ICMP,DIP,ISYLAB)
        Call DESYM_SONTO(DIP,SIZ(1),DIPs,ISYLAB)
        Call mma_deallocate(DIP)
        write(6,*) '  For istate ', ISTATE, ' jstate ',JSTATE
        write(6,*) '  Component ',ICMP
c Get complex matrices
        Call MMA_ALLOCATE(DIPsC,nb2,LABEL="DIPsC")
        do i=1,nb2
          TDMZZC(i)=cmplx(TDMZZ(di,i),TDMZZ(di+3,i),8)
          DIPsC(i)=cmplx(DIPs(i),zero,8)
        enddo
c TDM
        Call ZGEMM_('N','N',nb,nb,nb,(1.0D0,0.0D0),DIPsC,nb,
     &              TDMZZC,nb,(0.0D0,0.0D0),BUFF,nb)
C Trace the resulting matrix
        Transition_Dipole = cmplx(zero,zero,8)
        do i=1, nb
          Transition_Dipole = Transition_Dipole +
     &    BUFF((i-1)*nb+i)
        enddo
        ttdr(di)=real(Transition_Dipole)
        ttdi(di)=aimag(Transition_Dipole)
        write(6,*) '  Transition_Dipole :',Transition_Dipole
        write(6,*)
        Call MMA_DEALLOCATE(DIPs)
        Call MMA_DEALLOCATE(DIPsC)
      enddo
      Call MMA_DEALLOCATE(BUFF)
      Call MMA_DEALLOCATE(TDMZZC)
c      Endif
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c phase factor = cos phi + i * sin phi
c phi = 1/2 * arc tan (2*(x_r*x_i + y_r*y_i + z_r*z_i)/
c                        (x_i**2 + y_i**2 + z_i**2 -
c                         x_r**2 - y_r**2 - z_r**2))
c to minimize the imaginary part of total transition dipole
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      If(IFARGU) then
c check if minimum of maximum
        If((ttdi(1)**2+ttdi(2)**2+ttdi(3)**2
     &     -ttdr(1)**2-ttdr(2)**2-ttdr(3)**2).EQ.0.0D0) Then
          phi = 0.0D0
        Else
          phi = 0.5D0*atan(2*(ttdr(1)*ttdi(1)
     &                   +ttdr(2)*ttdi(2)
     &                   +ttdr(3)*ttdi(3))/
     &                   (ttdi(1)**2+ttdi(2)**2+ttdi(3)**2
     &                   -ttdr(1)**2-ttdr(2)**2-ttdr(3)**2))
        Endif
        sd = 2*cos(2*phi)*(ttdr(1)**2+ttdr(2)**2+ttdr(3)**2
     &                    -ttdi(1)**2-ttdi(2)**2-ttdi(3)**2)
     &      -4*sin(2*phi)*(ttdr(1)*ttdi(1)+ttdr(2)*ttdi(2)
     &                    +ttdr(3)*ttdi(3))
c make sure it's minimum
        If (sd.LT.Zero) phi=phi+pi/Two
c multipole phase factor with tdm as a whole
        write(6,*) 'Phase factor turned on with calculated'
        write(6,'(2X,A,F6.2)') "argument Phi: ", phi
        call mma_allocate(TMPR,nb2,Label='TMPR')
        call mma_allocate(TMPI,nb2,Label='TMPI')
        do i=1,nb2
          TMPR(i)=TDMZZ(3,i)*cos(phi)-TDMZZ(6,i)*sin(phi)
          TMPI(i)=TDMZZ(6,i)*cos(phi)+TDMZZ(3,i)*sin(phi)
        enddo
        TDMZZ(1,:) = TMPR(1:nb2)
        TDMZZ(2,:) = TMPR(1:nb2)
        TDMZZ(3,:) = TMPR(1:nb2)
        TDMZZ(4,:) = TMPI(1:nb2)
        TDMZZ(5,:) = TMPI(1:nb2)
        TDMZZ(6,:) = TMPI(1:nb2)
        call mma_deallocate(TMPI)
        call mma_deallocate(TMPR)
      EndIf
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

C Make NTO output file names without spin
      write(STATENAME,'(I3)') ISTATE
      write(STATENAMETMP,'(I3,a1,a)')
     &      JSTATE,'_',trim(adjustl(STATENAME))
      write(STATENAME,'(a)') trim(adjustl(STATENAMETMP))
Cc Everything is in C1 symmetry for now
C Do the Lowdin Orthogonalization assuming C1 symmetry
c SZZ  - AO Overlap integral
c SZZs - AO Overlap integral in square
c EIG  - AO Overlap eigenvalues
      call mma_allocate(SZZ,NBTRI,Label='SZZ')
      call mma_allocate(SZZs,nb2,Label='SZZs')
      call mma_allocate(EIG,nb,Label='EIG')
      SZZ(:)=0.0D0
      SZZs(:)=0.0D0
      EIG(:)=0.0D0
C AO OVERLAP MATRIX
      IRC=-1
c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      call RDONE(IRC,IOPT,LABEL,ICMP,SZZ,ISYLAB)
      IF (IRC.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
      Call DESYM_SONTO(SZZ,NBTRI,SZZs,ISYLAB)
      call mma_deallocate(SZZ)
**************
* For tests
**************
      call mma_allocate(TMP,nb2,Label='TMP')
      TMP(1:nb2) = TDMZZ(3,:)
      Call mma_allocate(BFF,nb2,Label='BFF')
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SZZs,nb,
     &             TMP,nb,0.0D0,BFF,nb)
C Trace the resulting matrix
      NumOfEc = Zero
      do i=0, nb-1
        do j=0, nb-1
          if(i.eq.j) then
            NumOfEc = NumOfEc + BFF(1+i*nb+j)
          endif
        enddo
      enddo
c      write(6,*) 'NumOfEc ',NumOfEc
      Call mma_deallocate(BFF)
**************

c DIAGONALIZE AO OVERLAP MATRIX
c Set LWORK=-1 to get the optimal scratch space RESI
c then let LWORK equal to length of scratch space
c free and reallocate memory for RESI using that length
      call mma_allocate(RESIX,1,Label='RESIX')
      LWORK=-1
      call DSYEV_('V','U',nb,SZZs,nb,EIG,
     &            RESIX,LWORK,INFO)
      LWORK=INT(RESIX(1))
      call mma_deallocate(RESIX)
      call mma_allocate(RESIX,LWORK,Label='RESIX')
c SZZs in as the AO overlap sqaure matrix
c out as the eigenvector matrix of SZZs
c with eigenvalues in EIG
      call DSYEV_('V','U',nb,SZZs,nb,EIG,
     &            RESIX,LWORK,INFO)
c Put EIG in sqrt and in diagonal in EIGM
      call mma_allocate(EIGM,nb2,Label='EIGM')
      do i=0,nb-1
        do j=0,nb-1
          If (i.eq.j) then
            EIGM(1+i*nb+j)=sqrt(EIG(1+i))
          Else
            EIGM(1+i*nb+j)=zero
          Endif
        enddo
      enddo
      call mma_deallocate(EIG)
      call mma_deallocate(RESIX)
c Get S^1/2 from S^1/2 = U S_diag^1/2 U^T
      call mma_allocate(SM,nb2,Label='SM')
      call DGEMM_('N','T',nb,nb,nb,1.0D0,EIGM,nb,
     &             SZZs,nb,0.0D0,TMP,nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SZZs,nb,
     &             TMP,nb,0.0D0,SM,nb)
      call mma_deallocate(SZZs)
      call mma_deallocate(EIGM)
c Get inverse of S^1/2 -> S^-1/2
c Before calling DGETRI, call DGETRF to factorize SM
c Set LWORK=-1 to get the optimal scratch space RESI
c then let LWORK equal to length of scratch space
c free and reallocate memory for RESI using that length
      call mma_allocate(SMI,nb2,Label='SMI')
      SMI(:)=SM(:)
      call mma_allocate(PIV,nb,Label='PIV')
      call DGETRF_(nb,nb,SMI,nb,PIV,INFO)
      call mma_allocate(RESIX,1,Label='RESIX')
      LWORK=-1
      call DGETRI_(nb,SMI,nb,PIV,RESIX,LWORK,INFO)
      LWORK=INT(RESIX(1))
      call mma_deallocate(RESIX)
      call mma_allocate(RESIX,LWORK,Label='RESIX')
      call DGETRI_(nb,SMI,nb,PIV,RESIX,LWORK,INFO)
      call mma_deallocate(PIV)
      call mma_deallocate(RESIX)

c Note: The density matrix should transform as S^1/2 D S^1/2
c Transform TDMZZ and TSDMZZ as S^1/2 T S^1/2
      Call MMA_ALLOCATE(TDMZZL,6,nb2,LABEL='LTDMZZL')
      Call MMA_ALLOCATE(TSDMZZL,6,nb2,LABEL='LTSDMZZL')
      call mma_allocate(TMPR,nb2,Label='TMPR')
c Real part of TDMZZ
      TMPR(1:nb2) = TDMZZ(3,:)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,TMPR,nb,
     &               SM,nb,0.0D0,TMP,nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SM,nb,
     &               TMP,nb,0.0D0,TMPR,nb)
      TDMZZL(3,:) = TMPR(1:nb2)
c Imaginary part of TDMZZ
      TMPR(1:nb2) = TDMZZ(6,:)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,TMPR,nb,
     &               SM,nb,0.0D0,TMP,nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SM,nb,
     &               TMP,nb,0.0D0,TMPR,nb)
      TDMZZL(6,:) = TMPR(1:nb2)
c Do for all components of TSDMZZ
      do i=1, 6
        TMPR(1:nb2) = TSDMZZ(i,:)
        call DGEMM_('N','N',nb,nb,nb,1.0D0,TMPR,nb,
     &               SM,nb,0.0D0,TMP,nb)
        call DGEMM_('N','N',nb,nb,nb,1.0D0,SM,nb,
     &               TMP,nb,0.0D0,TMPR,nb)
        TSDMZZL(i,:) = TMPR(1:nb2)
      enddo
      call mma_deallocate(TMPR)
C End of the Lowdin Orthogonalization

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c Do SVD for transition density matrix as a whole
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c In order to use zgesvd first we combine TDMZZL(3,:)
c and TDMZZL(6,:) as a complex matrix
      Call MMA_ALLOCATE(TDMZZLC,nb2,LABEL='LTDMZZLC')
      SumofTDMZZLC = 0.0D0
      do i=1, nb2
        TDMZZLC(i) = cmplx(TDMZZL(3,i),TDMZZL(6,i),8)
        SumofTDMZZLC = SumofTDMZZLC+abs(TDMZZLC(i))
      enddo
C Do SVD by using ZGESVD, see lapack for documentation
c A = U * SIGMA * V^dagger
C Get work space for U, SIGMA, and V^dagger, VH
      Call MMA_ALLOCATE(SVDU,nb2,LABEL='SVDU')
      Call MMA_ALLOCATE(SVDS,nb,LABEL='SVDS')
      Call MMA_ALLOCATE(SVDVH,nb2,LABEL='SVDVH')
      Call MMA_ALLOCATE(SIZC,1,LABEL='RESI')
      call mma_allocate(RESIR,5*nb,Label='RESIR')
c Set LWORK=-1 to get the optimal scratch space in SIZC
c then let LWORK equal to length of scratch space
c free and reallocate memory for SIZC using that length
      LWORK=-1
      Call ZGESVD_('A','A',NB,NB,TDMZZLC,NB,SVDS,
     &            SVDU,NB,SVDVH,NB,SIZC,
     &            LWORK,RESIR,INFO)
      LWORK=max(1,int(SIZC(1)))
      Call MMA_DEALLOCATE(SIZC)
      Call MMA_ALLOCATE(RESI,LWORK,LABEL='RESI')
c Do SVD for TDMZZLC
      call ZCOPY_(nb2,[(0.0d0,0.0d0)],0,SVDU,1)
      call DCOPY_(nb,[0.0d0],0,SVDS,1)
      call ZCOPY_(nb2,[(0.0d0,0.0d0)],0,SVDVH,1)
      call ZCOPY_(LWORK,[(0.0d0,0.0d0)],0,RESI,1)
      If(SumofTDMZZLC.GE.1.0D-20) Then
        Call ZGESVD_('A','A',nb,nb,TDMZZLC,nb,SVDS,
     &            SVDU,nb,SVDVH,nb,RESI,
     &            LWORK,RESIR,INFO)
      EndIf
      If(INFO.ne.zero) write(6,*) "SVD convergence issue"
c End testing SVD
c For partitioning properties in the NTO basis
c Partition of the MLTPL 1, dipole moment intergals
      Call mma_allocate(DIPs,nb2,Label='DIPs')
      Call MMA_ALLOCATE(DIPsC,nb2,LABEL="LDIPsC")
      Call MMA_ALLOCATE(BUFF1,nb2,LABEL='BUFF1')
      Call MMA_ALLOCATE(BUFF2,nb2,LABEL='BUFF2')
      Call MMA_ALLOCATE(YMAT,3,nb2,LABEL='YMAT')
      Call MMA_ALLOCATE(SumofYdiag,3,LABEL='SumofYdiag')
c The three components of dipole
      do di=1, 3
        LABEL='MLTPL  1'
        IRC = -1
        ICMP = di
        ISYLAB = 1
        IOPT = ibset(0,sOpSiz)
        CALL IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
        IOPT = 6
        Call mma_allocate(DIP,SIZ(1),Label='DIP')
        CALL RDONE(IRC,IOPT,LABEL,ICMP,DIP,ISYLAB)
        Call DESYM_SONTO(DIP,SIZ(1),DIPs,ISYLAB)
        Call mma_deallocate(DIP)
c Perform Lowdin orthogonalization on operator matrix
c They transform as S^-1/2 P S^-1/2
        Call DGEMM_('N','N',nb,nb,nb,1.0D0,DIPs,nb,
     &              SMI,nb,0.0D0,TMP,nb)
        Call DGEMM_('N','N',nb,nb,nb,1.0D0,SMI,nb,
     &              TMP,nb,0.0D0,DIPs,nb)
c ZGESVD destroys TDMZZLC after it finishes
c reconstruct TDMZZLC and DIPsC
        do i=1, nb2
          TDMZZLC(i) = cmplx(TDMZZL(3,i),TDMZZL(6,i),8)
          DIPsC(i)=cmplx(DIPs(i),zero,8)
        enddo
c Do U^H TDMZZLC DIP U = Y, Diagonal of Y contains the partition
        Call ZGEMM_('N','N',nb,nb,nb,(1.0D0,0.0D0),TDMZZLC(:),nb,
     &              DIPsC,nb,(0.0D0,0.0D0),BUFF1(:),nb)
        Call ZGEMM_('N','N',nb,nb,nb,(1.0D0,0.0D0),BUFF1(:),nb,
     &              SVDU(:),nb,(0.0D0,0.0D0),BUFF2(:),nb)
        Call ZGEMM_('C','N',nb,nb,nb,(1.0D0,0.0D0),SVDU(:),nb,
     &              BUFF2(:),nb,(0.0D0,0.0D0),BUFF1(:),nb)
        YMAT(di,:) = BUFF1(:)
        SumofYdiag(di) = cmplx(zero,zero,8)
        do i=1, nb
          SumofYdiag(di)=SumofYdiag(di) + YMAT(di,(i-1)*nb+i)
        enddo
      enddo
      Call mma_deallocate(DIPs)
      Call MMA_DEALLOCATE(BUFF1)
      Call MMA_DEALLOCATE(BUFF2)
      Call MMA_DEALLOCATE(DIPsC)
      Call MMA_DEALLOCATE(RESI)
      Call MMA_DEALLOCATE(RESIR)

c      Call MMA_DEALLOCATE(YMAT)

c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(1)
c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(2)
c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(3)
c End of partitioning properties

c But we still need to transform U and V back to original AO
c basis using U=S^{-1/2) U' and V=S^{-1/2} V'
c for V^t it is V'^t S^{-1/2} = V^t
      call mma_allocate(SVDUR,nb2,Label='SVDUR')
      call mma_allocate(SVDUI,nb2,Label='SVDUI')
      SVDUR(:)=0.0D0
      SVDUI(:)=0.0D0
      call mma_allocate(SVDVHR,nb2,Label='SVDVHR')
      call mma_allocate(SVDVHI,nb2,Label='SVDVHR')
      SVDVHR(:)=0.0D0
      SVDVHI(:)=0.0D0

      TMP(:)=0.0D0

      SVDUR(:)= real(SVDU(:))
      SVDUI(:)=aimag(SVDU(:))
      SVDVHR(:)= real(SVDVH(:))
      SVDVHI(:)=aimag(SVDVH(:))
c U
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SMI,nb,
     &             SVDUR,nb,0.0D0,TMP,nb)
      SVDUR(:)=TMP(:)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SMI,nb,
     &             SVDUI,nb,0.0D0,TMP,nb)
      SVDUI(:)=TMP(:)
c V^H
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SVDVHR,nb,
     &             SMI,nb,0.0D0,TMP,nb)
      SVDVHR(:)=TMP(:)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,SVDVHI,nb,
     &             SMI,nb,0.0D0,TMP,nb)
      SVDVHI(:)=TMP(:)
c V
      call mma_allocate(SVDVR,nb2,Label='SVDVR')
      call mma_allocate(SVDVI,nb2,Label='SVDVI')
      SVDVR(:)=0.0D0
      SVDVI(:)=0.0D0
      do i=0, nb-1
        do j=0, nb-1
          SVDVR(1+i*nb+j)=SVDVHR(1+j*nb+i)
c imaginary part takes a negative sign
          SVDVI(1+i*nb+j)=-1.D0*SVDVHI(1+j*nb+i)
        enddo
      enddo
      call mma_deallocate(SVDVHR)
      call mma_deallocate(SVDVHI)
c tests
      CALL ADD_INFO("LAMBDA",SVDS,5,4)

c singular values
      Sumofeigen=zero
      do i=0, nb-1
        Sumofeigen=Sumofeigen+SVDS(i+1)**2
      enddo

c Head of the output
      write(6,*)
      write(6,'(6X,A)') repeat('*',90)
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,29X,A31,28X,A)')
     & '*','Natural transition orbitals','*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,27X,A25,I2,A12,I2,20X,A)')
     &'*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A)') repeat('*',90)
      write(6,*)
c Start output singular value information for positive spin values
      write(6,'(6X,A)') repeat('=',90)
      eigen_print_limit=1.0D-8
      write(6,'(5X,A12,A12,A16,A51)')'EXCITATION','EIGENVALUE',
     &'EXCITATION',
     &'TRANSITION DIPOLE MOMENT'
      write(6,'(5X,A12,12X,A16,3A17)')'AMPLITUDE',
     &'CONTRIBUTION(%)',
     &'(1)','(2)','(3)'
      write(6,'(6X,A)') repeat('-',90)
      do i=0,nb-1
        IF(SVDS(i+1)**2.lt.eigen_print_limit)  EXIT
        write(6,'(4X,3X,F8.5,4X,F8.5,8X,F8.2,2X,
     &           3(F9.4,SP,F7.4,"i",SS))')
     &  SVDS(i+1),SVDS(i+1)**2,
     &  SVDS(i+1)**2/Sumofeigen*1.0D2,
     &  YMAT(1,i*nb+i+1),
     &  YMAT(2,i*nb+i+1),
     &  YMAT(3,i*nb+i+1)
      enddo
      write(6,'(6X,A,F8.5)')'SUM OF EIGENVALUES ',Sumofeigen
      write(6,'(6X,A24,15X,3(F9.4,SP,F7.4,"i",SS))')
     &                       'SUM OF TRANSITION DIPOLE',
     &                       SumofYdiag(1),
     &                       SumofYdiag(2),
     &                       SumofYdiag(3)
      write(6,'(6X,A)') repeat('=',90)
      write(6,*)
      write(6,*)
c Write NTOs to file in C1 symmetry
      Do i=1, nb
        SVDS(i)=SVDS(i)**2/Sumofeigen
      EndDo
      LU=50
      LU=ISFREEUNIT(LU)
      Note='*  Spin-orbit Natural Transition Orbitals'
c U real
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','PART','.','Re'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDUR,
     &           SVDS ,Dummy,iDummy,Note)
c U imaginary
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','PART','.','Im'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDUI,
     &           SVDS ,Dummy,iDummy,Note)
c V real
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Re'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDVR,
     &           SVDS ,Dummy,iDummy,Note)
c V imaginary
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Im'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],SVDVI,
     &           SVDS ,Dummy,iDummy,Note)
c End of output
      write(6,*)
      write(6,'(6X,A)') repeat('*',90)
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,28X,A34,25X,A)')
     & '*','End of natural transition orbitals','*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,27X,A25,I2,A12,I2,20X,A)')
     &'*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A)') repeat('*',90)
      write(6,*)

c Free up workspace
      Call MMA_DEALLOCATE(SVDU)
      Call MMA_DEALLOCATE(SVDS)
      Call MMA_DEALLOCATE(SVDVH)
      Call MMA_DEALLOCATE(TDMZZL)
      Call MMA_DEALLOCATE(TSDMZZL)
      Call MMA_DEALLOCATE(TDMZZLC)
      Call MMA_DEALLOCATE(YMAT)
      Call MMA_DEALLOCATE(SumofYdiag)
      Call MMA_DEALLOCATE(SM)
      Call MMA_DEALLOCATE(SMI)
      Call MMA_DEALLOCATE(TMP)
      Call MMA_DEALLOCATE(SVDUR)
      Call MMA_DEALLOCATE(SVDUI)
      Call MMA_DEALLOCATE(SVDVR)
      Call MMA_DEALLOCATE(SVDVI)

      END SUBROUTINE DO_AOTDMNTO
