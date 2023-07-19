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
      IMPLICIT None
#include "prgm.fh"
      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='DO_AOTDMNTO')
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
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
      Integer i,j,info,lwork,di,icmp,iopt,irc,isylab,LDIP,LDIPs,LEIG
      Integer LEIGM,LP,LRESI,LRESIR,LSM,LSMI,LSZZ,LSZZs,LTMP,LSVDUR
      Integer LSVDUI,LSVDVHR,LSVDVHI,LSVDVR,LSVDVI
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
      Integer LTMPR,LTMPI
c end
      CHARACTER*8 LABEL
      CHARACTER(LEN=7) STATENAME,STATENAMETMP
      CHARACTER(LEN=128) FNAME
      CHARACTER(LEN=72) NOTE
c      Logical TestPrint
c Test variables
      Integer LBUFF1
c End test variables
c      TestPrint=.TRUE.
c      TestPrint=.FALSE.
      Zero=0.0D0
      Two=2.0D0
      pi=DACOS(-1.0D0)

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
      Call DSCAL_(nb2,-1.0D0,TDMZZ(4,:),1)
      Call DSCAL_(nb2,-1.0D0,TDMZZ(5,:),1)
      Call DSCAL_(nb2,-1.0D0,TDMZZ(6,:),1)
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
        Call GETMEM('MSq','ALLO','REAL',LDIPs,nb2)
        LABEL(1:8)='MLTPL  1'
        IRC = -1
        ICMP = di
        ISYLAB = 1
        IOPT = ibset(0,sOpSiz) ! Only read the size of the array
        CALL IRDONE(IRC,IOPT,LABEL,ICMP,SIZ,ISYLAB)
        !no nuclear contrib, no origin of operator
        IOPT = ibset(ibset(0,sNoOri),sNoNuc)
        Call GETMEM('MLTPL  1','ALLO','REAL',LDIP,SIZ(1))
        CALL RDONE(IRC,IOPT,LABEL,ICMP,WORK(LDIP),ISYLAB)
        Call DESYM_SONTO(WORK(LDIP),SIZ(1),WORK(LDIPs),ISYLAB)
        Call GETMEM('MLTPL  1','FREE','REAL',LDIP,SIZ(1))
        write(6,*) '  For istate ', ISTATE, ' jstate ',JSTATE
        write(6,*) '  Component ',ICMP
c Get complex matrices
        Call MMA_ALLOCATE(DIPsC,nb2,LABEL="LDIPsC")
        do i=1,nb2
          TDMZZC(i)=cmplx(TDMZZ(di,i),TDMZZ(di+3,i),8)
          DIPsC(i)=cmplx(WORK(LDIPs+i-1),zero,8)
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
        Call GETMEM('MSq','FREE','REAL',LDIPs,nb2)
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
        call GETMEM('TMPR  ','ALLO','REAL',LTMPR,nb2)
        call GETMEM('TMPI  ','ALLO','REAL',LTMPI,nb2)
        do i=1,nb2
          WORK(LTMPR+i-1)=TDMZZ(3,i)*cos(phi)-TDMZZ(6,i)*sin(phi)
          WORK(LTMPI+i-1)=TDMZZ(6,i)*cos(phi)+TDMZZ(3,i)*sin(phi)
        enddo
        do i=1,3
          Call DCOPY_(nb2,WORK(LTMPR),1,TDMZZ(i,:),1)
          Call DCOPY_(nb2,WORK(LTMPI),1,TDMZZ(i+3,:),1)
        enddo
        call GETMEM('TMPI  ','FREE','REAL',LTMPI,nb2)
        call GETMEM('TMPR  ','FREE','REAL',LTMPR,nb2)
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
c LSZZ  - AO Overlap integral
c LSZZs - AO Overlap integral in square
c LEIG  - AO Overlap eigenvalues
      call GETMEM('SZZ   ','ALLO','REAL',LSZZ,NBTRI)
      call GETMEM('SZZs  ','ALLO','REAL',LSZZs,nb2)
      call GETMEM('EIG   ','ALLO','REAL',LEIG,nb)
      call DCOPY_(NBTRI,[0.0D00],0,WORK(LSZZ),1)
      call DCOPY_(nb2,[0.0D00],0,WORK(LSZZs),1)
      call DCOPY_(nb,[0.0D00],0,WORK(LEIG),1)
C AO OVERLAP MATRIX
      IRC=-1
c IOPT=6, origin and nuclear contrib not read
      IOPT=ibset(ibset(0,sNoOri),sNoNuc)
      ICMP=1
      ISYLAB=1
      LABEL='MLTPL  0'
      call RDONE(IRC,IOPT,LABEL,ICMP,WORK(LSZZ),ISYLAB)
      IF (IRC.NE.0) THEN
        WRITE(6,*)
        WRITE(6,*)'      *** ERROR IN SUBROUTINE  SONATORB ***'
        WRITE(6,*)'      OVERLAP INTEGRALS ARE NOT AVAILABLE'
        WRITE(6,*)
        CALL ABEND()
      ENDIF
      Call DESYM_SONTO(WORK(LSZZ),NBTRI,WORK(LSZZs),ISYLAB)
      call GETMEM('SZZ   ','FREE','REAL',LSZZ,NBTRI)
**************
* For tests
**************
      Call GETMEM('BUFF1','ALLO','REAL',LBUFF1,nb2)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSZZs),nb,
     &              TDMZZ(3,:),nb,0.0D0,WORK(LBUFF1),nb)
C Trace the resulting matrix
      NumOfEc = Zero
      do i=0, nb-1
        do j=0, nb-1
          if(i.eq.j) then
            NumOfEc = NumOfEc +
     &      WORK(LBUFF1+i*nb+j)
          endif
        enddo
      enddo
c      write(6,*) 'NumOfEc ',NumOfEc
      Call GETMEM('BUFF1','FREE','REAL',LBUFF1,nb2)
**************

c DIAGONALIZE AO OVERLAP MATRIX
c Set LWORK=-1 to get the optimal scratch space WORK(LRESI)
c then let LWORK equal to length of scratch space
c free and reallocate memory for LRESI using that length
      call GETMEM('RESI  ','ALLO','REAL',LRESI,1)
      LWORK=-1
      call DSYEV_('V','U',nb,WORK(LSZZs),nb,WORK(LEIG),
     &            WORK(LRESI),LWORK,INFO)
      LWORK=INT(WORK(LRESI))
      call GETMEM('RESI  ','FREE','REAL',LRESI,1)
      call GETMEM('RESI  ','ALLO','REAL',LRESI,LWORK)
c WORK(LSZZs) in as the AO overlap sqaure matrix
c out as the eigenvector matrix of WORK(LSZZs)
c with eigenvalues in WORK(LEIG)
      call DSYEV_('V','U',nb,WORK(LSZZs),nb,WORK(LEIG),
     &            WORK(LRESI),LWORK,INFO)
c Put WORK(LEIG) in sqrt and in diagonal in WORK(LEIGM)
      call GETMEM('EIGM  ','ALLO','REAL',LEIGM,nb2)
      do i=0,nb-1
        do j=0,nb-1
          If (i.eq.j) then
            WORK(LEIGM+i*nb+j)=sqrt(WORK(LEIG+i))
          Else
            WORK(LEIGM+i*nb+j)=zero
          Endif
        enddo
      enddo
      call GETMEM('EIG   ','FREE','REAL',LEIG,nb)
      call GETMEM('RESI  ','FREE','REAL',LRESI,LWORK)
c Get S^1/2 from S^1/2 = U S_diag^1/2 U^T
      call GETMEM('SsqrtM','ALLO','REAL',LSM,nb2)
      call GETMEM('TMP   ','ALLO','REAL',LTMP,nb2)
      call DGEMM_('N','T',nb,nb,nb,1.0D0,WORK(LEIGM),nb,
     &             WORK(LSZZs),nb,0.0D0,WORK(LTMP),nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSZZs),nb,
     &             WORK(LTMP),nb,0.0D0,WORK(LSM),nb)
      call GETMEM('TMP   ','FREE','REAL',LTMP,nb2)
      call GETMEM('SZZs  ','FREE','REAL',LSZZs,nb2)
      call GETMEM('EIGM  ','FREE','REAL',LEIGM,nb2)
c Get inverse of S^1/2 -> S^-1/2
c Before calling DGETRI, call DGETRF to factorize WORK(LSM)
c Set LWORK=-1 to get the optimal scratch space WORK(LRESI)
c then let LWORK equal to length of scratch space
c free and reallocate memory for LRESI using that length
      call GETMEM('SsqrtMI','ALLO','REAL',LSMI,nb2)
      do i=0,nb2-1
        WORK(LSMI+i)=WORK(LSM+i)
      enddo
      call GETMEM('PIV   ','ALLO','INTE',LP,nb)
      call DGETRF_(nb,nb,WORK(LSMI),nb,IWORK(LP),INFO)
      call GETMEM('RESI  ','ALLO','REAL',LRESI,1)
      LWORK=-1
      call DGETRI_(nb,WORK(LSMI),nb,IWORK(LP),WORK(LRESI),LWORK,INFO)
      LWORK=INT(WORK(LRESI))
      call GETMEM('RESI  ','FREE','REAL',LRESI,1)
      call GETMEM('RESI  ','ALLO','REAL',LRESI,LWORK)
      call DGETRI_(nb,WORK(LSMI),nb,IWORK(LP),WORK(LRESI),LWORK,INFO)
      call GETMEM('PIV   ','FREE','INTE',LP,nb)
      call GETMEM('RESI  ','FREE','REAL',LRESI,LWORK)

c Note: The density matrix should transform as S^1/2 D S^1/2
c Transform TDMZZ and TSDMZZ as S^1/2 T S^1/2
      Call GETMEM('TMP   ','ALLO','REAL',LTMP,nb2)
      Call MMA_ALLOCATE(TDMZZL,6,nb2,LABEL='LTDMZZL')
      Call MMA_ALLOCATE(TSDMZZL,6,nb2,LABEL='LTSDMZZL')
c Real part of TDMZZ
      call DGEMM_('N','N',nb,nb,nb,1.0D0,TDMZZ(3,:),nb,
     &               WORK(LSM),nb,0.0D0,WORK(LTMP),nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSM),nb,
     &               WORK(LTMP),nb,0.0D0,TDMZZL(3,:),nb)
c Imaginary part of TDMZZ
      call DGEMM_('N','N',nb,nb,nb,1.0D0,TDMZZ(6,:),nb,
     &               WORK(LSM),nb,0.0D0,WORK(LTMP),nb)
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSM),nb,
     &               WORK(LTMP),nb,0.0D0,TDMZZL(6,:),nb)
c Do for all components of TSDMZZ
      do i=1, 6
        call DGEMM_('N','N',nb,nb,nb,1.0D0,TSDMZZ(i,:),nb,
     &               WORK(LSM),nb,0.0D0,WORK(LTMP),nb)
        call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSM),nb,
     &               WORK(LTMP),nb,0.0D0,TSDMZZL(i,:),nb)
      enddo
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
c      call GETMEM('SVDRESI','ALLO','REAL',LRESI,1)
      call GETMEM('SVDRESIR','ALLO','REAL',LRESIR,5*nb)
c Set LWORK=-1 to get the optimal scratch space in SIZC
c then let LWORK equal to length of scratch space
c free and reallocate memory for SIZC using that length
      LWORK=-1
      Call ZGESVD_('A','A',NB,NB,TDMZZLC,NB,SVDS,
c     &            SVDU,NB,SVDVH,NB,WORK(LRESI),
     &            SVDU,NB,SVDVH,NB,SIZC,
     &            LWORK,WORK(LRESIR),INFO)
c      LWORK=max(1,int(WORK(LRESI)))
      LWORK=max(1,int(SIZC(1)))
      Call MMA_DEALLOCATE(SIZC)
c      call GETMEM('SVDRESI','FREE','INTE',LRESI,1)
      Call MMA_ALLOCATE(RESI,LWORK,LABEL='RESI')
c Do SVD for TDMZZLC
      call ZCOPY_(nb2,[(0.0d0,0.0d0)],0,SVDU,1)
      call DCOPY_(nb,[0.0d0],0,SVDS,1)
      call ZCOPY_(nb2,[(0.0d0,0.0d0)],0,SVDVH,1)
      call ZCOPY_(LWORK,[(0.0d0,0.0d0)],0,RESI,1)
      If(SumofTDMZZLC.GE.1.0D-20) Then
        Call ZGESVD_('A','A',nb,nb,TDMZZLC,nb,SVDS,
     &            SVDU,nb,SVDVH,nb,RESI,
     &            LWORK,WORK(LRESIR),INFO)
      EndIf
      If(INFO.ne.zero) write(6,*) "SVD convergence issue"
c End testing SVD
c For partitioning properties in the NTO basis
c Partition of the MLTPL 1, dipole moment intergals
      Call GETMEM('MSq','ALLO','REAL',LDIPs,nb2)
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
        Call GETMEM('MLTPL  1','ALLO','REAL',LDIP,SIZ(1))
        CALL RDONE(IRC,IOPT,LABEL,ICMP,WORK(LDIP),ISYLAB)
        Call DESYM_SONTO(WORK(LDIP),SIZ(1),WORK(LDIPs),ISYLAB)
        Call GETMEM('MLTPL  1','FREE','REAL',LDIP,SIZ(1))
c Perform Lowdin orthogonalization on operator matrix
c They transform as S^-1/2 P S^-1/2
        Call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LDIPs),nb,
     &              WORK(LSMI),nb,0.0D0,WORK(LTMP),nb)
        Call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSMI),nb,
     &              WORK(LTMP),nb,0.0D0,WORK(LDIPs),nb)
c ZGESVD destroys TDMZZLC after it finishes
c reconstruct TDMZZLC and DIPsC
        do i=1, nb2
          TDMZZLC(i) = cmplx(TDMZZL(3,i),TDMZZL(6,i),8)
          DIPsC(i)=cmplx(WORK(LDIPs+i-1),zero,8)
        enddo
c Do U^H TDMZZLC DIP U = Y, Diagonal of Y contains the partition
        Call ZGEMM_('N','N',nb,nb,nb,(1.0D0,0.0D0),TDMZZLC(:),nb,
     &              DIPsC,nb,(0.0D0,0.0D0),BUFF1(:),nb)
        Call ZGEMM_('N','N',nb,nb,nb,(1.0D0,0.0D0),BUFF1(:),nb,
     &              SVDU(:),nb,(0.0D0,0.0D0),BUFF2(:),nb)
        Call ZGEMM_('C','N',nb,nb,nb,(1.0D0,0.0D0),SVDU(:),nb,
     &              BUFF2(:),nb,(0.0D0,0.0D0),YMAT(di,:),nb)
        SumofYdiag(di) = cmplx(zero,zero,8)
        do i=1, nb
          SumofYdiag(di)=SumofYdiag(di) + YMAT(di,(i-1)*nb+i)
        enddo
      enddo
      Call GETMEM('MSq','FREE','REAL',LDIPs,nb2)
      Call MMA_DEALLOCATE(BUFF1)
      Call MMA_DEALLOCATE(BUFF2)
      Call MMA_DEALLOCATE(DIPsC)
      Call MMA_DEALLOCATE(RESI)
      call GETMEM('SVDRESIR','FREE','REAL',LRESIR,5*nb)

c      Call MMA_DEALLOCATE(YMAT)

c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(1)
c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(2)
c      write(6,'(F11.5,SP,F8.5,"i")') SumofYdiag(3)
c End of partitioning properties

c But we still need to transform U and V back to original AO
c basis using U=S^{-1/2) U' and V=S^{-1/2} V'
c for V^t it is V'^t S^{-1/2} = V^t
      call GETMEM('LSVDUR','ALLO','REAL',LSVDUR,nb2)
      call GETMEM('LSVDUI','ALLO','REAL',LSVDUI,nb2)
      call GETMEM('LSVDVHR','ALLO','REAL',LSVDVHR,nb2)
      call GETMEM('LSVDVHI','ALLO','REAL',LSVDVHI,nb2)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDUR),1)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDUI),1)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDVHR),1)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDVHI),1)
      call DCOPY_(nb2,[0.0D0],0,WORK(LTMP),1)
      do i=0, nb2-1
        WORK(LSVDUR+i)=real(SVDU(i+1))
        WORK(LSVDUI+i)=aimag(SVDU(i+1))
        WORK(LSVDVHR+i)=real(SVDVH(i+1))
        WORK(LSVDVHI+i)=aimag(SVDVH(i+1))
      enddo
c U
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSMI),nb,
     &             WORK(LSVDUR),nb,0.0D0,WORK(LTMP),nb)
      do i=0, nb2-1
        WORK(LSVDUR+i)=WORK(LTMP+i)
      enddo
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSMI),nb,
     &             WORK(LSVDUI),nb,0.0D0,WORK(LTMP),nb)
      do i=0, nb2-1
        WORK(LSVDUI+i)=WORK(LTMP+i)
      enddo
c V^H
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSVDVHR),nb,
     &             WORK(LSMI),nb,0.0D0,WORK(LTMP),nb)
      do i=0, nb2-1
        WORK(LSVDVHR+i)=WORK(LTMP+i)
      enddo
      call DGEMM_('N','N',nb,nb,nb,1.0D0,WORK(LSVDVHI),nb,
     &             WORK(LSMI),nb,0.0D0,WORK(LTMP),nb)
      do i=0, nb2-1
        WORK(LSVDVHI+i)=WORK(LTMP+i)
      enddo
c V
      call GETMEM('LSVDVR','ALLO','REAL',LSVDVR,nb2)
      call GETMEM('LSVDVI','ALLO','REAL',LSVDVI,nb2)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDVR),1)
      call DCOPY_(nb2,[0.0D0],0,WORK(LSVDVI),1)
      do i=0, nb-1
        do j=0, nb-1
          WORK(LSVDVR+i*nb+j)=WORK(LSVDVHR+j*nb+i)
c imaginary part takes a negative sign
          WORK(LSVDVI+i*nb+j)=-1.D0*WORK(LSVDVHI+j*nb+i)
        enddo
      enddo
      call GETMEM('LSVDVHR','FREE','REAL',LSVDVHR,nb2)
      call GETMEM('LSVDVHI','FREE','REAL',LSVDVHI,nb2)
c tests
      CALL ADD_INFO("LAMBDA",SVDS,5,4)

c singular values
      Sumofeigen=zero
      do i=0, nb-1
        Sumofeigen=Sumofeigen+SVDS(i+1)**2
      enddo

c Head of the output
      write(6,*)
      write(6,'(6X,90A1)') ('*',i=1,90)
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,29X,A31,28X,A)')
     & '*','Natural transition orbitals','*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,27X,A25,I2,A12,I2,20X,A)')
     &'*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,90A1)') ('*',i=1,90)
      write(6,*)
c Start output singular value information for positive spin values
      write(6,'(6X,90A1)') ('=',i=1,90)
      eigen_print_limit=1.0D-8
      write(6,'(5X,A12,A12,A16,A51)')'EXCITATION','EIGENVALUE',
     &'EXCITATION',
     &'TRANSITION DIPOLE MOMENT'
      write(6,'(5X,A12,12X,A16,3A17)')'AMPLITUDE',
     &'CONTRIBUTION(%)',
     &'(1)','(2)','(3)'
      write(6,'(6X,90A1)') ('-',i=1,90)
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
      write(6,'(6X,90A1)') ('=',i=1,90)
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
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],WORK(LSVDUR),
     &           SVDS ,Dummy,iDummy,Note)
c U imaginary
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','PART','.','Im'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],WORK(LSVDUI),
     &           SVDS ,Dummy,iDummy,Note)
c V real
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Re'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],WORK(LSVDVR),
     &           SVDS ,Dummy,iDummy,Note)
c V imaginary
      write(FNAME,'(6(a))')
     &      'NTORB.SO.',trim(adjustl(STATENAME)),'.','HOLE','.','Im'
      write(6,'(4(a))')
     & '      NATURAL TRANSITION ORBITALS FOR SPIN-ORBIT STATE ',
     & trim(STATENAME),
     & ' ARE WRITTEN ONTO FILE ',
     & FNAME
      call WRVEC(FNAME,LU,'CO',1,[NB],[NB],WORK(LSVDVI),
     &           SVDS ,Dummy,iDummy,Note)
c End of output
      write(6,*)
      write(6,'(6X,90A1)') ('*',i=1,90)
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,28X,A34,25X,A)')
     & '*','End of natural transition orbitals','*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,A,27X,A25,I2,A12,I2,20X,A)')
     &'*','Between spin-orbit state ',ISTATE,' and state ',JSTATE,'*'
      write(6,'(6X,A,88X,A)') '*','*'
      write(6,'(6X,90A1)') ('*',i=1,90)
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
      call GETMEM('SsqrtM','FREE','REAL',LSM,nb2)
      call GETMEM('SsqrtMI','FREE','REAL',LSMI,nb2)
      Call GETMEM('TMP   ','FREE','REAL',LTMP,nb2)
      call GETMEM('LSVDUR','FREE','REAL',LSVDUR,nb2)
      call GETMEM('LSVDUI','FREE','REAL',LSVDUI,nb2)
      call GETMEM('LSVDVR','FREE','REAL',LSVDVR,nb2)
      call GETMEM('LSVDVI','FREE','REAL',LSVDVI,nb2)

      END
