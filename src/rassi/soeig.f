************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      SUBROUTINE SOEIG(PROP,USOR,USOI,ENSOR,NSS,ENERGY)
      !> module dependencies
      use rassi_aux, only: ipglob
      use rassi_global_arrays, only: JBNUM
      use sorting, only : argsort
      use sorting_funcs, only : leq_r
#ifdef _HDF5_
      use Dens2HDF5
      use mh5, only: mh5_put_dset
      use RASSIWfn, only: wfn_SOS_CoefI, wfn_SOS_CoefR, wfn_SOS_Energy,
     &                    wfn_SOS_HSOI, wfn_SOS_HSOR,
     &                    wfn_SOS_VSOI, wfn_SOS_VSOR
      use Cntrl, only: IFTDM, IFTRD1, RHODYN
#endif
#ifdef _DMRG_
      use rasscf_global, only: doDMRG
      use qcmaquis_interface_cfg
#endif
      use Constants, only: auTocm, auToeV
      use stdalloc, only: mma_allocate, mma_deallocate
      use Cntrl, only: NSTATE, NPROP, NSOThr_PRT,
     &                 SOThr_PRT, EMIN, IFJ2, IFJZ, REDUCELOOP,
     &                 LOOPDIVIDE, ICOMP, MLTPLT, PNAME

      IMPLICIT NONE
#include "rassi.fh"

      INTEGER NSS
      REAL*8 PROP(NSTATE,NSTATE,NPROP)
      REAL*8 USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
      REAL*8 ENERGY(NSTATE)

      INTEGER I,N
      INTEGER ITOL,IDX
      INTEGER JOB
      INTEGER IPROP
      INTEGER IAMFIX,IAMFIY,IAMFIZ,IAMX,IAMY,IAMZ
      INTEGER ISS,JSS,ISTATE,JSTATE
      INTEGER MAGN
      INTEGER MPLET,MPLET1,MPLET2,MSPROJ,MSPROJ1,MSPROJ2

      REAL*8 AMFIX,AMFIY,AMFIZ
      REAL*8 CG0,CGM,CGP,CGX,CGY
      REAL*8 E,E0,E1,E2,E3,E_TMP,FACT,FRAC,EI,EPSH,EPSS,ERMS,V2SUM
      REAL*8 HSOI,HSOR,HSOTOT
      REAL*8 OMEGA
      REAL*8 S1,S2,SM1,SM2
      REAL*8 SOTHR_MIN
      REAL*8 X,X_THR,XJEFF
      REAL*8, ALLOCATABLE :: ESO(:), HAMSOR(:,:), HAMSOI(:,:)
#ifdef _DMRG_
      complex*16, allocatable :: hso_tmp(:,:)
      complex*16, allocatable :: ccwork(:)
      real*8    , allocatable :: rwork(:)
      integer                 :: lcwork, info
#endif
      Integer, Allocatable :: IndexE(:)

      REAL*8, EXTERNAL :: DCLEBS

      Logical lOMG, lJ2
      Integer, External :: cho_x_gettol
      LOGICAL :: debug_dmrg_rassi_code = .false.
      Integer, allocatable:: MAPST(:), MAPSP(:), MAPMS(:)
      Real*8, allocatable:: HTOTR(:,:), HTOTI(:,:)
      Real*8, allocatable:: LXI(:), LYI(:), LZI(:)
      Real*8, allocatable:: JXR(:), JYR(:), JZR(:)
      Real*8, allocatable:: JXI(:), JYI(:), JZI(:)
      Real*8, allocatable:: OMGR(:,:), OMGI(:,:)
      Real*8, allocatable:: J2R(:,:), J2I(:,:)


C CONSTANTS:
      lOMG=.False.
      lJ2 =.False.

C Identify AMFI and ANGMOM matrix elements:
      IAMFIX=0
      IAMFIY=0
      IAMFIZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:4).EQ.'AMFI') THEN
         IF(ICOMP(IPROP).EQ.1) IAMFIX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMFIY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMFIZ=IPROP
       ELSE IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO

C Mapping from spin states to spin-free state and to spin:
      CALL mma_allocate(MAPST,NSS,Label='MAPST')
      CALL mma_allocate(MAPSP,NSS,Label='MAPSP')
      CALL mma_allocate(MAPMS,NSS,Label='MAPMS')
      ISS=0
      DO ISTATE=1,NSTATE
       JOB=JBNUM(ISTATE)
       MPLET=MLTPLT(JOB)
       DO MSPROJ=-MPLET+1,MPLET-1,2
        ISS=ISS+1
        MAPST(ISS)=ISTATE
        MAPSP(ISS)=MPLET
        MAPMS(ISS)=MSPROJ
       END DO
      END DO
C Complex hamiltonian matrix elements over spin states:
      CALL mma_allocate(HTOTR,NSS,NSS,Label='HTOTR')
      CALL mma_allocate(HTOTI,NSS,NSS,Label='HTOTI')
      HTOTR(:,:)=0.0D0
      HTOTI(:,:)=0.0D0

      IF(IPGLOB.GE.1) THEN
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A,34X,A)')
     &      '*','       Spin-orbit section     ','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A)') repeat('*',100)
       WRITE(6,*)
      ENDIF

      if(debug_dmrg_rassi_code)then
        write(6,*) 'BLUBB BLUBB debug print of property matrix'
        do istate = 1, nstate
        do jstate = 1, nstate
        DO IPROP=1,NPROP
          if(abs(prop(istate,jstate,iprop)) > 1.0d-14)
     &    write(6,*) 'prop(',istate,',',jstate,',',iprop,') = ',
     &                prop(istate,jstate,iprop)
        end do
        end do
        end do
      end if

      DO ISS=1,NSS
        ISTATE=MAPST(ISS)
        MPLET1=MAPSP(ISS)
        MSPROJ1=MAPMS(ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=MAPST(JSS)
          MPLET2=MAPSP(JSS)
          MSPROJ2=MAPMS(JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
          AMFIX=0.0D0
          IF(IAMFIX.NE.0) AMFIX=PROP(ISTATE,JSTATE,IAMFIX)
          AMFIY=0.0D0
          IF(IAMFIY.NE.0) AMFIY=PROP(ISTATE,JSTATE,IAMFIY)
          AMFIZ=0.0D0
          IF(IAMFIZ.NE.0) AMFIZ=PROP(ISTATE,JSTATE,IAMFIZ)
* PAM07          HSCAL=0.0D0
* PAM07          IF(ISS.EQ.JSS) HSCAL=ENERGY(ISTATE)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX= SQRT(0.5D0)*(CGM-CGP)
          CGY=-SQRT(0.5D0)*(CGM+CGP)
C SPIN-ORBIT HAMILTONIAN MATRIX ELEMENTS:
C  according to expressions between eqs. (5) and (6)
C  in Malmqvist et all CPL 357 (2002) 230-240, but real and
C  imaginary parts are swapped here!!! more precisely the Hamiltonian
C  is multiplied by imaginary unit to keep its hermicity
          HSOR=CGY*AMFIY
          HSOI=CGX*AMFIX+CG0*AMFIZ
* PAM07: Delay addition of diagonal scalar part until later, see below:
          HTOTR(ISS,JSS)=HSOR
          HTOTI(ISS,JSS)=HSOI

        END DO
      END DO

* VKochetov 2021 put SOC matrix elements to hdf5:
#ifdef _HDF5_
      if (rhodyn) then
        call mh5_put_dset(wfn_sos_vsor, HTOTR)
        call mh5_put_dset(wfn_sos_vsoi, HTOTI)
      endif
#endif

* Perhaps write out large spin-orbit coupling elements:
      IF(NSOTHR_PRT.GT.0) THEN
* Prevent infinite loop below:
       SOTHR_MIN=MAX(SOTHR_PRT,1.0D-6)
* And work with quantities larger than 1:
       X_THR=SOTHR_PRT/SOTHR_MIN
  13   CONTINUE
       N=0
       DO ISS=1,NSS
        DO JSS=1,ISS
         HSOR=HTOTR(ISS,JSS)
         HSOI=HTOTI(ISS,JSS)
         HSOTOT=SQRT(HSOR**2+HSOI**2)
         IF(HSOTOT*auTocm.GT.SOTHR_PRT) N=N+1
        END DO
       END DO
       IF(N.LE.NSOTHR_PRT) GOTO 17
       X_THR=X_THR*1.2D0
       MAGN=INT(LOG10(X_THR))
       X=X_THR/DBLE(10**MAGN)
       X_THR=0.1D0*DBLE(NINT(10.0D0*X))*DBLE(10**(MAGN))
       SOTHR_PRT=X_THR*SOTHR_MIN
       GOTO 13

  17   CONTINUE

       IF (N.GT.0) THEN
       WRITE(6,*)
       WRITE(6,*)'Complex SO-Hamiltonian matrix elements over'
       WRITE(6,*)'spin components of spin-free eigenstates (SFS):'
       WRITE(6,'(1x,A,F10.3,A)')
     &              '(In cm-1. Print threshold: ',SOTHR_PRT,' cm-1)'
       WRITE(6,'(1X,10A7)')('-------',I=1,10)
       WRITE(6,*)
         WRITE(6,'(A)')'  I1  S1  MS1    I2  S2  MS2 '//
     &            '   Real part    Imag part      Absolute'
       DO  ISS=1,NSS
        ISTATE=MAPST(ISS)
        MPLET1=MAPSP(ISS)
        MSPROJ1=MAPMS(ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,ISS
         HSOR=HTOTR(ISS,JSS)
         HSOI=HTOTI(ISS,JSS)
         HSOTOT=SQRT(HSOR**2+HSOI**2)
         IF(HSOTOT*auTocm.GE.SOTHR_PRT) THEN
          JSTATE=MAPST(JSS)
          MPLET2=MAPSP(JSS)
          MSPROJ2=MAPMS(JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
         WRITE(6,'(1X,I5,F5.1,F5.1,I5,F5.1,F5.1,3F14.3)') ISS,S1,SM1,
     &         JSS,S2,SM2,HSOR*auTocm,HSOI*auTocm,HSOTOT*auTocm
         END IF
        END DO
       END DO
       WRITE(6,'(1X,10A7)')('-------',I=1,10)

       END IF
      ENDIF

* PAM07: Addition of scalar diagonal part was delayed until here, see above.
      DO ISS=1,NSS
        ISTATE=MAPST(ISS)
        HSOR=HTOTR(ISS,ISS)
        HTOTR(ISS,ISS)=HSOR+ENERGY(ISTATE)
      END DO


      IF(IPGLOB.GE.3) THEN
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)'Complex Hamiltonian matrix including SO-coupling'
       WRITE(6,*)'over spin components of spin-free eigenstates (SFS):'
       WRITE(6,'(1X,11A7)')('-------',I=1,11)
       CALL PRCHAM(NSS,HTOTR,HTOTI)
       WRITE(6,'(1X,11A7)')('-------',I=1,11)
      ENDIF
      ! save the Hamiltonian
      call mma_allocate(HAMSOR,NSS,NSS,'HAMSOR')
      call mma_allocate(HAMSOI,NSS,NSS,'HAMSOI')
      call dcopy_(NSS*NSS,[0.d0],0,HAMSOR,1)
      call dcopy_(NSS*NSS,[0.d0],0,HAMSOI,1)
      call dcopy_(NSS*NSS,HTOTR,1,HAMSOR,1)
      call dcopy_(NSS*NSS,HTOTI,1,HAMSOI,1)
      call put_darray('HAMSOR_SINGLE',HAMSOR,NSS*NSS)
      call put_darray('HAMSOI_SINGLE',HAMSOI,NSS*NSS)
#ifdef _HDF5_
      call mh5_put_dset(wfn_sos_hsor,HAMSOR)
      call mh5_put_dset(wfn_sos_hsoi,HAMSOI)
#endif

      !> use complex matrix diagonalization
#ifdef _DMRG_
      if(doDMRG)then
        call mma_allocate(hso_tmp,nss,nss)
        call mma_allocate(ccwork,(2*nss-1))
        call mma_allocate(rwork,(3*nss-2))
        hso_tmp = 0; ccwork = 0; rwork = 0

        DO jss = 1, nss
          DO iss = 1, nss
            hso_tmp(iss,jss) = cmplx(HTOTR(ISS,JSS),
     &                               HTOTI(ISS,JSS),
     &                               kind=8)
!         write(6,*) ' hso_tmp(',iss,',',jss,') = ',hso_tmp(iss,jss)
          END DO
        END DO

        lcwork = (2*nss-1); info = 0
        call zheev_('V','U',nss,hso_tmp,nss,ensor,ccwork,lcwork,
     &             rwork,info)

        if(info /= 0)then
          write(6,*) '* WARNING in rassi/soeig.f *'
          write(6,*) 'zheev did return with an error message,info=',info
        else
          write(6,*) 'zheev in rassi/soeig.f succeeded!'
        end if

!       write(6,*) 'eigenvalues of zheev, info=',info
!       do iss = 1, nss
!         write(6,*) 'ensor(',iss,') =',ensor(iss)
!       end do

        !> save eigenvectors
        DO jss = 1, nss
          DO iss = 1, nss
            usor(iss,jss) = dble(hso_tmp(iss,jss))
            usoi(iss,jss) = aimag(hso_tmp(iss,jss))
          END DO
        END DO

        !> sort eigenvalues in increasing sequence (using the same algorithm as zjac)
        call zorder(nss,nss,usor,usoi,ensor,0)

        CALL MMA_DEALLOCATE(hso_tmp)
        CALL MMA_DEALLOCATE(ccwork)
        CALL MMA_DEALLOCATE(rwork)
      else
#endif
        !> diagonalize H_SO and get array of eigenvalues/eigenvectors
        CALL ZJAC(NSS,HTOTR,HTOTI,NSS,USOR,USOI)
        DO ISS=1,NSS
         ENSOR(ISS)=HTOTR(ISS,ISS)
        END DO
#ifdef _DMRG_
      end if
#endif

!     write(6,*) 'eigenvectors of zheev/zjac (real)'
!     do iss = 1, nss
!       do jss = 1, nss
!         write(6,*) 'usor(',iss,',',jss,') =',usor(iss,jss)
!       end do
!     end do
!     write(6,*) 'eigenvectors of zheev/zjac (imag)'
!     do iss = 1, nss
!       do jss = 1, nss
!         write(6,*) 'usoi(',iss,',',jss,') =',usoi(iss,jss)
!       end do
!     end do

#ifdef _HDF5_
      call mh5_put_dset(wfn_sos_energy, ENSOR)
      call mh5_put_dset(wfn_sos_coefr,USOR)
      call mh5_put_dset(wfn_sos_coefi,USOI)
#endif
      !> free memory for H_SO - do not use it below!
      !> eigenvalues are stored in ENSOR!
      CALL mma_deallocate(HTOTR)
      CALL mma_deallocate(HTOTI)
C
C     BOR in Krapperup 070227
C     Compute J-values and Omega here instead of in subroutine PRPROP
C
      IAMFIX=0
      IAMFIY=0
      IAMFIZ=0
      IAMX=0
      IAMY=0
      IAMZ=0
      DO IPROP=1,NPROP
       IF(PNAME(IPROP)(1:4).EQ.'AMFI') THEN
         IF(ICOMP(IPROP).EQ.1) IAMFIX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMFIY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMFIZ=IPROP
       ELSE IF(PNAME(IPROP)(1:6).EQ.'ANGMOM') THEN
         IF(ICOMP(IPROP).EQ.1) IAMX=IPROP
         IF(ICOMP(IPROP).EQ.2) IAMY=IPROP
         IF(ICOMP(IPROP).EQ.3) IAMZ=IPROP
       END IF
      END DO

* The following matrix elements  require angular moment integrals:
C     IF(IAMX.eq.0 .or. IAMY.eq.0 .or. IAMZ.eq.0) GOTO 910
C Complex matrix elements of Jx, Jy, and/or Jz over spin states:
      CALL mma_allocate(LXI,NSS**2,Label='LXI')
      LXI(:)=0.0D0
      CALL mma_allocate(LYI,NSS**2,Label='LYI')
      LYI(:)=0.0D0
      CALL mma_allocate(LZI,NSS**2,Label='LZI')
      LZI(:)=0.0D0

      IF(IAMX.GT.0) CALL SMMAT(PROP,LXI,NSS,IAMX,0)
      IF(IAMY.GT.0) CALL SMMAT(PROP,LYI,NSS,IAMY,0)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,LZI,NSS,IAMZ,0)

      CALL mma_allocate(JXR,NSS**2,Label='JXR')
      CALL mma_allocate(JXI,NSS**2,Label='JIR')
      JXR(:)=0.0D0
      JXI(:)=0.0D0
      CALL mma_allocate(JYR,NSS**2,Label='JYR')
      CALL mma_allocate(JYI,NSS**2,Label='JYR')
      JYR(:)=0.0D0
      JYI(:)=0.0D0
      CALL mma_allocate(JZR,NSS**2,Label='JZR')
      CALL mma_allocate(JZI,NSS**2,Label='JZR')
      JZR(:)=0.0D0
      JZI(:)=0.0D0

      CALL SMMAT(PROP,JXR,NSS,0,1)
      CALL SMMAT(PROP,JYI,NSS,0,2)
      CALL SMMAT(PROP,JZR,NSS,0,3)

      CALL DAXPY_(NSS**2,1.0D0,LXI,1,JXI,1)
      CALL DAXPY_(NSS**2,1.0D0,LYI,1,JYI,1)
      CALL DAXPY_(NSS**2,1.0D0,LZI,1,JZI,1)

      CALL mma_deallocate(LXI)
      CALL mma_deallocate(LYI)
      CALL mma_deallocate(LZI)

      CALL mma_allocate(OMGR,NSS,NSS,Label='OMGR')
      CALL mma_allocate(OMGI,NSS,NSS,Label='OMGI')
      lOMG=.True.
      OMGR(:,:)=0.0D0
      OMGI(:,:)=0.0D0

      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JZR,NSS,
     &     JZR,NSS,0.0D0,OMGR,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,JZI,NSS,
     &     JZI,NSS,1.0D0,OMGR,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JZR,NSS,
     &     JZI,NSS,0.0D0,OMGI,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JZI,NSS,
     &     JZR,NSS,1.0D0,OMGI,NSS)

      CALL mma_deallocate(JZR)
      CALL mma_deallocate(JZI)

      lJ2 =.True.
      CALL mma_allocate(J2R,NSS,NSS,Label='J2R')
      CALL mma_allocate(J2I,NSS,NSS,Label='J2I')
      J2R(:,:)=0.0D0
      J2I(:,:)=0.0D0

      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JXR,NSS,
     &     JXR,NSS,1.0D0,J2R,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,JXI,NSS,
     &     JXI,NSS,1.0D0,J2R,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JYR,NSS,
     &     JYR,NSS,1.0D0,J2R,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,JYI,NSS,
     &     JYI,NSS,1.0D0,J2R,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JXR,NSS,
     &     JXI,NSS,1.0D0,J2I,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JXI,NSS,
     &     JXR,NSS,1.0D0,J2I,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JYR,NSS,
     &     JYI,NSS,1.0D0,J2I,NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,JYI,NSS,
     &     JYR,NSS,1.0D0,J2I,NSS)

      CALL mma_deallocate(JXR)
      CALL mma_deallocate(JXI)
      CALL mma_deallocate(JYR)
      CALL mma_deallocate(JYI)
      CALL ZTRNSF(NSS,USOR,USOI,OMGR,OMGI)
      CALL ZTRNSF(NSS,USOR,USOI,J2R,J2I)

* Jump here to skip computing omega and/or J:
C910  CONTINUE

      IF(IPGLOB.GE.1) THEN
       WRITE(6,*)
       WRITE(6,'(6X,A)')' Total energies including SO-coupling:'
       DO ISS=1,NSS
       E_tmp=ENSOR(ISS)+EMIN
       Call PrintResult(6, '(6x,A,I5,5X,A,F23.14)',
     &  'SO-RASSI State',ISS,'Total energy:',[E_tmp],1)
       END DO
      END IF

* Find E0=lowest energy, to use for printing table:
      IF(IPGLOB.GE.1) THEN
       E0=ENSOR(1)
       DO ISS=2,NSS
         E=ENSOR(ISS)
        IF(E.LT.E0) THEN
         E0=E
        END IF
       END DO
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)'  Eigenvalues of complex Hamiltonian:'
       WRITE(6,*)'  -----------------------------------'
       IF(EMIN.NE.0.0D0)
     &  WRITE(6,'(1X,A,F22.10,A1)')' (Shifted by EMIN (a.u.) =',EMIN,')'
       WRITE(6,*)
       if(ifj2.ne.0.and.ifjz.ne.0) then
        WRITE(6,*)'SO State       Relative EMIN(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)     J-value  Omega'
       else if(ifj2.ne.0.and.ifjz.eq.0) then
        WRITE(6,*)'SO State       Relative EMIN(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)     J-value'
       else if(ifj2.eq.0.and.ifjz.ne.0) then
        WRITE(6,*)'SO State       Relative EMIN(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)      Omega'
       else if(ifj2.eq.0.and.ifjz.eq.0) then
        WRITE(6,*)'SO State       Relative EMIN(au)   Rel lowest'//
     &          ' level(eV)    D:o, cm**(-1)'
       endif
       WRITE(6,*)
       E0=ENSOR(1)
       DO ISS=1,NSS
        E1=ENSOR(ISS)
        IF (E1.lt.E0) E0=E1
       END DO
       CALL MMA_ALLOCATE(ESO,NSS)
       DO ISS=1,NSS
        E1=ENSOR(ISS)
        E2=auToeV*(E1-E0)
        E3=auTocm*(E1-E0)
        IF (IFJ2.gt.0) THEN
         XJEFF=SQRT(0.25D0+J2R(ISS,ISS))-0.5D0
        END IF
        IF (IFJZ.gt.0) THEN
        OMEGA=SQRT(1.0D-12+OMGR(ISS,ISS))
        END IF

C Added by Ungur Liviu on 04.11.2009
C Saving the SO energies in ESO array.
        ESO(ISS)=E3
        if(ifj2.ne.0.and.ifjz.ne.0) then
          WRITE(6,'(1X,I5,7X,2(F18.10,2X),F18.4,4X,2(2X,F6.1))')
     &        ISS,E1,E2,E3,XJEFF,OMEGA
        else if(ifj2.ne.0.and.ifjz.eq.0) then
          WRITE(6,'(1X,I5,7X,2(F18.10,2X),F18.4,6X,F6.1)')
     &        ISS,E1,E2,E3,XJEFF
        else if(ifj2.eq.0.and.ifjz.ne.0) then
          WRITE(6,'(1X,I5,7X,2(F18.10,2X),F18.4,6X,F6.1)')
     &        ISS,E1,E2,E3,OMEGA
        else if(ifj2.eq.0.and.ifjz.eq.0) then
          WRITE(6,'(1X,I5,7X,2(F18.10,2X),F18.4)')
     &      ISS,E1,E2,E3
        endif
       ENDDO

C Added by Ungur Liviu on 04.11.2009
C Saving the ESO array in the RunFile.
       CALL Put_iscalar('NSS_SINGLE',NSS)
       CALL Put_dArray( 'ESO_SINGLE',ESO,NSS)
       CALL Put_dArray( 'ESO_LOW'   ,ENSOR+EMIN,NSS)
       CALL MMA_DEALLOCATE(ESO)
      END IF

      IF(lOMG) THEN
       CALL mma_deallocate(OMGR)
       CALL mma_deallocate(OMGI)
      END IF
      IF(lJ2) THEN
       CALL mma_deallocate(J2R)
       CALL mma_deallocate(J2I)
      END IF

C Put energy onto info file for automatic verification runs:
      EPSS=5.0D-11
      EPSH=MAX(5.0D-10,ABS(ENSOR(1)+EMIN)*EPSS)
      IDX=100
      DO ISS=1,NSS
       EI=(ENSOR(ISS)+EMIN)*EPSS
       V2SUM=0.0D0
       DO JSS=1,NSS
        V2SUM=V2SUM+USOR(JSS,ISS)**2+USOI(JSS,ISS)**2
       END DO
       ERMS=SQRT(EPSH**2+EI**2)*V2SUM
       IDX=MIN(IDX,INT(-LOG10(ERMS)))
      END DO
      iTol=cho_x_gettol(IDX) ! reset thr iff Cholesky
      Call Add_Info('ESO_LOW',ENSOR+EMIN,NSS,iTol)

      IF(IPGLOB.GE.3) THEN
       WRITE(6,*)
       WRITE(6,*)' Complex eigenvectors in basis of non-so eigenstates:'
       WRITE(6,*)'-----------------------------------------------------'
       WRITE(6,*)
       IF(IPGLOB.GE.4) THEN
         FRAC=0.0D0
       ELSE
         FRAC=0.25D0
         WRITE(6,*)'    (A selection of the largest components)'
       END IF
      END IF
       CALL PRCEVC(NSS,FRAC,ENSOR,MAPST,MAPSP,MAPMS,USOR,USOI)

C Update LoopDivide (SUBSets keyword)
C Assume the SO "ground states" are mostly formed by the SF "ground states"
      If (ReduceLoop) Then
        Call mma_Allocate(IndexE,nState,Label='IndexE')
        IndexE(:)=ArgSort(Energy, leq_r)
        n=0
        Do iState=1,LoopDivide
          Job=JbNum(IndexE(iState))
          n=n+Mltplt(Job)
        End Do
        LoopDivide=n
#ifdef _HDF5_
        If (IFTRD1.or.IFTDM)
     &    Call UpdateIdx(IndexE,nSS,USOR,USOI,MapSt)
#endif
        Call mma_deAllocate(IndexE)
      End If

      CALL mma_deallocate(MAPST)
      CALL mma_deallocate(MAPSP)
      CALL mma_deallocate(MAPMS)

      call mma_deallocate(HAMSOR)
      call mma_deallocate(HAMSOI)

      END SUBROUTINE SOEIG
