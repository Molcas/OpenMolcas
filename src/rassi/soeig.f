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
      use rassi_global_arrays, only: JBNUM
      use sorting, only : argsort
      use sorting_funcs, only : leq_r
#ifdef _HDF5_
      use Dens2HDF5
#endif
#ifdef _DMRG_
      use qcmaquis_interface_cfg
#endif
      IMPLICIT NONE
#include "prgm.fh"
#include "SysDef.fh"
#include "Molcas.fh"
#include "cntrl.fh"
#include "rassi.fh"
#include "symmul.fh"
#include "Files.fh"
#include "WrkSpc.fh"
#include "constants.fh"
#include "stdalloc.fh"
#include "rassiwfn.fh"

      CHARACTER*16 ROUTINE
      PARAMETER (ROUTINE='SOEIG')

      INTEGER NSS
      REAL*8 USOR(NSS,NSS),USOI(NSS,NSS),ENSOR(NSS)
      REAL*8 PROP(NSTATE,NSTATE,NPROP),ENERGY(NSTATE)

      INTEGER I,N
      INTEGER ITOL
      INTEGER JOB
      INTEGER IPROP
      INTEGER IAMFIX,IAMFIY,IAMFIZ,IAMX,IAMY,IAMZ
      INTEGER ISS,JSS,IJSS,ISTATE,JSTATE
      INTEGER LHTOTI,LHTOTR
      INTEGER LJ2I,LJ2R,LJXI,LJXR,LJYI,LJYR,LJZI,LJZR,LLXI,LLYI,LLZI
      INTEGER LMAPMS,LMAPSP,LMAPST,LOMGI,LOMGR
      INTEGER MAGN
      INTEGER MPLET,MPLET1,MPLET2,MSPROJ,MSPROJ1,MSPROJ2

      REAL*8 AU2EV,AU2CM
      REAL*8 AMFIX,AMFIY,AMFIZ
      REAL*8 CG0,CGM,CGP,CGX,CGY
      REAL*8 E,E0,E1,E2,E3,E_TMP,FACT,FRAC
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
      Integer  cho_x_gettol
      External cho_x_gettol
      LOGICAL :: debug_dmrg_rassi_code = .false.




C CONSTANTS:
      AU2EV=CONV_AU_TO_EV_
      AU2CM=CONV_AU_TO_CM1_
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
      CALL GETMEM('MAPST','ALLO','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','ALLO','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','ALLO','INTE',LMAPMS,NSS)
      ISS=0
      DO ISTATE=1,NSTATE
       JOB=JBNUM(ISTATE)
       MPLET=MLTPLT(JOB)
       DO MSPROJ=-MPLET+1,MPLET-1,2
        ISS=ISS+1
        IWORK(LMAPST-1+ISS)=ISTATE
        IWORK(LMAPSP-1+ISS)=MPLET
        IWORK(LMAPMS-1+ISS)=MSPROJ
       END DO
      END DO
C Complex hamiltonian matrix elements over spin states:
      CALL GETMEM('HTOTR','ALLO','REAL',LHTOTR,NSS**2)
      CALL GETMEM('HTOTI','ALLO','REAL',LHTOTI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LHTOTR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LHTOTI),1)

      IF(IPGLOB.GE.TERSE) THEN
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,A,34X,A,34X,A)')
     &      '*','       Spin-orbit section     ','*'
       WRITE(6,'(6X,A,98X,A)') '*','*'
       WRITE(6,'(6X,100A1)') ('*',i=1,100)
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
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,NSS
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
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
          IJSS=ISS+NSS*(JSS-1)
C WIGNER-ECKART THEOREM:
          FACT=1.0D0/SQRT(DBLE(MPLET1))
          IF(MPLET1.EQ.MPLET2-2) FACT=-FACT
          CGM=FACT*DCLEBS(S2,1.0D0,S1,SM2,-1.0D0,SM1)
          CG0=FACT*DCLEBS(S2,1.0D0,S1,SM2, 0.0D0,SM1)
          CGP=FACT*DCLEBS(S2,1.0D0,S1,SM2,+1.0D0,SM1)
          CGX= SQRT(0.5D0)*(CGM-CGP)
          CGY=-SQRT(0.5D0)*(CGM+CGP)
C SPIN-ORBIT HAMILTONIAN MATRIX ELEMENTS:
          HSOR=CGY*AMFIY
          HSOI=CGX*AMFIX+CG0*AMFIZ
* PAM07: Delay addition of diagonal scalar part until later, see below:
*          WORK(LHTOTR-1+IJSS)=HSCAL+HSOR
          WORK(LHTOTR-1+IJSS)=HSOR
          WORK(LHTOTI-1+IJSS)=HSOI

        END DO
      END DO

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
         IJSS=ISS+NSS*(JSS-1)
         HSOR=WORK(LHTOTR-1+IJSS)
         HSOI=WORK(LHTOTI-1+IJSS)
         HSOTOT=SQRT(HSOR**2+HSOI**2)
         IF(HSOTOT*AU2CM.GT.SOTHR_PRT) N=N+1
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
        ISTATE=IWORK(LMAPST-1+ISS)
        MPLET1=IWORK(LMAPSP-1+ISS)
        MSPROJ1=IWORK(LMAPMS-1+ISS)
        S1=0.5D0*DBLE(MPLET1-1)
        SM1=0.5D0*DBLE(MSPROJ1)
        DO JSS=1,ISS
         IJSS=ISS+NSS*(JSS-1)
         HSOR=WORK(LHTOTR-1+IJSS)
         HSOI=WORK(LHTOTI-1+IJSS)
         HSOTOT=SQRT(HSOR**2+HSOI**2)
         IF(HSOTOT*AU2CM.GE.SOTHR_PRT) THEN
          JSTATE=IWORK(LMAPST-1+JSS)
          MPLET2=IWORK(LMAPSP-1+JSS)
          MSPROJ2=IWORK(LMAPMS-1+JSS)
          S2=0.5D0*DBLE(MPLET2-1)
          SM2=0.5D0*DBLE(MSPROJ2)
         WRITE(6,'(1X,I5,F5.1,F5.1,I5,F5.1,F5.1,3F14.3)') ISS,S1,SM1,
     &         JSS,S2,SM2,HSOR*AU2CM,HSOI*AU2CM,HSOTOT*AU2CM
         END IF
        END DO
       END DO
       WRITE(6,'(1X,10A7)')('-------',I=1,10)

       END IF
      ENDIF

* PAM07: Addition of scalar diagonal part was delayed until here, see above.
      DO ISS=1,NSS
        ISTATE=IWORK(LMAPST-1+ISS)
        IJSS=ISS+NSS*(ISS-1)
        HSOR=WORK(LHTOTR-1+IJSS)
        WORK(LHTOTR-1+IJSS)=HSOR+ENERGY(ISTATE)
      END DO


      IF(IPGLOB.GE.VERBOSE) THEN
       WRITE(6,*)
       WRITE(6,*)
       WRITE(6,*)'Complex Hamiltonian matrix including SO-coupling'
       WRITE(6,*)'over spin components of spin-free eigenstates (SFS):'
       WRITE(6,'(1X,11A7)')('-------',I=1,11)
       CALL PRCHAM(NSS,WORK(LHTOTR),WORK(LHTOTI))
       WRITE(6,'(1X,11A7)')('-------',I=1,11)
      ENDIF
      ! save the Hamiltonian
      call mma_allocate(HAMSOR,NSS,NSS,'HAMSOR')
      call mma_allocate(HAMSOI,NSS,NSS,'HAMSOI')
      call dcopy_(NSS*NSS,[0.d0],0,HAMSOR,1)
      call dcopy_(NSS*NSS,[0.d0],0,HAMSOI,1)
      call dcopy_(NSS*NSS,WORK(LHTOTR),1,HAMSOR,1)
      call dcopy_(NSS*NSS,WORK(LHTOTI),1,HAMSOI,1)
      call put_darray('HAMSOR_SINGLE',HAMSOR,NSS*NSS)
      call put_darray('HAMSOI_SINGLE',HAMSOI,NSS*NSS)
#ifdef _HDF5_
      call mh5_put_dset_array_real(wfn_sos_hsor,HAMSOR)
      call mh5_put_dset_array_real(wfn_sos_hsoi,HAMSOI)
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
            hso_tmp(iss,jss) = dcmplx(WORK(LHTOTR-1+ISS+NSS*(JSS-1)),
     &                                WORK(LHTOTI-1+ISS+NSS*(JSS-1)))
!         write(6,*) ' hso_tmp(',iss,',',jss,') = ',hso_tmp(iss,jss)
          END DO
        END DO

        lcwork = (2*nss-1); info = 0
        call zheev('V','U',nss,hso_tmp,nss,ensor,ccwork,lcwork,
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
        CALL ZJAC(NSS,WORK(LHTOTR),WORK(LHTOTI),
     &          NSS,USOR,USOI)
        DO ISS=1,NSS
         ENSOR(ISS)=WORK(LHTOTR-1+ISS+NSS*(ISS-1))
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
      call mh5_put_dset_array_real(wfn_sos_coefr,USOR)
      call mh5_put_dset_array_real(wfn_sos_coefi,USOI)
#endif
      !> free memory for H_SO - do not use it below!
      !> eigenvalues are stored in ENSOR!
      CALL GETMEM('HTOTR','FREE','REAL',LHTOTR,NSS**2)
      CALL GETMEM('HTOTI','FREE','REAL',LHTOTI,NSS**2)
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
      CALL GETMEM('LXI','ALLO','REAL',LLXI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLXI),1)
      CALL GETMEM('LYI','ALLO','REAL',LLYI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLYI),1)
      CALL GETMEM('LZI','ALLO','REAL',LLZI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LLZI),1)
      IF(IAMX.GT.0) CALL SMMAT(PROP,WORK(LLXI),NSS,IAMX,0)
      IF(IAMY.GT.0) CALL SMMAT(PROP,WORK(LLYI),NSS,IAMY,0)
      IF(IAMZ.GT.0) CALL SMMAT(PROP,WORK(LLZI),NSS,IAMZ,0)

      CALL GETMEM('JXR','ALLO','REAL',LJXR,NSS**2)
      CALL GETMEM('JXI','ALLO','REAL',LJXI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJXR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJXI),1)
      CALL GETMEM('JYR','ALLO','REAL',LJYR,NSS**2)
      CALL GETMEM('JYI','ALLO','REAL',LJYI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJYR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJYI),1)
      CALL GETMEM('JZR','ALLO','REAL',LJZR,NSS**2)
      CALL GETMEM('JZI','ALLO','REAL',LJZI,NSS**2)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJZR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LJZI),1)

      CALL SMMAT(PROP,WORK(LJXR),NSS,0,1)
      CALL SMMAT(PROP,WORK(LJYI),NSS,0,2)
      CALL SMMAT(PROP,WORK(LJZR),NSS,0,3)

      CALL DAXPY_(NSS**2,1.0D0,WORK(LLXI),1,WORK(LJXI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLYI),1,WORK(LJYI),1)
      CALL DAXPY_(NSS**2,1.0D0,WORK(LLZI),1,WORK(LJZI),1)

      CALL GETMEM('LXI','FREE','REAL',LLXI,NSS**2)
      CALL GETMEM('LYI','FREE','REAL',LLYI,NSS**2)
      CALL GETMEM('LZI','FREE','REAL',LLZI,NSS**2)

      CALL GETMEM('OMGR','ALLO','REAL',LOMGR,NSS**2)
      CALL GETMEM('OMGI','ALLO','REAL',LOMGI,NSS**2)
      lOMG=.True.
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LOMGR),1)
      CALL DCOPY_(NSS**2,[0.0D0],0,WORK(LOMGI),1)

      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJZR),NSS,
     &     WORK(LJZR),NSS,0.0D0,WORK(LOMGR),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,WORK(LJZI),NSS,
     &     WORK(LJZI),NSS,1.0D0,WORK(LOMGR),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJZR),NSS,
     &     WORK(LJZI),NSS,0.0D0,WORK(LOMGI),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJZI),NSS,
     &     WORK(LJZR),NSS,1.0D0,WORK(LOMGI),NSS)

      CALL GETMEM('JZR','FREE','REAL',LJZR,NSS**2)
      CALL GETMEM('JZI','FREE','REAL',LJZI,NSS**2)

      lJ2 =.True.
      CALL GETMEM('J2R','ALLO','REAL',LJ2R,NSS**2)
      CALL GETMEM('J2I','ALLO','REAL',LJ2I,NSS**2)
      CALL DCOPY_(NSS**2,WORK(LOMGR),1,WORK(LJ2R),1)
      CALL DCOPY_(NSS**2,WORK(LOMGI),1,WORK(LJ2I),1)

      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJXR),NSS,
     &     WORK(LJXR),NSS,1.0D0,WORK(LJ2R),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,WORK(LJXI),NSS,
     &     WORK(LJXI),NSS,1.0D0,WORK(LJ2R),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJYR),NSS,
     &     WORK(LJYR),NSS,1.0D0,WORK(LJ2R),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS,-1.0D0,WORK(LJYI),NSS,
     &     WORK(LJYI),NSS,1.0D0,WORK(LJ2R),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJXR),NSS,
     &     WORK(LJXI),NSS,1.0D0,WORK(LJ2I),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJXI),NSS,
     &     WORK(LJXR),NSS,1.0D0,WORK(LJ2I),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJYR),NSS,
     &     WORK(LJYI),NSS,1.0D0,WORK(LJ2I),NSS)
      CALL DGEMM_('NSS','NSS',NSS,NSS,NSS, 1.0D0,WORK(LJYI),NSS,
     &     WORK(LJYR),NSS,1.0D0,WORK(LJ2I),NSS)

      CALL GETMEM('JXR','FREE','REAL',LJXR,NSS**2)
      CALL GETMEM('JXI','FREE','REAL',LJXI,NSS**2)
      CALL GETMEM('JYR','FREE','REAL',LJYR,NSS**2)
      CALL GETMEM('JYI','FREE','REAL',LJYI,NSS**2)
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LOMGR),WORK(LOMGI))
      CALL ZTRNSF(NSS,USOR,USOI,WORK(LJ2R),WORK(LJ2I))

* Jump here to skip computing omega and/or J:
C910  CONTINUE

      IF(IPGLOB.GE.TERSE) THEN
       WRITE(6,*)
       WRITE(6,'(6X,A)')' Total energies including SO-coupling:'
       DO ISS=1,NSS
       E_tmp=ENSOR(ISS)+EMIN
       Call PrintResult(6, '(6x,A,I5,5X,A,F16.8)',
     &  'SO-RASSI State',ISS,'Total energy:',[E_tmp],1)
       END DO
      END IF

* Find E0=lowest energy, to use for printing table:
      IF(IPGLOB.GE.TERSE) THEN
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
        E2=AU2EV*(E1-E0)
        E3=AU2CM*(E1-E0)
        IF (IFJ2.gt.0) THEN
         XJEFF=SQRT(0.25D0+WORK(LJ2R-1+ISS+NSS*(ISS-1)))-0.5D0
        END IF
        IF (IFJZ.gt.0) THEN
        OMEGA=SQRT(1.0D-12+WORK(LOMGR-1+ISS+NSS*(ISS-1)))
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
       CALL MMA_DEALLOCATE(ESO)
      END IF

      IF(lOMG) THEN
       CALL GETMEM('OMGR','FREE','REAL',LOMGR,NSS**2)
       CALL GETMEM('OMGI','FREE','REAL',LOMGI,NSS**2)
      END IF
      IF(lJ2) THEN
       CALL GETMEM('J2R','FREE','REAL',LJ2R,NSS**2)
       CALL GETMEM('J2I','FREE','REAL',LJ2I,NSS**2)
      END IF

C Put energy onto info file for automatic verification runs:
      iTol=cho_x_gettol(8) ! reset thr iff Cholesky
      Call Add_Info('ESO_LOW',ENSOR+EMIN,NSS,iTol)

      IF(IPGLOB.GE.VERBOSE) THEN
       WRITE(6,*)
       WRITE(6,*)' Complex eigenvectors in basis of non-so eigenstates:'
       WRITE(6,*)'-----------------------------------------------------'
       WRITE(6,*)
       IF(IPGLOB.GE.DEBUG) THEN
         FRAC=0.0D0
       ELSE
         FRAC=0.25D0
         WRITE(6,*)'    (A selection of the largest components)'
       END IF
      END IF
       CALL PRCEVC(NSS,FRAC,ENSOR,IWORK(LMAPST),IWORK(LMAPSP),
     &            IWORK(LMAPMS),USOR,USOI)

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
        If (IFTRD1.or.IFTRD2)
     &    Call UpdateIdx(IndexE,nSS,USOR,USOI,iWork(lMapSt))
#endif
        Call mma_deAllocate(IndexE)
      End If

      CALL GETMEM('MAPST','FREE','INTE',LMAPST,NSS)
      CALL GETMEM('MAPSP','FREE','INTE',LMAPSP,NSS)
      CALL GETMEM('MAPMS','FREE','INTE',LMAPMS,NSS)

      call mma_deallocate(HAMSOR)
      call mma_deallocate(HAMSOI)
      RETURN

      END
