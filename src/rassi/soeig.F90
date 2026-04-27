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

!ifdef _DEBUGPRINT_
subroutine SOEIG(PROP,USOR,USOI,ENSOR,NSS,ENERGY)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: JBNUM
use sorting, only: argsort
use sorting_funcs, only: leq_r
use Cntrl, only: EMIN, ICOMP, IFJ2, IFJZ, LOOPDIVIDE, MLTPLT, NPROP, NSOThr_PRT, NSTATE, PNAME, REDUCELOOP, SOThr_PRT
#ifdef _HDF5_
use Dens2HDF5, only: UpdateIdx
use mh5, only: mh5_put_dset
use RASSIWfn, only: wfn_SOS_CoefI, wfn_SOS_CoefR, wfn_SOS_Energy, wfn_SOS_HSOI, wfn_SOS_HSOR, wfn_SOS_VSOI, wfn_SOS_VSOR
use Cntrl, only: IFTDM, IFTRD1, RHODYN
#endif
#ifdef _DMRG_
use rasscf_global, only: doDMRG
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Ten, Half, Quart, auTocm, auToeV
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp), intent(in) :: PROP(NSTATE,NSTATE,NPROP)
integer(kind=iwp), intent(in) :: NSS
real(kind=wp), intent(inout) :: USOR(NSS,NSS), USOI(NSS,NSS), ENERGY(NSTATE)
real(kind=wp), intent(out) :: ENSOR(NSS)
integer(kind=iwp) :: IAMFIX, IAMFIY, IAMFIZ, IAMX, IAMY, IAMZ, IDX, IPROP, ISS, ISTATE, ITOL, JOB, JSS, JSTATE, MAGN, MPLET, &
                     MPLET1, MPLET2, MSPROJ, MSPROJ1, MSPROJ2, N
real(kind=wp) :: AMFIX, AMFIY, AMFIZ, CG0, CGM, CGP, CGX, CGY, E0, E1, E2, E3, E_TMP, EI, EPSH, EPSS, ERMS, FACT, FRAC, HSOI, &
                 HSOR, HSOTOT, OMEGA, S1, S2, SM1, SM2, SOTHR_MIN, X, X_THR, XJEFF
integer(kind=iwp), allocatable :: IndexE(:), MAPMS(:), MAPSP(:), MAPST(:)
real(kind=wp), allocatable :: ESO(:), HAMSOR(:,:), HTOTI(:,:), HTOTR(:,:), J2I(:,:), J2R(:,:), JXI(:), JXR(:), JYI(:), JYR(:), &
                              JZI(:), JZR(:), LXI(:), LYI(:), LZI(:), OMGI(:,:), OMGR(:,:)
#ifdef _DMRG_
integer(kind=iwp) :: info, lcwork
real(kind=wp), allocatable :: rwork(:)
complex(kind=wp), allocatable :: ccwork(:), hso_tmp(:,:)
#endif
integer(kind=iwp), external :: cho_x_gettol
real(kind=wp), external :: DCLEBS

! Identify AMFI and ANGMOM matrix elements:
IAMFIX = 0
IAMFIY = 0
IAMFIZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:4) == 'AMFI') then
    if (ICOMP(IPROP) == 1) IAMFIX = IPROP
    if (ICOMP(IPROP) == 2) IAMFIY = IPROP
    if (ICOMP(IPROP) == 3) IAMFIZ = IPROP
  else if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    if (ICOMP(IPROP) == 1) IAMX = IPROP
    if (ICOMP(IPROP) == 2) IAMY = IPROP
    if (ICOMP(IPROP) == 3) IAMZ = IPROP
  end if
end do

! Mapping from spin states to spin-free state and to spin:
call mma_allocate(MAPST,NSS,Label='MAPST')
call mma_allocate(MAPSP,NSS,Label='MAPSP')
call mma_allocate(MAPMS,NSS,Label='MAPMS')
ISS = 0
do ISTATE=1,NSTATE
  JOB = JBNUM(ISTATE)
  MPLET = MLTPLT(JOB)
  do MSPROJ=-MPLET+1,MPLET-1,2
    ISS = ISS+1
    MAPST(ISS) = ISTATE
    MAPSP(ISS) = MPLET
    MAPMS(ISS) = MSPROJ
  end do
end do
! Complex hamiltonian matrix elements over spin states:
call mma_allocate(HTOTR,NSS,NSS,Label='HTOTR')
call mma_allocate(HTOTI,NSS,NSS,Label='HTOTI')

if (IPGLOB >= 1) then
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,34X,A,34X,A)') '*','       Spin-orbit section     ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'BLUBB BLUBB debug print of property matrix'
do istate=1,nstate
  do jstate=1,nstate
    do IPROP=1,NPROP
      if (abs(prop(istate,jstate,iprop)) > 1.0e-14_wp) &
        write(u6,*) 'prop(',istate,',',jstate,',',iprop,') = ',prop(istate,jstate,iprop)
    end do
  end do
end do
#endif

do ISS=1,NSS
  ISTATE = MAPST(ISS)
  MPLET1 = MAPSP(ISS)
  MSPROJ1 = MAPMS(ISS)
  S1 = Half*real(MPLET1-1,kind=wp)
  SM1 = Half*real(MSPROJ1,kind=wp)
  do JSS=1,NSS
    JSTATE = MAPST(JSS)
    MPLET2 = MAPSP(JSS)
    MSPROJ2 = MAPMS(JSS)
    S2 = Half*real(MPLET2-1,kind=wp)
    SM2 = Half*real(MSPROJ2,kind=wp)
    AMFIX = Zero
    if (IAMFIX /= 0) AMFIX = PROP(ISTATE,JSTATE,IAMFIX)
    AMFIY = Zero
    if (IAMFIY /= 0) AMFIY = PROP(ISTATE,JSTATE,IAMFIY)
    AMFIZ = Zero
    if (IAMFIZ /= 0) AMFIZ = PROP(ISTATE,JSTATE,IAMFIZ)
    ! PAM07 HSCAL = Zero
    ! PAM07 if (ISS == JSS) HSCAL = ENERGY(ISTATE)
    ! WIGNER-ECKART THEOREM:
    FACT = One/sqrt(real(MPLET1,kind=wp))
    if (MPLET1 == MPLET2-2) FACT = -FACT
    CGM = FACT*DCLEBS(S2,One,S1,SM2,-One,SM1)
    CG0 = FACT*DCLEBS(S2,One,S1,SM2,Zero,SM1)
    CGP = FACT*DCLEBS(S2,One,S1,SM2,One,SM1)
    CGX = sqrt(Half)*(CGM-CGP)
    CGY = -sqrt(Half)*(CGM+CGP)
    ! SPIN-ORBIT HAMILTONIAN MATRIX ELEMENTS:
    !  according to expressions between eqs. (5) and (6)
    !  in Malmqvist et all CPL 357 (2002) 230-240, but real and
    !  imaginary parts are swapped here!!! more precisely the Hamiltonian
    !  is multiplied by imaginary unit to keep its hermicity
    HSOR = CGY*AMFIY
    HSOI = CGX*AMFIX+CG0*AMFIZ
    ! PAM07: Delay addition of diagonal scalar part until later, see below:
    HTOTR(ISS,JSS) = HSOR
    HTOTI(ISS,JSS) = HSOI

  end do
end do

! VKochetov 2021 put SOC matrix elements to hdf5:
#ifdef _HDF5_
if (rhodyn) then
  call mh5_put_dset(wfn_sos_vsor,HTOTR)
  call mh5_put_dset(wfn_sos_vsoi,HTOTI)
end if
#endif

! Perhaps write out large spin-orbit coupling elements:
if (NSOTHR_PRT > 0) then
  ! Prevent infinite loop below:
  SOTHR_MIN = max(SOTHR_PRT,1.0e-6_wp)
  ! And work with quantities larger than 1:
  X_THR = SOTHR_PRT/SOTHR_MIN
  do
    N = count(sqrt(HTOTR(:,:)**2+HTOTI(:,:)**2)*auTocm > SOTHR_PRT)
    if (N <= NSOTHR_PRT) exit
    X_THR = X_THR*1.2_wp
    MAGN = int(log10(X_THR))
    X = X_THR/real(10**MAGN,kind=wp)
    X_THR = 0.1_wp*real(nint(Ten*X),kind=wp)*real(10**MAGN,kind=wp)
    SOTHR_PRT = X_THR*SOTHR_MIN
  end do

  if (N > 0) then
    write(u6,*)
    write(u6,*) 'Complex SO-Hamiltonian matrix elements over'
    write(u6,*) 'spin components of spin-free eigenstates (SFS):'
    write(u6,'(1x,A,F10.3,A)') '(In cm-1. Print threshold: ',SOTHR_PRT,' cm-1)'
    write(u6,'(1X,A)') repeat('-',70)
    write(u6,*)
    write(u6,'(A)') '  I1  S1  MS1    I2  S2  MS2    Real part    Imag part      Absolute'
    do ISS=1,NSS
      ISTATE = MAPST(ISS)
      MPLET1 = MAPSP(ISS)
      MSPROJ1 = MAPMS(ISS)
      S1 = Half*real(MPLET1-1,kind=wp)
      SM1 = Half*real(MSPROJ1,kind=wp)
      do JSS=1,ISS
        HSOR = HTOTR(ISS,JSS)
        HSOI = HTOTI(ISS,JSS)
        HSOTOT = sqrt(HSOR**2+HSOI**2)
        if (HSOTOT*auTocm >= SOTHR_PRT) then
          JSTATE = MAPST(JSS)
          MPLET2 = MAPSP(JSS)
          MSPROJ2 = MAPMS(JSS)
          S2 = Half*real(MPLET2-1,kind=wp)
          SM2 = Half*real(MSPROJ2,kind=wp)
          write(u6,'(1X,I5,F5.1,F5.1,I5,F5.1,F5.1,3F14.3)') ISS,S1,SM1,JSS,S2,SM2,HSOR*auTocm,HSOI*auTocm,HSOTOT*auTocm
        end if
      end do
    end do
    write(u6,'(1X,A)') repeat('-',70)

  end if
end if

! PAM07: Addition of scalar diagonal part was delayed until here, see above.
do ISS=1,NSS
  ISTATE = MAPST(ISS)
  HTOTR(ISS,ISS) = HTOTR(ISS,ISS)+ENERGY(ISTATE)
end do

if (IPGLOB >= 3) then
  write(u6,*)
  write(u6,*)
  write(u6,*) 'Complex Hamiltonian matrix including SO-coupling'
  write(u6,*) 'over spin components of spin-free eigenstates (SFS):'
  write(u6,'(1X,A)') repeat('-',77)
  call PRCHAM(NSS,HTOTR,HTOTI)
  write(u6,'(1X,A)') repeat('-',77)
end if
! save the Hamiltonian
call mma_allocate(HAMSOR,NSS,NSS,'HAMSOR')
HAMSOR(:,:) = HTOTR(:,:)
call put_darray('HAMSOR_SINGLE',HTOTR,NSS**2)
call put_darray('HAMSOI_SINGLE',HTOTI,NSS**2)
#ifdef _HDF5_
! unshift the diagonal
do ISS=1,NSS
  HAMSOR(ISS,ISS) = HAMSOR(ISS,ISS)+EMIN
end do
call mh5_put_dset(wfn_sos_hsor,HAMSOR)
call mh5_put_dset(wfn_sos_hsoi,HTOTI)
#endif
call mma_deallocate(HAMSOR)

!> use complex matrix diagonalization
#ifdef _DMRG_
if (doDMRG) then
  call mma_allocate(hso_tmp,nss,nss)
  call mma_allocate(ccwork,(2*nss-1))
  call mma_allocate(rwork,(3*nss-2))
  hso_tmp = 0; ccwork = 0; rwork = 0

  hso_tmp(:,:) = cmplx(HTOTR(:,:),HTOTI(:,:),kind=wp)
  !do jss=1,nss
  !  do iss=1,nss
  !    write(u6,*) ' hso_tmp(',iss,',',jss,') = ',hso_tmp(iss,jss)
  !  end do
  !end do

  lcwork = (2*nss-1); info = 0
  call zheev_('V','U',nss,hso_tmp,nss,ensor,ccwork,lcwork,rwork,info)

  if (info /= 0) then
    write(u6,*) '* WARNING in rassi/soeig *'
    write(u6,*) 'zheev did return with an error message,info=',info
  else
    write(u6,*) 'zheev in rassi/soeig succeeded!'
  end if

  !write(u6,*) 'eigenvalues of zheev, info=',info
  !do iss=1,nss
  !  write(u6,*) 'ensor(',iss,') =',ensor(iss)
  !end do

  !> save eigenvectors
  usor(:,:) = real(hso_tmp(:,:))
  usoi(:,:) = aimag(hso_tmp(:,:))

  !> sort eigenvalues in increasing sequence (using the same algorithm as zjac)
  call zorder(nss,nss,usor,usoi,ensor,0)

  call MMA_DEALLOCATE(hso_tmp)
  call MMA_DEALLOCATE(ccwork)
  call MMA_DEALLOCATE(rwork)
else
#endif
  !> diagonalize H_SO and get array of eigenvalues/eigenvectors
  call ZJAC(NSS,HTOTR,HTOTI,NSS,USOR,USOI)
  do ISS=1,NSS
    ENSOR(ISS) = HTOTR(ISS,ISS)
  end do
#ifdef _DMRG_
end if
#endif

!write(u6,*) 'eigenvectors of zheev/zjac (real)'
!do iss=1,nss
!  do jss=1,nss
!    write(u6,*) 'usor(',iss,',',jss,') =',usor(iss,jss)
!  end do
!end do
!write(u6,*) 'eigenvectors of zheev/zjac (imag)'
!do iss=1,nss
!  do jss=1,nss
!    write(u6,*) 'usoi(',iss,',',jss,') =',usoi(iss,jss)
!  end do
!end do

#ifdef _HDF5_
call mh5_put_dset(wfn_sos_energy,ENSOR(:)+EMIN)
call mh5_put_dset(wfn_sos_coefr,USOR)
call mh5_put_dset(wfn_sos_coefi,USOI)
#endif
!> free memory for H_SO - do not use it below!
!> eigenvalues are stored in ENSOR!
call mma_deallocate(HTOTR)
call mma_deallocate(HTOTI)

! BOR in Krapperup 070227
! Compute J-values and Omega here instead of in subroutine PRPROP

IAMFIX = 0
IAMFIY = 0
IAMFIZ = 0
IAMX = 0
IAMY = 0
IAMZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:4) == 'AMFI') then
    if (ICOMP(IPROP) == 1) IAMFIX = IPROP
    if (ICOMP(IPROP) == 2) IAMFIY = IPROP
    if (ICOMP(IPROP) == 3) IAMFIZ = IPROP
  else if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    if (ICOMP(IPROP) == 1) IAMX = IPROP
    if (ICOMP(IPROP) == 2) IAMY = IPROP
    if (ICOMP(IPROP) == 3) IAMZ = IPROP
  end if
end do

! The following matrix elements  require angular moment integrals:
!if ((IAMX /= 0) .and. (IAMY /= 0) .and. (IAMZ /= 0)) then
! Complex matrix elements of Jx, Jy, and/or Jz over spin states:
call mma_allocate(LXI,NSS**2,Label='LXI')
LXI(:) = Zero
call mma_allocate(LYI,NSS**2,Label='LYI')
LYI(:) = Zero
call mma_allocate(LZI,NSS**2,Label='LZI')
LZI(:) = Zero

if (IAMX > 0) call SMMAT(PROP,LXI,NSS,IAMX,0)
if (IAMY > 0) call SMMAT(PROP,LYI,NSS,IAMY,0)
if (IAMZ > 0) call SMMAT(PROP,LZI,NSS,IAMZ,0)

call mma_allocate(JXR,NSS**2,Label='JXR')
call mma_allocate(JXI,NSS**2,Label='JIR')
JXR(:) = Zero
JXI(:) = Zero
call mma_allocate(JYR,NSS**2,Label='JYR')
call mma_allocate(JYI,NSS**2,Label='JYR')
JYR(:) = Zero
JYI(:) = Zero
call mma_allocate(JZR,NSS**2,Label='JZR')
call mma_allocate(JZI,NSS**2,Label='JZR')
JZR(:) = Zero
JZI(:) = Zero

call SMMAT(PROP,JXR,NSS,0,1)
call SMMAT(PROP,JYI,NSS,0,2)
call SMMAT(PROP,JZR,NSS,0,3)

JXI(:) = JXI(:)+LXI(:)
JYI(:) = JYI(:)+LYI(:)
JZI(:) = JZI(:)+LZI(:)

call mma_deallocate(LXI)
call mma_deallocate(LYI)
call mma_deallocate(LZI)

call mma_allocate(OMGR,NSS,NSS,Label='OMGR')
call mma_allocate(OMGI,NSS,NSS,Label='OMGI')

call DGEMM_('N','N',NSS,NSS,NSS,One,JZR,NSS,JZR,NSS,Zero,OMGR,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,-One,JZI,NSS,JZI,NSS,One,OMGR,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JZR,NSS,JZI,NSS,Zero,OMGI,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JZI,NSS,JZR,NSS,One,OMGI,NSS)

call mma_deallocate(JZR)
call mma_deallocate(JZI)

call mma_allocate(J2R,NSS,NSS,Label='J2R')
call mma_allocate(J2I,NSS,NSS,Label='J2I')

call DGEMM_('N','N',NSS,NSS,NSS,One,JXR,NSS,JXR,NSS,Zero,J2R,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,-One,JXI,NSS,JXI,NSS,One,J2R,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JYR,NSS,JYR,NSS,One,J2R,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,-One,JYI,NSS,JYI,NSS,One,J2R,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JXR,NSS,JXI,NSS,Zero,J2I,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JXI,NSS,JXR,NSS,One,J2I,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JYR,NSS,JYI,NSS,One,J2I,NSS)
call DGEMM_('N','N',NSS,NSS,NSS,One,JYI,NSS,JYR,NSS,One,J2I,NSS)

call mma_deallocate(JXR)
call mma_deallocate(JXI)
call mma_deallocate(JYR)
call mma_deallocate(JYI)
call ZTRNSF(NSS,USOR,USOI,OMGR,OMGI)
call ZTRNSF(NSS,USOR,USOI,J2R,J2I)

! Jump here to skip computing omega and/or J:
!end if

if (IPGLOB >= 1) then
  write(u6,*)
  write(u6,'(6X,A)') ' Total energies including SO-coupling:'
  do ISS=1,NSS
    E_tmp = ENSOR(ISS)+EMIN
    call PrintResult(u6,'(6x,A,I5,5X,A,F23.14)','SO-RASSI State',ISS,'Total energy:',[E_tmp],1)
  end do
end if

! Find E0=lowest energy, to use for printing table:
if (IPGLOB >= 1) then
  E0 = minval(ENSOR(:))
  write(u6,*)
  write(u6,*)
  write(u6,*) '  Eigenvalues of complex Hamiltonian:'
  write(u6,*) '  -----------------------------------'
  if (EMIN /= Zero) write(u6,'(1X,A,F22.10,A1)') ' (Shifted by EMIN (a.u.) =',EMIN,')'
  write(u6,*)
  if ((ifj2 /= 0) .and. (ifjz /= 0)) then
    write(u6,*) 'SO State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)     J-value  Omega'
  else if ((ifj2 /= 0) .and. (ifjz == 0)) then
    write(u6,*) 'SO State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)     J-value'
  else if ((ifj2 == 0) .and. (ifjz /= 0)) then
    write(u6,*) 'SO State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)      Omega'
  else if ((ifj2 == 0) .and. (ifjz == 0)) then
    write(u6,*) 'SO State       Relative EMIN(au)   Rel lowest level(eV)    D:o, cm**(-1)'
  end if
  write(u6,*)
  E0 = minval(ENSOR(:))
  call MMA_ALLOCATE(ESO,NSS)
  do ISS=1,NSS
    E1 = ENSOR(ISS)
    E2 = auToeV*(E1-E0)
    E3 = auTocm*(E1-E0)
    if (IFJ2 > 0) XJEFF = sqrt(Quart+J2R(ISS,ISS))-Half
    if (IFJZ > 0) OMEGA = sqrt(1.0e-12_wp+OMGR(ISS,ISS))

    ! Added by Ungur Liviu on 04.11.2009
    ! Saving the SO energies in ESO array.
    ESO(ISS) = E3
    if ((ifj2 /= 0) .and. (ifjz /= 0)) then
      write(u6,'(1X,I5,7X,2(F18.10,2X),F18.4,4X,2(2X,F6.1))') ISS,E1,E2,E3,XJEFF,OMEGA
    else if ((ifj2 /= 0) .and. (ifjz == 0)) then
      write(u6,'(1X,I5,7X,2(F18.10,2X),F18.4,6X,F6.1)') ISS,E1,E2,E3,XJEFF
    else if ((ifj2 == 0) .and. (ifjz /= 0)) then
      write(u6,'(1X,I5,7X,2(F18.10,2X),F18.4,6X,F6.1)') ISS,E1,E2,E3,OMEGA
    else if ((ifj2 == 0) .and. (ifjz == 0)) then
      write(u6,'(1X,I5,7X,2(F18.10,2X),F18.4)') ISS,E1,E2,E3
    end if
  end do

  ! Added by Ungur Liviu on 04.11.2009
  ! Saving the ESO array in the RunFile.
  call Put_iscalar('NSS_SINGLE',NSS)
  call Put_dArray('ESO_SINGLE',ESO,NSS)
  call Put_dArray('ESO_LOW',ENSOR(:)+EMIN,NSS)
  call MMA_DEALLOCATE(ESO)
end if

call mma_deallocate(OMGR)
call mma_deallocate(OMGI)
call mma_deallocate(J2R)
call mma_deallocate(J2I)

! Put energy onto info file for automatic verification runs:
EPSS = 5.0e-11_wp
EPSH = max(5.0e-10_wp,abs(ENSOR(1)+EMIN)*EPSS)
IDX = 100
do ISS=1,NSS
  EI = (ENSOR(ISS)+EMIN)*EPSS
  ERMS = sqrt(EPSH**2+EI**2)*sum(USOR(:,ISS)**2+USOI(:,ISS)**2)
  IDX = min(IDX,int(-log10(ERMS)))
end do
iTol = cho_x_gettol(IDX) ! reset thr iff Cholesky
call Add_Info('ESO_LOW',ENSOR+EMIN,NSS,iTol)

if (IPGLOB >= 3) then
  write(u6,*)
  write(u6,*) ' Complex eigenvectors in basis of non-so eigenstates:'
  write(u6,*) '-----------------------------------------------------'
  write(u6,*)
  if (IPGLOB >= 4) then
    FRAC = Zero
  else
    FRAC = Quart
    write(u6,*) '    (A selection of the largest components)'
  end if
end if
call PRCEVC(NSS,FRAC,ENSOR,MAPST,MAPSP,MAPMS,USOR,USOI)

! Update LoopDivide (SUBSets keyword)
! Assume the SO "ground states" are mostly formed by the SF "ground states"
if (ReduceLoop) then
  call mma_Allocate(IndexE,nState,Label='IndexE')
  IndexE(:) = ArgSort(Energy,leq_r)
  n = 0
  do iState=1,LoopDivide
    Job = JbNum(IndexE(iState))
    n = n+Mltplt(Job)
  end do
  LoopDivide = n
# ifdef _HDF5_
  if (IFTRD1 .or. IFTDM) call UpdateIdx(IndexE,nSS,USOR,USOI,MapSt)
# endif
  call mma_deAllocate(IndexE)
end if

call mma_deallocate(MAPST)
call mma_deallocate(MAPSP)
call mma_deallocate(MAPMS)

end subroutine SOEIG
