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

subroutine PRPROP(PROP,USOR,USOI,ENSOR,NSS,OVLP,ENERGY,JBNUM,EigVec)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: SODYSAMPS
use kVectors, only: k_Vector, nk_Vector
use Cntrl, only: BAngRes, BIncre, BStart, DIPR, Do_SK, Do_TMom, DoCD, DYSO, EPrThr, iComp, IfACAL, IfGCAL, IfGTCALSA, IfGTSHSA, &
                 IfMCal, IFSO, IfvanVleck, IfXCal, IPUSED, ISOCMP, LoopDivide, LoopMax, LPRPR, MLTPLT, MULTIP, NBSTep, NPROP, &
                 NSOPR, NSTATE, NTS, nTStep, OSThr_DIPR, OSthr_QIPR, PNAME, PNUC, PORIG, PRDIPCOM, PRMEE, PRMES, PRXVE, PTYPE, &
                 QIALL, QIPR, ReduceLoop, RSPR, RSThr, SOPRNM, SOPRTP, TIncre, TMaxs, TMins, Tolerance, TStart
#ifdef _HDF5_
use mh5, only: mh5_put_dset
use RASSIWfn, only: wfn_sfs_amfi, wfn_sfs_angmom, wfn_sfs_edipmom, wfn_sos_angmomi, wfn_sos_angmomr, wfn_sos_dys, &
                    wfn_sos_edipmomi, wfn_sos_edipmomr, wfn_sos_spini, wfn_sos_spinr
use Cntrl, only: RhoDyn
#endif
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Six, Nine, Ten, Half, Quart, OneHalf, Pi, cZero, auTocm, auToeV, auTofs, auTokJ, &
                     auToT, c_in_au, Debye, deg2rad, gElectron, kBoltzmann, mBohr, rNAVO
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: NSS, JBNUM(NSTATE)
real(kind=wp) :: PROP(NSTATE,NSTATE,NPROP), USOR(NSS,NSS), USOI(NSS,NSS), ENSOR(NSS), OVLP(NSTATE,NSTATE), ENERGY(NSTATE), &
                 EigVec(NSTATE,NSTATE)
integer(kind=iwp) :: I, I2Tot, I_Have_DL, I_Have_DV, i_print, I_Print_Header, iAMFIx, iAMFIy, iAMFIz, iAMx, iAMy, iAMz, IBStep, &
                     IC, ICMP, IEND, iERR, IfAnyD, IfAnyM, IfAnyO, IfAnyQ, IfAnyS, iFinal, IFUNCT, ijXYZ, iMLTPL, IPAM(3), &
                     IPAMFI(3), iPhi, iPhiStep, iPrDXY, iPrDXZ, iPrDYZ, iProp, ISO, iSOPr, ISS, ISTA, iStart, iState, IT, ITHE, &
                     ITStep, iVec, iXYZ, J, JC, JEND, JSO, JSS, JSTART, jState, jXYZ, K, KDGN, kXYZ, LMStep, MPLET, MPLET1, &
                     MPLET2, nCol, nMiss, NORIENT, nPhiStep, NPMSIZ, nTheStep, nVec, SECORD(4)
real(kind=wp) :: A, AFactor, AX, AY, AZ, B, BFinal, bPhiRes, Bx, By, Bz, c_1(3,3), c_2(3,3), Chi(3), COMPARE, CONTRIB, curit(3,3), &
                 D_MXI, D_MXR, D_MYI, D_MYR, D_MZI, D_MZR, D_XI, D_XR, D_YI, D_YR, D_ZI, D_ZR, DELTA, DIPSOM_SA, DLT, DLT_E, DLTT, &
                 DTENS(3,3), DTIJ, DX2, DXX2, DXXDYY, DXXDZZ, DXXXDX, DXXYDY, DXXZDZ, DXY2, DXYDZ, DXZ2, DXZDY, DY2, DysThr, &
                 DYXDZ, DYY2, DYYDZZ, DYYXDX, DYYYDY, DYYZDZ, DYZ2, DYZDX, DZ2, DZXDY, DZYDX, DZZ2, DZZXDX, DZZYDY, DZZZDZ, EDiff, &
                 EDIFF2, EDIFF3, EEX, EEY, EEZ, EVI(3), EVR(3), F, Fact, FACT0, FACTM, Factor, FACTP, FX, FXX, FXXFYY, FXXFZZ, &
                 FXXX, FXXY, FXXZ, FXY, FXZ, FY, FYX, FYY, FYYFZZ, FYYX, FYYY, FYYZ, FYZ, FZ, FZX, FZY, FZZ, FZZX, FZZY, FZZZ, G, &
                 GSEnergy, GSTENS(3,3), GTENS(3,3), GTij, GTOTAL(9), HZer, OSthr, OSThr2, p_Boltz, paramt(3,3), Phi, PLIMIT, PMAX, &
                 Q_XXI, Q_XXR, Q_XYI, Q_XYR, Q_XZI, Q_XZR, Q_YYI, Q_YYR, Q_YZI, Q_YZR, Q_ZZI, Q_ZZR, R, RKT, RMAGM(3), rMagm2, &
                 rMagMO, RPart, Rtensor(6), RXX, RXXY, RXXZ, RXY, RXYX, RXYY, RXYZ, RXZ, RXZX, RXZY, RXZZ, RYX, RYY, RYYX, RYYZ, &
                 RYZ, RYZX, RYZY, RYZZ, RZX, RZY, RZZ, RZZX, RZZY, S, S1, S2, SOSTERM(9), T, TFinal, THE, ThreEJ, TMPm(NTS), &
                 TMPMAT(3,3), TMPVEC(3,3), Zstat
complex(kind=wp) :: T0(3), TM1
logical(kind=iwp) :: IFAMFI, IFANGM, IFDIP1
integer(kind=iwp), allocatable :: PMAP(:)
real(kind=wp), allocatable :: AMFIINT(:,:,:), ANGMOME(:,:,:), chicuriT_tens(:,:,:), chiparamT_tens(:,:,:), chiT_tens(:,:,:), &
                              DL(:,:), DV(:,:), DXI(:,:), DXR(:,:), DXXXI(:,:), DXXXR(:,:), DXXYI(:,:), DXXYR(:,:), DXXZI(:,:), &
                              DXXZR(:,:), DYI(:,:), DYR(:,:), DYYXI(:,:), DYYXR(:,:), DYYYI(:,:), DYYYR(:,:), DYYZI(:,:), &
                              DYYZR(:,:), DZI(:,:), DZR(:,:), DZZXI(:,:), DZZXR(:,:), DZZYI(:,:), DZZYR(:,:), DZZZI(:,:), &
                              DZZZR(:,:), EDIP1MOM(:,:,:), ESO(:), LXI(:,:), LYI(:,:), LZI(:,:), MAGM(:), MDXI(:,:), MDXR(:,:), &
                              MDYI(:,:), MDYR(:,:), MDZI(:,:), MDZR(:,:), MQXYI(:,:), MQXYR(:,:), MQXZI(:,:), MQXZR(:,:), &
                              MQYXI(:,:), MQYXR(:,:), MQYZI(:,:), MQYZR(:,:), MQZXI(:,:), MQZXR(:,:), MQZYI(:,:), MQZYR(:,:), &
                              MXYZI(:,:,:), MXYZR(:,:,:), QXXI(:,:), QXXR(:,:), QXYI(:,:), QXYR(:,:), QXZI(:,:), QXZR(:,:), &
                              QYYI(:,:), QYYR(:,:), QYZI(:,:), QYZR(:,:), QZZI(:,:), QZZR(:,:), SOPRI(:,:), SOPRR(:,:), SXI(:,:), &
                              SXR(:,:), SXYI(:,:), SXYR(:,:), SXZI(:,:), SXZR(:,:), SYI(:,:), SYR(:,:), SYXI(:,:), SYXR(:,:), &
                              SYZI(:,:), SYZR(:,:), SZI(:,:), SZR(:,:), SZXI(:,:), SZXR(:,:), SZYI(:,:), SZYR(:,:), TOT2K(:,:), &
                              UZR(:,:), UZI(:,:), ZI(:,:), ZR(:,:), ZXYZI(:,:,:), ZXYZR(:,:,:)
#ifdef _HDF5_
real(kind=wp), allocatable :: TMP(:,:,:)
#endif
complex(kind=wp), allocatable :: DIPSOm(:,:,:), DIPSOn(:,:,:), GCONT(:,:), SPNSFS(:,:,:), Z(:,:), ZEKL(:,:,:,:)
logical(kind=iwp), allocatable :: ISGS(:)
character, parameter :: xyzchr(3) = ['x','y','z']
real(kind=wp), parameter :: AU2J = auTokJ*1.0e3_wp, AU2JTM = (AU2J/auToT)*rNAVO, AU2REDR = 200.0_wp*Debye, &
                            BOLTZ = kBoltzmann/AU2J, BOLTZ_K = kBoltzmann*auTocm/AU2J, &
                            coeff_chi = 0.1_wp*rNAVO/kBoltzmann*mBohr**2, FEGVAL = -gElectron, J2CM = auTocm/AU2J, &
                            ONEOVER10C = One/(Ten*c_in_au**2), ONEOVER30C = ONEOVER10C/Three, ONEOVER6C2 = One/(Six*c_in_au**2), &
                            ONEOVER9C2 = One/(Nine*c_in_au**2), Rmu0 = 4.0e-7_wp*Pi, THRSH = 1.0e-10_wp, Two3rds = Two/Three, &
                            TWOOVERM45C = -Two/(45.0_wp*c_in_au**2)

! D[XYZ][RI]           Dipole
! MD[XYZ][RI]          Magnetic-Dipole
! Q[XYZ][XYZ][RI]      Electric-Quadrupole
! MQ[XYZ][XYZ][RI]     Magnetic-Quadrupole
! D[XYZ][XYZ][XYZ][RI] Octupole
! S[XYZ][RI]           Spin-Magnetic-Dipole
! S[XYZ][XYZ][RI]      Spin-Magnetic-Quadrupole

if (IPGLOB <= 0) goto 400

!*****************************************************
! printout of properties over the spin-free states
!*****************************************************
if (PRXVE .or. PRMEE) then
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,34X,A,34X,A)') '*',' Spin-free properties section ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
  write(u6,*)
end if

! Did the user want printed expectation values?
if (PRXVE) then
  call CollapseOutput(1,'Expectation values')
  write(u6,*)
  write(u6,*) ' ============================================'
  write(u6,*) '  EXPECTATION VALUES OF 1-ELECTRON OPERATORS'
  write(u6,*) '  FOR THE SPIN-FREE EIGENSTATES:'
  write(u6,*) ' ============================================'
  write(u6,*) ' (note: negative sign used for electronic multipoles)'
  write(u6,*)
  NCOL = 4
  do IPROP=1,NPROP
    if (IPUSED(IPROP) == 0) goto 100

    ! Skip printing if all the diagonal values are very small
    !  (presumed zero for reasons of selection rules)
    PLIMIT = 1.0e-10_wp
    PMAX = ZERO

    do I=1,NSTATE
      PMAX = max(PMAX,abs(PROP(I,I,IPROP)+PNUC(IPROP)*OVLP(I,I)))
    end do
    if (PMAX < PLIMIT) goto 100

    do ISTA=1,NSTATE,NCOL
      IEND = min(NSTATE,ISTA+NCOL-1)
      write(u6,*)
      write(u6,'(1X,A,A8,A,I4)') 'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
      write(u6,'(1X,A,3(1X,ES16.9))') 'ORIGIN    :',(PORIG(I,IPROP),I=1,3)
      write(u6,'(1X,A,I8,4I17)') 'STATE     :',(I,I=ISTA,IEND)
      write(u6,*)
      write(u6,'(1X,A,4(1X,ES16.9))') 'ELECTRONIC:',(PROP(I,I,IPROP),I=ISTA,IEND)
      write(u6,'(1X,A,4(1X,ES16.9))') 'NUCLEAR   :',(PNUC(IPROP),I=ISTA,IEND)
      write(u6,'(1X,A,4(1X,ES16.9))') 'TOTAL     :',(PROP(I,I,IPROP)+PNUC(IPROP),I=ISTA,IEND)
      write(u6,*)
    end do
100 continue
  end do
  call CollapseOutput(0,'Expectation values')
  write(u6,*)
end if

! include nuclear contribution
do IPROP=1,NPROP
  do I=1,NSTATE
    PROP(I,I,IPROP) = PROP(I,I,IPROP)+PNUC(IPROP)
  end do
end do

! Did the user want printed matrix elements?
if (PRMEE) then
  call CollapseOutput(1,'Matrix elements')
  write(u6,*)
  write(u6,*) ' ========================================='
  write(u6,*) '  MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
  write(u6,*) '  FOR THE SPIN-FREE EIGENSTATES:'
  write(u6,*) ' ========================================='
  write(u6,*) ' (including nuclear contrib.)'
  write(u6,*)
  write(u6,*) ' SELECTED PROPERTIES:'
  do I=1,NPROP,5
    write(u6,'(1X,5(A8,1X,I2,4X))') (PNAME(IPROP),ICOMP(IPROP),IPROP=I,min(NPROP,I+4))
  end do

  NCOL = 4
  do IPROP=1,NPROP
    if (IPUSED(IPROP) == 0) goto 200
    write(u6,*)
    write(u6,'(1X,A,A8,A,I4)') 'PROPERTY: ',PNAME(IPROP),'   COMPONENT:',ICOMP(IPROP)
    write(u6,'(1X,A,3(1X,ES16.9))') 'ORIGIN: ',(PORIG(I,IPROP),I=1,3)
    do ISTA=1,NSTATE,NCOL
      IEND = min(NSTATE,ISTA+NCOL-1)
      write(u6,*)
      write(u6,'(1X,A,I8,3I17)') ' STATE   ',(I,I=ISTA,IEND)
      write(u6,*)
      do J=1,NSTATE
        write(u6,'(1X,I4,6X,4(1X,ES16.9))') J,(PROP(J,I,IPROP),I=ISTA,IEND)
      end do
    end do
200 continue
  end do

  call CollapseOutput(0,'Matrix elements')
  write(u6,*)
end if
! Added by Ungur Liviu on 04.11.2009.
! Addition of ANGMOM to Runfile.

call mma_allocate(ANGMOME,3,NSTATE,NSTATE,Label='ANGMOME')
call mma_allocate(EDIP1MOM,3,NSTATE,NSTATE,Label='EDIP1MOM')
call mma_allocate(AMFIINT,3,NSTATE,NSTATE,Label='AMFIINT')

IFANGM = .false.
IFDIP1 = .false.
IFAMFI = .false.
ANGMOME(:,:,:) = Zero
EDIP1MOM(:,:,:) = Zero
AMFIINT(:,:,:) = Zero
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    IFANGM = .true.
    do I=1,NSTATE
      do J=1,NSTATE
        ANGMOME(ICOMP(IPROP),I,J) = PROP(I,J,IPROP)
      end do
    end do
  end if
  ! add dipole moment integrals:
  if ((PNAME(IPROP) == 'MLTPL  1') .and. (PTYPE(IPROP) == 'HERMSING')) then
    IFDIP1 = .true.
    do I=1,NSTATE
      do J=1,NSTATE
        EDIP1MOM(ICOMP(IPROP),I,J) = PROP(I,J,IPROP)
      end do
    end do
  end if
  ! add spin-orbit AMFI integrals:
  if (PNAME(IPROP)(1:8) == 'AMFI    ') then
    IFAMFI = .true.
    do I=1,NSTATE
      do J=1,NSTATE
        AMFIINT(ICOMP(IPROP),I,J) = PROP(I,J,IPROP)
      end do
    end do
  end if
end do
if (IFANGM) call Put_dArray('ANGM_SINGLE',ANGMOME,3*NSTATE*NSTATE)
if (IFDIP1) call Put_dArray('DIP1_SINGLE',EDIP1MOM,3*NSTATE*NSTATE)
if (IFAMFI) call Put_dArray('AMFI_SINGLE',AMFIINT,3*NSTATE*NSTATE)
#ifdef _HDF5_
call mma_allocate(TMP,NSTATE,NSTATE,3,Label='TMP')
if (IFANGM) then
  do i=1,3
    TMP(:,:,i) = ANGMOME(i,:,:)
  end do
  call mh5_put_dset(wfn_sfs_angmom,TMP,[NSTATE,NSTATE,3],[0,0,0])
end if
if (IFDIP1) then
  do i=1,3
    TMP(:,:,i) = EDIP1MOM(i,:,:)
  end do
  call mh5_put_dset(wfn_sfs_edipmom,TMP,[NSTATE,NSTATE,3],[0,0,0])
end if
if (IFAMFI) then
  do i=1,3
    TMP(:,:,i) = AMFIINT(i,:,:)
  end do
  call mh5_put_dset(wfn_sfs_amfi,TMP,[NSTATE,NSTATE,3],[0,0,0])
end if
call mma_deallocate(TMP)
#endif
call mma_deallocate(ANGMOME)
call mma_deallocate(EDIP1MOM)
call mma_deallocate(AMFIINT)
!******************************************************
! printout of properties over the spin-orbit states
!******************************************************
! If PRPR requested, print the spin matrices
#ifdef _HDF5_
if (LPRPR .or. PRMES) then
#else
if (LPRPR) then
#endif
  call mma_Allocate(SOPRR,NSS,NSS,Label='SOPRR')
  call mma_Allocate(SOPRI,NSS,NSS,Label='SOPRI')
  do ICMP=1,3
    SOPRR(:,:) = Zero
    SOPRI(:,:) = Zero
    if (ICMP == 2) then
      call SMMAT(PROP,SOPRI,NSS,0,ICMP)
    else
      call SMMAT(PROP,SOPRR,NSS,0,ICMP)
    end if
    call ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
#   ifdef _HDF5_
    call mh5_put_dset(wfn_sos_spinr,SOPRR,[NSS,NSS,1],[0,0,ICMP-1])
    call mh5_put_dset(wfn_sos_spini,SOPRI,[NSS,NSS,1],[0,0,ICMP-1])
#   endif
    if (LPRPR) call PRCMAT3(NSS,SOPRR,SOPRI,ICMP)
  end do
  call mma_deallocate(SOPRR)
  call mma_deallocate(SOPRI)
end if

if (.not. IFSO) goto 300
NPMSIZ = NSOPR
if (NSOPR == 0) goto 300

if (PRMES) then
  ! match the SO property list to the SF property list
  call mma_allocate(PMAP,NPMSIZ,Label='PMap')
  NMISS = 0
  do ISOPR=1,NSOPR
    PMAP(ISOPR) = 0
    do IPROP=1,NPROP
      if ((PNAME(IPROP) == SOPRNM(ISOPR)) .and. (ICOMP(IPROP) == ISOCMP(ISOPR))) then
        PMAP(ISOPR) = IPROP
        goto 10
      end if
    end do
    NMISS = NMISS+1
10  continue
  end do

  ! check for inconsistencies
  if (NMISS > 0) then
    call WarningMessage(1,'Missing data integrals.')
    write(u6,*) 'WARNING: You have requested matrix elements'
    write(u6,*) 'over spin states of some operators. The present'
    write(u6,*) 'code uses matrix elements computed over spin-free'
    write(u6,*) 'states to compute those over spin states.'
    write(u6,*) 'Matrix elements of the following operator(s)'
    write(u6,*) 'were never computed and must be skipped.'
    write(u6,*) '   (If you need these properties, change the'
    write(u6,*) '    input to SEWARD and recompute.)'
    do ISOPR=1,NSOPR
      if (PMAP(ISOPR) == 0) write(u6,*) 'Property:',SOPRNM(ISOPR),'      Component:',ISOCMP(ISOPR)
    end do
    write(u6,*)
  end if
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,34X,A,34X,A)') '*','Spin-orbit properties section ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
  write(u6,*)

  call CollapseOutput(1,'Matrix elements over SO states')
  write(u6,*)
  write(u6,*) ' ========================================='
  write(u6,*) '  MATRIX ELEMENTS OF 1-ELECTRON OPERATORS'
  write(u6,*) '  FOR THE SPIN-ORBIT EIGENSTATES:'
  write(u6,*) ' ========================================='
  write(u6,*)
  write(u6,*) ' SELECTED PROPERTIES:'
  do I=1,NPROP,5
    write(u6,'(1X,5(A8,1X,I2,4X))') (SOPRNM(ISOPR),ISOCMP(ISOPR),ISOPR=I,min(NSOPR,I+4))
  end do

  ! Remove zeroes to make SOPRNM and ISOCMP lists contiguous. New NSOPR.
  ISOPR = 0
  do I=1,NSOPR
    IPROP = PMAP(I)
    if (IPROP > 0) then
      ISOPR = ISOPR+1
      SOPRNM(ISOPR) = SOPRNM(I)
      ISOCMP(ISOPR) = ISOCMP(I)
    end if
  end do
  call mma_deallocate(PMAP)
  NSOPR = ISOPR

  call mma_Allocate(SOPRR,NSS,NSS,Label='SOPRR')
  call mma_Allocate(SOPRI,NSS,NSS,Label='SOPRI')
  ! Print out the matrix elements:
  NCOL = 4
  do ISOPR=1,NSOPR
    write(u6,*)
    write(u6,'(1X,A,A8,A,I4)') 'PROPERTY: ',SOPRNM(ISOPR),'   COMPONENT:',ISOCMP(ISOPR)
    !IFG  should print the origin, but where is it stored (for SO properties)?
    SOPRR(:,:) = Zero
    SOPRI(:,:) = Zero

    call SMMAT(PROP,SOPRR,NSS,ISOPR,0)
    call ZTRNSF(NSS,USOR,USOI,SOPRR,SOPRI)
    call PRCMAT(NSS,SOPRR,SOPRI)
    ! prpr keyword: Print selected spin-orbit properties to ext. data files
    if (LPRPR .and. (SOPRNM(ISOPR)(1:5) == 'MLTPL')) then
      if (SOPRTP(ISOPR) == 'HERMSING') call PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
    else if (LPRPR .and. (SOPRNM(ISOPR)(1:6) == 'ANGMOM')) then
      if (SOPRTP(ISOPR) == 'ANTISING') call PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
    else if (LPRPR .and. (SOPRNM(ISOPR)(1:8) == 'VELOCITY')) then
      if (SOPRTP(ISOPR) == 'ANTISING') call PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
    else if (LPRPR .and. (SOPRNM(ISOPR)(1:5) == 'MLTPV')) then
      if (SOPRTP(ISOPR) == 'ANTISING') call PRCMAT2(ISOPR,NSS,SOPRR,SOPRI)
    end if
    ! prpr end

#   ifdef _HDF5_
    if (SOPRNM(ISOPR)(1:6) == 'ANGMOM') then
      call mh5_put_dset(wfn_sos_angmomr,SOPRR,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
      call mh5_put_dset(wfn_sos_angmomi,SOPRI,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
    end if

    if ((SOPRNM(ISOPR)(1:8) == 'MLTPL  1') .and. (SOPRTP(ISOPR) == 'HERMSING')) then
      call mh5_put_dset(wfn_sos_edipmomr,SOPRR,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
      call mh5_put_dset(wfn_sos_edipmomi,SOPRI,[NSS,NSS,1],[0,0,ISOCMP(ISOPR)-1])
    end if
#   endif

  end do
  call mma_deallocate(SOPRR)
  call mma_deallocate(SOPRI)
  call CollapseOutput(0,'Matrix elements over SO states')
  write(u6,*)

end if

300 continue

400 continue

!*****************************************************
! printout of special properties
!*****************************************************

! AFACTOR = 2*pi*e^2*E_h^2 / eps_0*m_e*c^3*h^2
! numerically: 2/c^3 (in a.u. of time ^ -1)
AFACTOR = Two/c_in_au**3/(auTofs*1.0e-15_wp)

if (IPGLOB >= 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A,34X,A,34X,A)') '*','  Special properties section  ','*'
  write(u6,'(6X,A,98X,A)') '*','*'
  write(u6,'(6X,A)') repeat('*',100)
  write(u6,*)
  write(u6,*)
end if

! Compute transition strengths for spin-orbit states:
if (.not. IFSO) goto 500

! Initial setup for both dipole, quadrupole etc. and exact operator

! printing threshold
!if (IPGLOB <= 2) then
!  OSTHR = 1.0e-8_wp ! first order
!  OSTHR2 = 1.0e-12_wp ! second order (weaker)
!else
!  OSTHR = Zero
!  OSTHR2 = Zero
!end if
OSTHR = 1.0e-5_wp
OSTHR2 = 1.0e-5_wp
if (DIPR) OSTHR = OSTHR_DIPR
if (DIPR) write(u6,30) 'Dipole printing threshold changed to ',OSTHR
! Again to avoid total negative transition strengths
if (QIPR) then
  OSTHR = OSTHR_QIPR
  write(u6,49) 'Printing threshold changed to ',OSTHR,'since quadrupole threshold is given'
  OSTHR2 = OSTHR_QIPR
  write(u6,30) 'Quadrupole printing threshold changed to ',OSTHR2
end if
if (QIALL) write(u6,*) ' Will write all quadrupole contributions'

! Rotatory strength threshold
if (RSPR) then
  write(u6,30) 'Rotatory strength printing threshold changed to ',RSTHR
else
  RSTHR = 1.0e-7_wp !Default
end if

! Reducing the loop over states - good for X-rays
! At the moment memory is not reduced

JEND = NSS
if (REDUCELOOP) then
  IEND = LOOPDIVIDE
  JSTART = LOOPDIVIDE+1
  if (LOOPMAX > 0) JEND = min(NSS,LOOPDIVIDE+LOOPMAX)
else
  IEND = NSS
  JSTART = 1
end if

if (IPGLOB >= 1) then

  ! Initialize arrays for indentifying problematic transitions
  ! These stores all dipole oscillator strengths in
  ! length and velocity gauge for a later comparison.

  call mma_allocate(DL,NSS,NSS,Label='DL')
  call mma_allocate(DV,NSS,NSS,Label='DV')
  DL(:,:) = Zero
  DV(:,:) = Zero
  I_HAVE_DL = 0
  I_HAVE_DV = 0

  ! Electric-Dipole Electric-Dipole transitions

  if (Do_SK) then
    nVec = nk_Vector
  else
    nVec = 1
  end if

  call Allocate_and_Load_electric_dipoles(IFANYD)

  if (IFANYD /= 0) then

    do iVec=1,nVec

      i_Print = 0

      do ISS=1,IEND
        do JSS=JSTART,JEND
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (abs(EDIFF) <= 1.0e-8_wp) cycle
          if (EDIFF > Zero) then
            T0(1) = cmplx(DXR(JSS,ISS),DXI(JSS,ISS),kind=wp)
            T0(2) = cmplx(DYR(JSS,ISS),DYI(JSS,ISS),kind=wp)
            T0(3) = cmplx(DZR(JSS,ISS),DZI(JSS,ISS),kind=wp)
            if (Do_SK) then
              TM1 = k_vector(1,iVec)*T0(1)+k_vector(2,iVec)*T0(2)+k_vector(3,iVec)*T0(3)
              T0(1) = T0(1)-TM1*k_vector(1,iVec)
              T0(2) = T0(2)-TM1*k_vector(2,iVec)
              T0(3) = T0(3)-TM1*k_vector(3,iVec)
            end if
            DX2 = abs(conjg(T0(1))*T0(1))
            DY2 = abs(conjg(T0(2))*T0(2))
            DZ2 = abs(conjg(T0(3))*T0(3))
            FX = Two3rds*EDIFF*(DX2)
            FY = Two3rds*EDIFF*(DY2)
            FZ = Two3rds*EDIFF*(DZ2)
            F = FX+FY+FZ
            AX = (AFACTOR*EDIFF**2)*FX
            AY = (AFACTOR*EDIFF**2)*FY
            AZ = (AFACTOR*EDIFF**2)*FZ
            A = (AFACTOR*EDIFF**2)*F
            ! Store dipole oscillator strength
            DL(JSS,ISS) = F
            if (abs(F) >= OSTHR) then
              if (i_Print == 0) then
                i_Print = 1

                ! Print full COMPLEX transition dipole moment vectors?
                ! J. Norell 7/5 2020
                if (PRDIPCOM) then

                  call CollapseOutput(1,'Complex transition dipole vectors (SO states):')
                  write(u6,'(3X,A)') '----------------------------------------'
                  if (OSTHR > Zero) then
                    write(u6,30) '   for osc. strength at least ',OSTHR
                    write(u6,*)
                  end if
                  write(u6,*) '     From   To','       Re(Dx)       Im(Dx)','       Re(Dy)       Im(Dy)', &
                              '       Re(Dz)       Im(Dz)'
                  write(u6,32)
                  goto 137 ! Skip past "regular" print
                end if
                ! END print COMPLEX vectors

                call CollapseOutput(1,'Dipole transition strengths (SO states):')
                write(u6,'(3X,A)') '----------------------------------------'
                if (OSTHR > Zero) then
                  write(u6,30) '   for osc. strength at least ',OSTHR
                  write(u6,*)
                end if
                if (Do_SK) then
                  write(u6,*)
                  write(u6,'(4x,a,3F10.6,a)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3),' (a.u.)'
                  write(u6,'(4x,a)') 'The light is assumed to be unpolarized.'
                  write(u6,*)
                end if
                write(u6,31) 'From','To','Osc. strength','Einstein coefficients Ax, Ay, Az (sec-1)   ','Total A (sec-1)'
                write(u6,32)
              end if

              ! Print full COMPLEX transition dipole moment vectors?
137           if (PRDIPCOM) then
                write(u6,'(5X,I5,I5,A,A,ES12.3,A,ES12.3,A,ES12.3,A,ES12.3,A,ES12.3,A,ES12.3)') ISS,JSS,'   ',' ',real(T0(1)),' ', &
                                                                                               aimag(T0(1)),' ',real(T0(2)),' ', &
                                                                                               aimag(T0(2)),' ',real(T0(3)),' ', &
                                                                                               aimag(T0(3))
              else
                write(u6,33) ISS,JSS,F,AX,AY,AZ,A ! "Regular" print instead
              end if
              ! END print COMPLEX vectors

            end if
            call Add_Info('TMS(SO,Len)',[F],1,6)
            if (PRDIPCOM) call Add_Info('TVC(SO,Len)',[DX2+DY2+DZ2],1,6)

          end if
        end do
      end do

      if (i_Print == 1) then
        write(u6,32)
        call CollapseOutput(0,'Dipole transition strengths (SO states):')
        write(u6,*)
      end if

    end do ! iVec

    I_HAVE_DL = 1
  end if

  call Deallocate_electric_dipoles()

  ! Now the same in velocity representation

  if (Do_SK) then
    nVec = nk_Vector
  else
    nVec = 1
  end if

  call Allocate_and_Load_velocities(IFANYD)

  if (IFANYD /= 0) then

    do iVec=1,nVec

      i_Print = 0

      do ISS=1,IEND
        do JSS=JSTART,JEND
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (abs(EDIFF) <= 1.0e-8_wp) cycle
          if (EDIFF > Zero) then
            T0(1) = cmplx(DXR(JSS,ISS),DXI(JSS,ISS),kind=wp)
            T0(2) = cmplx(DYR(JSS,ISS),DYI(JSS,ISS),kind=wp)
            T0(3) = cmplx(DZR(JSS,ISS),DZI(JSS,ISS),kind=wp)
            if (Do_SK) then
              TM1 = k_vector(1,iVec)*T0(1)+k_vector(2,iVec)*T0(2)+k_vector(3,iVec)*T0(3)
              T0(1) = T0(1)-TM1*k_vector(1,iVec)
              T0(2) = T0(2)-TM1*k_vector(2,iVec)
              T0(3) = T0(3)-TM1*k_vector(3,iVec)
            end if
            DX2 = abs(conjg(T0(1))*T0(1))
            DY2 = abs(conjg(T0(2))*T0(2))
            DZ2 = abs(conjg(T0(3))*T0(3))
            FX = Two3rds*(DX2)/EDIFF
            FY = Two3rds*(DY2)/EDIFF
            FZ = Two3rds*(DZ2)/EDIFF
            F = FX+FY+FZ
            AX = (AFACTOR*EDIFF**2)*FX
            AY = (AFACTOR*EDIFF**2)*FY
            AZ = (AFACTOR*EDIFF**2)*FZ
            A = (AFACTOR*EDIFF**2)*F
            ! Store dipole oscillator strength
            DV(JSS,ISS) = F
            if (abs(F) >= OSTHR) then
              if (i_Print == 0) then
                i_Print = 1
                call CollapseOutput(1,'Velocity transition strengths (SO states):')
                write(u6,'(3X,A)') '------------------------------------------'
                if (OSTHR > Zero) then
                  write(u6,30) '   for osc. strength at least ',OSTHR
                  write(u6,*)
                end if
                if (Do_SK) then
                  write(u6,*)
                  write(u6,'(4x,a,3F10.6,a)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3),' (a.u.)'
                  write(u6,'(4x,a)') 'The light is assumed to be unpolarized.'
                  write(u6,*)
                end if
                write(u6,31) 'From','To','Osc. strength','Einstein coefficients Ax, Ay, Az (sec-1)   ','Total A (sec-1)'
                write(u6,32)
              end if
              write(u6,33) ISS,JSS,F,AX,AY,AZ,A
            end if
            call Add_Info('TMS(SO,Vel)',[F],1,6)
          end if
        end do
      end do

      if (i_Print == 1) then
        write(u6,32)
        call CollapseOutput(0,'Velocity transition strengths (SO states):')
        write(u6,*)
      end if

    end do ! iVec

    I_HAVE_DV = 1
  end if

  call Deallocate_electric_dipoles()

  ! Compare oscillator strengths in length and velocity gauge
  ! All differences in oscillator strengths above the tolerance
  ! of 0.1 (10 percent) will be printed.

  if ((I_HAVE_DL == 1) .and. (I_HAVE_DV == 1)) then
    call CollapseOutput(1,'Length and velocity gauge comparison (SO states):')

    ! I guess that I have to explain it when I print a warning

    write(u6,*)
    write(u6,*) '--------------------------------------------------'
    write(u6,*) 'A comparison between the dipole oscillator strengths in'
    write(u6,*) 'length and velocity gauge will be performed'
    write(u6,*)
    write(u6,49) 'All dipole oscillator differences above the tolerance of ',TOLERANCE,' will be printed'
    write(u6,*)
    write(u6,*) 'Due to basis set deficiency these oscillator may be problematic'
    write(u6,*)
    write(u6,*) 'The tolerance is defined as ABS(1-O_l/O_v)'
    write(u6,*) 'O_l : dipole oscillator strength in length gauge'
    write(u6,*) 'O_p : dipole oscillator strength in velocity gauge'
    write(u6,*) '--------------------------------------------------'

    I_PRINT_HEADER = 0
    do I=1,IEND
      do J=JSTART,JEND
        EDIFF = ENSOR(J)-ENSOR(I)
        if (EDIFF < Zero) cycle
        COMPARE = Zero
        dlt = 1.0e-18_wp ! Add small value to avoid zero divide.
        if ((DL(J,I) >= OSTHR+dlt) .and. (DV(J,I) >= OSTHR+dlt)) then
          COMPARE = abs(1-DL(J,I)/DV(J,I))
        else if ((DL(J,I) >= OSTHR+dlt) .and. (DL(J,I) > Zero)) then
          COMPARE = -OneHalf
        else if ((DV(J,I) >= OSTHR+dlt) .and. (DV(J,I) > Zero)) then
          COMPARE = -2.5_wp
        end if
        if (abs(COMPARE) >= TOLERANCE) then
          I_PRINT_HEADER = I_PRINT_HEADER+1
          if (I_PRINT_HEADER == 1) then
            write(u6,*)
            write(u6,*) ' Problematic transitions have been found'
            write(u6,*)
            write(u6,39) 'From','To','Difference (%)','Osc. st. (len.)','Osc. st. (vel.)'
            write(u6,40)
          end if
          if (COMPARE >= Zero) then
            write(u6,38) I,J,COMPARE*100.0_wp,DL(J,I),DV(J,I)
          else if (COMPARE >= -Two) then
            write(u6,36) I,J,DL(J,I),'below threshold'
          else
            write(u6,37) I,J,'below threshold',DV(J,I)
          end if
        end if
      end do
    end do
    if (I_PRINT_HEADER == 0) then
      write(u6,*)
      write(u6,*) 'No problematic oscillator strengths above the tolerance ',TOLERANCE,' have been found'
      write(u6,*)
    else
      write(u6,40)
      write(u6,*)
      write(u6,*) 'Number of problematic transitions = ',I_PRINT_HEADER
      write(u6,*)
    end if
    call CollapseOutput(0,'Length and velocity gauge comparison (SO states):')
    write(u6,*)
  end if

  ! Free the memory

  call mma_deallocate(DL)
  call mma_deallocate(DV)

  ! We will first allocate a matrix for the total of the second order wave vector

  call mma_allocate(TOT2K,NSS,NSS,Label='TOT2K')
  TOT2K(:,:) = Zero

  ! Checking if all are in
  SECORD = 0

  ! Magnetic-Dipole - Magnetic-Dipole transitions and
  ! Spin-Magnetic-Dipole - Spin-Magnetic-Dipole transitions

  ! I will not separate these for SO states since there would then be
  ! M^2 + Ms^2 + 2*MMs (three terms to be programmed)
  ! M^2 and Ms^2 can be calculated separately but the cross term not directly

  ! Magnetic-Dipole
  call Allocate_and_Load_Magnetic_Dipoles(IFANYM)
  ! Spin-Magnetic-Dipole ---- notice the S
  call Allocate_and_Load_Spin_Magnetic_dipoles(IFANYS)

  if ((IFANYM /= 0) .or. (IFANYS /= 0)) then

    ! Only print the part calculated

    if (QIALL) then
      if ((IFANYM /= 0) .and. (IFANYS /= 0)) then
        call CollapseOutput(1,'Magnetic-dipole - magnetic-dipole and spin-magnetic-dipole - spin-magnetic-dipole transition '// &
                            'strengths (SO states):')
        write(u6,'(3X,A)') '---------------------------------------------------------------------------------------------'// &
                           '----------------------'
      else if ((IFANYM /= 0) .and. (IFANYS == 0)) then
        call CollapseOutput(1,'Magnetic-dipole - magnetic-dipole transition strengths (SO states):')
        write(u6,'(3X,A)') '-------------------------------------------------------------------'
      else if ((IFANYM == 0) .and. (IFANYS /= 0)) then
        call CollapseOutput(1,'Spin-magnetic-dipole - spin-magnetic-dipole transition strengths (SO states):')
        write(u6,'(3X,A)') '-----------------------------------------------------------------------------'
      end if
      if (OSTHR2 > Zero) then
        write(u6,30) '   for osc. strength at least ',OSTHR2
        write(u6,*)
      end if
      write(u6,31) 'From','To','Osc. strength'
      write(u6,35)
    end if

    ! Magnetic-Dipole

    ! Spin-Magnetic-Dipole

    g = FEGVAL
    do ISS=1,IEND
      do JSS=JSTART,NSS
        EDIFF = ENSOR(JSS)-ENSOR(ISS)
        if (abs(EDIFF) < 1.0e-8_wp) cycle
        if (EDIFF > Zero) then

          DX2 = (MDXI(JSS,ISS)+g*SXR(JSS,ISS))**2+(MDXR(JSS,ISS)-g*SXI(JSS,ISS))**2
          DY2 = (MDYI(JSS,ISS)+g*SYR(JSS,ISS))**2+(MDYR(JSS,ISS)-g*SYI(JSS,ISS))**2
          DZ2 = (MDZI(JSS,ISS)+g*SZR(JSS,ISS))**2+(MDZR(JSS,ISS)-g*SZI(JSS,ISS))**2

          F = (DX2+DY2+DZ2)*EDIFF*ONEOVER6C2
          ! Add it to the total
          TOT2K(JSS,ISS) = TOT2K(JSS,ISS)+F
          if (abs(F) >= OSTHR2) then
            if (QIALL) write(u6,33) ISS,JSS,F
          end if
        end if
      end do
    end do

    if (QIALL) then
      write(u6,35)
      if ((IFANYM /= 0) .and. (IFANYS /= 0)) then
        call CollapseOutput(0,'Magnetic-dipole - magnetic-dipole and spin-magnetic-dipole - spin-magnetic-dipole transition '// &
                            'strengths (SO states):')
        write(u6,*)
      else if ((IFANYM /= 0) .and. (IFANYS == 0)) then
        call CollapseOutput(0,'Magnetic-dipole - magnetic-dipole transition strengths (SO states):')
        write(u6,*)
      else if ((IFANYM == 0) .and. (IFANYS /= 0)) then
        call CollapseOutput(0,'Spin-magnetic-dipole - Spin-magnetic-dipole transition strengths (SO states):')
        write(u6,*)
      end if
    end if
    SECORD(1) = 1
  end if

  ! Magnetic-Dipole
  call Deallocate_Magnetic_dipoles()

  ! Spin-Magnetic-Dipole
  call Deallocate_Spin_Magnetic_dipoles()

  ! Electric-Quadrupole Electric-Quadrupole transitions

  call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)

  if (IFANYQ /= 0) then
    if (QIALL) then
      call CollapseOutput(1,'Quadrupole transition strengths (SO states):')
      write(u6,'(3X,A)') '--------------------------------------------'
      if (OSTHR2 > Zero) then
        write(u6,30) '   for osc. strength at least ',OSTHR2
        write(u6,*)
      end if
      write(u6,31) 'From','To','Osc. strength'
      write(u6,35)
    end if

    do ISS=1,IEND
      do JSS=JSTART,NSS
        EDIFF = ENSOR(JSS)-ENSOR(ISS)
        if (abs(EDIFF) < 1.0e-8_wp) cycle
        if (EDIFF > Zero) then

          ! D should be purely real since D is a real symmetric matrix

          EDIFF3 = EDIFF**3

          DXX2 = QXXR(JSS,ISS)**2+QXXI(JSS,ISS)**2
          DYY2 = QYYR(JSS,ISS)**2+QYYI(JSS,ISS)**2
          DZZ2 = QZZR(JSS,ISS)**2+QZZI(JSS,ISS)**2
          FXX = ONEOVER30C*EDIFF3*(DXX2)
          FYY = ONEOVER30C*EDIFF3*(DYY2)
          FZZ = ONEOVER30C*EDIFF3*(DZZ2)

          DXY2 = QXYR(JSS,ISS)**2+QXYI(JSS,ISS)**2
          DXZ2 = QXZR(JSS,ISS)**2+QXZI(JSS,ISS)**2
          DYZ2 = QYZR(JSS,ISS)**2+QYZI(JSS,ISS)**2
          FXY = ONEOVER10C*EDIFF3*(DXY2)
          FXZ = ONEOVER10C*EDIFF3*(DXZ2)
          FYZ = ONEOVER10C*EDIFF3*(DYZ2)

          DXXDYY = QXXR(JSS,ISS)*QYYR(JSS,ISS)+QXXI(JSS,ISS)*QYYI(JSS,ISS)
          DXXDZZ = QXXR(JSS,ISS)*QZZR(JSS,ISS)+QXXI(JSS,ISS)*QZZI(JSS,ISS)
          DYYDZZ = QYYR(JSS,ISS)*QZZR(JSS,ISS)+QYYI(JSS,ISS)*QZZI(JSS,ISS)
          FXXFYY = -ONEOVER30C*EDIFF3*(DXXDYY)
          FXXFZZ = -ONEOVER30C*EDIFF3*(DXXDZZ)
          FYYFZZ = -ONEOVER30C*EDIFF3*(DYYDZZ)

          F = FXX+FXY+FXZ+FYY+FYZ+FZZ+FXXFYY+FXXFZZ+FYYFZZ
          ! Add it to the total
          TOT2K(JSS,ISS) = TOT2K(JSS,ISS)+F

          if (abs(F) >= OSTHR2) then
            if (QIALL) write(u6,33) ISS,JSS,F
          end if
        end if
      end do
    end do

    if (QIALL) then
      write(u6,35)
      call CollapseOutput(0,'Quadrupole transition strengths (SO states):')
      write(u6,*)
    end if
    SECORD(2) = 1
  end if

  call Deallocate_Electric_Quadrupoles()

  ! Electric-Dipole Electric-Octupole transitions

  ! Octupole
  call Allocate_and_Load_Octupoles(IFANYO)
  ! Dipole
  call Allocate_and_Load_electric_dipoles(IFANYD)

  if ((IFANYD /= 0) .and. (IFANYO /= 0)) then
    if (QIALL) then
      call CollapseOutput(1,'Electric-dipole - electric-octupole transition strengths (SO states):')
      write(u6,'(3X,A)') '---------------------------------------------------------------------'
      if (OSTHR2 > Zero) then
        write(u6,30) '   for osc. strength at least ',OSTHR2
        write(u6,*)
      end if
      write(u6,31) 'From','To','Osc. strength'
      write(u6,35)
    end if

    do ISS=1,IEND
      do JSS=JSTART,NSS
        EDIFF = ENSOR(JSS)-ENSOR(ISS)
        if (abs(EDIFF) < 1.0e-8_wp) cycle
        if (EDIFF > Zero) then

          EDIFF3 = EDIFF**3

          DXXXDX = DXXXR(JSS,ISS)*DXR(JSS,ISS)+DXXXI(JSS,ISS)*DXI(JSS,ISS)
          DYYXDX = DYYXR(JSS,ISS)*DXR(JSS,ISS)+DYYXI(JSS,ISS)*DXI(JSS,ISS)
          DZZXDX = DZZXR(JSS,ISS)*DXR(JSS,ISS)+DZZXI(JSS,ISS)*DXI(JSS,ISS)
          FXXX = TWOOVERM45C*EDIFF3*(DXXXDX)
          FYYX = TWOOVERM45C*EDIFF3*(DYYXDX)
          FZZX = TWOOVERM45C*EDIFF3*(DZZXDX)

          DXXYDY = DXXYR(JSS,ISS)*DYR(JSS,ISS)+DXXYI(JSS,ISS)*DYI(JSS,ISS)
          DYYYDY = DYYYR(JSS,ISS)*DYR(JSS,ISS)+DYYYI(JSS,ISS)*DYI(JSS,ISS)
          DZZYDY = DZZYR(JSS,ISS)*DYR(JSS,ISS)+DZZYI(JSS,ISS)*DYI(JSS,ISS)
          FXXY = TWOOVERM45C*EDIFF3*(DXXYDY)
          FYYY = TWOOVERM45C*EDIFF3*(DYYYDY)
          FZZY = TWOOVERM45C*EDIFF3*(DZZYDY)

          DXXZDZ = DXXZR(JSS,ISS)*DZR(JSS,ISS)+DXXZI(JSS,ISS)*DZI(JSS,ISS)
          DYYZDZ = DYYZR(JSS,ISS)*DZR(JSS,ISS)+DYYZI(JSS,ISS)*DZI(JSS,ISS)
          DZZZDZ = DZZZR(JSS,ISS)*DZR(JSS,ISS)+DZZZI(JSS,ISS)*DZI(JSS,ISS)
          FXXZ = TWOOVERM45C*EDIFF3*(DXXZDZ)
          FYYZ = TWOOVERM45C*EDIFF3*(DYYZDZ)
          FZZZ = TWOOVERM45C*EDIFF3*(DZZZDZ)

          F = FXXX+FYYX+FZZX+FXXY+FYYY+FZZY+FXXZ+FYYZ+FZZZ
          ! Add it to the total
          TOT2K(JSS,ISS) = TOT2K(JSS,ISS)+F

          if (abs(F) >= OSTHR2) then
            if (QIALL) write(u6,33) ISS,JSS,F
          end if
        end if
      end do
    end do

    if (QIALL) then
      write(u6,35)
      call CollapseOutput(0,'Electric-dipole - electric-octupole transition strengths (SO states):')
      write(u6,*)
    end if
    SECORD(3) = 1
  end if

  call Deallocate_Octupoles()

  call Deallocate_electric_dipoles()

  ! Electric-Dipole - Magnetic-Quadrupole transitions and
  ! Electric-Dipole - Spin-Magnetic-Quadrupole transitions

  ! Again I will just include the spin-term so both terms are calculated
  ! (Can also be done separately)
  ! DM + DMs

  ! Magnetic-Quadrupole
  ! Spin-Magnetic-Quadrupole
  ! Spin-Magnetic-Quadrupole = M^s_ab = r_b * s_a

  ! Electric-Dipole
  call Allocate_and_Load_electric_dipoles(IFANYD)
  if (IFANYD /= 0) then
    ! Magnetic-Quadrupole
    call Allocate_and_Load_Magnetic_Quadrupoles(IFANYD)
    ! Spin-Magnetic-Quadrupole
    call Allocate_and_Load_Spin_Magnetic_Quadrupoles(IFANYS)

    if ((IFANYD /= 0) .or. (IFANYS /= 0)) then
      if (QIALL) then
        if ((IFANYD /= 0) .and. (IFANYS /= 0)) then
          call CollapseOutput(1,'Electric-dipole - magnetic-quadrupole and electric-dipole - spin-magnetic-quadrupole '// &
                              'transition strengths (SO states):')
          write(u6,'(3X,A)') '-------------------------------------------------------------------------------------'// &
                             '---------------------------------'
        else if ((IFANYD /= 0) .and. (IFANYS == 0)) then
          call CollapseOutput(1,'Electric-dipole - magnetic-quadrupole transition strengths (SO states):')
          write(u6,'(3X,A)') '-----------------------------------------------------------------------'
        else if ((IFANYD == 0) .and. (IFANYS /= 0)) then
          call CollapseOutput(1,'Electric-dipole - spin-magnetic-quadrupole transition strengths (SO states):')
          write(u6,'(3X,A)') '----------------------------------------------------------------------------'
        end if

        if (OSTHR2 > Zero) then
          write(u6,30) '   for osc. strength at least ',OSTHR2
          write(u6,*)
        end if
        write(u6,31) 'From','To','Osc. strength'
        write(u6,35)
      end if

      g = FEGVAL*Three/Two ! To remove the 2/3 factor in ONEOVER9C2
      g = g*Two ! Seem to be needed to agree with the exact term,
      ! needs to be looked further into!
      do ISS=1,IEND
        do JSS=JSTART,NSS
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (abs(EDIFF) < 1.0e-8_wp) cycle
          if (EDIFF > Zero) then

            EDIFF2 = EDIFF**2

            ! Since the Spin-Magnetic-Quadrupole is made from the multiplication of two complex integrals we have
            ! M^s = (a+ib)(c+id) = ac-bd + i(ad+bc) hence the long expressions below
            ! Also, since the magnetic quadrupole terms are real and the electric dipole are imaginary
            ! we multiply the real components of MQ with the imaginary of the dipole term, and vice versa.
            ! However, the spin y component is imaginary

            ! Magnetic-Quadrupole   Spin-Magnetic-Quadrupole
            ! Electric-Dipole
            DXYDZ = -((MQXYI(JSS,ISS)+g*SXYI(JSS,ISS))*DZI(JSS,ISS))+((MQXYR(JSS,ISS)+g*SXYR(JSS,ISS))*DZR(JSS,ISS))
            DYXDZ = -((MQYXI(JSS,ISS)+g*SYXR(JSS,ISS))*DZI(JSS,ISS))+((MQYXR(JSS,ISS)+g*SYXI(JSS,ISS))*DZR(JSS,ISS))
            FXY = ONEOVER9C2*EDIFF2*(DXYDZ)
            FYX = -ONEOVER9C2*EDIFF2*(DYXDZ)

            DZXDY = -((MQZXI(JSS,ISS)+g*SZXR(JSS,ISS))*DYI(JSS,ISS))+((MQZXR(JSS,ISS)+g*SZXI(JSS,ISS))*DYR(JSS,ISS))
            DXZDY = -((MQXZI(JSS,ISS)+g*SXZR(JSS,ISS))*DYI(JSS,ISS))+((MQXZR(JSS,ISS)+g*SXZI(JSS,ISS))*DYR(JSS,ISS))
            FZX = ONEOVER9C2*EDIFF2*(DZXDY)
            FXZ = -ONEOVER9C2*EDIFF2*(DXZDY)

            DYZDX = -((MQYZI(JSS,ISS)+g*SYZR(JSS,ISS))*DXI(JSS,ISS))+((MQYZR(JSS,ISS)+g*SYZI(JSS,ISS))*DXR(JSS,ISS))
            DZYDX = -((MQZYI(JSS,ISS)+g*SZYI(JSS,ISS))*DXI(JSS,ISS))+((MQZYR(JSS,ISS)+g*SZYR(JSS,ISS))*DXR(JSS,ISS))
            FYZ = ONEOVER9C2*EDIFF2*(DYZDX)
            FZY = -ONEOVER9C2*EDIFF2*(DZYDX)

            F = FYX+FXY+FZX+FXZ+FYZ+FZY
            ! Add it to the total
            TOT2K(JSS,ISS) = TOT2K(JSS,ISS)+F

            if (abs(F) >= OSTHR2) then
              if (QIALL) write(u6,33) ISS,JSS,F
            end if
          end if
        end do
      end do

      if (QIALL) then
        write(u6,35)
        if ((IFANYD /= 0) .and. (IFANYS /= 0)) then
          call CollapseOutput(0,'Electric-dipole - magnetic-quadrupole and electric-dipole - spin-magnetic-quadrupole '// &
                              'transition strengths (SO states):')
          write(u6,*)
        else if ((IFANYD /= 0) .and. (IFANYS == 0)) then
          call CollapseOutput(0,'Electric-dipole - magnetic-quadrupole transition strengths (SO states):')
          write(u6,*)
        else if ((IFANYD == 0) .and. (IFANYS /= 0)) then
          call CollapseOutput(0,'Electric-dipole - spin-magnetic-quadrupole transition strengths (SO states):')
          write(u6,*)
        end if
      end if
      SECORD(4) = 1
    end if

    ! Magnetic-Quadrupole
    call Deallocate_Magnetic_Quadrupoles()
    ! Spin-Magnetic-Quadrupole
    call Deallocate_Spin_Magnetic_Quadrupoles()
  end if
  ! Electric-Dipole
  call Deallocate_electric_dipoles()

  ! Now write out the total

  ! Add it to the total

  I2TOT = 0
  do I=1,4
    if (SECORD(I) == 1) I2TOT = I2TOT+1
  end do
  if (I2TOT >= 1) then
    if (SECORD(1) == 0) write(u6,*) 'Magnetic-dipole - magnetic-dipole not included'
    if (SECORD(2) == 0) write(u6,*) 'Electric-quadrupole - electric-quadrupole not included'
    if (SECORD(3) == 0) write(u6,*) 'Electric-dipole - electric-octupole not included'
    if (SECORD(4) == 0) write(u6,*) 'Electric-dipole - magnetic-quadrupole not included'
    i_Print = 0
    do ISS=1,IEND
      do JSS=JSTART,JEND
        EDIFF = ENSOR(JSS)-ENSOR(ISS)
        if (abs(EDIFF) < 1.0e-8_wp) cycle
        if (EDIFF > Zero) then

          F = TOT2K(JSS,ISS)
          if (abs(F) >= OSTHR2) then
            if (i_Print == 0) then
              i_Print = 1
              call CollapseOutput(1,'Second-order contribution to the transition strengths (SO states):')
              write(u6,'(3X,A)') '------------------------------------------------------------------'

              if (OSTHR2 > Zero) then
                write(u6,30) '   for osc. strength at least ',OSTHR2
                write(u6,*)
              end if
              write(u6,31) 'From','To','Osc. strength'
              write(u6,35)
            end if
            write(u6,33) ISS,JSS,F
            call Add_Info('TMS(SO,2nd)',[F],1,6)
          end if
        end if
      end do
    end do
    if (i_Print == 1) then
      write(u6,35)
      call CollapseOutput(0,'Second-order contribution to the transition strengths (SO states):')
      write(u6,*)
    end if
  end if
  ! release the memory again
  call mma_deallocate(TOT2K)

end if

if (DOCD) then
  ! Lasse 2019
  ! New CD here with electric dipole and magnetic-dipole - velocity gauge

  ! Electric dipole (linear momentum, p)
  call Allocate_and_Load_velocities(IFANYD)

  ! Magnetic-Dipole (angular momentum, l = r x p)
  call Allocate_and_Load_Magnetic_Dipoles(IFANYM)

  if ((IFANYD /= 0) .and. (IFANYM /= 0)) then

    ! Spin-Magnetic-Dipole
    call Allocate_and_Load_Spin_Magnetic_Dipoles(IFANYS)

    ! Electric quadrupole (r:p+p:r)

    call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)

    ! Only print the part calculated

    write(u6,*)
    call CollapseOutput(1,'Circular Dichroism - velocity gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (SO states):')
    write(u6,'(3X,A)') '-----------------------------------------------------------------------------------------------------'
    if (DO_SK) then
      write(u6,30) 'For red. rot. strength at least',RSTHR
    else
      write(u6,30) 'For isotropic red. rot. strength at least',RSTHR
    end if
    write(u6,*)

    if (Do_SK .and. (IFANYQ /= 0)) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      if (Do_SK .and. (IFANYQ /= 0)) then
        write(u6,*)
        write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
        write(u6,*)
        write(u6,31) 'From','To','Red. rot. str.'
      else
        write(u6,31) 'From','To','Red. rot. str.'
        if (IFANYQ /= 0) write(u6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
      end if
      write(u6,35)

      g = FEGVAL
      do ISS=1,IEND
        do JSS=JSTART,JEND
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (abs(EDIFF) < 1.0e-8_wp) cycle
          if (EDIFF > Zero) then

            ! These are all complex quantities, and their products are complex too,
            ! but eventually every piece will be:
            !   <I|a|J> <J|b|I> + <I|b|J> <J|b|I>
            ! for a and b Hermitian operators, so it will be reduced to:
            !   2*(Re(A_ij)*Re(B_ij)+Im(A_ij)*Im(B_ij))
            ! and the imaginary parts of the products can be ignored.

            ! Note p = -i*hbar*nabla
            D_XR = DXI(JSS,ISS)
            D_YR = DYI(JSS,ISS)
            D_ZR = DZI(JSS,ISS)
            D_XI = -DXR(JSS,ISS)
            D_YI = -DYR(JSS,ISS)
            D_ZI = -DZR(JSS,ISS)
            ! Note r x p = -i*hbar * (r x nabla)
            D_MXR = MDXI(JSS,ISS)+g*SXR(JSS,ISS)
            D_MYR = MDYI(JSS,ISS)+g*SYR(JSS,ISS)
            D_MZR = MDZI(JSS,ISS)+g*SZR(JSS,ISS)
            D_MXI = -MDXR(JSS,ISS)+g*SXI(JSS,ISS)
            D_MYI = -MDYR(JSS,ISS)+g*SYI(JSS,ISS)
            D_MZI = -MDZR(JSS,ISS)+g*SZI(JSS,ISS)

            ! R = 1/3 tr(Rtensor)
            RXX = D_XR*D_MXR+D_XI*D_MXI
            RYY = D_YR*D_MYR+D_YI*D_MYI
            RZZ = D_ZR*D_MZR+D_ZI*D_MZI
            if (abs(EDIFF) > 1.0e-8_wp) then
              R = Half/EDIFF*AU2REDR*(RXX+RYY+RZZ)
            else
              R = ZERO
            end if
            write(u6,43) '1/3 Tr(RTensor): ',R

            ! Compute full rotatory strength tensor
            ! (see Hansen and Bak, 10.1021/jp001899+)

            if (IFANYQ /= 0) then
              ! Note r:p+p:r = -i*hbar * (r:nabla+nabla:r)
              Q_XXR = QXXI(JSS,ISS)
              Q_XYR = QXYI(JSS,ISS)
              Q_XZR = QXZI(JSS,ISS)
              Q_YYR = QYYI(JSS,ISS)
              Q_YZR = QYZI(JSS,ISS)
              Q_ZZR = QZZI(JSS,ISS)
              Q_XXI = -QXXR(JSS,ISS)
              Q_XYI = -QXYR(JSS,ISS)
              Q_XZI = -QXZR(JSS,ISS)
              Q_YYI = -QYYR(JSS,ISS)
              Q_YZI = -QYZR(JSS,ISS)
              Q_ZZI = -QZZR(JSS,ISS)
              RXY = D_XR*D_MYR+D_XI*D_MYI
              RXZ = D_XR*D_MZR+D_XI*D_MZI
              RYX = D_YR*D_MXR+D_YI*D_MXI
              RYZ = D_YR*D_MZR+D_YI*D_MZI
              RZX = D_ZR*D_MXR+D_ZI*D_MXI
              RZY = D_ZR*D_MYR+D_ZI*D_MYI
              RXXY = Q_XXR*D_YR+Q_XXI*D_YI
              RXXZ = Q_XXR*D_ZR+Q_XXI*D_ZI
              RXYX = Q_XYR*D_XR+Q_XYI*D_XI
              RXYZ = Q_XYR*D_ZR+Q_XYI*D_ZI
              RXZX = Q_XZR*D_XR+Q_XZI*D_XI
              RXZY = Q_XZR*D_YR+Q_XZI*D_YI
              RXYY = Q_XYR*D_YR+Q_XYI*D_YI
              RYYX = Q_YYR*D_XR+Q_YYI*D_XI
              RYYZ = Q_YYR*D_ZR+Q_YYI*D_ZI
              RYZX = Q_YZR*D_XR+Q_YZI*D_XI
              RYZY = Q_YZR*D_YR+Q_YZI*D_YI
              RXZZ = Q_XZR*D_ZR+Q_XZI*D_ZI
              RYZZ = Q_YZR*D_ZR+Q_YZI*D_ZI
              RZZX = Q_ZZR*D_XR+Q_ZZI*D_XI
              RZZY = Q_ZZR*D_YR+Q_ZZI*D_YI
              ! xx, xy, xz, yy, yz, zz
              Rtensor(1) = 0.75_wp*(RYY+RZZ+(RXYZ-RXZY))
              Rtensor(2) = -0.375_wp*(RXY+RYX+(RXXZ+RYZY-RXZX-RYYZ))
              Rtensor(3) = -0.375_wp*(RXZ+RZX+(RXYX+RZZY-RXXY-RYZZ))
              Rtensor(4) = 0.75_wp*(RXX+RZZ+(RYZX-RXYZ))
              Rtensor(5) = -0.375_wp*(RYZ+RZY+(RYYX+RXZZ-RXYY-RZZX))
              Rtensor(6) = 0.75_wp*(RXX+RYY+(RXZY-RYZX))
              if (abs(EDIFF) > 1.0e-8_wp) then
                call DSCAL_(6,AU2REDR/EDIFF,Rtensor,1)
              else
                Rtensor(:) = ZERO
              end if
              if (Do_SK) then
                ! k^T R k
                R = k_vector(1,iVec)**2*Rtensor(1)+k_vector(2,iVec)**2*Rtensor(4)+k_vector(3,iVec)**2*Rtensor(6)+ &
                    Two*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+Two*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+ &
                    Two*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
              else
                write(u6,43) 'tensor: ',Rtensor(:)
              end if
            end if

            if (abs(R) > RSTHR) write(u6,33) ISS,JSS,R

            call Add_Info('CD_V(SO)',[R],1,6)
          end if
        end do
      end do

      write(u6,35)
    end do

    call Deallocate_Spin_Magnetic_Dipoles()

    call Deallocate_Electric_Quadrupoles()

    call CollapseOutput(0,'Circular Dichroism - velocity gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (SO states):')
  end if

  call Deallocate_electric_dipoles()

  call Deallocate_magnetic_dipoles()

  ! Lasse 2019
  ! New CD here with electric dipole and magnetic-dipole - mixed gauge

  ! Electric dipole (r)
  call Allocate_and_Load_electric_dipoles(IFANYD)
  ! Magnetic-Dipole (angular momentum, l = r x p)
  call Allocate_and_Load_Magnetic_dipoles(IFANYM)

  if ((IFANYD /= 0) .and. (IFANYM /= 0)) then

    ! Spin-Magnetic-Dipole
    call Allocate_and_Load_Spin_Magnetic_Dipoles(IFANYS)

    ! Electric quadrupole (r:r)
    call Allocate_and_Load_Electric_Quadrupoles(IFANYQ)

    ! Only print the part calculated

    write(u6,*)
    call CollapseOutput(1,'Circular Dichroism - mixed gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (SO states):')
    write(u6,'(3X,A)') '--------------------------------------------------------------------------------------------------'
    write(u6,*)
    write(u6,*) ' WARNING WARNING WARNING !!!'
    write(u6,*)
    write(u6,*) ' Circular Dichroism in the mixed gauge'
    write(u6,*) ' is NOT origin independent - check your results'
    if (DO_SK) then
      write(u6,30) 'For red. rot. strength at least',RSTHR
    else
      write(u6,30) 'For isotropic red. rot. strength at least',RSTHR
    end if
    write(u6,*)

    if (Do_SK .and. (IFANYQ /= 0)) then
      nVec = nk_Vector
    else
      nVec = 1
    end if

    do iVec=1,nVec

      if (Do_SK .and. (IFANYQ /= 0)) then
        write(u6,*)
        write(u6,'(4x,a,3F10.6)') 'Direction of the k-vector: ',(k_vector(k,iVec),k=1,3)
        write(u6,*)
        write(u6,31) 'From','To','Red. rot. str.'
      else
        write(u6,31) 'From','To','Red. rot. str.'
        if (IFANYQ /= 0) write(u6,44) 'Rxx','Rxy','Rxz','Ryy','Ryz','Rzz'
      end if
      write(u6,35)

      g = FEGVAL
      do ISS=1,IEND
        do JSS=JSTART,JEND
          EDIFF = ENSOR(JSS)-ENSOR(ISS)
          if (abs(EDIFF) < 1.0e-8_wp) cycle
          if (EDIFF > Zero) then

            ! These are all complex quantities, and their products are complex too,
            ! but eventually every piece will be:
            !   <I|a|J> <J|b|I> + <I|b|J> <J|b|I>
            ! for a and b Hermitian operators, so it will be reduced to:
            !   2*(Re(A_ij)*Re(B_ij)+Im(A_ij)*Im(B_ij))
            ! and the imaginary parts of the products can be ignored.

            D_XR = DXR(JSS,ISS)
            D_YR = DYR(JSS,ISS)
            D_ZR = DZR(JSS,ISS)
            D_XI = DXI(JSS,ISS)
            D_YI = DYI(JSS,ISS)
            D_ZI = DZI(JSS,ISS)
            ! Note r x p = -i*hbar * (r x nabla),
            ! but we will need i * (r x p) = hbar * r x nabla
            D_MXI = MDXI(JSS,ISS)+g*SXR(JSS,ISS)
            D_MYI = MDYI(JSS,ISS)+g*SYR(JSS,ISS)
            D_MZI = MDZI(JSS,ISS)+g*SZR(JSS,ISS)
            D_MXR = MDXR(JSS,ISS)+g*SXI(JSS,ISS)
            D_MYR = MDYR(JSS,ISS)+g*SYI(JSS,ISS)
            D_MZR = MDZR(JSS,ISS)+g*SZI(JSS,ISS)

            ! R = 1/3 tr(Rtensor)
            RXX = D_XR*D_MXR+D_XI*D_MXI
            RYY = D_YR*D_MYR+D_YI*D_MYI
            RZZ = D_ZR*D_MZR+D_ZI*D_MZI
            R = Half*AU2REDR*(RXX+RYY+RZZ)

            ! Compute full rotatory strength tensor
            ! (see Hansen and Bak, 10.1021/jp001899+)

            if (IFANYQ /= 0) then
              Q_XXR = QXXR(JSS,ISS)
              Q_XYR = QXYR(JSS,ISS)
              Q_XZR = QXZR(JSS,ISS)
              Q_YYR = QYYR(JSS,ISS)
              Q_YZR = QYZR(JSS,ISS)
              Q_ZZR = QZZR(JSS,ISS)
              Q_XXI = QXXI(JSS,ISS)
              Q_XYI = QXYI(JSS,ISS)
              Q_XZI = QXZI(JSS,ISS)
              Q_YYI = QYYI(JSS,ISS)
              Q_YZI = QYZI(JSS,ISS)
              Q_ZZI = QZZI(JSS,ISS)
              RXY = D_XR*D_MYR+D_XI*D_MYI
              RXZ = D_XR*D_MZR+D_XI*D_MZI
              RYX = D_YR*D_MXR+D_YI*D_MXI
              RYZ = D_YR*D_MZR+D_YI*D_MZI
              RZX = D_ZR*D_MXR+D_ZI*D_MXI
              RZY = D_ZR*D_MYR+D_ZI*D_MYI
              RXXY = Q_XXR*D_YR+Q_XXI*D_YI
              RXXZ = Q_XXR*D_ZR+Q_XXI*D_ZI
              RXYX = Q_XYR*D_XR+Q_XYI*D_XI
              RXYZ = Q_XYR*D_ZR+Q_XYI*D_ZI
              RXZX = Q_XZR*D_XR+Q_XZI*D_XI
              RXZY = Q_XZR*D_YR+Q_XZI*D_YI
              RXYY = Q_XYR*D_YR+Q_XYI*D_YI
              RYYX = Q_YYR*D_XR+Q_YYI*D_XI
              RYYZ = Q_YYR*D_ZR+Q_YYI*D_ZI
              RYZX = Q_YZR*D_XR+Q_YZI*D_XI
              RYZY = Q_YZR*D_YR+Q_YZI*D_YI
              RXZZ = Q_XZR*D_ZR+Q_XZI*D_ZI
              RYZZ = Q_YZR*D_ZR+Q_YZI*D_ZI
              RZZX = Q_ZZR*D_XR+Q_ZZI*D_XI
              RZZY = Q_ZZR*D_YR+Q_ZZI*D_YI
              ! xx, xy, xz, yy, yz, zz
              Rtensor(1) = 0.75_wp*(RYY+RZZ+EDIFF*(RXYZ-RXZY))
              Rtensor(2) = -0.375_wp*(RXY+RYX+EDIFF*(RXXZ+RYZY-RXZX-RYYZ))
              Rtensor(3) = -0.375_wp*(RXZ+RZX+EDIFF*(RXYX+RZZY-RXXY-RYZZ))
              Rtensor(4) = 0.75_wp*(RXX+RZZ+EDIFF*(RYZX-RXYZ))
              Rtensor(5) = -0.375_wp*(RYZ+RZY+EDIFF*(RYYX+RXZZ-RXYY-RZZX))
              Rtensor(6) = 0.75_wp*(RXX+RYY+EDIFF*(RXZY-RYZX))
              call DSCAL_(6,AU2REDR,Rtensor,1)
              if (Do_SK) then
                ! k^T R k
                R = k_vector(1,iVec)**2*Rtensor(1)+k_vector(2,iVec)**2*Rtensor(4)+k_vector(3,iVec)**2*Rtensor(6)+ &
                    Two*k_vector(1,iVec)*k_vector(2,iVec)*Rtensor(2)+Two*k_vector(1,iVec)*k_vector(3,iVec)*Rtensor(3)+ &
                    Two*k_vector(2,iVec)*k_vector(3,iVec)*Rtensor(5)
              else
                write(u6,43) 'tensor: ',Rtensor(:)
              end if
            end if

            if (abs(R) > RSTHR) write(u6,33) ISS,JSS,R

            call Add_Info('CD_M(SO)',[R],1,6)
          end if
        end do
      end do
      write(u6,35)
    end do

    call Deallocate_Spin_Magnetic_Dipoles()

    call Deallocate_Electric_Quadrupoles()

    call CollapseOutput(0,'Circular Dichroism - mixed gauge Electric-Dipole - Magnetic-Dipole rotatory strengths (SO states):')
  end if

  call Deallocate_electric_dipoles()
  call Deallocate_Magnetic_Dipoles()
end if
! CD end

! +++ J. Norell 19/7 - 2018
! Dyson amplitudes for (1-electron) ionization transitions
if (DYSO) then
  call Add_Info('SODYSAMPS',SODYSAMPS,NSS*NSS,4)
  DYSTHR = 1.0e-5_wp
  write(u6,*)
  call CollapseOutput(1,'Dyson amplitudes (SO states):')
  write(u6,'(3X,A)') '-----------------------------------------------'
  if (DYSTHR > Zero) then
    write(u6,*) 'for Dyson intensities at least',DYSTHR
    write(u6,*)
  end if
  write(u6,*) '       From      To        BE (eV)       Dyson intensity'
  write(u6,'(3X,A)') '-----------------------------------------------------------------------------------------'
  do I=1,NSS
    do J=1,NSS
      F = SODYSAMPS(I,J)*SODYSAMPS(I,J)
      EDIFF = auToeV*(ENSOR(J)-ENSOR(I))
      if (abs(EDIFF) < 1.0e-8_wp) cycle
      if (F > 1.0e-5_wp) then
        if (EDIFF > Zero) write(u6,'(A,I8,I8,F15.3,ES22.5)') '    ',I,J,EDIFF,F
      end if
    end do ! J
  end do ! I
  write(u6,*)
  write(u6,*)
  call CollapseOutput(0,'Dyson amplitudes (SO states):')
  write(u6,*)
  ! VKochetov 2021 put SO Dyson amplitudes to hdf5
# ifdef _HDF5_
  if (rhodyn) call mh5_put_dset(wfn_sos_dys,SODYSAMPS)
# endif
end if
! +++ J. Norell

!***********************************************************************
!                                                                      *
!     Start of section for transition moments using the exact operator *
!     for the vector potential.                                        *
!                                                                      *
!***********************************************************************

if (Do_TMOM) call PRPROP_TM_Exact(PROP,USOR,USOI,ENSOR,NSS,JBNUM,EigVec)

500 continue

! CALCULATION OF THE D-TENSOR (experimental)
! IFDCAL to implement keyword that will activate computation
! of d-tensor
!if (.not. IFDCAL) goto 600
goto 600

write(u6,*)
write(u6,*) '  D-Matrix'
write(u6,*) '  ========================================='
write(u6,*) '  calculated using 2nd order perturbation'
write(u6,*) '  > any spin degeneracy, no spatial degeneracy'
write(u6,*) '  > weak spin-orbit coupling'
write(u6,*)

! SVC 2006: Compute D-tensor through second order perturbation theory.
! no orbitally-degenerate groundstates!
IAMFIX = 0
IAMFIY = 0
IAMFIZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:4) == 'AMFI') then
    if (ICOMP(IPROP) == 1) IAMFIX = IPROP
    if (ICOMP(IPROP) == 2) IAMFIY = IPROP
    if (ICOMP(IPROP) == 3) IAMFIZ = IPROP
  end if
end do
IPAMFI(1) = IAMFIX
IPAMFI(2) = IAMFIY
IPAMFI(3) = IAMFIZ
! initialisations
do IXYZ=1,3
  do JXYZ=1,3
    DTENS(IXYZ,JXYZ) = Zero
  end do
end do

! loop over all excited states, different factors will arise depending
! on the difference in spin between ground and excited states.
ISTATE = 1
MPLET1 = MLTPLT(JBNUM(ISTATE))
S1 = Half*real(MPLET1-1,kind=wp)
FACT0 = THREEJ(S1,One,S1,S1,Zero,-S1)*THREEJ(S1,One,S1,S1,Zero,-S1)/(S1*S1)
FACTP = THREEJ(S1+One,One,S1,S1+One,-One,-S1)*THREEJ(S1,One,S1+One,S1,One,-(S1+One))/((S1+One)*(Two*S1+One))
FACTM = THREEJ(S1-One,One,S1,S1-One,One,-S1)*THREEJ(S1,One,S1-One,S1,-One,-(S1-One))/(S1*(Two*S1-One))
!write(u6,*)
!write(u6,*) 'S1 ', S1
!write(u6,*) 'FACT0 ', FACT0
!write(u6,*) 'FACTP ', FACTP
!write(u6,*) 'FACTM ', FACTM
!write(u6,*)
do IXYZ=1,3
  do JXYZ=1,3
    DTIJ = Zero
    do JSTATE=2,NSTATE
      MPLET2 = MLTPLT(JBNUM(JSTATE))
      S2 = Half*real(MPLET2-1,kind=wp)
      DELTA = ENERGY(JSTATE)-ENERGY(ISTATE)
      if (DELTA < 1.0e-5_wp) goto 600
      CONTRIB = Zero
      CONTRIB = PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*PROP(JSTATE,ISTATE,IPAMFI(JXYZ))/DELTA
      !write(u6,*) 'ISTATE, JSTATE, IXYZ, JXYZ, DELTA ',ISTATE,JSTATE,IXYZ,JXYZ,DELTA
      !write(u6,*) 'CONTRIB ', CONTRIB
      !write(u6,*) 'PROP(ISTATE,JSTATE,IXYZ) ',PROP(ISTATE,JSTATE,IPAMFI(IXYZ))
      !write(u6,*) 'PROP(JSTATE,ISTATE,JXYZ) ',PROP(JSTATE,ISTATE,IPAMFI(JXYZ))
      if (S2 == S1) then
        DTIJ = DTIJ+FACT0*CONTRIB
      else if (S2 == S1+One) then
        DTIJ = DTIJ+FACTP*CONTRIB
      else if (S2 == S1-One) then
        DTIJ = DTIJ+FACTM*CONTRIB
      end if
    end do
    DTENS(IXYZ,JXYZ) = DTIJ
  end do
end do

! diagonalisation of the D-tensor matrix
do I=1,3
  EVR(I) = Zero
  EVI(I) = Zero
end do
do IXYZ=1,3
  do JXYZ=1,3
    TMPMAT(IXYZ,JXYZ) = DTENS(IXYZ,JXYZ)
  end do
end do
call unitmat(TMPVEC,3)
call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

! D-factor printout
! D = D_zz - 1/2 * (D_xx + D_yy)
! E = 1/2 * (D_xx - D_yy)
write(u6,*)
write(u6,*) 'The D matrix and eigenvalues:'
write(u6,*)
write(u6,'(2x,2x,2x,3(5x,a2,5x),4x,4x,2x,10x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('D_',IXYZ,IXYZ=1,3)
write(u6,*)
do IXYZ=1,3
  write(u6,'(2x,a2,2x,3(1x,f10.8,1x),4x,a2,i1,a1,2x,f10.8,2x,a2,2x,3(1x,f8.4,1x),3x,a2,i1,a1,2x,f8.3,2x,a5)') &
    xyzchr(IXYZ),(DTENS(IXYZ,JXYZ),JXYZ=1,3),'D_',IXYZ,':',EVR(IXYZ),xyzchr(IXYZ),(TMPVEC(IXYZ,JXYZ),JXYZ=1,3),'D_',IXYZ,':', &
    EVR(IXYZ)*auTocm,'cm^-1'
end do

600 continue

! CALCULATION OF THE G-TENSOR
! IFGCAL is set by keyword EPRG
if (.not. IFGCAL) goto 800
! PAM 2005 Experimental: Compute g-tensor through 2-nd order
! perturbation approach, mixed ang mom / spin orbit terms
! Declarations for gtensor(3,3) and some other odds and ends
! for nice output have been added to declaration head above.

if (.not. IFSO) then
  write(u6,*) 'keyword SPIN needed together with EPRG'
  write(u6,*)
  goto 800
end if

write(u6,*)
write(u6,*) '  g-Matrix Approach I'
write(u6,*) '  ========================================='
write(u6,*) '  calculated using 2nd order perturbation'
write(u6,*) '  > any spin degeneracy, no spatial degeneracy'
write(u6,*) '  > weak spin-orbit coupling'
write(u6,*)

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
IPAMFI(1) = IAMFIX
IPAMFI(2) = IAMFIY
IPAMFI(3) = IAMFIZ
IPAM(1) = IAMX
IPAM(2) = IAMY
IPAM(3) = IAMZ

! start loop over the states ISTATE:
ISTATE = 1
do while (((ENERGY(min(ISTATE,NSTATE))-ENERGY(1)) <= EPRTHR) .and. (ISTATE <= NSTATE))

  do IXYZ=1,3
    do JXYZ=1,3
      GTENS(IXYZ,JXYZ) = Zero
    end do
  end do

  MPLET = MLTPLT(JBNUM(ISTATE))
  S = Half*real(MPLET-1,kind=wp)

  write(u6,*)
  write(u6,'(3x,A6,I4,3x,A4,F4.1,3x,A4,F18.8)') 'STATE ',ISTATE,'S = ',S,'E = ',ENERGY(ISTATE)
  write(u6,'(3x,A46)') '----------------------------------------------'

  if (MPLET /= 1) then
    FACTOR = One/sqrt(S*(S+One)*(Two*S+One))
  else
    goto 690
  end if

  ! print separate contributions if verbose
  if (IPGLOB >= 3) then
    write(u6,*)
    write(u6,*) 'contributions from the SOS expansion to delta(g_pq) in *ppt* (p,q=x,y,z)'
    write(u6,*)
    write(u6,'(2x,a8,2x,9(4x,a2,3x))') ' states ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
    write(u6,*)
    do JSTATE=1,NSTATE
      if (JSTATE /= ISTATE) then
        DELTA = ENERGY(JSTATE)-ENERGY(ISTATE)
        if (abs(DELTA) < 1.0e-4_wp) then
          write(u6,'(1x,i3,2x,i3,3x,A20,1x,A20,F18.8)') ISTATE,JSTATE,'possible degeneracy,','energy difference = ',DELTA
          goto 610
        end if
        do IXYZ=1,3
          do JXYZ=1,3
            CONTRIB = PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*PROP(ISTATE,JSTATE,IPAM(JXYZ))
            CONTRIB = CONTRIB/DELTA
            SOSTERM(3*(IXYZ-1)+JXYZ) = -Two*FACTOR*CONTRIB
          end do
        end do
        write(u6,'(1x,i3,2x,i3,3x,9(f8.3,1x))') ISTATE,JSTATE,(SOSTERM(I)*1.0e3_wp,I=1,9)
      end if
610   continue
    end do
  end if

  ! calculate sum-over-states for each g_pq (p,q = x,y,z)
  do IXYZ=1,3
    do JXYZ=1,3
      GTIJ = Zero
      DELTA = Zero
      CONTRIB = Zero
      do JSTATE=1,NSTATE
        if (JSTATE /= ISTATE) then
          DELTA = ENERGY(JSTATE)-ENERGY(ISTATE)
          ! SVC 2008: no good criterium for spatial degeneracy, use rasscf energies ?
          if (abs(DELTA) < 1.0e-4_wp) then
            write(u6,*)
            write(u6,*) 'SPATIALLY DEGENERATE STATE: sum-over-states not applicable'
            !write(u6,*) '> lower the degeneracy treshold if this is not a spatially degenerate state'
            write(u6,*)
            goto 690
          end if
          CONTRIB = PROP(ISTATE,JSTATE,IPAMFI(IXYZ))*PROP(ISTATE,JSTATE,IPAM(JXYZ))
          CONTRIB = CONTRIB/DELTA
          GTIJ = GTIJ+CONTRIB
        end if
      end do
      GTENS(IXYZ,JXYZ) = -Two*FACTOR*GTIJ
    end do
  end do

  ! put g_e on the diagonal
  do IXYZ=1,3
    do JXYZ=IXYZ,3
      if (IXYZ == JXYZ) GTENS(IXYZ,JXYZ) = GTENS(IXYZ,JXYZ)+FEGVAL
    end do
  end do

  ! determine symmetric G = gg+ tensor, this is what can be measured
  ! experimentally, and store as GSTENS
  do IXYZ=1,3
    do JXYZ=1,3
      GSTENS(IXYZ,JXYZ) = Zero
      do KXYZ=1,3
        GSTENS(IXYZ,JXYZ) = GSTENS(IXYZ,JXYZ)+GTENS(IXYZ,KXYZ)*GTENS(JXYZ,KXYZ)
      end do
    end do
  end do

  ! determine the eigenvalues of the g matrix
  do I=1,3
    EVR(I) = Zero
    EVI(I) = Zero
  end do
  ! XEIGEN alters the input matrix! copy GTENS to TMPMAT
  do IXYZ=1,3
    do JXYZ=1,3
      TMPMAT(IXYZ,JXYZ) = GTENS(IXYZ,JXYZ)
    end do
  end do
  call unitmat(TMPVEC,3)

  IERR = 0
  call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)
  if (IERR /= 0) then
    write(u6,*) 'Error: xEigen returned IERR = ',IERR
    return
  end if

  write(u6,*)
  write(u6,*) 'The g matrix and eigenvalues:'
  write(u6,*)
  write(u6,'(6x,3(5x,a2,5x))') (xyzchr(IXYZ),IXYZ=1,3)
  write(u6,*)
  do IXYZ=1,3
    write(u6,'(2x,a2,2x,3(1x,f10.8,1x),4x,a2,i1,a1,2x,f8.4,3x,a8,i1,a2,2x,f10.3,2x,a3)') xyzchr(IXYZ),(GTENS(IXYZ,JXYZ),JXYZ=1,3), &
                                                                                         'g_',IXYZ,':',EVR(IXYZ),'delta(g_', &
                                                                                         IXYZ,'):',(EVR(IXYZ)-FEGVAL)*1.0e3_wp,'ppt'
  end do

  !write(u6,'(6x,3(5x,i1,4x))') (IXYZ,IXYZ=1,3)

  ! determine the eigenvalues of the G = gg* matrix
  do I=1,3
    EVR(I) = Zero
    EVI(I) = Zero
  end do
  do IXYZ=1,3
    do JXYZ=1,3
      TMPMAT(IXYZ,JXYZ) = GSTENS(IXYZ,JXYZ)
    end do
  end do
  call unitmat(TMPVEC,3)

  IERR = 0
  call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)
  if (IERR /= 0) then
    write(u6,*) 'Error: xEigen returned IERR = ',IERR
    return
  end if

  ! reconstruct g_s from the square root of the eigenvalues
  ! and the eigenvectors of G = gg+ by back transformation
  do IXYZ=1,3
    do JXYZ=1,3
      GSTENS(IXYZ,JXYZ) = Zero
      do KXYZ=1,3
        GSTENS(IXYZ,JXYZ) = GSTENS(IXYZ,JXYZ)+TMPVEC(IXYZ,KXYZ)*sqrt(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
      end do
    end do
  end do

  write(u6,*)
  write(u6,*) 'The symmetric g matrix is the actual experimentally determined g matrix.'
  write(u6,*) 'The sign of the eigenvalues is undetermined (assumed positive).'
  write(u6,*)
  write(u6,'(2x,2x,2x,3(5x,a2,5x),4x,4x,2x,10x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('g_',IXYZ,IXYZ=1,3)
  write(u6,*)
  do IXYZ=1,3
    write(u6,'(2x,a2,2x,3(1x,f10.8,1x),4x,a2,i1,a1,2x,f10.6,2x,a2,2x,3(1x,f8.4,1x),3x,a8,i1,a2,2x,f8.3,2x,a3)') &
      xyzchr(IXYZ),(GSTENS(IXYZ,JXYZ),JXYZ=1,3),'g_',IXYZ,':',sqrt(EVR(IXYZ)),xyzchr(IXYZ),(TMPVEC(IXYZ,JXYZ),JXYZ=1,3), &
      'delta(g_',IXYZ,'):',(sqrt(EVR(IXYZ))-FEGVAL)*1.0e3_wp,'ppt'
  end do
  do I=1,3
    EVR(I) = sqrt(EVR(I))-FEGVAL
  end do
  call Add_Info('EPRGVAL',EVR,3,6)

690 continue

  ISTATE = ISTATE+1

  ! end long loop over states ISTATE
end do

! SVC alternative approach for the g-tensor:
! using first order degenerate perturbation theory

if (IFVANVLECK) then
  write(u6,*)
  write(u6,*) '  VAN VLECK Tensor and g-Matrix Approach II'
  write(u6,*) '  ========================================='
  write(u6,*) '  1st order degenerate perturbation theory'
  write(u6,*) '  within isolated kramers doublets.'
  write(u6,*) '  > spatial degeneracy'
  write(u6,*) '  > strong spin-orbit coupling'
  write(u6,*)
else
  write(u6,*)
  write(u6,*) '  g-Matrix Approach II'
  write(u6,*) '  ========================================='
  write(u6,*) '  1st order degenerate perturbation theory'
  write(u6,*) '  within isolated kramers doublets.'
  write(u6,*) '  > spatial degeneracy'
  write(u6,*) '  > strong spin-orbit coupling'
  write(u6,*)
end if

IAMX = 0
IAMY = 0
IAMZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    !write(u6,*) '3****ANGMOM rassi/prprop'
    if (ICOMP(IPROP) == 1) IAMX = IPROP
    if (ICOMP(IPROP) == 2) IAMY = IPROP
    if (ICOMP(IPROP) == 3) IAMZ = IPROP
  end if
end do

call mma_allocate(LXI,NSS,NSS,Label='LXI')
LXI(:,:) = Zero
call mma_allocate(LYI,NSS,NSS,Label='LYI')
LYI(:,:) = Zero
call mma_allocate(LZI,NSS,NSS,Label='LZI')
LZI(:,:) = Zero

if (IAMX > 0) call SMMAT(PROP,LXI,NSS,IAMX,0)
if (IAMY > 0) call SMMAT(PROP,LYI,NSS,IAMY,0)
if (IAMZ > 0) call SMMAT(PROP,LZI,NSS,IAMZ,0)

call mma_allocate(ZXYZR,NSS,NSS,3,Label='ZXYZR')
call mma_allocate(ZXYZI,NSS,NSS,3,Label='ZXYZI')
ZXYZR(:,:,:) = Zero
ZXYZI(:,:,:) = Zero

call SMMAT(PROP,ZXYZR(:,:,1),NSS,0,1)
call SMMAT(PROP,ZXYZI(:,:,2),NSS,0,2)
call SMMAT(PROP,ZXYZR(:,:,3),NSS,0,3)

call DSCAL_(NSS**2,FEGVAL,ZXYZR(:,:,1),1)
call DSCAL_(NSS**2,FEGVAL,ZXYZI(:,:,2),1)
call DSCAL_(NSS**2,FEGVAL,ZXYZR(:,:,3),1)

call DAXPY_(NSS**2,One,LXI,1,ZXYZI(:,:,1),1)
call DAXPY_(NSS**2,One,LYI,1,ZXYZI(:,:,2),1)
call DAXPY_(NSS**2,One,LZI,1,ZXYZI(:,:,3),1)

call mma_deallocate(LXI)
call mma_deallocate(LYI)
call mma_deallocate(LZI)

! SVC 20090926 Experimental
! Add analysis of different contributions

! Establish which spin components of SFS belong to the ground state
call mma_allocate(ISGS,NSS,Label='ISGS')
do I=1,NSS
  ISGS(I) = .false.
end do

GSENERGY = ENERGY(1)
do ISTATE=2,NSTATE
  if (ENERGY(ISTATE) < GSENERGY) GSENERGY = ENERGY(ISTATE)
end do

IMLTPL = 1
do ISTATE=1,NSTATE
  if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then
    do I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
      ISGS(I) = .true.
    end do
  else
    do I=IMLTPL,IMLTPL-1+MLTPLT(JBNUM(ISTATE))
      ISGS(I) = .false.
    end do
  end if
  IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
end do

! Analyze the different contributions to the GS Kramers doublet
! Zeeman matrix elements.  There are 4 different ME's: <1|Ze|1>,
! <1|Ze|2>, <2|Ze|1>, and <2|Ze|2>, stored in ZEKL.  Contributions
! of SFS i,j to SOS k,l (k,l=1,2): <k|Ze|l> = Sum(i,j) U(i,k)*
! <i|Ze|j> U(j,l).  This sum is decomposed into parts belonging to
! each SFS state i as follows:
! -> GS's contain only MEs with themselves and other GS's
! -> ES's contain MEs with themselves, the GS's (2x) and other ES's
!    The ME's with the GS's are counted twice as they do not belong to
!    any GS's (they contain only ME's within their own GS group)
!    The contributions with other ES's are split between the ES's,
!    counting them double (<i|Ze|j> and <j|Ze|i>) and divide by two later.

call mma_allocate(ZEKL,2,2,3,NSTATE,Label='ZEKL')

IMLTPL = 1
ZEKL(:,:,:,:) = cZero

do ISTATE=1,NSTATE

  ISTART = IMLTPL
  IFINAL = IMLTPL-1+MLTPLT(JBNUM(ISTATE))

  if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

    ! Contribution of the GS spin components
    do IXYZ=1,3
      do ISS=ISTART,IFINAL
        do JSS=1,NSS
          if (ISGS(JSS)) then
            call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
            call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
            !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
          end if
        end do
      end do
    end do

  else

    ! Contributions of the ES spin components
    do IXYZ=1,3
      do ISS=ISTART,IFINAL
        do JSS=1,NSS
          call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
          call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
          if (ISGS(JSS)) then
            call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,ISS,JSS)
            call ZECON(NSTATE,NSS,USOR,USOI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ),ZEKL,IXYZ,ISTATE,JSS,ISS)
          end if
          !write(u6,710) 'ZEKL',ISTATE,IXYZ,ISS,JSS,ZEKL(:,:,IXYZ,ISTATE)
        end do
      end do
    end do

  end if

  !do IXYZ=1,3
  !  write(u6,720) 'ZEKL',IXYZ,ISTATE,ZEKL(:,:,IXYZ,ISTATE)
  !end do

  IMLTPL = IMLTPL+MLTPLT(JBNUM(ISTATE))
end do

call mma_deallocate(ISGS)

! We now have decomposed the <k|Ze|l> into terms belonging to either
! a GS or an ES for each k,l=1,2 and p=x,y,z stored in ZEKL(k,l,p,SFS)
! Now, these new decomposed terms of the ZEKL ME's are combined to
! form the G tensor.  Consider e.g. that <k|Ze|l> is decomposed into
! <k|GS|l> + <k|ES1|l> + <k|ES2|l>, then the contributions to G are given as:
! -> G_pq/2 = <k|Ze_p|l> <l|Ze_q|k>
!         = (<k|GS_p|l> + <k|ES1_p|l> + <k|ES2_p|l>)
!         * (<l|GS_q|k> + <l|ES1_q|k> + <l|ES2_q|k>)

! from GS: (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_p|l>/2 * <l|GS_q|k>/2)/2
! from ES1: 2*((<k|ES1_p|l>/2 * <l|GS_q|k>/2)/2 + (<k|GS_q|l>/2 * <l|ES1_p|k>/2)/2)
!           + (<k|ES1_p|l>/2 * <l|ES1_q|k>/2)/2 + (<k|ES1_q|l>/2 * <l|ES2_p|k>/2)/2
!           + (<k|ES1_p|l>/2 * <l|ES2_q|k>/2)/2 + (<k|ES2_q|l>/2 * <l|ES1_p|k>/2)/2
! In the end, the outer division by 2 cancels on both sides, and the
! inner divisions by two combine to a division by 4.

call mma_allocate(GCONT,9,NSTATE,Label='GCONT')

GCONT(:,:) = cZero
GTOTAL(:) = Zero

do ISTATE=1,NSTATE
  do IXYZ=1,3
    do JXYZ=1,3
      IJXYZ = 3*(IXYZ-1)+JXYZ

      if (abs(ENERGY(ISTATE)-GSENERGY) < 1.0e-6_wp) then

        ! Contributions for the GS's
        do JSTATE=1,NSTATE
          if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) then
            do I=1,2
              do J=1,2
                GCONT(IJXYZ,ISTATE) = GCONT(IJXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE))*Quart+ &
                                      (ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
              end do
            end do
          end if
        end do

      else

        ! Contributions for the ES's
        do JSTATE=1,NSTATE
          do I=1,2
            do J=1,2
              GCONT(IJXYZ,ISTATE) = GCONT(IJXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE))*Quart+ &
                                    (ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
              if (abs(ENERGY(JSTATE)-GSENERGY) < 1.0e-6_wp) &
                GCONT(IJXYZ,ISTATE) = GCONT(IJXYZ,ISTATE)+(ZEKL(I,J,IXYZ,ISTATE)*ZEKL(J,I,JXYZ,JSTATE))*Quart+ &
                                      (ZEKL(I,J,IXYZ,JSTATE)*ZEKL(J,I,JXYZ,ISTATE))*Quart
            end do
          end do
        end do

      end if

    end do
  end do

  do IJXYZ=1,9
    GTOTAL(IJXYZ) = GTOTAL(IJXYZ)+real(GCONT(IJXYZ,ISTATE))
  end do
end do

call mma_deallocate(ZEKL)

call mma_allocate(DIPSOm,3,NSS,NSS,Label='DIPSOm')
call mma_allocate(DIPSOn,3,NSS,NSS,Label='DIPSOn')
call mma_allocate(ESO,NSS,Label='ESO')
call mma_allocate(Z,NSS,NSS,Label='Z')

DIPSOm(:,:,:) = cZero
DIPSOn(:,:,:) = cZero

! Continue original calculation of G tensor (=gg^*)
call get_dArray('ESO_SINGLE',ESO,NSS)
call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,1),ZXYZI(:,:,1))
call MULMAT(NSS,ZXYZR(:,:,1),ZXYZI(:,:,1),eex,Z)
do ISS=1,NSS
  do JSS=1,NSS
    DIPSOm(1,ISS,JSS) = Half*Z(ISS,JSS)
    DIPSOn(1,ISS,JSS) = -Z(ISS,JSS)
  end do
end do
call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,2),ZXYZI(:,:,2))
call MULMAT(NSS,ZXYZR(:,:,2),ZXYZI(:,:,2),eey,Z)
do ISS=1,NSS
  do JSS=1,NSS
    DIPSOm(2,ISS,JSS) = Half*Z(ISS,JSS)
    DIPSOn(2,ISS,JSS) = -Z(ISS,JSS)
  end do
end do
call ZTRNSF(NSS,USOR,USOI,ZXYZR(:,:,3),ZXYZI(:,:,3))
call MULMAT(NSS,ZXYZR(:,:,3),ZXYZI,eez,Z)
do ISS=1,NSS
  do JSS=1,NSS
    DIPSOm(3,ISS,JSS) = Half*Z(ISS,JSS)
    DIPSOn(3,ISS,JSS) = -Z(ISS,JSS)
  end do
end do
write(u6,*) ''

call mma_deallocate(Z)

if (IFVANVLECK) then

  call mma_allocate(chiT_tens,NTS,3,3,Label='chiT_tens')
  call mma_allocate(chicuriT_tens,NTS,3,3,Label='chicuriT_tens')
  call mma_allocate(chiparamT_tens,NTS,3,3,Label='chiparamT_tens')

  iT = 0
  do iT=1,NTS
    do ic=1,3
      do jc=1,3
        chiT_tens(iT,ic,jc) = Zero
        chicuriT_tens(iT,ic,jc) = Zero
        chiparamT_tens(iT,ic,jc) = Zero
      end do
    end do
  end do
  iT = 0
  do iT=1,NTS
    if (iT == 1) then
      TMPm(iT) = TMINS+1.0e-4_wp
    else
      DLTT = (TMAXS-TMINS)/real(NTS-1,kind=wp)
      TMPm(iT) = TMINS+DLTT*real(iT-1,kind=wp)
    end if
    Zstat = Zero
    do Iss=1,Nss
      p_Boltz = exp(-ESO(Iss)/Boltz_k/TMPm(iT))
      Zstat = Zstat+p_Boltz
      do IC=1,3
        do JC=1,3
          c_2(IC,JC) = Zero
          curit(IC,JC) = Zero
          paramt(IC,JC) = Zero
        end do
      end do
      do Jss=1,Nss
        dlt_E = Eso(Iss)-Eso(Jss)
        do IC=1,3
          do JC=1,3
            c_1(IC,JC) = Zero
          end do
        end do
        do ic=1,3
          do jc=1,3
            c_1(ic,jc) = real(DIPSOn(ic,Iss,Jss)*conjg(DIPSOn(jc,Iss,Jss)))
            if (abs(dlt_E) < 10.97_wp) then ! what is this number?
              c_2(ic,jc) = c_2(ic,jc)+c_1(ic,jc)
              curit(ic,jc) = curit(ic,jc)+c_1(ic,jc)
              !paramt(ic,jc) = paramt(ic,jc)+Zero*c_1(ic,jc)
            else
              c_2(ic,jc) = c_2(ic,jc)-Two*Boltz_k*TMPm(iT)*c_1(ic,jc)/dlt_E
              !curit(ic,jc) = curit(ic,jc)-Zero*(Two*Boltz_k*TMPm(iT)*c_1(ic,jc)/dlt_E)
              paramt(ic,jc) = paramt(ic,jc)-Two*Boltz_k*TMPm(iT)*c_1(ic,jc)/dlt_E
            end if
          end do
        end do
      end do !Jss
      do ic=1,3
        do jc=1,3
          chiT_tens(iT,ic,jc) = chiT_tens(iT,ic,jc)+p_Boltz*c_2(ic,jc)
          chicuriT_tens(iT,ic,jc) = chicuriT_tens(iT,ic,jc)+p_Boltz*curit(ic,jc)
          chiparamT_tens(iT,ic,jc) = chiparamT_tens(iT,ic,jc)+p_Boltz*paramt(ic,jc)
        end do
      end do
    end do !Iss
    !Zstat1m(iT) = Zstat
    do ic=1,3
      do jc=1,3
        chiT_tens(iT,ic,jc) = coeff_chi*(chiT_tens(iT,ic,jc)/Zstat)
        chicuriT_tens(iT,ic,jc) = coeff_chi*(chicuriT_tens(iT,ic,jc)/Zstat)
        chiparamT_tens(iT,ic,jc) = coeff_chi*(chiparamT_tens(iT,ic,jc)/Zstat)
      end do
    end do
  end do ! iT

  write(u6,'(/)')
  write(u6,'(A)') repeat('-',120)
  write(u6,'(30X,A)') 'VAN VLECK SUSCEPTIBILITY TENSOR  (cm3*K/mol)'
  write(u6,'(A)') repeat('-',120)
  write(u6,*)
  !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
  write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
  write(u6,*)
  do iT=1,NTS
    write(u6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),((chiT_tens(iT,ic,jc),jc=1,3),ic=1,3)
  end do

  write(u6,'(/)')
  write(u6,'(A)') repeat('-',120)
  write(u6,'(30X,A)') 'Curie contrib. to VAN VLECK TENSOR  (cm3*K/mol)'
  write(u6,'(A)') repeat('-',120)
  write(u6,*)
  !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
  write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
  write(u6,*)
  do iT=1,NTS
    write(u6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),((chicuriT_tens(iT,ic,jc),jc=1,3),ic=1,3)
  end do
  write(u6,'(/)')
  write(u6,'(A)') repeat('-',120)
  write(u6,'(30X,A)') 'Parama. contrib. to VAN VLECK TENSOR  (cm3*K/mol)'
  write(u6,'(A)') repeat('-',120)
  write(u6,*)
  !write(u6,'(8X,A,9(7X,A))') 'T','(1,1)','(1,2)','(1,3)','(2,1)','(2,2)','(2,3)','(3,1)','(3,2)','(3,3)'
  write(u6,'(6X,A,8X,9(A,9X))') 'T(K)','xx','xy','xz','yx','yy','yz','zx','zy','zz'
  write(u6,*)
  do iT=1,NTS
    write(u6,'(4X,F6.1,3X,11(F9.4,2X),F8.4)') TMPm(iT),((chiparamT_tens(iT,ic,jc),jc=1,3),ic=1,3)
  end do
  write(u6,*)
  write(u6,*)
  write(u6,*) '  g-Matrix'
  write(u6,*) '  =========='

  !do I=1,3
  !  do J=1,3
  !    do iT=1,NT
  !      chiT_tens(iT,I,J) = Zero
  !    end do
  !  end do
  !end do
  call mma_deallocate(chiT_tens)
  call mma_deallocate(chicuriT_tens)
  call mma_deallocate(chiparamT_tens)
end if ! IFVANVLECK

call mma_allocate(SPNSFS,3,NSS,NSS,Label='SPNSFS')

ISS = 1
do while ((ISS <= NSS) .and. (ENSOR(min(ISS,NSS))-ENSOR(1) <= EPRTHR))

  do IXYZ=1,3
    do JXYZ=1,3
      GTENS(IXYZ,JXYZ) = Zero
    end do
  end do

  KDGN = 1
  do JSS=ISS+1,NSS
    EDIFF = ENSOR(JSS)-ENSOR(ISS)
    if (IFGTCALSA .and. IFGTSHSA) then
      KDGN = MULTIP
      !write(u6,*) 'KDGN=',KDGN
    else if (abs(EDIFF) < 1.0e-6_wp) then
      KDGN = KDGN+1
    end if
  end do

  write(u6,*)
  do I=1,KDGN
    write(u6,'(3x,A9,I4,3x,A4,F18.8)') 'SO-STATE ',ISS-1+I,'E = ',ENSOR(ISS-1+I)
  end do
  write(u6,'(3x,A46)') '----------------------------------------------'

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  if (.not. IFGTCALSA) goto 450
  if (ISS == 1) IFUNCT = 0
  call SINANI(KDGN,IFUNCT,NSS,DIPSOn,SPNSFS,DIPSOm_SA)
  IFUNCT = IFUNCT+KDGN
450 continue
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

  if (KDGN /= 2) then
    write(u6,*) 'no twofold degeneracy'
    goto 780
  end if

  if ((ISS == 1) .and. (KDGN == 2) .and. (IPGLOB >= 3)) then
    write(u6,*) 'Experimental: SFS contributions to G=gg+'
    write(u6,*)
    write(u6,'(a6,9(5x,a2,5x))') 'state ','xx','xy','xz','yx','yy','yz','zx','zy','zz'
    write(u6,*)
    do ISTATE=1,NSTATE
      write(u6,'(2x,I2,2x,9(F12.6))') ISTATE,(real(GCONT(IJXYZ,ISTATE)),IJXYZ=1,9)
    end do

    write(u6,*)
    write(u6,'(A6,9(F12.6))') 'total ',(GTOTAL(IJXYZ),IJXYZ=1,9)
  end if

  JSS = ISS+1

  do IXYZ=1,3
    do JXYZ=1,3
      GTIJ = Zero
      CONTRIB = Zero
      do ISO=ISS,JSS
        do JSO=ISS,JSS
          CONTRIB = ZXYZR(ISO,JSO,IXYZ)*ZXYZR(JSO,ISO,JXYZ)-ZXYZI(ISO,JSO,IXYZ)*ZXYZI(JSO,ISO,JXYZ)
          GTIJ = GTIJ+CONTRIB
        end do
      end do
      GTENS(IXYZ,JXYZ) = Two*GTIJ
    end do
  end do

  if (IPGLOB > 3) then
    write(u6,*) 'G tensor = gg+'
    write(u6,*)
    write(u6,'(6x,3(6x,a2,4x))') (xyzchr(IXYZ),IXYZ=1,3)
    do IXYZ=1,3
      write(u6,'(2x,a2,2x,3(1x,f18.8,1x))') xyzchr(IXYZ),(GTENS(IXYZ,JXYZ),JXYZ=1,3)
    end do
  end if

  do I=1,3
    EVR(I) = Zero
    EVI(I) = Zero
  end do
  do IXYZ=1,3
    do JXYZ=1,3
      TMPMAT(IXYZ,JXYZ) = GTENS(IXYZ,JXYZ)
    end do
  end do
  call unitmat(TMPVEC,3)

  call XEIGEN(1,3,3,TMPMAT,EVR,EVI,TMPVEC,IERR)

  ! construct g_s matrix from G by back-transormation of the
  ! square root of the G eigenvalues
  do IXYZ=1,3
    do JXYZ=1,3
      GTENS(IXYZ,JXYZ) = Zero
      do KXYZ=1,3
        GTENS(IXYZ,JXYZ) = GTENS(IXYZ,JXYZ)+TMPVEC(IXYZ,KXYZ)*sqrt(EVR(KXYZ))*TMPVEC(JXYZ,KXYZ)
      end do
    end do
  end do

  write(u6,'(6x,3(5x,a2,5x),4x,4x,2x,8x,2x,2x,2x,3(4x,a2,i1,3x))') (xyzchr(IXYZ),IXYZ=1,3),('g_',IXYZ,IXYZ=1,3)
  write(u6,*)
  do IXYZ=1,3
    write(u6,'(2x,a2,2x,3(1x,f10.6,1x),4x,a2,i1,a1,2x,f8.4,2x,a2,2x,3(1x,f8.4,1x))') xyzchr(IXYZ), &
                                                                                     (GTENS(IXYZ,JXYZ),JXYZ=1,3),'g_',IXYZ, &
                                                                                     ':',sqrt(EVR(IXYZ)),xyzchr(IXYZ), &
                                                                                     (TMPVEC(IXYZ,JXYZ),JXYZ=1,3)
  end do

780 continue

  ISS = ISS+KDGN

end do

call mma_deallocate(ZXYZR)
call mma_deallocate(ZXYZI)
call mma_deallocate(DIPSOn)
call mma_deallocate(GCONT)
call mma_deallocate(SPNSFS)

800 continue

!*****************************************************
!* Experimental hyperfine tensor stuff starts here
!*****************************************************

! Skip if not a hyperfine calculation
if (IFACAL) call HFCTS(PROP,USOR,USOI,ENSOR,NSS,ENERGY,JBNUM,DIPSOM,ESO,XYZCHR,BOLTZ_K)

!*****************************************************
!* Experimental hyperfine tensor stuff ends here
!*****************************************************

call mma_deallocate(DIPSOm,safe='*')
call mma_deallocate(ESO,safe='*')

! SVC20080312 calculation of magnetization

if (.not. IFXCAL) goto 900

if (.not. IFSO) then
  write(u6,*) 'keyword SPIN needed with MAGN'
  write(u6,*)
  goto 900
end if

write(u6,*)
write(u6,*) '  ========================================='
write(u6,*) '  Magnetization and Magnetic Susceptibility'
write(u6,*) '  ========================================='
write(u6,*)

! initialization same as G-tensor, construct L+gS matrix elements
IAMX = 0
IAMY = 0
IAMZ = 0
do IPROP=1,NPROP
  if (PNAME(IPROP)(1:6) == 'ANGMOM') then
    !write(u6,*) '4****ANGMOM rassi/prprop'
    if (ICOMP(IPROP) == 1) IAMX = IPROP
    if (ICOMP(IPROP) == 2) IAMY = IPROP
    if (ICOMP(IPROP) == 3) IAMZ = IPROP
  end if
end do

call mma_allocate(LXI,NSS,NSS,Label='LXI')
LXI(:,:) = Zero
call mma_allocate(LYI,NSS,NSS,Label='LYI')
LYI(:,:) = Zero
call mma_allocate(LZI,NSS,NSS,Label='LZI')
LZI(:,:) = Zero

if (IAMX > 0) call SMMAT(PROP,LXI,NSS,IAMX,0)
if (IAMY > 0) call SMMAT(PROP,LYI,NSS,IAMY,0)
if (IAMZ > 0) call SMMAT(PROP,LZI,NSS,IAMZ,0)

call mma_allocate(MXYZR,NSS,NSS,3,Label='MXYZR')
call mma_allocate(MXYZI,NSS,NSS,3,Label='MXYZI')
MXYZR(:,:,:) = Zero
MXYZI(:,:,:) = Zero

call SMMAT(PROP,MXYZR(:,:,1),NSS,0,1)
call SMMAT(PROP,MXYZI(:,:,2),NSS,0,2)
call SMMAT(PROP,MXYZR(:,:,3),NSS,0,3)

call DSCAL_(NSS**2,FEGVAL,MXYZR(:,:,1),1)
call DSCAL_(NSS**2,FEGVAL,MXYZI(:,:,2),1)
call DSCAL_(NSS**2,FEGVAL,MXYZR(:,:,3),1)

call DAXPY_(NSS**2,One,LXI,1,MXYZI(:,:,1),1)
call DAXPY_(NSS**2,One,LYI,1,MXYZI(:,:,2),1)
call DAXPY_(NSS**2,One,LZI,1,MXYZI(:,:,3),1)

call mma_deallocate(LXI)
call mma_deallocate(LYI)
call mma_deallocate(LZI)

call ZTRNSF(NSS,USOR,USOI,MXYZR(:,:,1),MXYZI(:,:,1))
call ZTRNSF(NSS,USOR,USOI,MXYZR(:,:,2),MXYZI(:,:,2))
call ZTRNSF(NSS,USOR,USOI,MXYZR(:,:,3),MXYZI(:,:,3))

call mma_allocate(ZXYZR,NSS,NSS,3,Label='ZXYZR')
call mma_allocate(ZXYZI,NSS,NSS,3,Label='ZXYZI')
ZXYZR(:,:,:) = Zero
ZXYZR(:,:,:) = Zero

call mma_allocate(ZR,NSS,NSS,Label='ZR')
call mma_allocate(ZI,NSS,NSS,Label='ZI')
call mma_allocate(UZR,NSS,NSS,Label='UZR')
call mma_allocate(UZI,NSS,NSS,Label='UZI')

BFINAL = BSTART+(NBSTEP-1)*BINCRE
TFINAL = TSTART+(NTSTEP-1)*TINCRE

write(u6,*) 'Magnetic flux density range (T):'
write(u6,'(2x,f6.2,a3,f6.2,a4,i4,a6)') BSTART,' - ',BFINAL,' in ',NBSTEP,' steps'
write(u6,*)
write(u6,*) 'Temperature range (K):'
write(u6,'(2x,f6.2,a3,f6.2,a4,i4,a6)') TSTART,' - ',TFINAL,' in ',NTSTEP,' steps'

call mma_allocate(MAGM,9*NBSTEP*NTSTEP,Label='MAGM')

LMSTEP = 0

do IXYZ=1,3

  write(u6,*)
  write(u6,'(3x,a1,3x,8(1x,a12,1x))') 'T','    B'//xyzchr(IXYZ)//' (T)  ','   M (J/T)  ','  Mx (J/T)  ','  My (J/T)  ', &
                                      '  Mz (J/T)  ','Xx'//xyzchr(IXYZ)//' (m3/mol)','Xy'//xyzchr(IXYZ)//' (m3/mol)', &
                                      'Xz'//xyzchr(IXYZ)//' (m3/mol)'
  write(u6,*)

  do IBSTEP=1,NBSTEP
    B = BSTART+BINCRE*(IBSTEP-1)
    ZR(:,:) = Zero
    ZI(:,:) = Zero
    call DAXPY_(NSS**2,Half*B/auToT,MXYZR(:,:,IXYZ),1,ZR,1)
    call DAXPY_(NSS**2,Half*B/auToT,MXYZI(:,:,IXYZ),1,ZI,1)
    do ISS=1,NSS
      HZER = ZR(ISS,ISS)
      ZR(ISS,ISS) = HZER+ENSOR(ISS)
    end do
    call DCOPY_(NSS**2,[Zero],0,UZR,1)
    call DCOPY_(NSS,[One],0,UZR,NSS+1)
    call DCOPY_(NSS**2,[Zero],0,UZI,1)
    call ZJAC(NSS,ZR,ZI,NSS,UZR,UZI)
    do JXYZ=1,3
      call DCOPY_(NSS**2,MXYZR(:,:,JXYZ),1,ZXYZR(:,:,JXYZ),1)
      call DCOPY_(NSS**2,MXYZI(:,:,JXYZ),1,ZXYZI(:,:,JXYZ),1)
      call DSCAL_(NSS**2,-Half,ZXYZR(:,:,JXYZ),1)
      call DSCAL_(NSS**2,-Half,ZXYZI(:,:,JXYZ),1)
      call ZTRNSF(NSS,UZR,UZI,ZXYZR(:,:,JXYZ),ZXYZI(:,:,JXYZ))
    end do
    do ITSTEP=1,NTSTEP
      T = TSTART+TINCRE*(ITSTEP-1)
      RkT = T*BOLTZ
      RMAGM(1) = Zero
      RMAGM(2) = Zero
      RMAGM(3) = Zero
      RPART = Zero
      if (IPGLOB > 2) then
        write(u6,*)
        write(u6,'(2x,a14,3(4x,a4,4x),2x,a6)') 'Energy (cm^-1)','mu_x','mu_y','mu_z','weight'
        write(u6,*)
      end if
      do ISS=1,NSS
        DELTA = ZR(ISS,ISS)-ZR(1,1)
        FACT = exp(-DELTA/RkT)
        RMAGM(1) = RMAGM(1)+ZXYZR(ISS,ISS,1)*FACT
        RMAGM(2) = RMAGM(2)+ZXYZR(ISS,ISS,2)*FACT
        RMAGM(3) = RMAGM(3)+ZXYZR(ISS,ISS,3)*FACT
        RPART = RPART+FACT
        if (IPGLOB > 2) &
          write(u6,'(2x,f14.3,3(1x,f10.6,1x),2x,f6.3)') (ZR(ISS,ISS)-ZR(1,1))*auTocm,ZXYZR(ISS,ISS,1),ZXYZR(ISS,ISS,2), &
                                                        ZXYZR(ISS,ISS,3),FACT
      end do
      if (IPGLOB > 2) write(u6,*)
      RMAGM(1) = (RMAGM(1)/RPART)*AU2JTM
      RMAGM(2) = (RMAGM(2)/RPART)*AU2JTM
      RMAGM(3) = (RMAGM(3)/RPART)*AU2JTM
      RMAGM2 = RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
      RMAGMO = sqrt(RMAGM2)
      do JXYZ=1,3
        LMSTEP = LMSTEP+1
        MAGM(LMSTEP) = RMAGM(JXYZ)
        if (IBSTEP > 1) then
          Chi(JXYZ) = RMAGM(JXYZ)-MAGM(LMSTEP-3*NTSTEP)
          Chi(JXYZ) = Chi(JXYZ)*Rmu0/BINCRE
        end if
      end do
      if (IBSTEP == 1) then
        write(u6,'(1x,f6.2,5(1x,es12.5,1x))') T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
      else
        write(u6,'(1x,f6.2,8(1x,es12.5,1x))') T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3),Chi(1),Chi(2),Chi(3)
      end if
    end do
  end do
end do

call mma_deallocate(MAGM)

write(u6,*)

! powder magnetization, useful in nonlinear cases

if (.not. IFMCAL) goto 810

write(u6,*)
write(u6,*) 'Powder Magnetization'
write(u6,*)
write(u6,'(3x,a1,3x,5(1x,a12,1x))') 'T','    B  (T)  ','   M (J/T)  ','  Mx (J/T)  ','  My (J/T)  ','  Mz (J/T)  '

call mma_allocate(MAGM,3*NBSTEP*NTSTEP,Label='MAGM')
MAGM(:) = Zero

NPHISTEP = int(360.0_wp/BANGRES)
NTHESTEP = int(180.0_wp/BANGRES)

! scale number of points on phi via sin(theta)
NORIENT = 0
do ITHE=1,NTHESTEP+1
  THE = BANGRES*(ITHE-1)*deg2rad
  IPHISTEP = int((NPHISTEP-1)*sin(THE)+1)
  BPHIRES = 360/IPHISTEP
  do IPHI=1,IPHISTEP
    PHI = BPHIRES*(IPHI-1)*deg2rad

    NORIENT = NORIENT+1

    LMSTEP = 0
    !write(u6,*)
    !write(u6,'(1x,2(A6,I4))') ' ITHE ',ITHE,' IPHI',IPHI
    !write(u6,'(6(5x,A4,5x))') ' B  ','THE ','PHI ',' Mx ',' My ',' Mz '
    do IBSTEP=1,NBSTEP
      B = BSTART+BINCRE*(IBSTEP-1)
      BX = B*sin(THE)*cos(PHI)
      BY = B*sin(THE)*sin(PHI)
      BZ = B*cos(THE)
      call DCOPY_(NSS**2,[Zero],0,ZR,1)
      call DAXPY_(NSS**2,Half*BX/auToT,MXYZR(:,:,1),1,ZR,1)
      call DAXPY_(NSS**2,Half*BY/auToT,MXYZR(:,:,2),1,ZR,1)
      call DAXPY_(NSS**2,Half*BZ/auToT,MXYZR(:,:,3),1,ZR,1)
      call DCOPY_(NSS**2,[Zero],0,ZI,1)
      call DAXPY_(NSS**2,Half*BX/auToT,MXYZI(:,:,1),1,ZI,1)
      call DAXPY_(NSS**2,Half*BY/auToT,MXYZI(:,:,2),1,ZI,1)
      call DAXPY_(NSS**2,Half*BZ/auToT,MXYZI(:,:,3),1,ZI,1)
      do ISS=1,NSS
        HZER = ZR(ISS,ISS)
        ZR(ISS,ISS) = HZER+ENSOR(ISS)
      end do
      call DCOPY_(NSS**2,[Zero],0,UZR,1)
      call DCOPY_(NSS,[One],0,UZR,NSS+1)
      call DCOPY_(NSS**2,[Zero],0,UZI,1)
      call ZJAC(NSS,ZR,ZI,NSS,UZR,UZI)
      do IXYZ=1,3
        call DCOPY_(NSS**2,MXYZR(:,:,IXYZ),1,ZXYZR(:,:,IXYZ),1)
        call DCOPY_(NSS**2,MXYZI(:,:,IXYZ),1,ZXYZI(:,:,IXYZ),1)
        call DSCAL_(NSS**2,-Half,ZXYZR(:,:,IXYZ),1)
        call DSCAL_(NSS**2,-Half,ZXYZI(:,:,IXYZ),1)
        call ZTRNSF(NSS,UZR,UZI,ZXYZR(:,:,IXYZ),ZXYZI(:,:,IXYZ))
      end do
      do ITSTEP=1,NTSTEP
        T = TSTART+TINCRE*(ITSTEP-1)
        RkT = T*BOLTZ
        RMAGM(1) = Zero
        RMAGM(2) = Zero
        RMAGM(3) = Zero
        RPART = Zero
        do ISS=1,NSS
          DELTA = ZR(ISS,ISS)-ZR(1,1)
          FACT = exp(-DELTA/RkT)
          RMAGM(1) = RMAGM(1)+ZXYZR(ISS,ISS,1)*FACT
          RMAGM(2) = RMAGM(2)+ZXYZR(ISS,ISS,2)*FACT
          RMAGM(3) = RMAGM(3)+ZXYZR(ISS,ISS,3)*FACT
          RPART = RPART+FACT
        end do
        RMAGM(1) = (RMAGM(1)/RPART)*AU2JTM
        RMAGM(2) = (RMAGM(2)/RPART)*AU2JTM
        RMAGM(3) = (RMAGM(3)/RPART)*AU2JTM
        !write(u6,'(6(1x,es12.5,1x))') B,THE,PHI,RMAGM(1),RMAGM(2),RMAGM(3)
        ! backtransformation in two steps, -phi and -theta
        A = RMAGM(1)
        B = RMAGM(2)
        RMAGM(1) = A*cos(PHI)+B*sin(PHI)
        RMAGM(2) = B*cos(PHI)-A*sin(PHI)
        A = RMAGM(1)
        B = RMAGM(3)
        RMAGM(1) = A*cos(THE)-B*sin(THE)
        RMAGM(3) = B*cos(THE)+A*sin(THE)
        do IXYZ=1,3
          LMSTEP = LMSTEP+1
          MAGM(LMSTEP) = MAGM(LMSTEP)+RMAGM(IXYZ)
        end do
      end do
    end do

  end do
end do

write(u6,*)
LMSTEP = 0
do IBSTEP=1,NBSTEP
  B = BSTART+BINCRE*(IBSTEP-1)
  do ITSTEP=1,NTSTEP
    T = TSTART+TINCRE*(ITSTEP-1)
    do IXYZ=1,3
      LMSTEP = LMSTEP+1
      RMAGM(IXYZ) = MAGM(LMSTEP)/NORIENT
    end do
    RMAGM2 = RMAGM(1)*RMAGM(1)+RMAGM(2)*RMAGM(2)+RMAGM(3)*RMAGM(3)
    RMAGMO = sqrt(RMAGM2)
    write(u6,'(1x,f6.2,5(1x,es12.5,1x))') T,B,RMAGMO,RMAGM(1),RMAGM(2),RMAGM(3)
  end do
end do

call mma_deallocate(MAGM)

810 continue

write(u6,*)

call mma_deallocate(ZR)
call mma_deallocate(ZI)
call mma_deallocate(UZR)
call mma_deallocate(UZI)

call mma_deallocate(MXYZR)
call mma_deallocate(MXYZI)
call mma_deallocate(ZXYZR)
call mma_deallocate(ZXYZI)

900 continue

30 format(5X,A,1X,ES15.8)
31 format(5X,2(1X,A4),6X,A15,1X,A47,1X,A15)
32 format(5X,95('-'))
33 format(5X,2(1X,I4),5X,5(1X,ES15.8))
35 format(5X,31('-'))
36 format(5X,2(1X,I4),6X,15('-'),1X,ES15.8,1X,A15)
37 format(5X,2(1X,I4),6X,15('-'),1X,A15,1X,ES15.8)
38 format(5X,2(1X,I4),6X,F15.6,4(1X,ES15.8))
39 format(5X,2(1X,A4),6X,A15,1X,A15,1X,A15)
40 format(5X,63('-'))
43 format(12X,A8,6(1X,ES15.6))
44 format(20X,6(1X,A15))
49 format(5X,A,1X,ES15.8,1X,A)
!710 format(A4,4I4,4(2X,'('F12.8','F12.8')'))
!720 format(A4,2I4,4(2X,'('F12.8','F12.8')'))

contains

subroutine Allocate_and_Load_electric_dipoles(IFANY)

  integer ISOPR
  integer IPRDX, IPRDY, IPRDZ, IFANY

  IPRDX = 0
  IPRDY = 0
  IPRDZ = 0
  IFANY = 0
  do ISOPR=1,NSOPR
    if ((SOPRNM(ISOPR) == 'MLTPL  1') .and. (SOPRTP(ISOPR) == 'HERMSING')) then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRDX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRDY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRDZ = ISOPR
    end if
  end do
  call mma_allocate(DXR,NSS,NSS,Label='DXR')
  call mma_allocate(DXI,NSS,NSS,Label='DXI')
  call mma_allocate(DYR,NSS,NSS,Label='DYR')
  call mma_allocate(DYI,NSS,NSS,Label='DYI')
  call mma_allocate(DZR,NSS,NSS,Label='DZR')
  call mma_allocate(DZI,NSS,NSS,Label='DZI')
  DXR(:,:) = Zero
  DXI(:,:) = Zero
  DYR(:,:) = Zero
  DYI(:,:) = Zero
  DZR(:,:) = Zero
  DZI(:,:) = Zero
  if (IPRDX > 0) then
    call SMMAT(PROP,DXR,NSS,IPRDX,0)
    call ZTRNSF(NSS,USOR,USOI,DXR,DXI)
  end if
  if (IPRDY > 0) then
    call SMMAT(PROP,DYR,NSS,IPRDY,0)
    call ZTRNSF(NSS,USOR,USOI,DYR,DYI)
  end if
  if (IPRDZ > 0) then
    call SMMAT(PROP,DZR,NSS,IPRDZ,0)
    call ZTRNSF(NSS,USOR,USOI,DZR,DZI)
  end if

end subroutine Allocate_and_Load_electric_dipoles

subroutine Allocate_and_Load_velocities(IFANY)

  integer ISOPR
  integer IPRDX, IPRDY, IPRDZ, IFANY

  IPRDX = 0
  IPRDY = 0
  IPRDZ = 0
  IFANY = 0
  do ISOPR=1,NSOPR
    if (SOPRNM(ISOPR) == 'VELOCITY') then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRDX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRDY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRDZ = ISOPR
    end if
  end do
  call mma_allocate(DXR,NSS,NSS,Label='DXR')
  call mma_allocate(DXI,NSS,NSS,Label='DXI')
  call mma_allocate(DYR,NSS,NSS,Label='DYR')
  call mma_allocate(DYI,NSS,NSS,Label='DYI')
  call mma_allocate(DZR,NSS,NSS,Label='DZR')
  call mma_allocate(DZI,NSS,NSS,Label='DZI')
  DXR(:,:) = Zero
  DXI(:,:) = Zero
  DYR(:,:) = Zero
  DYI(:,:) = Zero
  DZR(:,:) = Zero
  DZI(:,:) = Zero
  if (IPRDX > 0) then
    call SMMAT(PROP,DXR,NSS,IPRDX,0)
    call ZTRNSF(NSS,USOR,USOI,DXR,DXI)
  end if
  if (IPRDY > 0) then
    call SMMAT(PROP,DYR,NSS,IPRDY,0)
    call ZTRNSF(NSS,USOR,USOI,DYR,DYI)
  end if
  if (IPRDZ > 0) then
    call SMMAT(PROP,DZR,NSS,IPRDZ,0)
    call ZTRNSF(NSS,USOR,USOI,DZR,DZI)
  end if

end subroutine Allocate_and_Load_velocities

subroutine Deallocate_electric_dipoles()

  call mma_deallocate(DXR)
  call mma_deallocate(DXI)
  call mma_deallocate(DYR)
  call mma_deallocate(DYI)
  call mma_deallocate(DZR)
  call mma_deallocate(DZI)

end subroutine Deallocate_electric_dipoles

subroutine Allocate_and_Load_magnetic_dipoles(IFANY)

  integer ISOPR
  integer IPRMDX, IPRMDY, IPRMDZ, IFANY

  IPRMDX = 0
  IPRMDY = 0
  IPRMDZ = 0

  IFANY = 0
  do ISOPR=1,NSOPR
    if (SOPRNM(ISOPR) == 'ANGMOM') then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRMDX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRMDY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRMDZ = ISOPR
    end if
  end do
  call mma_allocate(MDXR,NSS,NSS,Label='MDXR')
  call mma_allocate(MDXI,NSS,NSS,Label='MDXI')
  call mma_allocate(MDYR,NSS,NSS,Label='MDYR')
  call mma_allocate(MDYI,NSS,NSS,Label='MDYI')
  call mma_allocate(MDZR,NSS,NSS,Label='MDZR')
  call mma_allocate(MDZI,NSS,NSS,Label='MDZI')
  MDXR(:,:) = Zero
  MDXI(:,:) = Zero
  MDYR(:,:) = Zero
  MDYI(:,:) = Zero
  MDZR(:,:) = Zero
  MDZI(:,:) = Zero
  if (IPRMDX > 0) then
    call SMMAT(PROP,MDXR,NSS,IPRMDX,0)
    call ZTRNSF(NSS,USOR,USOI,MDXR,MDXI)
  end if
  if (IPRMDY > 0) then
    call SMMAT(PROP,MDYR,NSS,IPRMDY,0)
    call ZTRNSF(NSS,USOR,USOI,MDYR,MDYI)
  end if
  if (IPRMDZ > 0) then
    call SMMAT(PROP,MDZR,NSS,IPRMDZ,0)
    call ZTRNSF(NSS,USOR,USOI,MDZR,MDZI)
  end if

end subroutine Allocate_and_Load_magnetic_dipoles

subroutine Deallocate_magnetic_dipoles()

  call mma_deallocate(MDXR)
  call mma_deallocate(MDXI)
  call mma_deallocate(MDYR)
  call mma_deallocate(MDYI)
  call mma_deallocate(MDZR)
  call mma_deallocate(MDZI)

end subroutine Deallocate_magnetic_dipoles

subroutine Allocate_and_Load_Spin_Magnetic_dipoles(IFANY)

  integer ISOPR
  integer IPRSX, IPRSY, IPRSZ, IFANY

  IPRSX = 0
  IPRSY = 0
  IPRSZ = 0

  IFANY = 0
  do ISOPR=1,NSOPR
    if ((SOPRNM(ISOPR) == 'MLTPL  0') .and. (SOPRTP(ISOPR) == 'ANTITRIP')) then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRSX = ISOPR
      if (ISOCMP(ISOPR) == 1) IPRSY = ISOPR
      if (ISOCMP(ISOPR) == 1) IPRSZ = ISOPR
    end if
  end do
  call mma_allocate(SXR,NSS,NSS,Label='SXR')
  call mma_allocate(SXI,NSS,NSS,Label='SXI')
  call mma_allocate(SYR,NSS,NSS,Label='SYR')
  call mma_allocate(SYI,NSS,NSS,Label='SYI')
  call mma_allocate(SZR,NSS,NSS,Label='SZR')
  call mma_allocate(SZI,NSS,NSS,Label='SZI')
  SXR(:,:) = Zero
  SXI(:,:) = Zero
  SYR(:,:) = Zero
  SYI(:,:) = Zero
  SZR(:,:) = Zero
  SZI(:,:) = Zero
  if (IPRSX > 0) then
    call SMMAT(PROP,SXR,NSS,IPRSX,1)
    call ZTRNSF(NSS,USOR,USOI,SXR,SXI)
  end if
  if (IPRSY > 0) then
    call SMMAT(PROP,SYR,NSS,IPRSY,2)
    call ZTRNSF(NSS,USOR,USOI,SYR,SYI)
  end if
  if (IPRSZ > 0) then
    call SMMAT(PROP,SZR,NSS,IPRSZ,3)
    call ZTRNSF(NSS,USOR,USOI,SZR,SZI)
  end if

end subroutine Allocate_and_Load_Spin_Magnetic_dipoles

subroutine Deallocate_Spin_Magnetic_dipoles()

  call mma_deallocate(SXR)
  call mma_deallocate(SXI)
  call mma_deallocate(SYR)
  call mma_deallocate(SYI)
  call mma_deallocate(SZR)
  call mma_deallocate(SZI)

end subroutine Deallocate_Spin_Magnetic_dipoles

subroutine Allocate_and_Load_Spin_Magnetic_Quadrupoles(IFANY)

  integer ISOPR
  integer IPRSXY, IPRSXZ, IPRSYX, IPRSYZ, IPRSZX, IPRSZY, IFANY

  IPRSXY = 0
  IPRSXZ = 0

  IPRSYX = 0
  IPRSYZ = 0

  IPRSZX = 0
  IPRSZY = 0
  IFANY = 0
  do ISOPR=1,NSOPR
    if ((SOPRNM(ISOPR) == 'MLTPL  1') .and. (SOPRTP(ISOPR) == 'ANTITRIP')) then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRSXY = ISOPR
      if (ISOCMP(ISOPR) == 1) IPRSXZ = ISOPR

      if (ISOCMP(ISOPR) == 2) IPRSYX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRSYZ = ISOPR

      if (ISOCMP(ISOPR) == 3) IPRSZX = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRSZY = ISOPR

    end if
  end do
  call mma_allocate(SZXR,NSS,NSS,Label='SZXR')
  call mma_allocate(SZXI,NSS,NSS,Label='SZXI')
  SZXR(:,:) = Zero
  SZXI(:,:) = Zero
  call mma_allocate(SXZR,NSS,NSS,Label='SXZR')
  call mma_allocate(SXZI,NSS,NSS,Label='SXZI')
  SXZR(:,:) = Zero
  SXZI(:,:) = Zero

  call mma_allocate(SXYR,NSS,NSS,Label='SXYR')
  call mma_allocate(SXYI,NSS,NSS,Label='SXYI')
  SXYR(:,:) = Zero
  SXYI(:,:) = Zero
  call mma_allocate(SYXR,NSS,NSS,Label='SYXR')
  call mma_allocate(SYXI,NSS,NSS,Label='SYXI')
  SYXR(:,:) = Zero
  SYXI(:,:) = Zero

  call mma_allocate(SYZR,NSS,NSS,Label='SYZR')
  call mma_allocate(SYZI,NSS,NSS,Label='SYZI')
  SYZR(:,:) = Zero
  SYZI(:,:) = Zero
  call mma_allocate(SZYR,NSS,NSS,Label='SZYR')
  call mma_allocate(SZYI,NSS,NSS,Label='SZYI')
  SZYR(:,:) = Zero
  SZYI(:,:) = Zero
  if (IPRSXY > 0) then
    call SMMAT(PROP,SXYR,NSS,IPRSXY,2)
    call ZTRNSF(NSS,USOR,USOI,SXYR,SXYI)
  end if
  if (IPRSYX > 0) then
    call SMMAT(PROP,SYXR,NSS,IPRSYX,1)
    call ZTRNSF(NSS,USOR,USOI,SYXR,SYXI)
  end if

  if (IPRSXZ > 0) then
    call SMMAT(PROP,SXZR,NSS,IPRSXZ,3)
    call ZTRNSF(NSS,USOR,USOI,SXZR,SXZI)
  end if
  if (IPRSZX > 0) then
    call SMMAT(PROP,SZXR,NSS,IPRSZX,1)
    call ZTRNSF(NSS,USOR,USOI,SZXR,SZXI)
  end if

  if (IPRSYZ > 0) then
    call SMMAT(PROP,SYZR,NSS,IPRSYZ,3)
    call ZTRNSF(NSS,USOR,USOI,SYZR,SYZI)
  end if
  if (IPRSZY > 0) then
    call SMMAT(PROP,SZYR,NSS,IPRSZY,2)
    call ZTRNSF(NSS,USOR,USOI,SZYR,SZYI)
  end if

end subroutine Allocate_and_Load_Spin_Magnetic_Quadrupoles

subroutine Deallocate_Spin_Magnetic_Quadrupoles()

  call mma_deallocate(SXYR)
  call mma_deallocate(SXYI)
  call mma_deallocate(SYXR)
  call mma_deallocate(SYXI)

  call mma_deallocate(SYZR)
  call mma_deallocate(SYZI)
  call mma_deallocate(SZYR)
  call mma_deallocate(SZYI)

  call mma_deallocate(SZXR)
  call mma_deallocate(SZXI)
  call mma_deallocate(SXZR)
  call mma_deallocate(SXZI)

end subroutine Deallocate_Spin_Magnetic_Quadrupoles

subroutine Allocate_and_Load_Electric_Quadrupoles(IFANY)

  integer ISOPR
  integer IPRDXX, IPRDXY, IPRDXZ, IPRDYY, IPRDYZ, IPRDZZ, IFANY

  IPRDXX = 0
  IPRDXY = 0
  IPRDXZ = 0
  IPRDYY = 0
  IPRDYZ = 0
  IPRDZZ = 0

  IFANY = 0
  do ISOPR=1,NSOPR
    if (SOPRNM(ISOPR) == 'MLTPL  2') then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRDXX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRDXY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRDXZ = ISOPR
      if (ISOCMP(ISOPR) == 4) IPRDYY = ISOPR
      if (ISOCMP(ISOPR) == 5) IPRDYZ = ISOPR
      if (ISOCMP(ISOPR) == 6) IPRDZZ = ISOPR
    end if
  end do
  call mma_allocate(QXXR,NSS,NSS,Label='QXXR')
  call mma_allocate(QXXI,NSS,NSS,Label='DXXI')
  call mma_allocate(QXYR,NSS,NSS,Label='DXYR')
  call mma_allocate(QXYI,NSS,NSS,Label='DXYI')
  call mma_allocate(QXZR,NSS,NSS,Label='DXZR')
  call mma_allocate(QXZI,NSS,NSS,Label='DXZI')
  call mma_allocate(QYYR,NSS,NSS,Label='DYYR')
  call mma_allocate(QYYI,NSS,NSS,Label='DYYI')
  call mma_allocate(QYZR,NSS,NSS,Label='DYZR')
  call mma_allocate(QYZI,NSS,NSS,Label='DYZI')
  call mma_allocate(QZZR,NSS,NSS,Label='DZZR')
  call mma_allocate(QZZI,NSS,NSS,Label='DZZI')
  QXXR(:,:) = Zero
  QXXI(:,:) = Zero
  QXYR(:,:) = Zero
  QXYI(:,:) = Zero
  QXZR(:,:) = Zero
  QXZI(:,:) = Zero
  QYYR(:,:) = Zero
  QYYI(:,:) = Zero
  QYZR(:,:) = Zero
  QYZI(:,:) = Zero
  QZZR(:,:) = Zero
  QZZI(:,:) = Zero
  if (IPRDXX > 0) then
    call SMMAT(PROP,QXXR,NSS,IPRDXX,0)
    call ZTRNSF(NSS,USOR,USOI,QXXR,QXXI)
  end if
  if (IPRDXY > 0) then
    call SMMAT(PROP,QXYR,NSS,IPRDXY,0)
    call ZTRNSF(NSS,USOR,USOI,QXYR,QXYI)
  end if
  if (IPRDXZ > 0) then
    call SMMAT(PROP,QXZR,NSS,IPRDXZ,0)
    call ZTRNSF(NSS,USOR,USOI,QXZR,QXZI)
  end if
  if (IPRDYY > 0) then
    call SMMAT(PROP,QYYR,NSS,IPRDYY,0)
    call ZTRNSF(NSS,USOR,USOI,QYYR,QYYI)
  end if
  if (IPRDYZ > 0) then
    call SMMAT(PROP,QYZR,NSS,IPRDYZ,0)
    call ZTRNSF(NSS,USOR,USOI,QYZR,QYZI)
  end if
  if (IPRDZZ > 0) then
    call SMMAT(PROP,QZZR,NSS,IPRDZZ,0)
    call ZTRNSF(NSS,USOR,USOI,QZZR,QZZI)
  end if

end subroutine Allocate_and_Load_Electric_Quadrupoles

subroutine Deallocate_Electric_Quadrupoles()

  call mma_deallocate(QXXR)
  call mma_deallocate(QXXI)
  call mma_deallocate(QXYR)
  call mma_deallocate(QXYI)
  call mma_deallocate(QXZR)
  call mma_deallocate(QXZI)
  call mma_deallocate(QYYR)
  call mma_deallocate(QYYI)
  call mma_deallocate(QYZR)
  call mma_deallocate(QYZI)
  call mma_deallocate(QZZR)
  call mma_deallocate(QZZI)

end subroutine Deallocate_Electric_Quadrupoles

subroutine Allocate_and_Load_Magnetic_Quadrupoles(IFANY)

  integer ISOPR
  integer IPRDZX, IPRDYX, IPRDZY, IFANY

  IPRDXY = 0
  IPRDXZ = 0
  IPRDYX = 0
  IPRDYZ = 0
  IPRDZX = 0
  IPRDZY = 0

  IFANY = 0
  do ISOPR=1,NSOPR
    if (SOPRNM(ISOPR) == 'OMQ') then
      IFANY = 1
      if (ISOCMP(ISOPR) == 2) IPRDXY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRDXZ = ISOPR

      if (ISOCMP(ISOPR) == 4) IPRDYX = ISOPR
      if (ISOCMP(ISOPR) == 6) IPRDYZ = ISOPR

      if (ISOCMP(ISOPR) == 7) IPRDZX = ISOPR
      if (ISOCMP(ISOPR) == 8) IPRDZY = ISOPR
    end if
  end do

  call mma_allocate(MQXYR,NSS,NSS,Label='MQXYR')
  call mma_allocate(MQXYI,NSS,NSS,Label='MQXYI')
  call mma_allocate(MQYXR,NSS,NSS,Label='MQYXR')
  call mma_allocate(MQYXI,NSS,NSS,Label='MQYXI')
  call mma_allocate(MQXZR,NSS,NSS,Label='MQXZR')
  call mma_allocate(MQXZI,NSS,NSS,Label='MQXZI')
  call mma_allocate(MQZXR,NSS,NSS,Label='MQZXR')
  call mma_allocate(MQZXI,NSS,NSS,Label='MQZXI')
  call mma_allocate(MQYZR,NSS,NSS,Label='MQYZR')
  call mma_allocate(MQYZI,NSS,NSS,Label='MQYZI')
  call mma_allocate(MQZYR,NSS,NSS,Label='MQZYR')
  call mma_allocate(MQZYI,NSS,NSS,Label='MQZYI')
  MQXYR(:,:) = Zero
  MQXYI(:,:) = Zero
  MQYXR(:,:) = Zero
  MQYXI(:,:) = Zero
  MQXZR(:,:) = Zero
  MQXZI(:,:) = Zero
  MQZXR(:,:) = Zero
  MQZXI(:,:) = Zero
  MQYZR(:,:) = Zero
  MQYZI(:,:) = Zero
  MQZYR(:,:) = Zero
  MQZYI(:,:) = Zero
  if (IPRDXY > 0) then
    call SMMAT(PROP,MQXYR,NSS,IPRDXY,0)
    call ZTRNSF(NSS,USOR,USOI,MQXYR,MQXYI)
  end if
  if (IPRDYX > 0) then
    call SMMAT(PROP,MQYXR,NSS,IPRDYX,0)
    call ZTRNSF(NSS,USOR,USOI,MQYXR,MQYXI)
  end if

  if (IPRDXZ > 0) then
    call SMMAT(PROP,MQXZR,NSS,IPRDXZ,0)
    call ZTRNSF(NSS,USOR,USOI,MQXZR,MQXZI)
  end if
  if (IPRDZX > 0) then
    call SMMAT(PROP,MQZXR,NSS,IPRDZX,0)
    call ZTRNSF(NSS,USOR,USOI,MQZXR,MQZXI)
  end if

  if (IPRDYZ > 0) then
    call SMMAT(PROP,MQYZR,NSS,IPRDYZ,0)
    call ZTRNSF(NSS,USOR,USOI,MQYZR,MQYZI)
  end if
  if (IPRDZY > 0) then
    call SMMAT(PROP,MQZYR,NSS,IPRDZY,0)
    call ZTRNSF(NSS,USOR,USOI,MQZYR,MQZYI)
  end if

end subroutine Allocate_and_Load_Magnetic_Quadrupoles

subroutine Deallocate_Magnetic_Quadrupoles()

  call mma_deallocate(MQXYR)
  call mma_deallocate(MQXYI)
  call mma_deallocate(MQYXR)
  call mma_deallocate(MQYXI)
  call mma_deallocate(MQXZR)
  call mma_deallocate(MQXZI)
  call mma_deallocate(MQZXR)
  call mma_deallocate(MQZXI)
  call mma_deallocate(MQYZR)
  call mma_deallocate(MQYZI)
  call mma_deallocate(MQZYR)
  call mma_deallocate(MQZYI)

end subroutine Deallocate_Magnetic_Quadrupoles

subroutine Allocate_and_Load_Octupoles(IFANY)

  integer ISOPR
  integer IPRDZZX, IPRDZZY, IPRDZZZ, IPRDXXX, IPRDXXY, IPRDXXZ, IPRDYYX, IPRDYYY, IPRDYYZ, IFANY

  ! This is a real symmetric rank 3 tensor so only 10 and not 27 is needed
  ! The order which comes in
  IPRDXXX = 0 !
  IPRDXXY = 0 !
  IPRDXXZ = 0 !

  !IPRDXYX = 0
  !IPRDXYY = 0 ! YYX These are the same due to symmetry
  !IPRDXYZ = 0 ! Not present

  !IPRDXZX = 0
  !IPRDXZY = 0
  !IPRDXZZ = 0 ! ZZX

  !IPRDYXX = 0
  !IPRDYXY = 0
  !IPRDYXZ = 0

  IPRDYYX = 0 ! Taking the XYY order
  IPRDYYY = 0 !
  IPRDYYZ = 0 !

  !IPRDYZX = 0
  !IPRDYZY = 0
  !IPRDYZZ = 0 ! ZZY

  !IPRDZXX = 0
  !IPRDZXY = 0
  !IPRDZXZ = 0

  !IPRDZYX = 0
  !IPRDZYY = 0
  !IPRDZYZ = 0

  IPRDZZX = 0 ! Taking order from XZZ
  IPRDZZY = 0 ! Taking order from YZZ
  IPRDZZZ = 0 !

  IFANY = 0
  do ISOPR=1,NSOPR
    if (SOPRNM(ISOPR) == 'MLTPL  3') then
      IFANY = 1
      if (ISOCMP(ISOPR) == 1) IPRDXXX = ISOPR
      if (ISOCMP(ISOPR) == 2) IPRDXXY = ISOPR
      if (ISOCMP(ISOPR) == 3) IPRDXXZ = ISOPR
      if (ISOCMP(ISOPR) == 4) IPRDYYX = ISOPR ! Changed from XYY
      !if (ISOCMP(ISOPR) == 5) IPRDXYZ = ISOPR
      if (ISOCMP(ISOPR) == 6) IPRDZZX = ISOPR ! Changed from XZZ
      if (ISOCMP(ISOPR) == 7) IPRDYYY = ISOPR
      if (ISOCMP(ISOPR) == 8) IPRDYYZ = ISOPR
      if (ISOCMP(ISOPR) == 9) IPRDZZY = ISOPR ! Changed from YZZ
      if (ISOCMP(ISOPR) == 10) IPRDZZZ = ISOPR

    end if
  end do
  call mma_allocate(DXXXR,NSS,NSS,Label='DXXXR')
  call mma_allocate(DXXXI,NSS,NSS,Label='DXXXI')
  call mma_allocate(DXXYR,NSS,NSS,Label='DXXYR')
  call mma_allocate(DXXYI,NSS,NSS,Label='DXXYI')
  call mma_allocate(DXXZR,NSS,NSS,Label='DXXZR')
  call mma_allocate(DXXZI,NSS,NSS,Label='DXXZI')
  DXXXR(:,:) = Zero
  DXXXI(:,:) = Zero
  DXXYR(:,:) = Zero
  DXXYI(:,:) = Zero
  DXXZR(:,:) = Zero
  DXXZI(:,:) = Zero
  call mma_allocate(DYYXR,NSS,NSS,Label='DYYXR')
  call mma_allocate(DYYXI,NSS,NSS,Label='DYYXI')
  call mma_allocate(DYYYR,NSS,NSS,Label='DYYYR')
  call mma_allocate(DYYYI,NSS,NSS,Label='DYYYI')
  call mma_allocate(DYYZR,NSS,NSS,Label='DYYZR')
  call mma_allocate(DYYZI,NSS,NSS,Label='DYYZI')
  DYYXR(:,:) = Zero
  DYYXI(:,:) = Zero
  DYYYR(:,:) = Zero
  DYYYI(:,:) = Zero
  DYYZR(:,:) = Zero
  DYYZI(:,:) = Zero
  call mma_allocate(DZZXR,NSS,NSS,Label='DZZXR')
  call mma_allocate(DZZXI,NSS,NSS,Label='DZZXI')
  call mma_allocate(DZZYR,NSS,NSS,Label='DZZYR')
  call mma_allocate(DZZYI,NSS,NSS,Label='DZZYI')
  call mma_allocate(DZZZR,NSS,NSS,Label='DZZZR')
  call mma_allocate(DZZZI,NSS,NSS,Label='DZZZI')
  DZZXR(:,:) = Zero
  DZZXI(:,:) = Zero
  DZZYR(:,:) = Zero
  DZZYI(:,:) = Zero
  DZZZR(:,:) = Zero
  DZZZI(:,:) = Zero
  if (IPRDXXX > 0) then
    call SMMAT(PROP,DXXXR,NSS,IPRDXXX,0)
    call ZTRNSF(NSS,USOR,USOI,DXXXR,DXXXI)
  end if
  if (IPRDXXY > 0) then
    call SMMAT(PROP,DXXYR,NSS,IPRDXXY,0)
    call ZTRNSF(NSS,USOR,USOI,DXXYR,DXXYI)
  end if
  if (IPRDXXZ > 0) then
    call SMMAT(PROP,DXXZR,NSS,IPRDXXZ,0)
    call ZTRNSF(NSS,USOR,USOI,DXXZR,DXXZI)
  end if

  if (IPRDYYX > 0) then
    call SMMAT(PROP,DYYXR,NSS,IPRDYYX,0)
    call ZTRNSF(NSS,USOR,USOI,DYYXR,DYYXI)
  end if
  if (IPRDYYY > 0) then
    call SMMAT(PROP,DYYYR,NSS,IPRDYYY,0)
    call ZTRNSF(NSS,USOR,USOI,DYYYR,DYYYI)
  end if
  if (IPRDYYZ > 0) then
    call SMMAT(PROP,DYYZR,NSS,IPRDYYZ,0)
    call ZTRNSF(NSS,USOR,USOI,DYYZR,DYYZI)
  end if

  if (IPRDZZX > 0) then
    call SMMAT(PROP,DZZXR,NSS,IPRDZZX,0)
    call ZTRNSF(NSS,USOR,USOI,DZZXR,DZZXI)
  end if
  if (IPRDZZY > 0) then
    call SMMAT(PROP,DZZYR,NSS,IPRDZZY,0)
    call ZTRNSF(NSS,USOR,USOI,DZZYR,DZZYI)
  end if
  if (IPRDZZZ > 0) then
    call SMMAT(PROP,DZZZR,NSS,IPRDZZZ,0)
    call ZTRNSF(NSS,USOR,USOI,DZZZR,DZZZI)
  end if

end subroutine Allocate_and_Load_Octupoles

subroutine deallocate_Octupoles()

  call mma_deallocate(DXXXR)
  call mma_deallocate(DXXXI)
  call mma_deallocate(DXXYR)
  call mma_deallocate(DXXYI)
  call mma_deallocate(DXXZR)
  call mma_deallocate(DXXZI)
  call mma_deallocate(DYYXR)
  call mma_deallocate(DYYXI)
  call mma_deallocate(DYYYR)
  call mma_deallocate(DYYYI)
  call mma_deallocate(DYYZR)
  call mma_deallocate(DYYZI)
  call mma_deallocate(DZZXR)
  call mma_deallocate(DZZXI)
  call mma_deallocate(DZZYR)
  call mma_deallocate(DZZYI)
  call mma_deallocate(DZZZR)
  call mma_deallocate(DZZZI)

end subroutine deallocate_Octupoles

end subroutine PRPROP
