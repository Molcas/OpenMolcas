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

subroutine SINGLE_ANISO2(nH,nT,nTempMagn,nDir,nDirZee,nss,nstate,nMult,input_file_name,Ifrestart,IReturn,GRAD)

use Constants, only: Zero, Two, cZero
use Definitions, only: wp, u6

implicit none
#include "cntrl_sa.fh"
#include "stdalloc.fh"
integer :: mem, RtoB, CtoB, ItoB
integer, intent(in) :: nss, nstate
integer :: i, j, iReturn
integer :: idim, Ifunct, imltpl, iprint
integer :: axisoption
integer :: ndimcf, ldimcf
integer :: input_to_read, nK, MG, nm
integer :: encut_definition
integer :: ncut, nDirTot
integer :: nBlock
integer :: imanIfold
integer :: i1, i2, ldim2, lDIM
integer :: AngPoints
integer :: nlanth
real(kind=8) :: zmagn(3,3)
real(kind=8) :: HMIN, HMAX
real(kind=8) :: cryst(6), coord(3)
real(kind=8) :: encut_rate, em
real(kind=8) :: zJ, thrs
integer, allocatable :: multiplicity(:)
!---g-tens-------------------
integer, intent(in) :: nMult
integer, allocatable :: ndim(:)
real(kind=8), allocatable :: gtens(:,:), maxes(:,:,:)
!---M-------------------
integer :: nH, nTempMagn
!---MVEC and ZEEM-------------------
integer :: nDir, nDirZee
integer, allocatable :: LuZee(:)
real(kind=8), allocatable :: dir_weight(:,:)
real(kind=8), allocatable :: dirX(:), dirY(:), dirZ(:)
!---XT-------------------
integer :: nT
real(kind=8), allocatable :: Texp(:)
real(kind=8), allocatable :: chit_exp(:)
real(kind=8) :: xfield
real(kind=8) :: tmin, tmax
!---Oscillator strength----------
!real(kind=8) :: F, Fx, Fy, Fz, AT, Ax, Ay, Az, AF, dnrm, dE
real(kind=8) :: H_torq, T_torq
!----BIG ARRAYS------------------
real(kind=8), allocatable :: eso(:), eso_au(:)
real(kind=8), allocatable :: esfs(:), esfs_au(:)
real(kind=8), allocatable :: t(:)
real(kind=8), allocatable :: XTexp(:)
real(kind=8), allocatable :: XT_no_field(:)
real(kind=8), allocatable :: hexp(:)
real(kind=8), allocatable :: magn_exp(:,:)
real(kind=8), allocatable :: angmom(:,:,:)
real(kind=8), allocatable :: eDmom(:,:,:)
real(kind=8), allocatable :: amfi(:,:,:)
real(kind=8), allocatable :: TempMagn(:)
complex(kind=8), allocatable :: MM(:,:,:), MS(:,:,:), HSO(:,:), ML(:,:,:), DM(:,:,:), U(:,:)
character(len=180), intent(in) :: input_file_name
logical :: poly_file
logical :: ifrestart
logical :: Do_structure_abc
logical :: compute_magnetization
logical :: m_paranoid
logical :: compute_Mdir_vector
logical :: zeeman_energy
logical :: compute_torque
logical :: smagn
logical :: compute_cf
logical :: compute_g_tensors
logical :: compute_barrier
logical :: tinput, hinput
logical :: GRAD
logical :: DoPlot
logical :: DBG
integer :: l
integer :: nss2, nstate2

DBG = .false.

!-----------------------------------------------------------------------
! Allocate memory for all arrays:
!-----------------------------------------------------------------------
if (dbg) then
  write(u6,*) 'S_A2:: Memory Allocation Parameters'
  write(u6,*) 'S_A2:: nss             =',nss
  write(u6,*) 'S_A2:: nstate          =',nstate
  write(u6,*) 'S_A2:: nH              =',nH
  write(u6,*) 'S_A2:: nT              =',nT
  write(u6,*) 'S_A2:: nTempMagn       =',nTempMagn
  write(u6,*) 'S_A2:: nMult           =',nMult
  write(u6,*) 'S_A2:: nDir            =',nDir
  write(u6,*) 'S_A2:: nDirZee         =',nDirZee
  write(u6,*) 'S_A2:: input_file_name =',input_file_name
  write(u6,*) 'S_A2:: GRAD            =',GRAD
end if

mem = 0
RtoB = 8
CtoB = 16
ItoB = 8

! spin free energies
call mma_allocate(esfs,nstate,'esfs')
call mma_allocate(esfs_au,nstate,'esfs_au')
call dcopy_(nstate,[Zero],0,ESFS,1)
call dcopy_(nstate,[Zero],0,ESFS_AU,1)
mem = mem+2*nstate*RtoB
! angular momentum
call mma_allocate(ANGMOM,3,nstate,nstate,'angmom')
call dcopy_(3*nstate*nstate,[Zero],0,ANGMOM,1)
mem = mem+3*nstate*nstate*RtoB
! electric dipole moment
call mma_allocate(EDMOM,3,nstate,nstate,'edmom')
call dcopy_(3*nstate*nstate,[Zero],0,EDMOM,1)
mem = mem+3*nstate*nstate*RtoB
! amfi integrals
call mma_allocate(AMFI,3,nstate,nstate,'amfi')
call dcopy_(3*nstate*nstate,[Zero],0,AMFI,1)
mem = mem+3*nstate*nstate*RtoB
! multiplicity of each state
call mma_allocate(multiplicity,nstate,'multiplicity')
call icopy(nstate,[0],0,multiplicity,1)
mem = mem+nstate*ItoB
! allocated memory counter
if (dbg) write(u6,'(A,I16)') 'mem 1 =',mem
! spin orbit energies
call mma_allocate(eso,nss,'eso')
call mma_allocate(eso_au,nss,'eso_au')
call dcopy_(nss,[Zero],0,eso,1)
call dcopy_(nss,[Zero],0,eso_au,1)
mem = mem+2*nss*RtoB
! spin orbit eigenstates
call mma_allocate(U,nss,nss,'U')
call zcopy_(nss*nss,[cZero],0,U,1)
mem = mem+nss*nss*CtoB
! spin orbit hamiltonian
call mma_allocate(HSO,nss,nss,'HSO')
call zcopy_(nss*nss,[cZero],0,HSO,1)
mem = mem+nss*nss*CtoB
! magnetic moment
call mma_allocate(MM,3,nss,nss,'MM')
call zcopy_(3*nss*nss,[cZero],0,MM,1)
mem = mem+3*nss*nss*CtoB
! spin moment
call mma_allocate(MS,3,nss,nss,'MS')
call zcopy_(3*nss*nss,[cZero],0,MS,1)
mem = mem+3*nss*nss*CtoB
! orbital mooment
call mma_allocate(ML,3,nss,nss,'ML')
call zcopy_(3*nss*nss,[cZero],0,ML,1)
mem = mem+3*nss*nss*CtoB
! electric dipole moment
call mma_allocate(DM,3,nss,nss,'DM')
call zcopy_(3*nss*nss,[cZero],0,DM,1)
mem = mem+3*nss*nss*CtoB
! allocated memory counter
if (dbg) write(u6,'(A,I16)') 'mem 2 =',mem
! experimental magnetic field points
call mma_allocate(Hexp,nH,'Hexp')
call dcopy_(nH,[Zero],0,Hexp,1)
mem = mem+nH*RtoB
! experiemental magnetization
call mma_allocate(magn_exp,nH,nTempMagn,'magn_exp')
call dcopy_(nH*nTempMagn,[Zero],0,magn_exp,1)
mem = mem+nH*nTempMagn*RtoB
! temperature points for magnetization
call mma_allocate(TempMagn,nTempMagn,'TempMagn')
call dcopy_(nTempMagn,[Zero],0,TempMagn,1)
mem = mem+nTempMagn*RtoB
! allocated memory counter
!if (dbg) write(u6,'(A,I16)') 'mem 3 =',mem
!call mma_allocate(Hexp,nH,'Hexp')
!mem = mem+nH*RtoB
!call dcopy_(nH,[Zero],0,Hexp,1)
!call mma_allocate(magn_exp,nH,0,'magn_exp')
!call mma_allocate(TempMagn,0,'TempMagn')
! allocated memory counter
!if (dbg) write(u6,'(A,I16)') 'mem 3 =',mem
!call mma_allocate(Hexp,0,'Hexp')
!call mma_allocate(magn_exp,0,nTempMagn,'magn_exp')
!call mma_allocate(TempMagn,nTempMagn,'TempMagn')
!call dcopy_(nTempMagn,[Zero],0,TempMagn,1)
!mem = mem+nTempMagn*RtoB
! allocated memory counter
!if (dbg) write(u6,'(A,I16)') 'mem 3 =',mem
!call mma_allocate(Hexp,0,'Hexp')
!call mma_allocate(magn_exp,0,0,'magn_exp')
!call mma_allocate(TempMagn,0,'TempMagn')

! dimensions of pseudospins
call mma_allocate(ndim,nMult,'ndim')
call icopy(nMult,[0],0,ndim,1)
mem = mem+nMult*ItoB
! temperature points for magnetization
call mma_allocate(gtens,nMult,3,'gtens')
call dcopy_(3*nMult,[Zero],0,gtens,1)
mem = mem+3*nMult*RtoB
! temperature points for magnetization
call mma_allocate(maxes,nMult,3,3,'maxes')
call dcopy_(3*3*nMult,[Zero],0,maxes,1)
mem = mem+3*3*nMult*RtoB
! allocated memory counter
if (dbg) write(u6,'(A,I16)') 'mem 4 =',mem
call mma_allocate(T,(nT+nTempMagn),'Temperature')
call dcopy_((nT+nTempMagn),[Zero],0,T,1)
mem = mem+(nT+nTempMagn)*RtoB
call mma_allocate(XTexp,(nT+nTempMagn),'XTexp')
call dcopy_((nT+nTempMagn),[Zero],0,XTexp,1)
mem = mem+(nT+nTempMagn)*RtoB
call mma_allocate(XT_no_field,(nT+nTempMagn),'XT_no_field')
call dcopy_((nT+nTempMagn),[Zero],0,XT_no_field,1)
mem = mem+(nT+nTempMagn)*RtoB
! allocated memory counter
if (dbg) write(u6,'(A,I16)') 'mem 5 =',mem

! unit numbers for the files with Zeeman energies
call mma_allocate(LuZee,nDirZee,'LUZee')
call icopy(nDirZee,[0],0,LuZee,1)
mem = mem+nDirZee*ItoB
! directions for applied field for Zeeman states
call mma_allocate(dir_weight,nDirZee,3,'dir_weight')
call dcopy_(3*nDirZee,[Zero],0,dir_weight,1)
mem = mem+3*nDirZee*RtoB
! allocated memory counter
!if (dbg) write(u6,'(A,I16)') 'mem 6 =',mem
!call mma_allocate(LuZee,0,'LUZee')
!call mma_allocate(dir_weight,0,3,'dir_weight')

! magnetization vectors
call mma_allocate(dirX,nDir,'dirX')
call mma_allocate(dirY,nDir,'dirY')
call mma_allocate(dirZ,nDir,'dirZ')
mem = mem+3*nDir*RtoB
call dcopy_(nDir,[Zero],0,dirX,1)
call dcopy_(nDir,[Zero],0,dirY,1)
call dcopy_(nDir,[Zero],0,dirZ,1)
! allocated memory counter
!if (dbg) Write(u6,'(A,I16)') 'mem 7 =',mem
!call mma_allocate(dirX,0,'dirX')
!call mma_allocate(dirY,0,'dirY')
!call mma_allocate(dirZ,0,'dirZ')
! T expeirimental given by user in the input
call mma_allocate(Texp,nT,'Texp')
call dcopy_(nT,[Zero],0,Texp,1)
mem = mem+nT*RtoB
! XT expeirimental given by user in the input
call mma_allocate(chit_exp,nT,'chit_exp')
call dcopy_(nT,[Zero],0,chit_exp,1)
mem = mem+nT*RtoB
! allocated memory counter
if (dbg) write(u6,'(A,I16)') 'mem 8 =',mem

write(u6,'(A,I16,A)') 'The code allocated initially:',mem,' bytes of memory for this run.'
call xFlush(u6)
!-----------------------------------------------------------------------
IReturn = 0
IPRINT = 2
lDIM = 1
iDIM = 1
NM = 0
EM = Zero
POLY_FILE = .false.
compute_CF = .false.
axisoption = 1
nDIMcf = 1
lDIMcf = 1
H_torq = 0.1_wp ! in tesla
T_torq = Two    ! in K
! read the input
if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter readin_single'
call readin_single(iprint,nmult,ndim,ldim,ndimcf,ldimcf,nlanth,axisoption,poly_file,Ifrestart,input_to_read,nk,mg,zmagn, &
                   Do_structure_abc,cryst,coord,encut_definition,compute_g_tensors,compute_CF,nDirTot,nss,nstate, &
                   compute_magnetization,compute_torque,smagn,tinput,hinput,compute_Mdir_vector,zeeman_energy,LUZee,doplot, &
                   encut_rate,ncut,nTempMagn,TempMagn,m_paranoid,compute_barrier,nBlock,AngPoints,input_file_name,nT,nH,texp, &
                   chit_exp,zJ,hexp,magn_exp,hmin,hmax,nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,thrs,H_torq,T_torq)
if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit readin_single'

if (ifrestart) then
  ! if restart, fetch "big data" from the input file
  if (input_to_read == 1) then
    call read_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)

  else if (input_to_read == 2) then
    ! get the information from formatted aniso.input file:
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter read_formatted_aniso'
    if (DBG) write(u6,*) 'SA:',input_file_name
    nss2 = nss
    nstate2 = nstate
    call read_formatted_aniso(input_file_name,nss2,nstate2,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit  read_formatted_aniso'

  else if (input_to_read == 3) then
    ! get the information from RASSI-HDF5 file:
    if (DBG) write(u6,*) 'SA:',input_file_name
#   ifdef _HDF5_
    call read_hdf5_all(input_file_name,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)
#   else
    call WarningMessage(2,'File '//trim(input_file_name)//' cannot be opened. Molcas was compiled without HDF5 option.')
    call Quit_OnUserError()
#   endif

  else if (input_to_read == 4) then
    if (DBG) write(u6,*) 'SA:',input_file_name
    ! get the information from formatted aniso.input file:
    nss2 = nss
    nstate2 = nstate
    call read_formatted_aniso_old(input_file_name,nss2,nstate2,multiplicity,eso,MM,MS,ML)

  else if (input_to_read == 6) then !  using DATA keyword

    if (DBG) write(u6,*) 'SA:',input_file_name
    ! get the information from formatted new aniso-data file:
    nss2 = nss
    nstate2 = nstate
    call read_formatted_new_aniso(input_file_name,nss2,nstate2,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO,eso_au, &
                                  esfs_au)
  end if ! input_to_read

else
  ! ifrestart = .false., i.e. usual S-A calculation
  call fetch_data_RunFile_all(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO,eso_au,esfs_au)
  write(u6,'(A)') 'AFTER fetch_data_RunFile_all'
  call xFlush(u6)
  if (DBG) then
    write(u6,'(A)') 'SA: ANGMOM(x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(ANGMOM(l,i,j),l=1,3)
      end do
    end do
    write(u6,'(/)')
    write(u6,'(A)') 'SA: EDMOM(x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(EDMOM(l,i,j),l=1,3)
      end do
    end do
    write(u6,'(/)')
    write(u6,'(A)') 'SA: AMFI(x,y,z)'
    do i=1,nstate
      do j=1,nstate
        write(u6,'(2i4,A,3ES24.14)') i,j,' |',(AMFI(l,i,j),l=1,3)
      end do
    end do
    write(u6,'(/)')
    write(u6,'(A)') 'SA: HSO(i,j)'
    do i=1,nss
      do j=1,nss
        write(u6,'(2i4,A,4ES24.14)') i,j,' |',HSO(i,j),HSO(j,i)
      end do
    end do
  end if ! DBG

end if ! Ifrestart

! print some input data in the beginning of the output:
! so that the user knows which ws the input ...
write(u6,'(A)') 'LOW-LYING SPIN-ORBIT ENERGIES:'
do i=1,nss
  write(u6,'(A,I4,A,F25.14)') 'ENERGY OF THE SPIN-ORBIT STATE (',i,') =',ESO(i)
end do
write(u6,'(A)') 'LOW-LYING SPIN-FREE ENERGIES:'
do i=1,nstate
  write(u6,'(A,I4,A,F25.14)') 'ENERGY OF THE SPIN-FREE STATE  (',i,') =',ESFS(i)
end do
lDIM2 = 1 ! the same as the defult for lDIM
do i=2,nstate
  if (abs(ESFS(i)-ESFS(i-1)) < 20.0_wp) then
    lDIM2 = lDIM2+1
  else
    goto 107
  end if
end do
107 continue
! set the lDIM:
if (lDIM2 > lDIM) lDIM = lDIM2

!----- input processing finished -----
! save some important data, regardless of the following execution
! -- binary $Project.aniso
call write_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)
! ASCII -- anisoinput:
call write_formatted_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)

if (.not. ifrestart) then
  ! ASCII -- new_aniso file format:
  call write_new_formatted_aniso(nss,nstate,multiplicity,eso_au,esfs_au,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)
end if

!----- compute various properties ----------|
! calculation of magnetic Hamiltonians:
!IDIM = 0
IFUNCT = 0
if (compute_g_tensors .and. (nMult > 0)) then
  do IMLTPL=1,NMULT
    IReturn = 0
    if (ndim(imltpl) == 1) then
      write(u6,'(5X,A,I2,A)') 'THE DIMENSION OF THE ',IMLTPL,' MULTIPLET IS 1.'
      write(u6,'(5X,A)') 'THERE IS NO G TENSOR FOR EFFECTIVE S = 0.'
      go to 10
    end if

    i1 = 1+Ifunct
    i2 = ndim(imltpl)+Ifunct

    if (i2 <= nss) then
      if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter g_high',IMLTPL
      call g_high(eso(i1:i2),GRAD,MS(:,i1:i2,i1:i2),MM(:,i1:i2,i1:i2),imltpl,ndim(imltpl),Do_structure_abc,cryst,coord, &
                  gtens(imltpl,:),maxes(imltpl,:,:),iprint)
      if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit g_high',IMLTPL
      call Add_Info('GTENS_MAIN',gtens(imltpl,:),3,4)
    end if
10  continue

    IFUNCT = IFUNCT+NDIM(IMLTPL)
    if (IReturn /= 0) call ABEnd()
  end do
end if

!----------------------------------------------------------------------|
! >> AB INITIO CRYSTAL-FIELD <<

if (compute_CF .and. (nDIMCF > 0)) then
  if ((axisoption == 1) .and. (nMult > 0) .and. (nDIM(1) > 1)) then
    iDIM = NDIM(1)
    lDIM = NDIM(1)
  else
    iDIM = nDIMCF
    lDIM = lDIMCF
  end if

  ! compute the CF of the ground |J,MJ> multiplet:
  if ((nDIMCF > 1) .and. (nDIMCF <= nss)) then
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter CF',nDIMCF
    call CRYSTALFIELD(ESO(1:nDIMCF),MM(:,1:nDIMCF,1:nDIMCF),MS(:,1:nDIMCF,1:nDIMCF),nDIMcf,iDIM,nlanth,zmagn,axisoption,GRAD,iPrint)
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit CF',nDIMCF
  end if
end if

if (DBG) then
  write(u6,*) 'SINGLE_ANISO2:: nlanth=',nlanth
  write(u6,*) 'SINGLE_ANISO2:: lDIMCF=',lDIMCF
  write(u6,*) 'SINGLE_ANISO2:: nstate=',nstate
end if

if (compute_CF .and. (lDIMCF > 0) .and. (lDIMCF <= nstate)) then
  if ((axisoption == 1) .and. (nMult > 0) .and. (nDIM(1) > 1)) then
    iDIM = NDIM(1)
    lDIM = NDIM(1)
  else
    iDIM = nDIMCF
    lDIM = lDIMCF
  end if
  ! compute the CF of the ground |L,ML> term:
  if (.not. (ifrestart .and. (input_to_read == 4))) then
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter t-CF',lDIMCF
    call termCF(angmom(:,1:lDIMCF,1:lDIMCF),AMFI(:,1:lDIMCF,1:lDIMCF),esfs(1:lDIMCF),lDIMCF,lDIM,zmagn,axisoption,nlanth,iPrint)
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit t-CF',lDIMCF
  end if !ifrestart
end if !compute_CF

!----------------------------------------------------------------------|
! >> AB INITIO BLOCKING BARRIER FOR SMMs <<
if (compute_barrier) then
  if (nBlock /= 0) then
    imanifold = 1

    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter barrier',nBlock
    call BARRIER(nBlock,MM(:,1:nBlock,1:nBlock),eso(1:nBlock),imanIfold,nMult,nDim,doplot,iPrint)
    if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit barrier',nBlock

  else
    write(u6,'(A)') 'nBlock parameter is not defined. '
    write(u6,'(A)') 'Did you specify the MLTP keyword in the input?'
    write(u6,'(A)') 'If the problem persists, please, submit a bug report.'
  end if
end if

!----------------------------------------------------------------------|

call set_T(nT,nTempMagn,TINPUT,TempMagn,Tmin,Tmax,chit_exp,Texp,T,XTexp)

if (compute_torque .or. compute_magnetization .or. (Xfield /= Zero)) then
  ! set the NM- number of states to be exactly diagonalized
  ! in Zeeman Interaction
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter set_nm'
  call set_nm(nss,ncut,encut_definition,nk,mg,nTempMagn,hmax,ESO,encut_rate,TempMagn,nM,EM,dbg)
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  exit set_nm'

end if

!----------------------------------------------------------------------|
! >> MAGNETIC SUSCEPTIBILITY <<
iReturn = 0
if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter SUSCEPTIBILITY'
call SUSCEPTIBILITY(NSS,ESO,MS,MM,nT,nTempMagn,T,tmin,tmax,XTexp,zJ,tinput,XT_no_field,doplot,iPrint,mem)
if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit SUSCEPTIBILITY'

if (Xfield /= Zero) then
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter XT_dMoverdH_single'
  ! nM = nss
  call XT_dMoverdH_single(nss,nTempMagn,nT,nss,Tmin,Tmax,XTexp,ESO,T,zJ,Xfield,EM,MM,MS,XT_no_field,tinput,smagn,mem,DoPlot)
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit XT_dMoverdH_single'
end if
!----------------------------------------------------------------------|
! >> TORQUE <<
if (compute_torque) then

  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter TORQUE'
  ! nM = nss
  call torque(Nss,Nss,AngPoints,EM,eso,mm,ms,zJ,thrs,mem,m_paranoid,smagn,H_torq,T_torq,zmagn,dbg)
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit TORQUE'

end if
!----------------------------------------------------------------------|
! >> MAGNETIZATION <<

if (compute_magnetization) then

  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Enter magnetization'
  call magnetization(nss,nM,nTempMagn,nDirTot,nDir,nDirZee,nH,iPrint,LUZee,mem,compute_Mdir_vector,zeeman_energy,hinput, &
                     m_paranoid,smagn,doplot,TempMagn,eso,dirX,dirY,dirZ,dir_weight,hexp,magn_exp,zJ,hmin,hmax,EM,thrs,mm,ms,dbg)
  if (DBG) write(u6,*) 'SINGLE_ANISO2::  Exit magnetization'

else
  write(u6,*)
  write(u6,'(5X,A)') 'ON USER REQUEST, THE MAGNETIZATION VECTOR AND THE MEAN MAGNETIZATION WAS NOT CALCULATED'
end if

!-----------------------------------------------------------------------
! Deallocate memory for all arrays:
!-----------------------------------------------------------------------
call mma_deallocate(esfs)
call mma_deallocate(esfs_au)
call mma_deallocate(ANGMOM)
call mma_deallocate(EDMOM)
call mma_deallocate(AMFI)
call mma_deallocate(multiplicity)

call mma_deallocate(eso)
call mma_deallocate(eso_au)
call mma_deallocate(U)
call mma_deallocate(HSO)
call mma_deallocate(MM)
call mma_deallocate(MS)
call mma_deallocate(ML)
call mma_deallocate(DM)

call mma_deallocate(Hexp)
call mma_deallocate(magn_exp)
call mma_deallocate(TempMagn)

call mma_deallocate(ndim)
call mma_deallocate(gtens)
call mma_deallocate(maxes)
call mma_deallocate(T)
call mma_deallocate(XTexp)
call mma_deallocate(XT_no_field)

call mma_deallocate(LuZee)
call mma_deallocate(dir_weight)

call mma_deallocate(dirX)
call mma_deallocate(dirY)
call mma_deallocate(dirZ)

call mma_deallocate(Texp)
call mma_deallocate(chit_exp)

return

end subroutine SINGLE_ANISO2
