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

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, cZero
use Definitions, only: wp, iwp, u6, CtoB, ItoB, RtoB

implicit none
integer(kind=iwp), intent(inout) :: nH, nT, nTempMagn, nDir, nDirZee, nMult
integer(kind=iwp), intent(in) :: nss, nstate
character(len=180), intent(inout) :: input_file_name
logical(kind=iwp), intent(out) :: ifrestart
integer(kind=iwp), intent(out) :: iReturn
logical(kind=iwp), intent(in) :: GRAD
integer(kind=iwp) :: AngPoints, axisoption, d, encut_definition, i, i1, i2, Ifunct, imanifold, imltpl, input_to_read, iprint, &
                     lDIM, ldimcf, mem, MG, nBlock, ncut, ndimcf, nDirTot, ngrid, nK, nlanth, nm, nss2, nstate2, nsymm
real(kind=wp) :: coord(3), cryst(6), em, encut_rate, gtens(3), H_torq, HMAX, HMIN, maxes(3,3), T_torq, thrs, tmax, tmin, xfield, &
                 zJ, zmagn(3,3)
logical(kind=iwp) :: compute_barrier, compute_cf, compute_g_tensors, compute_magnetization, compute_Mdir_vector, compute_torque, &
                     Do_structure_abc, DoPlot, hinput, m_paranoid, poly_file, smagn, tinput, zeeman_energy
integer(kind=iwp), allocatable :: LuZee(:), multiplicity(:), ndim(:)
real(kind=wp), allocatable :: amfi(:,:,:), amfi_tmp(:,:,:), angmom(:,:,:), angmom_tmp(:,:,:), chit_exp(:), dir_weight(:,:), &
                              dirX(:), dirY(:), dirZ(:), eDmom(:,:,:), esfs(:), esfs_au(:), eso(:), eso_au(:), hexp(:), &
                              magn_exp(:,:), t(:), TempMagn(:), Texp(:), XT_no_field(:), XTexp(:)
complex(kind=wp), allocatable :: DM(:,:,:), HSO(:,:), ML(:,:,:), MM(:,:,:), MM_TMP(:,:,:), MS(:,:,:), MS_TMP(:,:,:), U(:,:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: j, l
#  define _DBG_ .true.
#else
#  define _DBG_ .false.
#endif
logical(kind=iwp), parameter :: dbg = _DBG_

!-----------------------------------------------------------------------
! Allocate memory for all arrays:
!-----------------------------------------------------------------------
#ifdef _DEBUGPRINT_
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
#endif

mem = 0

! spin free energies
call mma_allocate(esfs,nstate,'esfs')
call mma_allocate(esfs_au,nstate,'esfs_au')
mem = mem+size(esfs)*RtoB
mem = mem+size(esfs_au)*RtoB
! angular momentum
call mma_allocate(ANGMOM,3,nstate,nstate,'angmom')
mem = mem+size(ANGMOM)*RtoB
! electric dipole moment
call mma_allocate(EDMOM,3,nstate,nstate,'edmom')
mem = mem+size(EDMOM)*RtoB
! amfi integrals
call mma_allocate(AMFI,3,nstate,nstate,'amfi')
mem = mem+size(AMFI)*RtoB
! multiplicity of each state
call mma_allocate(multiplicity,nstate,'multiplicity')
mem = mem+size(multiplicity)*ItoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 1 =',mem
#endif
! spin orbit energies
call mma_allocate(eso,nss,'eso')
call mma_allocate(eso_au,nss,'eso_au')
mem = mem+size(eso)*RtoB
mem = mem+size(eso_au)*RtoB
! spin orbit eigenstates
call mma_allocate(U,nss,nss,'U')
mem = mem+size(U)*CtoB
! spin orbit hamiltonian
call mma_allocate(HSO,nss,nss,'HSO')
mem = mem+size(HSO)*CtoB
! magnetic moment
call mma_allocate(MM,3,nss,nss,'MM')
mem = mem+size(MM)*CtoB
! spin moment
call mma_allocate(MS,3,nss,nss,'MS')
mem = mem+size(MS)*CtoB
! orbital mooment
call mma_allocate(ML,3,nss,nss,'ML')
mem = mem+size(ML)*CtoB
! electric dipole moment
call mma_allocate(DM,3,nss,nss,'DM')
mem = mem+size(DM)*CtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 2 =',mem
#endif
! experimental magnetic field points
call mma_allocate(Hexp,nH,'Hexp')
mem = mem+size(Hexp)*RtoB
! experiemental magnetization
call mma_allocate(magn_exp,nH,nTempMagn,'magn_exp')
mem = mem+size(magn_exp)*RtoB
! temperature points for magnetization
call mma_allocate(TempMagn,nTempMagn,'TempMagn')
mem = mem+size(TempMagn)*RtoB

! dimensions of pseudospins
call mma_allocate(ndim,nMult,'ndim')
mem = mem+size(ndim)*ItoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 4 =',mem
#endif
call mma_allocate(T,nT+nTempMagn,'Temperature')
mem = mem+size(T)*RtoB
call mma_allocate(XTexp,nT+nTempMagn,'XTexp')
mem = mem+size(XTexp)*RtoB
call mma_allocate(XT_no_field,nT+nTempMagn,'XT_no_field')
mem = mem+size(XT_no_field)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 5 =',mem
#endif

! unit numbers for the files with Zeeman energies
call mma_allocate(LuZee,nDirZee,'LUZee')
mem = mem+size(LuZee)*ItoB
! directions for applied field for Zeeman states
call mma_allocate(dir_weight,nDirZee,3,'dir_weight')
mem = mem+size(dir_weight)*RtoB

! magnetization vectors
call mma_allocate(dirX,nDir,'dirX')
call mma_allocate(dirY,nDir,'dirY')
call mma_allocate(dirZ,nDir,'dirZ')
mem = mem+size(dirX)*RtoB
mem = mem+size(dirY)*RtoB
mem = mem+size(dirZ)*RtoB
! T experimental given by user in the input
call mma_allocate(Texp,nT,'Texp')
mem = mem+size(Texp)*RtoB
! XT experimental given by user in the input
call mma_allocate(chit_exp,nT,'chit_exp')
mem = mem+size(chit_exp)*RtoB
! allocated memory counter
#ifdef _DEBUGPRINT_
write(u6,'(A,I16)') 'mem 8 =',mem
#endif

write(u6,'(A,I16,A)') 'The code allocated initially:',mem,' bytes of memory for this run.'
call xFlush(u6)
!-----------------------------------------------------------------------
IReturn = 0
NM = 0
EM = Zero
! read the input
#ifdef _DEBUGPRINT_
write(u6,*) 'SINGLE_ANISO2::  Enter readin_single'
#endif
call readin_single(iprint,nmult,ndim,ndimcf,ldimcf,nlanth,axisoption,poly_file,Ifrestart,input_to_read,nk,mg,zmagn, &
                   Do_structure_abc,cryst,coord,encut_definition,compute_g_tensors,compute_CF,nDirTot,nss,nstate, &
                   compute_magnetization,compute_torque,smagn,tinput,hinput,compute_Mdir_vector,zeeman_energy,LUZee,doplot, &
                   encut_rate,ncut,nTempMagn,TempMagn,m_paranoid,compute_barrier,nBlock,AngPoints,input_file_name,nT,nH,texp, &
                   chit_exp,zJ,hexp,magn_exp,hmin,hmax,nDir,nDirZee,dirX,dirY,dirZ,dir_weight,xfield,tmin,tmax,thrs,H_torq,T_torq, &
                   nsymm,ngrid)
#ifdef _DEBUGPRINT_
write(u6,*) 'SINGLE_ANISO2::  Exit readin_single'
#endif

if (ifrestart) then
  ! if restart, fetch "big data" from the input file
  if (input_to_read == 1) then
    call read_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)
    esfs_au(:) = Zero
    eso_au(:) = Zero

  else if (input_to_read == 2) then
    ! get the information from formatted aniso.input file:
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Enter read_formatted_aniso'
    write(u6,*) 'SA:',input_file_name
#   endif
    nss2 = nss
    nstate2 = nstate
    call read_formatted_aniso(input_file_name,nss2,nstate2,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)
    esfs_au(:) = Zero
    eso_au(:) = Zero
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Exit  read_formatted_aniso'
#   endif

  else if (input_to_read == 3) then
    ! get the information from RASSI-HDF5 file:
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SA:',input_file_name
#   endif
#   ifdef _HDF5_
    call read_hdf5_all(input_file_name,nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO)
#   else
    call WarningMessage(2,'File '//trim(input_file_name)//' cannot be opened. Molcas was compiled without HDF5 option.')
    call Quit_OnUserError()
#   endif
    esfs_au(:) = Zero
    eso_au(:) = Zero

  else if (input_to_read == 4) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SA:',input_file_name
#   endif
    ! get the information from formatted aniso.input file:
    nss2 = nss
    nstate2 = nstate
    call read_formatted_aniso_old(input_file_name,nss2,nstate2,multiplicity,eso,MM,MS,ML)
    esfs(:) = Zero
    esfs_au(:) = Zero
    eso_au(:) = Zero
    ANGMOM(:,:,:) = Zero
    EDMOM(:,:,:) = Zero
    AMFI(:,:,:) = Zero
    U(:,:) = cZero
    HSO(:,:) = cZero
    DM(:,:,:) = cZero

  else if (input_to_read == 6) then !  using DATA keyword

#   ifdef _DEBUGPRINT_
    write(u6,*) 'SA:',input_file_name
#   endif
    ! get the information from formatted new aniso-data file:
    nss2 = nss
    nstate2 = nstate
    call read_formatted_new_aniso(input_file_name,nss2,nstate2,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO,eso_au, &
                                  esfs_au)
  end if ! input_to_read

else
  ! ifrestart == .false., i.e. usual S-A calculation
  call fetch_data_RunFile_all(nss,nstate,multiplicity,eso,esfs,U,MM,MS,ML,DM,ANGMOM,EDMOM,AMFI,HSO,eso_au,esfs_au)
  write(u6,'(A)') 'AFTER fetch_data_RunFile_all'
  call xFlush(u6)
# ifdef _DEBUGPRINT_
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
# endif

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

!----- input processing finished -----
! save some important data, regardless of the following execution
! -- binary $Project.aniso
call write_binary_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)
! ASCII -- anisoinput:
call write_formatted_aniso(nss,nstate,multiplicity,eso,esfs,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)

! ASCII -- new_aniso file format:
if (.not. ifrestart) call write_new_formatted_aniso(nss,nstate,multiplicity,eso_au,esfs_au,U,MM,MS,DM,ANGMOM,EDMOM,AMFI,HSO)

!----- compute various properties ----------|
! calculation of magnetic Hamiltonians:
IFUNCT = 0
if (compute_g_tensors .and. (nMult > 0)) then
  do IMLTPL=1,NMULT
    if (ndim(imltpl) == 1) then
      write(u6,'(5X,A,I2,A)') 'THE DIMENSION OF THE ',IMLTPL,' MULTIPLET IS 1.'
      write(u6,'(5X,A)') 'THERE IS NO G TENSOR FOR EFFECTIVE S = 0.'
    else
      i1 = 1+Ifunct
      i2 = ndim(imltpl)+Ifunct

      if (i2 <= nss) then
#       ifdef _DEBUGPRINT_
        write(u6,*) 'SINGLE_ANISO2::  Enter g_high',IMLTPL
#       endif
        call mma_allocate(MS_TMP,3,ndim(imltpl),ndim(imltpl),label='MS_TMP')
        call mma_allocate(MM_TMP,3,ndim(imltpl),ndim(imltpl),label='MM_TMP')
        MS_TMP(:,:,:) = MS(:,i1:i2,i1:i2)
        MM_TMP(:,:,:) = MM(:,i1:i2,i1:i2)
        call g_high(eso(i1:i2),GRAD,MS_TMP,MM_TMP,imltpl,ndim(imltpl),Do_structure_abc,cryst,coord,gtens,maxes,iprint)
        call mma_deallocate(MS_TMP)
        call mma_deallocate(MM_TMP)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'SINGLE_ANISO2::  Exit g_high',IMLTPL
#       endif
        call Add_Info('GTENS_MAIN',gtens,3,4)
      end if
    end if

    IFUNCT = IFUNCT+NDIM(IMLTPL)
  end do
end if

!----------------------------------------------------------------------|
! >> AB INITIO CRYSTAL-FIELD <<

if (compute_CF .and. (nDIMCF > 0)) then
  if ((axisoption == 1) .and. (nMult > 0) .and. (nDIM(1) > 1)) then
    d = NDIM(1)
  else
    d = nDIMCF
  end if

  ! compute the CF of the ground |J,MJ> multiplet:
  if ((nDIMCF > 1) .and. (nDIMCF <= nss)) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Enter CF',nDIMCF
#   endif
    call mma_allocate(MS_TMP,3,nDIMCF,nDIMCF,label='MS_TMP')
    call mma_allocate(MM_TMP,3,nDIMCF,nDIMCF,label='MM_TMP')
    MS_TMP(:,:,:) = MS(:,1:nDIMCF,1:nDIMCF)
    MM_TMP(:,:,:) = MM(:,1:nDIMCF,1:nDIMCF)
    call CRYSTALFIELD(ESO(1:nDIMCF),MM_TMP,MS_TMP,nDIMcf,d,nlanth,zmagn,axisoption,GRAD,iPrint)
    call mma_deallocate(MS_TMP)
    call mma_deallocate(MM_TMP)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Exit CF',nDIMCF
#   endif
  end if
end if

#ifdef _DEBUGPRINT_
write(u6,*) 'SINGLE_ANISO2:: nlanth=',nlanth
write(u6,*) 'SINGLE_ANISO2:: lDIMCF=',lDIMCF
write(u6,*) 'SINGLE_ANISO2:: nstate=',nstate
#endif

if (compute_CF .and. (lDIMCF > 0) .and. (lDIMCF <= nstate)) then
  if ((axisoption == 1) .and. (nMult > 0) .and. (nDIM(1) > 1)) then
    lDIM = NDIM(1)
  else
    lDIM = lDIMCF
  end if
  ! compute the CF of the ground |L,ML> term:
  if (.not. (ifrestart .and. (input_to_read == 4))) then
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Enter t-CF',lDIMCF
#   endif
    call mma_allocate(amfi_tmp,3,lDIMCF,lDIMCF,label='amfi_tmp')
    call mma_allocate(angmom_tmp,3,lDIMCF,lDIMCF,label='angmom_tmp')
    amfi_tmp(:,:,:) = AMFI(:,1:lDIMCF,1:lDIMCF)
    angmom_tmp(:,:,:) = angmom(:,1:lDIMCF,1:lDIMCF)
    call termCF(angmom_tmp,AMFI_tmp,esfs(1:lDIMCF),lDIMCF,lDIM,zmagn,axisoption,nlanth,iPrint)
    call mma_deallocate(amfi_tmp)
    call mma_deallocate(angmom_tmp)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Exit t-CF',lDIMCF
#   endif
  end if !ifrestart
end if !compute_CF

!----------------------------------------------------------------------|
! >> AB INITIO BLOCKING BARRIER FOR SMMs <<
if (compute_barrier) then
  if (nBlock /= 0) then
    imanifold = 1

#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Enter barrier',nBlock
#   endif
    call mma_allocate(MM_TMP,3,nBlock,nBlock,label='MM_TMP')
    MM_TMP(:,:,:) = MM(:,1:nBlock,1:nBlock)
    call BARRIER(nBlock,MM_TMP,eso(1:nBlock),imanifold,nMult,nDim,doplot,iPrint)
    call mma_deallocate(MM_TMP)
#   ifdef _DEBUGPRINT_
    write(u6,*) 'SINGLE_ANISO2::  Exit barrier',nBlock
#   endif

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
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Enter set_nm'
# endif
  call set_nm(nss,ncut,encut_definition,nk,mg,nTempMagn,hmax,ESO,encut_rate,TempMagn,nM,EM,dbg)
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  exit set_nm'
# endif

end if

!----------------------------------------------------------------------|
! >> MAGNETIC SUSCEPTIBILITY <<
#ifdef _DEBUGPRINT_
write(u6,*) 'SINGLE_ANISO2::  Enter SUSCEPTIBILITY'
#endif
call SUSCEPTIBILITY(NSS,ESO,MS,MM,nT,nTempMagn,T,tmin,tmax,XTexp,zJ,tinput,XT_no_field,doplot,iPrint,mem)
#ifdef _DEBUGPRINT_
write(u6,*) 'SINGLE_ANISO2::  Exit SUSCEPTIBILITY'
#endif

if (Xfield /= Zero) then
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Enter XT_dMoverdH_single'
# endif
  ! nM = nss
  call XT_dMoverdH_single(nss,nTempMagn,nT,nss,Tmin,Tmax,XTexp,ESO,T,zJ,Xfield,MM,MS,XT_no_field,tinput,smagn,mem,DoPlot)
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Exit XT_dMoverdH_single'
# endif
end if
!----------------------------------------------------------------------|
! >> TORQUE <<
if (compute_torque) then

# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Enter TORQUE'
# endif
  ! nM = nss
  call torque(Nss,Nss,AngPoints,EM,eso,mm,ms,zJ,thrs,mem,m_paranoid,smagn,H_torq,T_torq,zmagn,dbg)
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Exit TORQUE'
# endif

end if
!----------------------------------------------------------------------|
! >> MAGNETIZATION <<

if (compute_magnetization) then

# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Enter magnetization'
# endif
  call magnetization(nss,nM,nTempMagn,nDirTot,nDir,nDirZee,nH,iPrint,LUZee,mem,nsymm,ngrid,compute_Mdir_vector,zeeman_energy, &
                     hinput,m_paranoid,smagn,doplot,TempMagn,eso,dirX,dirY,dirZ,dir_weight,hexp,magn_exp,zJ,hmin,hmax,EM,thrs,mm, &
                     ms,dbg)
# ifdef _DEBUGPRINT_
  write(u6,*) 'SINGLE_ANISO2::  Exit magnetization'
# endif

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
