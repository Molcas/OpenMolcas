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
subroutine exctdm(SIJ,TRAD,TDMAB,iRC,CMO1,CMO2,TDMZZ,TRASD,TSDMAB,TSDMZZ,istate,jstate)

use Index_Functions, only: iTri, nTri_Elem
use Basis_Info, only: nBas
use Data_structures, only: Allocate_DT, Deallocate_DT, DSBA_Type
use frenkel_global_vars, only: doexch, eNucB, iTyp, labb, VNucB
use Cntrl, only: MLTPLT, NSTATE
use Symmetry_Info, only: nIrrep
use rassi_data, only: CHFRACMEM, NASHT, NBASF, NCMO, NTDMAB, NTDMZZ
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Half
use Definitions, only: wp, iwp
#ifdef _DEBUGPRINT_
use Definitions, only: u6
#endif

implicit none
real(kind=wp), intent(in) :: SIJ, TRAD(NASHT,NASHT), CMO1(NCMO), CMO2(NCMO), TRASD(NASHT,NASHT)
real(kind=wp), intent(out) :: TDMAB(NTDMAB), TDMZZ(NTDMZZ), TSDMAB(NTDMAB), TSDMZZ(NTDMZZ)
integer(kind=iwp), intent(out) :: iRC
integer(kind=iwp), intent(in) :: istate, jstate
integer(kind=iwp) :: I, INTEG, IPNB, J, LuT, m(1), n(1), nbas_A(1), nbas_B(1), nbas_tot(1), NNLTD, run
type(DSBA_Type) :: DLT(1), Salpha(1), Sbeta(1), SDLT(1)
character(len=13) :: filnam
real(kind=wp), allocatable :: STDMZZ_mtx(:,:), STDMZZ_new(:), TDMZZ_mtx(:,:), TDMZZ_new(:)
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: a, K
#endif
integer(kind=iwp), external :: isFreeUnit
real(kind=wp), external :: DDot_

!write(u6,*) 'MLTPL', MLTPLT(1)
! Write number of states in files, since they can differ on monomers
write(filnam,'(A,I1)') 'states_',iTyp
LuT = isFreeUnit(10)
call molcas_open(LuT,filnam)
write(LuT,*) NSTATE
close(LuT)

! get number of basis functions for different monomers

if (LABB) then
  nbas_B = NBASF(1)
# ifdef _DEBUGPRINT_
  write(u6,*) 'basis functions of mon B',nbas_B
# endif
  call NameRun('AUXRFIL1')
  call get_iArray('nBas',nBas,nIrrep)
  nbas_tot = nBas(0)
  call NameRun('#Pop')    ! switch back to old RUNFILE
# ifdef _DEBUGPRINT_
  write(u6,*) 'total num of basis functions',nbas_tot
# endif
  nbas_A = nbas_tot-nbas_B
# ifdef _DEBUGPRINT_
  write(u6,*) 'basis functions of mon A',nbas_A
# endif
end if

if (.not. LABB) then
  nbas_A = NBASF(1)
# ifdef _DEBUGPRINT_
  write(u6,*) 'basis functions of mon A',nbas_A
# endif
  call NameRun('AUXRFIL1')
  call get_iArray('nBas',nBas,nIrrep)
  nbas_tot = nBas(0)
  call NameRun('#Pop')    ! switch back to old RUNFILE
# ifdef _DEBUGPRINT_
  write(u6,*) 'total num of basis functions',nbas_tot
# endif
  nbas_B = nbas_tot-nbas_A
# ifdef _DEBUGPRINT_
  write(u6,*) 'basis functions of mon B',nbas_B
# endif
end if
! get the dimensions of the packed lower triangular TDM
NNLTD = nTri_Elem(nbas_tot(1))
m = nbas_tot
n = nbas_tot

! regular-TDM
call MKTDAB(SIJ,TRAD,TDMAB,iRC)
! transform to AO basis
call MKTDZZ(CMO1,CMO2,TDMAB,TDMZZ,iRC)
call Allocate_DT(DLT(1),n,m,nIrrep,aCase='TRI')

if (doexch) then
  ! for the exchange part we need the spin density
  ! spin-TDM
  call MKTDAB(Zero,TRASD,TSDMAB,iRC)
  ! transform to AO basis
  call MKTDZZ(CMO1,CMO2,TSDMAB,TSDMZZ,iRC)
  ! alpha part of the spin density
  if (MLTPLT(1) == 1) then
    call Allocate_DT(Salpha(1),nbas_tot,nbas_tot,nIrrep)
  else
    call Allocate_DT(Salpha(1),nbas_tot,nbas_tot,nIrrep)
    call Allocate_DT(Sbeta(1),nbas_tot,nbas_tot,nIrrep)
    ! allocate matrix for the spin density difference
    call Allocate_DT(SDLT(1),nbas_tot,nbas_tot,nIrrep)
    SDLT(1)%A00(:) = Zero
  end if
end if

! case where MON B of the exciton calculation was calculated first
if (labb) then
  call mma_allocate(TDMZZ_mtx,nbas_tot(1),nbas_tot(1))
  call mma_allocate(TDMZZ_new,nbas_tot(1)*nbas_tot(1))
  TDMZZ_mtx(:,:) = Zero
  TDMZZ_new(:) = Zero
  ! fill a matrix in combined basis with the TDM of MON B
  ! here, one shifts the entries of row and column of the matrix
  ! by the number of basis functions of MON A
  ! the fast running index is the column, because TDMZZ is filled in
  ! column-major order
  run = 0
  do i=1,nbas_B(1)
    do j=1,nbas_B(1)
      run = run+1
      TDMZZ_mtx(j+nbas_A(1),i+nbas_A(1)) = TDMZZ(run)
    end do
  end do
  ! fill the TDM in common basis in an array
  ! in column-major order
  run = 0
  do i=1,nbas_tot(1)
    do j=1,nbas_tot(1)
      run = run+1
      TDMZZ_new(run) = TDMZZ_mtx(j,i)
    end do
  end do
  ! fill DLT (=TDM in common basis in packed, lower trangular storage)
  do i=1,nbas_tot(1)
    do j=1,i-1
      DLT(1)%A00(iTri(i,j)) = TDMZZ_new(nbas_tot(1)*(j-1)+i)+TDMZZ_new(nbas_tot(1)*(i-1)+j)
    end do
    DLT(1)%A00(nTri_Elem(j)) = TDMZZ_new(nbas_tot(1)*(j-1)+j)
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'TDMZZvector:'
  k = nbas_B(1)*nbas_B(1)
  do i=1,k
    write(u6,*) TDMZZ(i)
  end do

  write(u6,*) 'TDMZZ:'
  k = int(sqrt(real(nTDMZZ,kind=wp)),kind=iwp)
  do i=1,k
    write(u6,'(100ES18.8)') (TDMZZ(i+(j-1)*k),j=1,k)
  end do

  write(u6,*) 'TDMZZ_new vector:'
  k = nbas_tot(1)*nbas_tot(1)
  do i=1,k
    write(u6,*) TDMZZ_new(i)
  end do

  write(u6,*) 'TDMZZ_new:'
  k = nbas_tot(1)
  do i=1,k
    write(u6,'(1000ES18.8)') (TDMZZ_new(i+(j-1)*k),j=1,k)
  end do
  write(u6,*) 'TDMZZ_mtx:'
  do i=1,k
    write(u6,'(1000ES18.8)') (TDMZZ_mtx(i,j),j=1,k)
  end do

  write(u6,*) 'DLT:'
  write(u6,*) DLT(1)%A00(:)
  write(u6,*) 'DLT as mtx:'
  k = 0
  do i=1,nbas_tot(1)
    write(u6,'(1000ES18.8)') (DLT(1)%SB(1)%A1(j+k),j=1,i)
    k = k+i
  end do
# endif
  ! spin part if KCOUL is requested in excitonic section
  if (DoExch) then
    call mma_allocate(STDMZZ_mtx,nbas_tot(1),nbas_tot(1))
    call mma_allocate(STDMZZ_new,nbas_tot(1)*nbas_tot(1))
    STDMZZ_mtx(:,:) = Zero
    STDMZZ_new(:) = Zero
    run = 0
    do i=1,nbas_B(1)
      do j=1,nbas_B(1)
        run = run+1
        STDMZZ_mtx(j+nbas_A(1),i+nbas_A(1)) = TSDMZZ(run)
      end do
    end do
    run = 0
    do i=1,nbas_tot(1)
      do j=1,nbas_tot(1)
        run = run+1
        STDMZZ_new(run) = STDMZZ_mtx(j,i)
      end do
    end do
    ! fill SDLT (spin density difference)
    if (MLTPLT(1) > 1) then
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          SDLT(1)%SB(1)%A1(run) = STDMZZ_mtx(i,j)
        end do
      end do
    end if
    ! fill Salpha from matrix, due to the row major data type
    run = 0
    do i=1,nbas_tot(1)
      do j=1,nbas_tot(1)
        run = run+1
        Salpha(1)%SB(1)%A1(run) = TDMZZ_mtx(i,j)
      end do
    end do

    if (MLTPLT(1) > 1) then
      ! Salpha(1) = TM = GAA+GBB
      ! Salpha(1) becomes TM + SDLT = GAA+GBB+GAA-GBB=2GAA
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          Salpha(1)%SB(1)%A1(run) = Salpha(1)%SB(1)%A1(run)+SDLT(1)%SB(1)%A1(run)
        end do
      end do
      ! Salpha(1) becomes GAA
      integ = int(nbas_tot(1)*nbas_tot(1),kind=iwp)
      Salpha(1)%A00(1:integ) = Half*Salpha(1)%A00(1:integ)
      ! Sbeta(1) becomes GBB
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          Sbeta(1)%SB(1)%A1(run) = TDMZZ_new(run)-Salpha(1)%SB(1)%A1(run)
        end do
      end do
    end if

    call mma_deallocate(STDMZZ_new)
    call mma_deallocate(STDMZZ_mtx)
    ! end Exch statement
  end if

  call mma_deallocate(TDMZZ_new)
  call mma_deallocate(TDMZZ_mtx)
  ! end LABB statement
end if

! case where MON A was calculated first
if (.not. LABB) then
  call mma_allocate(TDMZZ_mtx,nbas_tot(1),nbas_tot(1))
  call mma_allocate(TDMZZ_new,nbas_tot(1)*nbas_tot(1))
  TDMZZ_mtx(:,:) = Zero
  TDMZZ_new(:) = Zero
  ! fill matrix with the TDM in common basis (column-major order)
  run = 0
  do i=1,nbas_A(1)
    do j=1,nbas_A(1)
      run = run+1
      TDMZZ_mtx(j,i) = TDMZZ(run)
    end do
  end do
  ! fill TDM back as array
  run = 0
  do i=1,nbas_tot(1)
    do j=1,nbas_tot(1)
      run = run+1
      TDMZZ_new(run) = TDMZZ_mtx(j,i)
    end do
  end do

  ! fill DLT in packed lower triangular storage
  do i=1,nbas_tot(1)
    do j=1,i-1
      DLT(1)%A00(iTri(i,j)) = TDMZZ_new(nbas_tot(1)*(j-1)+i)+TDMZZ_new(nbas_tot(1)*(i-1)+j)
    end do
    DLT(1)%A00(nTri_Elem(j)) = TDMZZ_new(nbas_tot(1)*(j-1)+j)
  end do

# ifdef _DEBUGPRINT_
  write(u6,*) 'TDMZZ vector:'
  a = nbas_A(1)*nbas_A(1)
  do i=1,a
    write(u6,*) TDMZZ(i)
  end do
  write(u6,*) 'TDMZZ:'
  a = int(sqrt(real(nTDMZZ,kind=wp)),kind=iwp)
  do i=1,a
    write(u6,'(100ES18.8)') (TDMZZ(i+(j-1)*a),j=1,a)
  end do

  write(u6,*) 'TDMZZ_new vector:'
  a = nbas_tot(1)*nbas_tot(1)
  do i=1,a
    write(u6,*) TDMZZ_new(i)
  end do

  write(u6,*) 'TDMZZ_new:'
  a = nbas_tot(1)
  do i=1,a
    write(u6,'(1000ES18.8)') (TDMZZ_new(i+(j-1)*a),j=1,a)
  end do

  write(u6,*) 'TDMZZ_mtx:'
  do i=1,a
    write(u6,'(1000ES18.8)') (TDMZZ_mtx(i,j),j=1,a)
  end do

  write(u6,*) 'DLT:'
  write(u6,*) DLT(1)%A00(:)
# endif

  ! spin part
  if (DoExch) then
    call mma_allocate(STDMZZ_mtx,nbas_tot(1),nbas_tot(1))
    call mma_allocate(STDMZZ_new,nbas_tot(1)*nbas_tot(1))
    STDMZZ_mtx(:,:) = Zero
    STDMZZ_new(:) = Zero
    run = 0
    do i=1,nbas_A(1)
      do j=1,nbas_A(1)
        run = run+1
        STDMZZ_mtx(i,j) = TSDMZZ(run)
      end do
    end do
    ! fill back as array
    run = 0
    do i=1,nbas_tot(1)
      do j=1,nbas_tot(1)
        run = run+1
        STDMZZ_new(run) = STDMZZ_mtx(i,j)
      end do
    end do
    ! fill SDLT (spin density difference)
    if (MLTPLT(1) > 1) then
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          SDLT(1)%SB(1)%A1(run) = STDMZZ_mtx(i,j)
        end do
      end do
    end if

    ! fill Salpha from matrix, since here the row major data type
    run = 0
    do i=1,nbas_tot(1)
      do j=1,nbas_tot(1)
        run = run+1
        Salpha(1)%SB(1)%A1(run) = TDMZZ_mtx(i,j)
      end do
    end do

    if (MLTPLT(1) > 1) then
      ! Salpha(1) = TM = GAA+GBB
      ! Salpha(1) becomes TM + SDLT = GAA+GBB+GAA-GBB=2GAA
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          Salpha(1)%SB(1)%A1(run) = Salpha(1)%SB(1)%A1(run)+SDLT(1)%SB(1)%A1(run)
        end do
      end do
      ! Salpha(1) becomes GAA
      integ = int(nbas_tot(1)*nbas_tot(1),kind=iwp)
      Salpha(1)%A00(1:integ) = Half*Salpha(1)%A00(1:integ)
      ! Sbeta(1) becomes GBB
      run = 0
      do i=1,nbas_tot(1)
        do j=1,nbas_tot(1)
          run = run+1
          Sbeta(1)%SB(1)%A1(run) = TDMZZ_new(run)-Salpha(1)%SB(1)%A1(run)
        end do
      end do
    end if
    !Salpha(1)%A00(1:integ) = Salpha(1)%A00(1:integ)+SDLT(1)%A00(1:integ)
    call mma_deallocate(STDMZZ_new)
    call mma_deallocate(STDMZZ_mtx)
    ! end EXCH
  end if

  call mma_deallocate(TDMZZ_new)
  call mma_deallocate(TDMZZ_mtx)
  ! end filling MON A
end if
ipNB = nTri_Elem(ISTATE-1)+JSTATE

! calculate rho-nuc interaction (TDM with core potential)
eNucB(ipNB) = ddot_(NNLTD,DLT(1)%A0,1,VNucB(1),1)

! change the RunFile to the one
! from bsse calculation to get the cholesky vectors
call NameRun('AUXRFIL1')
call CHO_X_INIT(irc,ChFracMem)
call CHO_TRDENS(irc,DLT,Salpha,istate,jstate,iTyp,DoExch,LABB)
call NameRun('#Pop')    ! switch back to old RUNFILE
call deallocate_DT(DLT(1))
if (DOEXCH) then
  if (MLTPLT(1) == 1) then
    call deallocate_DT(Salpha(1))
  else
    call deallocate_DT(Salpha(1))
    call deallocate_DT(Sbeta(1))
    call deallocate_DT(SDLT(1))
  end if
end if

end subroutine exctdm
