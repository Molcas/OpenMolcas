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
! Copyright (C) 2025, Yoshio Nishimoto                                 *
!***********************************************************************
!
! Dynamically weighted solvation (reaction field)
!
! Could be used for dynamically weighted MCSCF (DW-MCSCF), and this was an initial plan...
! The density matrix used for generating the reaction field is constructed by a weighted averaged
! of state-specific density matrices. The weight is determined with the electronic energy and
! distribution function. This is similar to the well-known DW-MCSCF.
!
! At present, only the Boltzmann distribution function can be used.
! Note also that the electronic gradient and Hessian of the weight function is ignored.
!
! (hidden) options in &GATEWAY RF-Input
!
! DWSOlvation
!   DWZeta in the code. This regulates the distribution function.
! FIXRfroot
!   Generate fixed weight reaction fields; it is not a dynamically weighted, of course.
!   This is the only option written in the manual among the three.
! DWTYpe
!   Specifies the distribution function (see DWSol_func). Only DWTYPE = 1 is working.
!   I think we just need to update DWSol_func and DWSol_Der (zeroth- and first-order derivatives
!   wrt the energy)

module DWSol

use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp, u6

implicit none

type type_DW
  logical(kind=iwp) :: do_DW = .false.
  integer(kind=iwp) :: DWType = -99
  integer(kind=iwp) :: DWRoot = 0
  real(kind=wp) :: DWZeta = Zero
end type type_DW

!! DWSolv: dynamically weighted solvation
!! DWSCF : dynamically weighted MCSCF
type(type_DW), target :: DWSolv, DWSCF

integer(kind=iwp) :: IPCMROOT = 0, nRoots = 0

real(kind=wp), allocatable :: W_SOLV(:)

contains

!-----------------------------------------------------------------------

! dummy subroutine for initialization things
subroutine DWSol_dummy()

  return

end subroutine DWSol_dummy

!-----------------------------------------------------------------------

! should be called somewhere in RASSCF or MCLR after nRoots and
! IPCMROOT are determined
subroutine DWSol_init(IPCMROOT_,nRoots_,NonEq_)

  use rctfld_module, only: iSLPar, rSlPar
  use input_mclr, only: Weight_MCLR => Weight
  use rasscf_global, only: Weight_RASSCF => Weight
  use UnixInfo, only: ProgName

  integer(kind=iwp), intent(in) :: IPCMROOT_, nRoots_
  logical(kind=iwp), intent(in) :: NonEq_
  integer(kind=iwp) :: i, istate, jstate, mData
  real(kind=wp) :: wgt
  logical(kind=iwp) :: Found

  IPCMROOT = IPCMROOT_
  nRoots = nRoots_

  DWSolv%DWType = ISlPar(17)
  DWSolv%DWRoot = IPCMROOT
  DWSolv%DWZeta = RSlPar(53)

  if (NonEq_) DWSolv%DWZeta = Zero

  call mma_allocate(W_SOLV,nRoots,Label='W_SOLV')
  W_SOLV(:) = Zero

  call qpg_dArray('SolventWeight',Found,mData)
  if (Found .and. (mData /= 0)) then
    if (mData /= nRoots) call warningMessage(1,'The number of roots specified by FIXRFROOT/DWROOT and CIROOT may be inconsistent.')
    call Get_dArray('SolventWeight',W_SOLV,min(mData,nRoots))
    DWSolv%DWZeta = -12345.0_wp
  else
    if (DWSolv%DWZeta == Zero) then
      if (IPCMROOT > 0) then
        W_SOLV(IPCMROOT) = One
      else
        if (ProgName(1:6) == 'rasscf') W_SOLV(1:nRoots) = Weight_RASSCF(1:nRoots)
        if (ProgName(1:4) == 'mclr') W_SOLV(1:nRoots) = Weight_MCLR(1:nRoots)
      end if
    else if (DWSolv%DWZeta < Zero) then
      !! fixed averaging
      call DWSol_fixed(istate,jstate)
      do i=1,nRoots
        if ((i == istate) .or. (i == jstate)) then
          wgt = Half
        else
          wgt = Zero
        end if
        W_SOLV(i) = wgt
      end do
    else if (DWSolv%DWZeta > Zero) then
      DWSolv%do_DW = .true.
    end if
  end if

end subroutine DWSol_init

!-----------------------------------------------------------------------

! should be called somewhere in RASSCF or MCLR after nRoots and
! IPCMROOT are determined
subroutine DWSCF_init(mode,nRoots_)

  integer(kind=iwp), intent(in) :: mode, nRoots_
  integer(kind=iwp) :: idum
  real(kind=wp) :: rdum

  nRoots = nRoots_

  if (mode == 1) then
    !! call from RASSCF, nothing is done
  else if (mode == 2) then
    !! call from MCLR, read from RunFile(?)
    call Get_iScalar('DWTypeSCF',idum)
    call Get_dScalar('DWZetaSCF',rdum)
    DWSCF%DWRoot = ibits(idum,16,16)
    DWSCF%DWType = ibits(idum,0,16)
    DWSCF%DWZeta = rdum
  end if

  DWSCF%do_DW = .false.
  if (DWSCF%DWROOT /= 0) then
    write(u6,*) 'DW-MCSCF is enabled'
    DWSCF%do_DW = .true.
    if (DWSCF%DWType == -99) then
      write(u6,*) 'DWTYpe has not been specified.'
      write(u6,*) 'The default value DWTYpe = 0 will be used.'
      DWSCF%DWType = 0
    end if
    if (DWSCF%DWZeta == Zero) then
      select case (DWSCF%DWType)
        case (-1,0,1)
          DWSCF%DWZeta = One
      end select
    end if
  end if

end subroutine DWSCF_init

!-----------------------------------------------------------------------

subroutine DWSol_DWRO(LuInput,nRoots_,i_all)

  integer(kind=iwp), intent(in) :: LuInput, i_all
  integer(kind=iwp), intent(inout) :: nRoots_
  integer(kind=iwp) :: i, iError, iSum
  character(len=180) :: KWord
  integer(kind=iwp), allocatable :: IROOT(:), Temp1(:)
  real(kind=wp), allocatable :: W_local(:)
  !integer(kind=iwp), external :: nToken
  character(len=180), external :: Get_Ln

  ! Determine the roots to be used for solvation by reading RFROOT
  !
  !i_all = 0
  !call Get_I1(1,nRoots_)
  !call Get_I1(2,lRoots_) !! this is not used
  !if (nToken(KWord_) > 3) call Get_I1(3,i_all)

  call mma_allocate(W_local,nRoots_,Label='W_local')
  call mma_allocate(IROOT,nRoots_,Label='IROOT')

  W_local(1:nRoots_) = Zero
  IROOT(1:nRoots_) = 0

  if (i_all == 1) then
    do i=1,nRoots_
      iroot(i) = i
      W_local(i) = One/real(nRoots_,kind=wp)
    end do
  else
    KWord = Get_Ln(LuInput)
    read(KWord,*,iostat=iError) (IROOT(I),I=1,nRoots_)
    if (iError /= 0) call Error(iError)
    if (nRoots_ == 1) then
      W_local(1) = One
    else
      KWord = Get_Ln(LuInput)

      call mma_allocate(Temp1,nRoots_,Label='Temp1')
      read(KWord,*,iostat=iError) (Temp1(I),I=1,nRoots_)
      if (iError /= 0) call Error(iError)

      iSum = 0
      do i=1,nRoots_
        iSum = iSum+Temp1(i)
      end do
      do i=1,nRoots_
        if (Temp1(i) == 0) then
          W_local(i) = Zero
        else
          W_local(i) = real(Temp1(i),kind=wp)/real(iSum,kind=wp)
        end if
      end do
      call mma_deallocate(Temp1)
    end if
  end if

  call Put_dArray('SolventWeight',W_local,nRoots_)

  ! Set IPCMROOT
  if (i_all == 1) then
    nRoots_ = 0
  else
    nRoots_ = maxloc(W_local(1:nRoots_),1)
  end if

  call mma_deallocate(W_local)
  call mma_deallocate(IROOT)

contains

  subroutine Error(rc)

    integer(kind=iwp), intent(in) :: rc

    if (rc > 0) then
      write(u6,*) 'DWSol_DWRO: Error while reading input'
    else
      write(u6,*) 'DWSol_DWRORdInp: Premature end of input file'
    end if
    write(u6,'(A,A)') 'Last command:',KWord
    call Abend()

  end subroutine Error

end subroutine DWSol_DWRO

!-----------------------------------------------------------------------

subroutine DWSol_final()

  call mma_deallocate(W_SOLV)

end subroutine DWSol_final

!-----------------------------------------------------------------------

subroutine DWSCF_final()

  integer(kind=iwp) :: idum
  real(kind=wp) :: rdum

  if (DWSCF%do_DW) then
    !! Save some parameters for MCLR
    idum = DWSCF%DWType
    call mvbits(DWSCF%DWRoot,0,16,idum,16)
    rdum = DWSCF%DWZeta
  else
    idum = 0
    rdum = Zero
  end if

  call Put_iScalar('DWTypeSCF',idum)
  call Put_dScalar('DWZetaSCF',rdum)

end subroutine DWSCF_final

!-----------------------------------------------------------------------

subroutine DWSol_wgt(mode,ENER,weight)

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: ENER(*)
  real(kind=wp), intent(out), optional :: weight(*)
  integer(kind=iwp) :: i
  real(kind=wp) :: Ealpha, Ebeta, Egamma, wgt, Wtot, xi_ab, xi_ag
  type(type_DW), pointer :: DWlocal
  real(kind=wp), allocatable :: W_local(:)

  select case (mode)
    case (1)
      !! DW-MCSCF
      DWlocal => DWSCF
    case (2)
      !! DW solvation
      DWlocal => DWSolv
    !case (3)
    !  !! DW CASPT2?
    case default
      nullify(DWlocal)
  end select

  if (DWlocal%DWZeta <= Zero) return

  !! dynamical weighting
  call mma_allocate(W_local,nRoots,Label='W_local')

  !! Normalization
  Wtot = Zero
  Ealpha = ENER(DWlocal%DWRoot)
  do i=1,nRoots
    Egamma = ENER(i)
    xi_ag = Ealpha-Egamma
    wgt = DWSol_func(xi_ag,DWlocal)
    Wtot = Wtot+wgt
  end do

  W_local = Zero
  do i=1,nRoots
    Ebeta = ENER(i)
    xi_ab = Ealpha-Ebeta
    wgt = DWSol_func(xi_ab,DWlocal)
    W_local(i) = wgt/Wtot
  end do

  select case (mode)
    case (1)
      weight(1:nRoots) = W_local(1:nRoots)
    case (2)
      W_SOLV(1:nRoots) = W_local(1:nRoots)
  end select

  call mma_deallocate(W_local)

end subroutine DWSol_wgt

!-----------------------------------------------------------------------

function DWSol_func(xx,DWlocal)

  real(kind=wp) :: DWSol_func
  real(kind=wp), intent(in) :: xx
  type(type_DW), intent(in) :: DWlocal

  DWSol_func = Zero
  select case (DWlocal%DWType)
    case (-1)
      !! the oldest (?) DW-CASSCF (J. Chem. Phys. 120, 7281 (2004)
      DWSol_func = cosh(xx/DWlocal%DWZeta)**2
    case (0)
      !! new (?) DW-CASSCF: J. Chem. Phys. 141, 171102 (2014)
      if (xx <= Zero) then
        DWSol_func = One/nRoots
      else if (xx <= DWlocal%DWZeta) then
        DWSol_func = One-Three*(xx/DWlocal%DWZeta)**2+Two*(xx/DWlocal%DWZeta)**3
      else
        DWSol_func = Zero
      end if
    case (1)
      !! Boltzmann factor (XDW-CASPT2 with DWType = 1)
      DWSol_func = exp(-DWlocal%DWZeta*xx**2)
    case default
      write(u6,*) 'Unrecognized DWType in DWSol...'
      call abend()
  end select

  return

end function DWSol_func

!-----------------------------------------------------------------------

subroutine DWSol_fixed(istate,jstate)

  integer(kind=iwp), intent(out) :: istate, jstate
  integer(kind=iwp) :: DWZetaInt

  DWZetaInt = int(DWSolv%DWZeta,kind=iwp)
  if (real(DWZetaInt,kind=wp) /= DWSolv%DWZeta) DWZetaInt = 0
  select case (DWZetaInt)
    case (-2)
      istate = 1
      jstate = 2
    case (-6)
      istate = 2
      jstate = 3
    case (-12)
      istate = 3
      jstate = 4
    case (-20)
      istate = 4
      jstate = 5
    case (-30)
      istate = 5
      jstate = 6
    case (-42)
      istate = 6
      jstate = 7
    case (-56)
      istate = 7
      jstate = 8
    case (-72)
      istate = 8
      jstate = 9
    case (-90)
      istate = 9
      jstate = 10
    case default
      write(u6,*) 'Unrecognized negative DWZeta ...'
      istate = 0
      jstate = 0
  end select

end subroutine DWSol_fixed

!-----------------------------------------------------------------------

subroutine DWSol_der(mode,DEROMG,DERHII,ENER,weight)

  integer(kind=iwp), intent(in) :: mode
  real(kind=wp), intent(in) :: DEROMG(:), ENER(:)
  !! partial derivative (pseudo-density, more precisely) of H_{II}
  real(kind=wp), intent(inout) :: DERHII(:)
  real(kind=wp), intent(in), optional :: weight(1:nRoots)
  integer(kind=iwp) :: i, j, k
  real(kind=wp) :: Ealpha, Ebeta, Egamma, Scal
  type(type_DW), pointer :: DWlocal
  real(kind=wp), allocatable :: W_local(:)

  select case (mode)
    case (1)
      !! DW-MCSCF
      DWlocal => DWSCF
    case (2)
      !! DW solvation
      DWlocal => DWSolv
    case default
      nullify(DWlocal)
  end select

  DERHII(1:nRoots) = Zero
  if (DWlocal%DWZeta <= Zero) return

  call mma_allocate(W_local,nRoots,Label='W_local')

  select case (mode)
    case (1)
      !! DW-MCSCF
      W_local(1:nRoots) = weight(1:nRoots)
    case (2)
      !! DW solvation
      W_local(1:nRoots) = W_SOLV(1:nRoots)
  end select

  i = DWlocal%DWRoot
  Ealpha = ENER(i)

  if (DWlocal%DWType == 1) then
    do j=1,nRoots
      Ebeta = ENER(j)
      Scal = -Two*DWlocal%DWZeta*W_local(j)*(Ealpha-Ebeta)*DEROMG(j)
      DERHII(i) = DERHII(i)-Scal
      DERHII(j) = DERHII(j)+Scal

      do k=1,nRoots
        Egamma = ENER(k)
        Scal = Two*DWlocal%DWZeta*W_local(j)*W_local(k)*(Ealpha-Egamma)*DEROMG(j)
        DERHII(i) = DERHII(i)-Scal
        DERHII(k) = DERHII(k)+Scal
      end do
    end do
  end if

  call mma_deallocate(W_local)

end subroutine DWSol_Der

end module DWSol
