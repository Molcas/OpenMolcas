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
! Copyright (C) 2026, Roland Lindh                                     *
!***********************************************************************

!> Run-time verification options for CI developments.
!>
!> Read once from the environment variable MOLCAS_CI_VERIFY, a
!> comma-separated list:
!>
!>   (unset) or NO  production code path: Lucia, and Faroald only with the
!>                  RASSCF keyword FARO; no SGUGA cross-checks or sanity
!>                  checks
!>   YES    coupling-coefficient test of <I|Epq|J> and <I|EpqErs|J>
!>          (SG_Test_EpqErs_Ref)
!>   EPQ    coupling-coefficient test of <I|Epq|J> only
!>   SIGMA  cross-check the CI back end against SGUGA: sigma vectors, CI
!>          transformation and density matrices (the SGUGA sigma vector is
!>          then used in CI_Lib_Sigma)
!>   FARO   use Faroald automatically for C1 CASSCF (compared with Lucia)
!>   HDIAG  give the CI solver the exact diagonal <I|H|I> from SGUGA instead
!>          of the configuration average (preconditioner and selection of the
!>          explicit subspace; an experiment, not included in ALL)
!>   SGEX   explicit Hamiltonian of the CI start from SGUGA, as with the
!>          RASSCF keyword SGEX (not included in ALL)
!>   ALL    YES, SIGMA and FARO
!>
!> CI_Verify('ANY') is true if any of YES, EPQ, SIGMA or FARO is requested;
!> it switches on the internal sanity checks of the SGUGA and SD/CSF codes.
!>
!> This module sits in sguga_util so that the SGUGA code itself can use it;
!> CI_Interfaces re-exports CI_Verify.
module CI_Verify_Options

use Definitions, only: iwp, u6

implicit none
private

logical(kind=iwp) :: Verify_Parsed = .false., Verify_Yes = .false., Verify_Epq = .false., Verify_Sigma = .false., &
                     Verify_Faro = .false., Verify_HDiag = .false., Verify_SGEx = .false.

public :: CI_Verify

contains

!> Is the verification option Opt (YES, EPQ, SIGMA, FARO, HDIAG, SGEX or ANY) active?
function CI_Verify(Opt) result(Active)

  character(len=*), intent(in) :: Opt
  logical(kind=iwp) :: Active
  character(len=128) :: Env
  character(len=16) :: Tok
  integer(kind=iwp) :: iPos, iComma

  if (.not. Verify_Parsed) then
    Env = ' '
    call get_environment_variable('MOLCAS_CI_VERIFY',Env)
    call UpCase(Env)
    iPos = 1
    do while (iPos <= len_trim(Env))
      iComma = index(Env(iPos:),',')
      if (iComma == 0) then
        Tok = adjustl(Env(iPos:))
        iPos = len_trim(Env)+1
      else
        Tok = adjustl(Env(iPos:iPos+iComma-2))
        iPos = iPos+iComma
      end if
      select case (trim(Tok))
        case ('YES')
          Verify_Yes = .true.
        case ('EPQ')
          Verify_Epq = .true.
        case ('SIGMA')
          Verify_Sigma = .true.
        case ('FARO')
          Verify_Faro = .true.
        case ('HDIAG')
          Verify_HDiag = .true.
        case ('SGEX')
          Verify_SGEx = .true.
        case ('ALL')
          Verify_Yes = .true.
          Verify_Sigma = .true.
          Verify_Faro = .true.
        case ('NO','')
        case default
          write(u6,'(3A)') 'Illegal value in MOLCAS_CI_VERIFY: "',trim(Tok),'"'
          write(u6,'(A)') 'Valid values (comma-separated): NO, YES, EPQ, SIGMA, FARO, HDIAG, SGEX, ALL'
          call Abend()
      end select
    end do
    Verify_Parsed = .true.
  end if

  select case (Opt)
    case ('YES')
      Active = Verify_Yes
    case ('EPQ')
      Active = Verify_Epq
    case ('SIGMA')
      Active = Verify_Sigma
    case ('FARO')
      Active = Verify_Faro
    case ('HDIAG')
      Active = Verify_HDiag
    case ('SGEX')
      Active = Verify_SGEx
    case ('ANY')
      Active = Verify_Yes .or. Verify_Epq .or. Verify_Sigma .or. Verify_Faro
    case default
      write(u6,'(3A)') 'CI_Verify called with illegal option "',Opt,'"'
      call Abend()
      Active = .false.
  end select

end function CI_Verify

end module CI_Verify_Options
