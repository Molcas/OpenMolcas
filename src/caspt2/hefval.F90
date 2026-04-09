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

subroutine hefval(ist, jst, dvalue)
  use definitions, only: wp, iwp, u6
  use caspt2_global, only: luciex, idtcex, iPrGlb
  use printLevel, only: debug
  use eqsolv, only: iVecC, iVecW
  use stdalloc, only: mma_allocate, mma_deallocate
  use caspt2_module, only: iSCF, nAshT, nConf, nState, STSym
#if defined _DMRG_
  use caspt2_module, only: DMRG
#endif
      use pt2_guga, only: MxCI
  implicit none
  ! Apart from input call parameters, we need two vectors stored on
  ! LUSOLV. Vector nr IVECC (presently=2) contains the contravariant
  ! elements of the solution to the CASPT2 equations.
  ! IVECW is the number (presently=6) of the vector on LUSOLV
  ! where a contravariant representation of the CASPT2 Right-Hand Side
  ! vector is stored. This depends on the MOs used, but is actually
  ! the same for all the root states.


  integer(kind=iwp), intent(in) :: ist, jst
  real(kind=wp), intent(out) :: dvalue

  integer(kind=iwp) :: i
  integer(kind=iwp) :: ntg1, ntg2, ntg3
  integer(kind=iwp) :: idci
  real(kind=wp) :: ovl, dummy(1)
  real(kind=wp), allocatable :: tg1(:), tg2(:), tg3(:), ci1(:), ci2(:)

  ! We evaluate the effective Hamiltonian matrix element in two steps.

  ntg1 = nasht**2
  ntg2 = nasht**4
  ntg3 = (ntg1*(ntg1+1)*(ntg1+2))/6
  ! Note: Need proper allocation even if unused, since allocated
  ! arrays are in subroutine parameter lists of MKTG3, HCOUP.
  ntg1 = max(1, ntg1)
  ntg2 = max(1, ntg2)
  ntg3 = max(1, ntg3)
  call mma_allocate(tg1, ntg1, label='TG1')
  call mma_allocate(tg2, ntg2, label='TG2')
  call mma_allocate(tg3, ntg3, label='TG3')
  tg1(1) = 0.0e0_wp
  tg2(1) = 0.0e0_wp
  tg3(1) = 0.0e0_wp

#if defined _DMRG_
  if (DMRG) then
    call mktg3qcm(stsym, stsym, ist, jst, ovl, tg1, tg2, ntg3, tg3)
  else if (.not. DMRG) then
#endif
    call mma_allocate(ci1, mxci, label='CI1')
    call mma_allocate(ci2, mxci, label='CI2')
    if (iscf == 0) then
      ! Read root vectors nr. IST and JST from LUCI.
      idci = idtcex
      do i = 1, nstate
        if (i == ist) then
          call ddafile(luciex, 2, ci1, nconf, idci)
          if (i == jst) then
            call dcopy_(nconf, ci1, 1, ci2, 1)
          end if
        else if (i == jst) then
          call ddafile(luciex, 2, ci2, nconf, idci)
        else
          call ddafile(luciex, 0, dummy, nconf, idci)
        end if
      end do
    end if

    if (iPrGlb >= debug) then
      write(u6,*) "=== VANILLA: Building TRANSITION-RDM === "
      write(u6,*) "between bra", ist, " and  ket", jst
    end if
    call mktg3(stsym, stsym, ci1, ci2, ovl, tg1, tg2, ntg3, tg3)
    call mma_deallocate(ci1)
    call mma_deallocate(ci2)
#if defined _DMRG_
  end if
#endif

  call hcoup(ivecw, ivecc, ovl, tg1, tg2, tg3, dvalue)

  call mma_deallocate(tg1)
  call mma_deallocate(tg2)
  call mma_deallocate(tg3)

end subroutine hefval
