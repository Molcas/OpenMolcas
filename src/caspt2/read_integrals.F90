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
! Copyright (C) 2025, Stefano Battaglia                                *
!***********************************************************************

<<<<<<< HEAD
subroutine read_integrals()

use Symmetry_Info, only: Mul
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6
use caspt2_global, only: FIMO, iPrGlb
use printLevel, only: debug
use rasscf_global, only: Emy
#ifdef _DMRG_
use qcmaquis_interface, only: qcmaquis_interface_update_integrals_C, &
                              qcmaquis_interface_optimize, &
                              qcmaquis_interface_set_state, &
                              qcmaquis_interface_remove_param
      use caspt2_module, only: nState
#endif
use iso_c_binding, only: c_int
      use caspt2_module, only: nAshT, nSym, nFro, nIsh, nAsh, nOsh, nOrb

implicit none


real(kind=wp), allocatable :: ERI(:,:), SCR(:,:)
integer(kind=iwp) :: tSym, t_nFro, t_nIsh, t_nAsh, t_nOsh, t_nOrb, &
                     uSym, u_nFro, u_nIsh, u_nAsh, u_nOsh, u_nOrb, &
                     tuSym, &
                     vSym, v_nIsh, v_nAsh, v_nOsh, v, &
                     xSym, x_nIsh, x_nAsh, x_nOsh, x, t, u, n

! Arrays that get passed to QCMaquis
integer(c_int), dimension(:), allocatable :: indices
real(kind=wp), dimension(:), allocatable :: values
! calculate the size of the arrays
integer(kind=iwp) :: arr_size = 0
integer(kind=iwp) :: max_index2
! Offset for values and indices array
integer(kind=iwp) :: offset
real(kind=wp) , parameter                :: threshold = 1.0e-16_wp
integer(kind=iwp) NACPAR

write(u6,*) "=== QCM: Rotating Orbitals to SS === "
=======
#include "compiler_features.h"
#ifdef _DMRG_

subroutine read_integrals()

use, intrinsic :: iso_c_binding, only: c_int
use Index_Functions, only: iTri, nTri_Elem
use Symmetry_Info, only: Mul
use printLevel, only: debug
use rasscf_global, only: Emy
use qcmaquis_interface, only: qcmaquis_interface_optimize, qcmaquis_interface_remove_param, qcmaquis_interface_set_state, &
                              qcmaquis_interface_update_integrals_C
use caspt2_global, only: FIMO, iPrGlb
use general_data, only: NASH
use caspt2_module, only: nAshT, nFro, nIsh, nOrb, nOsh, nState, nSym
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp) :: arr_size, max_index2, n, NACPAR, offset, t, t_nAsh, t_nFro, t_nIsh, t_nOrb, t_nOsh, tSym, tuSym, u, u_nAsh, &
                     u_nFro, u_nIsh, u_nOrb, u_nOsh, uSym, v, v_nAsh, v_nIsh, v_nOsh, vSym, x, x_nAsh, x_nIsh, x_nOsh, xSym
integer(c_int), allocatable :: indices(:)
real(kind=wp), allocatable :: ERI(:,:), SCR(:,:), values(:)
real(kind=wp), parameter :: threshold = 1.0e-16_wp

! indices, values: Arrays that get passed to QCMaquis
! arr_size, max_index2: calculate the size of the arrays
! offset: Offset for values and indices array

write(u6,*) '=== QCM: Rotating Orbitals to SS === '
>>>>>>> upstream-openmolcas/master

if (iPrGlb >= debug) then
  write(u6,*) 'ERI in MO-basis'
  write(u6,*) '---------------'
end if

<<<<<<< HEAD
NACPAR = (nAshT*(nAshT + 1))/2

! calculate the maximum size of the integral array, can be less b/c of symmetry
max_index2 = NACPAR*(NACPAR+1)/2

arr_size = 1 & ! core energy
         + max_index2 & ! one-electron integrals
         + max_index2*(max_index2+1)/2 ! two-electron integrals
call mma_allocate(indices, 4*arr_size)
indices(:) = 0
call mma_allocate(values, arr_size)
=======
NACPAR = nTri_Elem(nAshT)

! calculate the maximum size of the integral array, can be less b/c of symmetry
max_index2 = nTri_Elem(NACPAR)

! core energy + one-electron integrals + two-electron integrals
arr_size = 1+max_index2+nTri_Elem(max_index2)
call mma_allocate(indices,4*arr_size)
indices(:) = 0
call mma_allocate(values,arr_size)
>>>>>>> upstream-openmolcas/master
values(:) = Zero

! Two Body
offset = 1
do tSym=1,nSym
  t_nFro = nFro(tSym)
  t_nIsh = nIsh(tSym)
  t_nAsh = nAsh(tSym)
  t_nOsh = nOsh(tSym)
  t_nOrb = nOrb(tSym)

  do uSym=1,tSym
    u_nFro = nFro(uSym)
    u_nIsh = nIsh(uSym)
    u_nAsh = nAsh(uSym)
    u_nOsh = nOsh(uSym)
    u_nOrb = nOrb(uSym)
    tuSym = Mul(tSym,uSym)

<<<<<<< HEAD
    call mma_allocate(ERI, t_nOrb, u_nOrb, Label='ERI')
    call mma_allocate(SCR, t_nOrb, u_nOrb, Label='SCR')
=======
    call mma_allocate(ERI,t_nOrb,u_nOrb,Label='ERI')
    call mma_allocate(SCR,t_nOrb,u_nOrb,Label='SCR')
>>>>>>> upstream-openmolcas/master
    ERI(:,:) = Zero
    SCR(:,:) = Zero

    if (iPrGlb >= debug) then
<<<<<<< HEAD
      write(u6,*) ' tSym = ', tSym, ' uSym = ', uSym
      write(u6,*) ' t_nFro = ', t_nFro, ' u_nFro = ', u_nFro
      write(u6,*) ' t_nIsh = ', t_nIsh, ' u_nIsh = ', u_nIsh
      write(u6,*) ' t_nAsh = ', t_nAsh, ' u_nAsh = ', u_nAsh
      write(u6,*) ' t_nOsh = ', t_nOsh, ' u_nOsh = ', u_nOsh
      write(u6,*) ' t_nOrb = ', t_nOrb, ' u_nOrb = ', u_nOrb
=======
      write(u6,*) ' tSym = ',tSym,' uSym = ',uSym
      write(u6,*) ' t_nFro = ',t_nFro,' u_nFro = ',u_nFro
      write(u6,*) ' t_nIsh = ',t_nIsh,' u_nIsh = ',u_nIsh
      write(u6,*) ' t_nAsh = ',t_nAsh,' u_nAsh = ',u_nAsh
      write(u6,*) ' t_nOsh = ',t_nOsh,' u_nOsh = ',u_nOsh
      write(u6,*) ' t_nOrb = ',t_nOrb,' u_nOrb = ',u_nOrb
>>>>>>> upstream-openmolcas/master
    end if

    do vSym=1,nSym
      v_nIsh = nIsh(vSym)
      v_nAsh = nAsh(vSym)
      v_nOsh = nOsh(vSym)
      xSym = Mul(tuSym,vSym)

      ! only totally symmetric integrals are different from zero
      if (xSym <= vSym) then
        x_nIsh = nIsh(xSym)
        x_nAsh = nAsh(xSym)
        x_nOsh = nOsh(xSym)

<<<<<<< HEAD
        write(u6,*) ' vSym = ', vSym, ' xSym = ', xSym
        write(u6,*) ' v_nIsh = ', v_nIsh, ' x_nIsh = ', x_nIsh
        write(u6,*) ' v_nAsh = ', v_nAsh, ' x_nAsh = ', x_nAsh
        write(u6,*) ' v_nOsh = ', v_nOsh, ' x_nOsh = ', x_nOsh

        ! loop over active orbital indices only
        do v = v_nIsh+1, v_nOsh
          do x = v, x_nOsh
            call coul(tSym, uSym, vSym, xSym, v, x, ERI, SCR)
            do t = v, t_nOsh
              do u = t, u_nOsh
                if (v .eq. t .and. u < x) then
                  cycle
                end if
                if (abs(ERI(t,u)) < threshold) then
                  cycle
                end if
                values(offset) = ERI(t,u)
                indices(4*(offset-1)+1:4*(offset-1)+4) = &
                  (/ int(v - v_nIsh, kind=c_int), int(x - x_nIsh, kind=c_int), &
                     int(t - t_nIsh, kind=c_int), int(u - u_nIsh, kind=c_int) /)
                offset = offset + 1
=======
        write(u6,*) ' vSym = ',vSym,' xSym = ',xSym
        write(u6,*) ' v_nIsh = ',v_nIsh,' x_nIsh = ',x_nIsh
        write(u6,*) ' v_nAsh = ',v_nAsh,' x_nAsh = ',x_nAsh
        write(u6,*) ' v_nOsh = ',v_nOsh,' x_nOsh = ',x_nOsh

        ! loop over active orbital indices only
        do v=v_nIsh+1,v_nOsh
          do x=v,x_nOsh
            call coul(tSym,uSym,vSym,xSym,v,x,ERI,SCR)
            do t=v,t_nOsh
              do u=t,u_nOsh
                if ((v == t) .and. (u < x)) cycle
                if (abs(ERI(t,u)) < threshold) cycle
                values(offset) = ERI(t,u)
                indices(4*(offset-1)+1:4*(offset-1)+4) = [int(v-v_nIsh,kind=c_int),int(x-x_nIsh,kind=c_int), &
                                                          int(t-t_nIsh,kind=c_int),int(u-u_nIsh,kind=c_int)]
                offset = offset+1
>>>>>>> upstream-openmolcas/master
              end do
            end do
          end do
        end do
      end if
      call mma_deallocate(ERI)
      call mma_deallocate(SCR)
    end do
  end do
end do

! One Body
do vSym=1,nSym
  v_nIsh = nIsh(vSym)
  v_nAsh = nAsh(vSym)
  v_nOsh = nOsh(vSym)
<<<<<<< HEAD
  do v = v_nIsh+1, v_nOsh
    do x = x_nIsh+1, v
      t = max(v, x)
      u = min(v, x)
      n = u + t * (t - 1) / 2
      if (abs(FIMO(n)) < threshold) then
        cycle
      end if
      values(offset) = FIMO(n)
      indices(4*(offset-1)+1:4*(offset-1)+4) = &
      (/ int(t- v_nIsh, kind=c_int), int(u - v_nIsh, kind=c_int), 0_c_int, 0_c_int /)
      offset = offset + 1
    END DO
  END DO
END DO

! Core Energy
offset = offset + 1
values(offset) = EMY
indices(4*(offset-1)+1:4*(offset-1)+4) = &
(/ 0_c_int, 0_c_int, 0_c_int, 0_c_int /)

! Rotate MPS wavefunction to new orbitals
#ifdef _DMRG_
call qcmaquis_interface_update_integrals_C(indices, values, int(offset, kind=c_int))
=======
  do v=v_nIsh+1,v_nOsh
    do x=x_nIsh+1,v
      n = iTri(v,x)
      t = max(v,x)
      u = min(v,x)
      if (abs(FIMO(n)) < threshold) cycle
      values(offset) = FIMO(n)
      indices(4*(offset-1)+1:4*(offset-1)+4) = [int(t-v_nIsh,kind=c_int),int(u-v_nIsh,kind=c_int),0_c_int,0_c_int]
      offset = offset+1
    end do
  end do
end do

! Core Energy
offset = offset+1
values(offset) = EMY
indices(4*(offset-1)+1:4*(offset-1)+4) = 0_c_int

! Rotate MPS wavefunction to new orbitals
call qcmaquis_interface_update_integrals_C(indices,values,int(offset,kind=c_int))
>>>>>>> upstream-openmolcas/master
do n=1,NSTATE
  call qcmaquis_interface_remove_param('MEASURE[trans1rdm]')
  call qcmaquis_interface_remove_param('MEASURE[trans2rdm]')
  call qcmaquis_interface_remove_param('MEASURE[trans3rdm]')

<<<<<<< HEAD
  call qcmaquis_interface_set_state(int(n-1, c_int))
  call qcmaquis_interface_optimize()
end do
#endif
=======
  call qcmaquis_interface_set_state(int(n-1,c_int))
  call qcmaquis_interface_optimize()
end do
>>>>>>> upstream-openmolcas/master

call mma_deallocate(indices)
call mma_deallocate(values)

end subroutine read_integrals
<<<<<<< HEAD
=======

#elif ! defined (EMPTY_FILES)

! Some compilers do not like empty files
#include "macros.fh"
dummy_empty_procedure(read_integrals)

#endif
>>>>>>> upstream-openmolcas/master
