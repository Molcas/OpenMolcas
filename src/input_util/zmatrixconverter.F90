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
! Copyright (C) 2006, Giovanni Ghigo                                   *
!***********************************************************************
!  ZMatrixConverter
!
!> @brief
!>   Routine for reading seward input in Z-Matrix format
!> @author Giovanni Ghigo
!>
!> @details
!> The input for seward is read and a string vector is generate (\p STDINP).
!> This vector contains a standard seward input and it will be read later as usual
!> by a modified copy of the section ``BASI`` code already present in ::RdCtl_Seward.
!> This new code is in the ::StdSewInput routine.
!> Only the standard basis present in the ``$MOLCAS/basis_library`` are allowed.
!>
!> @param[in]     LuRd    Input file unit number
!> @param[in]     LuWr    Output file unit number
!> @param[in]     mxAtom  Parameter
!> @param[out]    STDINP  String vector of seward standard input
!> @param[out]    lSTDINP Length of String vector \p STDINP
!> @param[in]     iglobal
!> @param[in,out] nxbas
!> @param[in]     xb_label
!> @param[in]     xb_bas
!> @param[out]    iErr    Error flag
!***********************************************************************

subroutine ZMatrixConverter(LuRd,LuWr,mxAtom,STDINP,lSTDINP,iglobal,nxbas,xb_label,xb_bas,iErr)
!***********************************************************************
! Author: Giovanni Ghigo                                               *
!         Torino (Italy)  October-November 2006                        *
!                                                                      *
! This is an adaptation of Program ZMatrixConverter                    *
! A converter of Z-Matrix in cartesian coordinates in MolCAS format.   *
! Version 1.0                                                          *
! The input for seward is read and a string vector is generate. This   *
! vector is a standard seward input and it will be read later as usual *
! by a copy on the code already present in RdCtl_Seward.               *
!***********************************************************************

use ZMatConv_Mod, only: BasAva, Base, BasReq, Coords, iZmat, MaxAtoms, NAT, Symbols, Zmat
use isotopes, only: MaxAtomNum
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, deg2rad
use Definitions, only: wp, iwp

implicit none
integer(kind=iwp), intent(in) :: LuRd, LuWr, mxAtom, iglobal
character(len=180), intent(out) :: STDINP(mxAtom*2)
integer(kind=iwp), intent(out) :: lSTDINP, iErr
integer(kind=iwp), intent(inout) :: nxbas
character(len=*), intent(in) :: xb_label(*), xb_bas(*)
integer(kind=iwp) :: i, iAtom, iSTDINP, j, k, nAtoms, NATprev, nBase, nBasis, nXAtoms
logical(kind=iwp) :: IfTest
real(kind=wp) :: r
character(len=180) :: aDebug
character(len=12) :: Angstring

#ifdef _DEBUGPRINT_
IfTest = .true.
#else
IfTest = .false.
#endif

!  ***  X dummy atoms (NA = 0 )
!  ***  Z ghost atoms (NA =-1 )
!  ***  nAskAtoms == -1  =>  Seward ZMAT input
!  ***  nAskAtoms /= -1  =>  GateWay ZMAT input

! nAtoms : nr. of atoms passed to SEWARD (includes X dummy atoms).
! nXAtoms: nr. of ghost Z atoms (not passed to SEWARD but resumed by OutZMat in SLAPAF).
! nBase  : number of BasisSets found in input.
! Base(i): BasisSet for atom with Atomic Number -i-.
! BasAva(i) & BasReq(i): Logical to check BasisSet-consistency.
! Coords(_,i): X, Y, Z, coordinates (in Angstrom) for atom -i-.

nAtoms = 0
nXAtoms = 0
nBase = 0
lSTDINP = 0
nBasis = 0
iErr = 0
Angstring = '  / Angstrom'

call mma_allocate(Base,MaxAtomNum,label='Base')
call mma_allocate(BasAva,MaxAtomNum,label='BasAva')
call mma_allocate(BasReq,MaxAtomNum,label='BasReq')
call mma_allocate(NAT,MaxAtoms,label='NAT')
call mma_allocate(Symbols,MaxAtoms,label='Symbols')
call mma_allocate(iZmat,3,MaxAtoms,label='iZmat')
call mma_allocate(Zmat,3,MaxAtoms,label='Zmat')
Base(:) = ''
BasAva(:) = .false.
BasReq(:) = .false.
NAT(:) = 0
Symbols(:) = ''
iZmat(:,:) = 0
Zmat(:,:) = Zero

! Reading input
call BasisReader(LuWr,nBase,iglobal,nxbas,xb_label,xb_bas,iErr)
if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'ZMatrixConverter - From BasisReader :'
  write(LuWr,*) '                   nBase=',nBase
  do i=1,size(Base)
    if (BasAva(i)) write(LuWr,'(I23,3X,A)') i,Base(i)
  end do
  write(LuWr,*)
end if
if (iErr /= 0) then
  write(LuWr,*) ' ERROR: Wrong input in Bases Set definition !'
  return
end if
call ZMatReader(LuRd,LuWr,nAtoms,nXAtoms,nBasis,-1,iErr)
if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'ZMatrixConverter - From ZMatReader :'
  write(LuWr,*) '                   nAtoms=',nAtoms,', nXAtoms=',nXAtoms,', Tot=',nAtoms+nXAtoms
  write(LuWr,*) ' Label NA   i    bond        j    angle       k    dihedral'
  do i=1,nAtoms+nXAtoms
    write(LuWr,'(1X,A,I3,3(1X,I4,1X,F11.6))') Symbols(i),NAT(i),(iZmat(j,i),Zmat(j,i),j=1,3)
  end do
  write(LuWr,*)
end if
if (iErr /= 0) then
  write(LuWr,*) ' ERROR: Wrong input in Z-Matrix definition !'
  return
end if

! Some checks
if (nBase == 0) then
  iErr = 1
  write(LuWr,*) 'ERROR: No basis set specified !'
  return
end if
if (nAtoms == 0) then
  iErr = 1
  write(LuWr,*) 'ERROR: No atom coordinates specified !'
  return
end if
if (nBase < nBasis) then
  iErr = 1
  write(LuWr,*) 'ERROR: Wrong number of basis sets !'
  write(LuWr,*) '       Available=',nBase,'  Required=',nBasis
  return
end if
call BasisConsistency(LuWr,iErr)
if (iErr /= 0) then
  iErr = 1
  write(LuWr,*) 'ERROR: Basis set inconsistency !'
  return
end if

call Put_iScalar('N ZMAT',nAtoms+nXAtoms)
call Put_cArray('Symbol ZMAT',Symbols(1),(nAtoms+nXAtoms)*len(Symbols))
call Put_iArray('Index ZMAT',iZmat,(nAtoms+nXAtoms)*3)
call Put_iArray('NAT ZMAT',NAT,nAtoms+nXAtoms)

call mma_allocate(Coords,3,nAtoms+nXAtoms,label='Coords')
Coords(:,:) = Zero

! Calculate coordinates
! Atom #1
if (nAtoms+nXAtoms > 1) then
  ! Atom #2
  Coords(3,2) = Zmat(1,2)  ! Z(2)=R
end if
if (nAtoms+nXAtoms > 2) then
  ! Atom #3
  if (iZmat(1,3) == 1) then
    Coords(1,3) = Zmat(1,3)*sin(Zmat(2,3)*deg2rad) ! X(2)=R sin(A)
    Coords(3,3) = Zmat(1,3)*cos(Zmat(2,3)*deg2rad) ! Z(3)=R cos(A)
  else
    Coords(1,3) = Zmat(1,3)*sin(Zmat(2,3)*deg2rad)
    Coords(3,3) = Coords(3,2)-Zmat(1,3)*cos(Zmat(2,3)*deg2rad)
  end if
end if
if (nAtoms+nXAtoms > 3) then
  ! Atom #4 ->
  do iAtom=4,nAtoms+nXAtoms
    call ZMatConv(LuWr,iAtom,iErr)
  end do
  if (iErr /= 0) return

  if (IfTest) then
    write(LuWr,*)
    write(LuWr,*) '------------------------------------------------'
    write(LuWr,*) 'ZMatrixConverter - XYZCoords (Angstroms) :'
    do i=1,nAtoms+nXAtoms
      write(LuWr,99) i,NAT(i),(Coords(j,i),j=1,3)
    end do
    write(LuWr,*)
  end if

  ! Check for superposed atoms
  do i=1,nAtoms+nXAtoms
    if (NAT(i) > 0) then
      do j=i+1,nAtoms+nXAtoms
        if (NAT(j) > 0) then
          r = Zero
          do k=1,3
            r = r+(Coords(k,i)-Coords(k,j))**2
          end do
          if (r < 1.0e-4_wp) then
            iErr = 1
            write(LuWr,*) ' ERROR: Superimposed atoms: ',i,j,'  r=',sqrt(r)
            return
          end if
        end if
      end do
    end if
  end do
end if

! Writing
if (NAT(1) == -1) then
  NATprev = -1
else
  NATprev = -9999
end if
iSTDINP = 1
do i=1,nAtoms+nXAtoms
  if (NAT(i) == -1) cycle
  if (NAT(i) /= NATprev) then
    if ((i /= 1) .and. (NATprev /= -1)) then
      write(STDINP(iSTDINP),'(A)') 'End of basis'
      iSTDINP = iSTDINP+1
    end if
    write(STDINP(iSTDINP),'(A)') 'Basis set'
    iSTDINP = iSTDINP+1
    if (NAT(i) > 0) then
      write(STDINP(iSTDINP),'(A)') Base(NAT(i))
    else
      write(STDINP(iSTDINP),'(A)') 'X..... / InLine '
      iSTDINP = iSTDINP+1
      write(STDINP(iSTDINP),'(A)') '0.   0'
      iSTDINP = iSTDINP+1
      write(STDINP(iSTDINP),'(A)') '0    0'
    end if
    iSTDINP = iSTDINP+1
    NATprev = NAT(i)
  end if
  write(STDINP(iSTDINP),'(A5,3F16.10,A)') Symbols(i),(Coords(j,i),j=1,3),Angstring
  iSTDINP = iSTDINP+1
end do
write(STDINP(iSTDINP),'(A)') 'End of basis'
lSTDINP = iSTDINP
if (IfTest) then
  write(LuWr,*)
  write(LuWr,*) '------------------------------------------------'
  write(LuWr,*) 'ZMatrixConverter - The input passed to SEWARD : '
  do i=1,iSTDINP
    aDebug = STDINP(i)
    write(LuWr,*) aDebug
  end do
  write(LuWr,*)
end if

call mma_deallocate(Base)
call mma_deallocate(BasAva)
call mma_deallocate(BasReq)
call mma_deallocate(NAT)
call mma_deallocate(Symbols)
call mma_deallocate(iZmat)
call mma_deallocate(Zmat)
call mma_deallocate(Coords)

return

99 format(I3,1X,I3,1X,3(F12.6))

end subroutine ZMatrixConverter
