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

subroutine Bond_Tester(Coor,nAtoms,iTab,nx,ny,nz,ix,iy,iz,iAtom,iRow,iANr,iTabBonds,nBonds,nBondMax,iTabAtoms,nMax,ThrB,ThrB_vdW)

use Slapaf_Info, only: Covalent_Bond, ddV_Schlegel, iOptC, vdW_Bond
#ifdef _DEBUGPRINT_
use Slapaf_Info, only: BondType
#endif
use ddvdt, only: aAV, alpha_vdW, r_ref_vdW, rAV
use Constants, only: Zero, One, Two
#define _OLD_CODE_
#ifndef _OLD_CODE_
use Constants, only: Four, Pi
#endif
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: nAtoms, nx, ny, nz, nMax, iTab(0:nMax,nx,ny,nz), ix, iy, iz, iAtom, iANr(nAtoms), nBondMax
real(kind=wp), intent(in) :: Coor(3,nAtoms), ThrB, ThrB_vdW
integer(kind=iwp), intent(out) :: iRow
integer(kind=iwp), intent(inout) :: iTabBonds(3,nBondMax), nBonds, iTabAtoms(2,0:nMax,nAtoms)
integer(kind=iwp) :: i, Ir, ivdW, jAtom, jRow, nn, nNeighbor, Nr, nVal_i, nVal_j
real(kind=wp) :: alpha, r0, r0_vdw, Rab, RabCov, rij2, test, test_vdw, x, y, z
#ifndef _OLD_CODE_
integer(kind=iwp) :: iBond, iType, kAtom
real(kind=wp) :: A(3), ANorm, B(3), BNorm, CosPhi, Phi
#endif
logical(kind=iwp) :: Help
integer(kind=iwp), external :: iTabRow
real(kind=wp), external :: CovRad

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
! Check box indices for consistency

if ((ix < 1) .or. (ix > nx)) return
if ((iy < 1) .or. (iy > ny)) return
if ((iz < 1) .or. (iz > nz)) return
Nr = iTab(0,ix,iy,iz) ! nr of atoms in the box.
if (Nr == 0) return

iRow = iTabRow(iANr(iAtom))
nVal_i = 0
nn = iTabAtoms(1,0,iAtom)
do i=1,nn
  if (iTabBonds(3,iTabAtoms(2,i,iAtom)) == Covalent_Bond) nVal_i = nVal_i+1
end do
#ifdef _DEBUGPRINT_
write(u6,*)
write(u6,*) 'Bond_Tester: iAtom,=',iAtom
write(u6,*) '                Box(ix,iy,iz)=(',ix,iy,iz,')'
write(u6,*) '                        nVal_i=',nVal_i
write(u6,*)
#endif

! Loop over all atoms in the box

if (ThrB < ThrB_vdw) ivdW = vdW_Bond
box: do Ir=1,Nr
  jAtom = iTab(Ir,ix,iy,iz)
  !if (iAtom <= jAtom) cycle box
  if (iAtom >= jAtom) cycle box
  jRow = iTabRow(iANr(jAtom))
  Help = (iRow > 3) .or. (jRow > 3)
# ifdef _DEBUGPRINT_
  write(u6,*) ' jAtom, iAnr(jAtom)=',jAtom,iAnr(jAtom)
  write(u6,*) 'Help=',Help
# endif

  x = Coor(1,iAtom)-Coor(1,jAtom)
  y = Coor(2,iAtom)-Coor(2,jAtom)
  z = Coor(3,iAtom)-Coor(3,jAtom)
# ifndef _OLD_CODE_
  A(:) = Coor(:,iAtom)-Coor(:,jAtom)
  ANorm = dot_product(A,A)
# endif
  rij2 = x**2+y**2+z**2
  r0 = rAv(iRow,jRow)
  alpha = aAv(iRow,jRow)

  ! Test if we have a bond iAtom-jAtom

  if (ddV_Schlegel .or. Help) then
    Rab = sqrt(rij2)
    RabCov = CovRad(iANr(iAtom))+CovRad(iANr(jAtom))
#   ifdef _DEBUGPRINT_
    write(u6,*) 'Rab=',Rab
    write(u6,*) CovRad(iANr(iAtom)),CovRad(iANr(jAtom))
#   endif
    if (Rab <= 1.25_wp*RabCov) then

      ! covalent bond

      test = One
      test_vdW = Zero

      ! Skip if we are looking for vdW bonds

      if (ThrB > ThrB_vdW) cycle box
    else if ((Rab > 1.25_wp*RabCov) .and. (Rab <= Two*RabCov)) then

      ! vdW's bond

      test = Zero
      test_vdW = ThrB_VdW
    else

      ! No Bond!

      cycle box
    end if

  else

    test = exp(alpha*(r0**2-rij2))
    if (btest(iOptC,11)) then
      r0_vdW = r_ref_vdW(iRow,jRow)
      test_vdW = exp(-alpha_vdW*(r0_vdW-sqrt(rij2))**2)
    else
      r0_vdW = Zero
      test_vdW = Zero
    end if
#   ifdef _DEBUGPRINT_
    write(u6,*)
    write(u6,*) 'Bond_Tester: iAtom,jAtom=',iAtom,jAtom
    write(u6,*) 'Bond_Tester: iRow, jRow =',iRow,jRow
    write(u6,*) 'Bond_Tester: test=',test,ThrB
    write(u6,*) 'Bond_Tester: test_vdW=',test_vdW,ThrB_vdW
    write(u6,*) 'Bond_Tester: r0_vdW=',r0_vdW,sqrt(rij2)
#   endif

    ! If the valence force constant small but not too small
    ! denote the bond as an vdW's bond.
    Test_vdW = max(Test_vdW,Test)

    ! If already valence bond skip if also vdW bond. We picked
    ! up this bond before!

    if ((test >= ThrB) .and. (Test_vdW >= ThrB_vdW)) cycle box

    ! If none skip

    if ((test < ThrB) .and. (test_vdW < ThrB_vdW)) cycle box

    ! Some logic to see if vdw bond should be included.
    ! Hydrogen-hydrogen is always included.

    if ((iANr(iAtom) /= 1) .or. (iANr(jAtom) /= 1)) then

      ! Skip if any of the atoms has more than 6 valence bonds
      ! and the other at least 1 valence bond.
      nVal_j = 0
      nn = iTabAtoms(1,0,jAtom)
      do i=1,nn
        if (iTabBonds(3,iTabAtoms(2,i,jAtom)) == Covalent_Bond) nVal_j = nVal_j+1
      end do
#     ifdef _DEBUGPRINT_
      write(u6,*) 'nVal_j=',nVal_j
#     endif
      if (((nVal_i >= 6) .and. (nVal_j >= 1)) .or. ((nVal_j >= 6) .and. (nVal_i >= 1))) cycle box
    end if
#   ifndef _OLD_CODE_

    ! We need to exclude vdW bonds if there is an atom close to
    ! being in between the two atoms being considered and
    ! forming a covalent bond.

    if (test < ThrB) then ! only in case of vdW bond

      ! Loop over all covalently bonded neighbors of atom iAtom

#     ifdef _DEBUGPRINT_
      write(u6,*) 'Test validity of vdW bond'
#     endif
      do i=1,iTabAtoms(1,0,iAtom)
        kAtom = iTabAtoms(1,i,iAtom)
        iBond = iTabAtoms(2,i,iAtom)
        iType = iTabBonds(3,iBond)
        if (iType /= Covalent_Bond) cycle
        B(:) = Coor(:,iAtom)-Coor(:,kAtom)
        BNorm = dot_product(B,B)
        CosPhi = dot_product(A,B)/sqrt(ANorm*BNorm)
        Phi = acos(CosPhi)
#       ifdef _DEBUGPRINT_
        write(u6,*) 'kAtom,Phi=',kAtom,Phi
#       endif
        if (abs(Phi) < Pi/Four) cycle box
      end do
    end if
#endif

  end if

  if (nBonds+1 > nBondMax) then
    write(u6,*) 'Bond_Tester: nBonds+1 > nBondMax'
    write(u6,*) 'nBonds+1=',nBonds+1
    write(u6,*) 'nBondMax=',nBondMax
    call Abend()
  end if

  nBonds = nBonds+1
  iTabBonds(1,nBonds) = iAtom
  iTabBonds(2,nBonds) = jAtom

  if (test >= ThrB) then
    ivdW = Covalent_Bond
  else if (test_vdW >= ThrB_vdW) then
    ivdW = vdW_Bond
  else
    write(u6,*) 'Bond_Tester: Illegal operation'
    call Abend()
    ivdW = 99
  end if
  iTabBonds(3,nBonds) = ivdW
# ifdef _DEBUGPRINT_
  write(u6,*)
  write(u6,*) 'Add bond to the list.'
  write(u6,*) 'Atom pair:',iAtom,jAtom
  write(u6,*) 'Bond type:',Bondtype(min(3,ivdW))
  write(u6,*)
# endif

  nNeighbor = iTabAtoms(1,0,iAtom)
  if (nNeighbor+1 > nMax) then
    write(u6,*) 'Bond_Tester(1): nNeighbor+1 > nMax'
    write(u6,*) 'iAtom=',iAtom
    write(u6,*) 'nNeighbor=',nNeighbor
    write(u6,*) 'nMax=',nMax
    call Abend()
  end if

  ! Update the neighbor lists of atoms iAtom and jAtom
  nNeighbor = nNeighbor+1
  iTabAtoms(1,0,iAtom) = nNeighbor
  iTabAtoms(1,nNeighbor,iAtom) = jAtom
  iTabAtoms(2,nNeighbor,iAtom) = nBonds

  nNeighbor = iTabAtoms(1,0,jAtom)
  if (nNeighbor+1 > nMax) then
    write(u6,*) 'Bond_Tester(2): nNeighbor+1 > nMax'
    write(u6,*) 'jAtom=',jAtom
    write(u6,*) 'nNeighbor=',nNeighbor
    write(u6,*) 'nMax=',nMax
    call Abend()
  end if
  nNeighbor = nNeighbor+1
  iTabAtoms(1,0,jAtom) = nNeighbor
  iTabAtoms(1,nNeighbor,jAtom) = iAtom
  iTabAtoms(2,nNeighbor,jAtom) = nBonds

end do box

return

end subroutine Bond_Tester
