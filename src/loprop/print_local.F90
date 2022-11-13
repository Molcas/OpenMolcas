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

subroutine Print_Local(rMP,nij,nElem,Coor,nAtoms,C_o_C,Q_Nuc,lMax,Lbl_Center,rMPq,EC,Polar,NoField,Temp,xrMP,xxRMP,xnrMP,iANr, &
                       nOcOb,Energy_Without_FFPT,Ene_Occ,MpProp_Level,Bond_Threshold,ChPolBB,LIonize)

use Real_Spherical, only: Sphere, Sphere_Free
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One, Two, Three, Half, auToeV
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp), intent(in) :: nij, nElem, nAtoms, lMax, iANr(nAtoms), nOcOb, MpProp_Level
real(kind=wp), intent(inout) :: rMP(nij,nElem), Temp(nij)
real(kind=wp), intent(in) :: Coor(3,nAtoms), C_o_C(3), Q_Nuc(nAtoms), rMPq(nElem), EC(3,nij), Polar(6,nij), Energy_Without_FFPT, &
                             Ene_Occ(nOcOb), Bond_Threshold, ChPolBB(6,nij)
character(len=LenIn), intent(in) :: Lbl_Center(nAtoms)
logical(kind=iwp), intent(in) :: NoField, LIonize
real(kind=wp), intent(out) :: xrMP(nij,nElem), xxrMP(nij,nElem), xnrMP(nij,nElem)
integer(kind=iwp) :: iAtom, i_Dim, iElem, iEnd, ij, iPol, iPrint, iStrt, ix, iy, iz, jAtom, l, nReal_Centers
real(kind=wp) :: A(3), Charge_center, CRN, CRX, CRY, CRZ, CX(3), Dip_Tot, DR, LA, LI, LI_MIN, LI_TOT, Pol_Tot, Polar_M(6), &
                 rMP_Tot, rMP_Tot_Electronic, rMP_Tot_Nuclear, rms, Sites, TP, TP_Q
character(len=80) :: Banner_Line(2)
logical(kind=iwp) :: Center_OK, Check_Bond, get_BasisType
real(kind=wp), allocatable :: Scratch_New(:), Scratch_Org(:)
real(kind=wp), external :: DDot_

Sites = Zero
LI_TOT = Zero
LI_MIN = huge(LI_MIN)
!                                                                      *
!***********************************************************************
!                                                                      *
Center_OK = .true.
call Sphere(lMax)
!                                                                      *
!***********************************************************************
!                                                                      *
! Make copy of rMP so we keep a copy of the
! electronic contributions to the multipoles.

i_Dim = nij*nElem
call dCopy_(i_Dim,rMP,1,xrMP,1)
!                                                                      *
!***********************************************************************
!                                                                      *
! Put info to check file.

!call RecPrt('Print_Local: rMP',' ',rMP,nij,nElem)
do iElem=1,nElem
  call Add_Info('LoProp MP',rMP(1,iElem),nij,4)
end do
if (.not. NoField) then
  do iPol=1,6
    call dcopy_(nij,Polar(iPol,1),6,Temp,1)
    call Add_Info('LoProp a',Temp,nAtoms*(nAtoms+1)/2,4)
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
call Set_Binom()
!                                                                      *
!***********************************************************************
!                                                                      *
! Set up array for nuclear contributions.
!
call dCopy_(nij*nElem,[Zero],0,xnrMP,1)
ij = 0
do iAtom=1,nAtoms
  do jAtom=1,iAtom
    ij = ij+1
    if (iAtom == jAtom) then
      xnrMP(ij,1) = Q_nuc(iAtom)
      CX(1) = Coor(1,iAtom)
      CX(2) = Coor(2,iAtom)
      CX(3) = Coor(3,iAtom)
      call ReExpand(xnrMP,nij,nElem,CX,EC(1,ij),ij,lMax)
    end if
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
! Loop over all domains

write(u6,*)
Banner_Line(1) = 'The Localized properties'
call Banner(Banner_Line(1),1,80)
write(u6,*)

nReal_Centers = 0
ij = 0
do iAtom=1,nAtoms
  do jAtom=1,iAtom
    ij = ij+1
    if (iAtom == jAtom) then
      !--------Atomic domain
      Charge_center = Q_nuc(iAtom)
      CX(1) = Coor(1,iAtom)
      CX(2) = Coor(2,iAtom)
      CX(3) = Coor(3,iAtom)
      Center_OK = .true.
    else
      !--------Bond domain
      Charge_center = Zero
      CX(1) = (Coor(1,iAtom)+Coor(1,jAtom))*Half
      CX(2) = (Coor(2,iAtom)+Coor(2,jAtom))*Half
      CX(3) = (Coor(3,iAtom)+Coor(3,jAtom))*Half
      Center_OK = Check_Bond(Coor(1,iAtom),Coor(1,jAtom),iANr(iAtom),iANr(jAtom),Bond_Threshold)
    end if
    call dcopy_(3,EC(1,ij),1,A,1)

    !------ Print out values for each domain.

    if (Center_OK) then
      nReal_Centers = nReal_Centers+1
      write(u6,*)
      write(u6,*)
      if (iAtom == jAtom) then
        write(u6,*) '===================='
        write(u6,'(A,A)') ' ATOMIC DOMAIN: ',Lbl_Center(iAtom)
        write(u6,*) '===================='
      else
        write(u6,*) '========================='
        write(u6,'(A,A,A,A)') ' BOND DOMAIN: ',Lbl_Center(iAtom),',',Lbl_Center(jAtom)
        write(u6,*) '========================='
      end if
      write(u6,'(A,3F12.8,A)') ' Domain center:  :',CX,' / bohr'
      write(u6,'(A,3F12.8,A)') ' Expansion center:',A,' / bohr'
      write(u6,'(A, F12.8  )') ' Total charge    :',Charge_center+rMP(ij,1)
      !*
      if (LIonize) then
        if (.not. NoField) then
          TP = (Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
          TP_Q = Charge_center+rMP(ij,1)
          Sites = Sites+One
          if (TP > Zero) then
            LI = 0.623_wp*(TP_Q+1.742_wp)*(One/TP)**(One/Three)
            LA = LI*(TP_Q+0.258_wp)/(TP_Q+1.742_wp)
            write(u6,*)
            write(u6,'(A)') 'Local Ionization Energy'
            write(u6,'(10F12.8)') LI
            write(u6,*)
            if (LI < LI_MIN) LI_MIN = LI
            LI_TOT = LI_TOT+LI
            write(u6,'(A)') 'Local Electron Affinity'
            write(u6,'(10F12.8)') LA
            write(u6,'(A)') 'Local Electronegativity'
            write(u6,'(10F12.8)') (LI+LA)/Two
            write(u6,'(A)') 'Local Chemical Hardness'
            write(u6,'(10F12.8)') (LI-LA)/Two

          else
            write(u6,'(A)') 'Negative isotropic polarizability!'
            write(u6,'(A)') 'The local ionization energy will'
            write(u6,'(A)') 'not be computed for this site'
          end if
        end if
      end if
      !*
      write(u6,*)
      write(u6,'(A)') ' Electronic multipole moments:'
      iStrt = 1

      do l=0,lMax
        iEnd = iStrt+(l+1)*(l+2)/2-1
        if (l == 0) write(u6,'(A)') 'Electronic Charge'
        if (l == 1) write(u6,'(A)') 'Electronic Dipole'
        if (l == 2) write(u6,'(A)') 'Electronic Quadrupole'
        if (l == 3) write(u6,'(A)') 'Electronic Octupole'
        if (l == 4) write(u6,'(A)') 'Electronic Hexadecapole'
        if (l >= 5) write(u6,'(A,I2)') 'Electronic multipole moment with l = ',l
        write(u6,'(10F12.8)') (rMP(ij,iElem),iElem=iStrt,iEnd)
        !****************************************************
        if ((iAtom == jAtom) .and. (l >= 1)) then
          write(u6,'(A)') '... with nuclear contribution'
          write(u6,'(10F12.8)') (rMP(ij,iElem)+xnrMP(ij,iElem),iElem=iStrt,iEnd)
        end if
        !****************************************************
        iStrt = iEnd+1
        write(u6,*)
      end do
      if (lMax >= 1) then
        Dip_Tot = sqrt(rMP(ij,2)**2+rMP(ij,3)**2+rMP(ij,4)**2)
        write(u6,*)
        write(u6,'(A,F12.8)') 'Dipole magnitude:',Dip_Tot
        write(u6,*)
      end if
      if (.not. NoField) then
        write(u6,*)
        call TriPrt('Symmetrized Local Polarizability Tensor',' ',Polar(1,ij),3)
        Pol_Tot = (Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
        write(u6,*)
        write(u6,'(A,F12.8)') 'Isotropic Polarizability:',Pol_Tot
        write(u6,*)
      end if
    end if ! Center_OK

    !debug
    !call mma_allocate(EVec,9,label='EVec')
    !call mma_allocate(EVal,6,label='EVal')
    !call unitmat(EVec,3)
    !call dcopy_(6,Polar(1,ij),1,EVal,1)
    !call Jacob(EVal,EVec,3,3)
    !call TriPrt('EVal',' ',EVal,3)
    !call RecPrt('EVec',' ',EVec,3,3)
    !call mma_deallocate(EVal)
    !call mma_deallocate(EVec)
    !debug
    ! Reexpand all multipole moments to the center of mass

    call ReExpand(rMP,nij,nElem,A,C_o_C,ij,lMax)
    !?? call ReExpand(xrMP,nij,nElem,A,C_o_C,ij,lMax)
    call ReExpand(xnrMP,nij,nElem,EC(1,ij),C_o_C,ij,lMax)

  end do
end do
! Write the charge capacitances for the bonds
if (.not. NoField) then
  write(u6,*) '=== Charge capacitance for bonds ==='
  ij = 0
  do jAtom=1,nAtoms
    do iAtom=jAtom+1,nAtoms
      ij = iAtom*(iAtom-1)/2+jAtom
      if (iAtom == jAtom) then
      else
        CX(1) = (Coor(1,iAtom)-Coor(1,jAtom))
        CX(2) = (Coor(2,iAtom)-Coor(2,jAtom))
        CX(3) = (Coor(3,iAtom)-Coor(3,jAtom))
        CRX = ChPolBB(1,ij)*CX(1)+ChPolBB(2,ij)*CX(2)+ChPolBB(4,ij)*CX(3)

        CRY = ChPolBB(2,ij)*CX(1)+ChPolBB(3,ij)*CX(2)+ChPolBB(5,ij)*CX(3)

        CRZ = ChPolBB(4,ij)*CX(1)+ChPolBB(5,ij)*CX(2)+ChPolBB(6,ij)*CX(3)

        CRN = sqrt(CRX*CRX+CRY*CRY+CRZ*CRZ)

        DR = sqrt(CX(1)*CX(1)+CX(2)*CX(2)+CX(3)*CX(3))
        if (0.001_wp < (CRN/(DR*DR*DR))) then
          write(u6,' (A, A, F12.8)') Lbl_Center(iAtom),Lbl_Center(jAtom),CRN/(DR*DR*DR)
        end if
      end if
    end do
  end do
  write(u6,*) '=== =========================== ==='
  write(u6,*)
end if
!*** Local ionizations are implemented by A. Holt
if (LIonize) then
  if (.not. NoField) then
    write(u6,'(A)') '===   Local Ionization energy   ==='
    write(u6,'(A)') 'The local ionization energies are'
    write(u6,'(A)') 'computed using the expression for atoms'
    write(u6,'(A)') 'found in: J.Phys.Chem. 1996, 100,4828'
    write(u6,'(A)')
    write(u6,'(A)') 'Average local ionization energy'
    write(u6,'(10F12.8)') LI_TOT/Sites*auToeV
    write(u6,'(A)') 'Lowest local ionization energy (eV)'
    write(u6,'(10F12.8)') LI_MIN*auToeV
    write(u6,'(A)') 'HOMO energy (absolute value, eV)'
    write(u6,'(10F12.8)') abs(Ene_Occ(nOcOb))*auToeV
    write(u6,'(A)') '=== =========================== ==='
    write(u6,'(A)')
  end if
end if

! Write out the molecular properties

write(u6,*)
write(u6,*)
write(u6,*)
Banner_Line(1) = 'The Molecular Multipole Moments'
call Banner(Banner_Line(1),1,80)
write(u6,'(A,3F12.8,A)') ' Expansion center:',C_o_C,' / bohr'
write(u6,*)
write(u6,*)
iElem = 0
do l=0,lMax
  write(u6,*)
  write(u6,'(A,I1)') 'l=',l
  write(u6,*)
  write(u6,'(A)') 'xyz    Nuclear        Electronic     Molecular   '
  write(u6,*)
  do ix=l,0,-1
    do iy=l-ix,0,-1
      iz = l-ix-iy
      iElem = iElem+1
      rMP_Tot_Electronic = DDot_(nij,[One],0,rMP(1,iElem),1)
      rMP_Tot_Nuclear = rMPq(iElem)
      rMP_Tot = rMP_Tot_Nuclear+rMP_Tot_Electronic
      write(u6,'(3I1,3F16.8)') ix,iy,iz,rMP_Tot_Nuclear,rMP_Tot_Electronic,rMP_Tot
    end do
  end do
end do
if (.not. NoField) then
  do iPol=1,6
    Polar_M(iPol) = DDot_(nij,[One],0,Polar(iPol,1),6)
  end do
  call TriPrt('Molecular Polarizability Tensor',' ',Polar_M,3)
  !debug
  !call mma_allocate(EVec,9,label='EVec')
  !call mma_allocate(EVal,6,label='EVal')
  !call dcopy_(9,Zero,0,EVec,1)
  !call dcopy_(3,One,0,EVec,4)
  !call dcopy_(6,Polar_M,1,EVal,1)
  !call Jacob(EVal,EVec,3,3)
  !call TriPrt('EVal',' ',EVal,3)
  !call RecPrt('EVec',' ',EVec,3,3)
  !call mma_deallocate(EVal)
  !call mma_deallocate(EVec)
  !debug
end if

!------ Generate mpprop file

call Print_MPPROP(rMP,xrMP,xnrMP,nij,nElem,lMax,EC,Polar,Lbl_Center,nAtoms,iANr,NoField,C_o_C,Coor,nOcOb,Energy_Without_FFPT, &
                  Ene_Occ,MpProp_Level,Bond_Threshold,nReal_Centers)

call mma_allocate(Scratch_New,nij*(2*lMax+1),label='Scratch_New')
call mma_allocate(Scratch_Org,nij*(2*lMax+1),label='Scratch_Org')

!------ Compare the molecular multipole moments to the ones arrising from truncation
!       xrMP contains the original multipole moments (electr. + nuclear)
!       xxrMP contains the approximate multipole moments.

call dCopy_(i_Dim,rMP,1,xrMP,1)
call daxpy_(i_Dim,One,xnrMP,1,xrMP,1)
call dCopy_(i_Dim,xrMP,1,xxrMP,1)

iPrint = 1
do l=lMax-1,0,-1
  call CutOff_Error(l,lMax,xrMP,xxrMP,nij,EC,C_o_C,nElem,Scratch_New,Scratch_Org,nAtoms,iPrint,rms)
end do

call mma_deallocate(Scratch_New)
call mma_deallocate(Scratch_Org)
call Sphere_Free()

!-- Print warning if non-ANO basis sets have been used

if (.not. get_BasisType('ANO')) then
  write(u6,*)
  write(u6,*) 'WARNING: The calculation were performed with at '
  write(u6,*) '         least one non-ANO basis! The results'
  write(u6,*) '         might therefore be erroneous.'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Print_Local
