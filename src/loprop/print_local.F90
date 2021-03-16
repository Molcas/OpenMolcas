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
                       nOcOb,Energy_Without_FFPT,ip_Ene_Occ,MpProp_Level,Bond_Threshold,XHole,XHoleLoc,D2,ChPol,ChPolBB,LIonize)

use Real_Spherical

implicit real*8(a-h,o-z)
#include "Molcas.fh"
#include "real.fh"
#include "WrkSpc.fh"
#include "stdalloc.fh"
#include "status.fh"
real*8 rMP(nij,nElem), Coor(3,nAtoms), A(3), C_o_C(3), Q_Nuc(nAtoms), rMPq(nElem), EC(3,nij), Polar(6,nij), Polar_M(6), Temp(nij), &
       CX(3), xrMP(nij,nElem), xxrMP(nij,nElem), xnrMP(nij,nElem), XHoleLoc(nij), ChPol(6,nij), ChPolBB(6,nij)
integer iANr(nAtoms)
character*(LENIN) Lbl_Center(nAtoms)
character*80 Banner_Line(2)
logical NoField, Center_OK, Check_Bond, get_BasisType, XHole
logical LIonize

real*8 Sites
real*8 TP
real*8 TP_Q
real*8 LI, LA
real*8 LI_TOT
real*8 LI_MIN

Sites = 0.0
LI_TOT = 0.0
LI_MIN = 99.0
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

iDim = nij*nElem
call dCopy_(iDim,rMP,1,xrMP,1)
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

write(6,*)
Banner_Line(1) = 'The Localized properties'
call Banner(Banner_Line(1),1,80)
write(6,*)

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
      write(6,*)
      write(6,*)
      if (iAtom == jAtom) then
        write(6,*) '===================='
        write(6,'(A,A)') ' ATOMIC DOMAIN: ',Lbl_Center(iAtom)
        write(6,*) '===================='
      else
        write(6,*) '========================='
        write(6,'(A,A,A,A)') ' BOND DOMAIN: ',Lbl_Center(iAtom),',',Lbl_Center(jAtom)
        write(6,*) '========================='
      end if
      write(6,'(A,3F12.8,A)') ' Domain center:  :',CX,' / bohr'
      write(6,'(A,3F12.8,A)') ' Expansion center:',A,' / bohr'
      write(6,'(A, F12.8  )') ' Total charge    :',Charge_center+rMP(ij,1)
      !*
      if (LIonize) then
        if (.not. NoField) then
          TP = (Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
          TP_Q = Charge_center+rMP(ij,1)
          Sites = Sites+1.0
          if (TP > 0.0) then
            LI = 0.623*(TP_Q+1.742)*(1/TP)**(1.0/3.0)
            LA = LI*(TP_Q+0.258)/(TP_Q+1.742)
            write(6,*)
            write(6,'(A)') 'Local Ionization Energy'
            write(6,'(10F12.8)') LI
            write(6,*)
            if (LI < LI_MIN) LI_MIN = LI
            LI_TOT = LI_TOT+LI
            write(6,'(A)') 'Local Electron Affinity'
            write(6,'(10F12.8)') LA
            write(6,'(A)') 'Local Electronegativity'
            write(6,'(10F12.8)') (LI+LA)/2.0
            write(6,'(A)') 'Local Chemical Hardness'
            write(6,'(10F12.8)') (LI-LA)/2.0

          else
            write(6,'(A)') 'Negative isotropic polarizability!'
            write(6,'(A)') 'The local ionization energy will'
            write(6,'(A)') 'not be computed for this site'
          end if
        end if
      end if
      !*
      write(6,*)
      write(6,'(A)') ' Electronic multipole moments:'
      iStrt = 1

      do l=0,lMax
        iEnd = iStrt+(l+1)*(l+2)/2-1
        if (l == 0) write(6,'(A)') 'Electronic Charge'
        if (l == 1) write(6,'(A)') 'Electronic Dipole'
        if (l == 2) write(6,'(A)') 'Electronic Quadrupole'
        if (l == 3) write(6,'(A)') 'Electronic Octupole'
        if (l == 4) write(6,'(A)') 'Electronic Hexadecapole'
        if (l >= 5) write(6,'(A,I2)') 'Electronic multipole moment with l = ',l
        write(6,'(10F12.8)') (rMP(ij,iElem),iElem=iStrt,iEnd)
        !****************************************************
        if (iAtom == jAtom .and. l >= 1) then
          write(6,'(A)') '... with nuclear contribution'
          write(6,'(10F12.8)') (rMP(ij,iElem)+xnrMP(ij,iElem),iElem=iStrt,iEnd)
        end if
        !****************************************************
        iStrt = iEnd+1
        write(6,*)
      end do
      if (lMax >= 1) then
        Dip_Tot = sqrt(rMP(ij,2)**2+rMP(ij,3)**2+rMP(ij,4)**2)
        write(6,*)
        write(6,'(A,F12.8)') 'Dipole magnitude:',Dip_Tot
        write(6,*)
      end if
      if (.not. NoField) then
        write(6,*)
        call TriPrt('Symmetrized Local Polarizability Tensor',' ',Polar(1,ij),3)
        Pol_Tot = (Polar(1,ij)+Polar(3,ij)+Polar(6,ij))/Three
        write(6,*)
        write(6,'(A,F12.8)') 'Isotropic Polarizability:',Pol_Tot
        write(6,*)
      end if
      if (XHole) then
        x2Dip = XHoleLoc(ij)
        write(6,'(A,F12.8)') 'Exchange hole second-moment:',x2Dip
        write(6,*)
      end if
    end if ! Center_OK

    !debug
    !call Allocate_Work(ip_EVec,9)
    !call Allocate_Work(ip_EVal,6)
    !call dcopy_(9,[Zero],0,Work(ip_EVec),1)
    !call dcopy_(3,[One],0,Work(ip_EVec),4)
    !call dcopy_(6,Polar(1,ij),1,Work(ip_EVal),1)
    !call Jacob(Work(ip_EVal),Work(ip_EVec),3,3)
    !call TriPrt('EVal',' ',Work(ip_EVal),3)
    !call RecPrt('EVec',' ',Work(ip_EVec),3,3)
    !call Free_Work(ip_EVal)
    !call Free_Work(ip_EVec)
    !debug
    ! Reexpand all multipole moments to the center of mass

    call ReExpand(rMP,nij,nElem,A,C_o_C,ij,lMax)
    !?? call ReExpand(xrMP,nij,nElem,A,C_o_C,ij,lMax)
    call ReExpand(xnrMP,nij,nElem,EC(1,ij),C_o_C,ij,lMax)

  end do
end do
! Write the charge capacitances for the bonds
if (.not. NoField) then
  write(6,*) '=== Charge capacitance for bonds ==='
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
        if (0.001 < (CRN/(DR*DR*DR))) then
          write(6,' (A, A, F12.8)') Lbl_Center(iAtom),Lbl_Center(jAtom),CRN/(DR*DR*DR)
        end if
      end if
    end do
  end do
  write(6,*) '=== =========================== ==='
  write(6,*)
end if
!*** Local ionizations are implemented by A. Holt
if (LIonize) then
  if (.not. NoField) then
    write(6,'(A)') '===   Local Ionization energy   ==='
    write(6,'(A)') 'The local ionization energies are'
    write(6,'(A)') 'computed using the expression for atoms'
    write(6,'(A)') 'found in; J.Phys.Chem. 1996, 100,4828'
    write(6,'(A)')
    write(6,'(A)') 'Average local ionization energy'
    write(6,'(10F12.8)') LI_TOT/Sites*27.2114
    write(6,'(A)') 'Lowest local ionization energy (eV)'
    write(6,'(10F12.8)') LI_MIN*27.2114
    write(6,'(A)') 'HOMO energy (absolute value, eV)'
    write(6,'(10F12.8)') abs(Work(ip_Ene_Occ+nOcOb-1))*27.2114
    write(6,'(A)') '=== =========================== ==='
    write(6,'(A)')
  end if
end if

! Write out the molecular properties

write(6,*)
write(6,*)
write(6,*)
Banner_Line(1) = 'The Molecular Multipole Moments'
call Banner(Banner_Line(1),1,80)
write(6,'(A,3F12.8,A)') ' Expansion center:',C_o_C,' / bohr'
write(6,*)
write(6,*)
iElem = 0
do l=0,lMax
  write(6,*)
  write(6,'(A,I1)') 'l=',l
  write(6,*)
  write(6,'(A)') 'xyz    Nuclear        Electronic     Molecular   '
  write(6,*)
  do ix=l,0,-1
    do iy=l-ix,0,-1
      iz = l-ix-iy
      iElem = iElem+1
      rMP_Tot_Electronic = DDot_(nij,[One],0,rMP(1,iElem),1)
      rMP_Tot_Nuclear = rMPq(iElem)
      rMP_Tot = rMP_Tot_Nuclear+rMP_Tot_Electronic
      write(6,'(3I1,3F16.8)') ix,iy,iz,rMP_Tot_Nuclear,rMP_Tot_Electronic,rMP_Tot
    end do
  end do
end do
if (.not. NoField) then
  do iPol=1,6
    Polar_M(iPol) = DDot_(nij,[One],0,Polar(iPol,1),6)
  end do
  call TriPrt('Molecular Polarizability Tensor',' ',Polar_M,3)
  !debug
  !call Allocate_Work(ip_EVec,9)
  !call Allocate_Work(ip_EVal,6)
  !call dcopy_(9,Zero,0,Work(ip_EVec),1)
  !call dcopy_(3,One,0,Work(ip_EVec),4)
  !call dcopy_(6,Polar_M,1,Work(ip_EVal),1)
  !call Jacob(Work(ip_EVal),Work(ip_EVec),3,3)
  !call TriPrt('EVal',' ',Work(ip_EVal),3)
  !call RecPrt('EVec',' ',Work(ip_EVec),3,3)
  !call Free_Work(ip_EVal)
  !call Free_Work(ip_EVec)
  !debug
end if
if (XHole) then
  write(6,*)
  write(6,'(A,F12.8)') 'Molecular exchange hole second moment:',D2
  write(6,*)
end if

!------ Generate mpprop file

call Print_MPPROP(rMP,xrMP,xnrMP,nij,nElem,lMax,EC,Polar,Lbl_Center,nAtoms,iANr,NoField,C_o_C,Coor,nOcOb,Energy_Without_FFPT, &
                  ip_Ene_Occ,MpProp_Level,Bond_Threshold,nReal_Centers)

call Allocate_Work(iScratch_New,nij*(2*lMax+1))
call Allocate_Work(iScratch_Org,nij*(2*lMax+1))

!------ Compare the molecular multipole moments to the ones arrising from truncation
!       xrMP contains the original multipole moments (electr. + nuclear)
!       xxrMP contains the approximate multipole moments.

call dCopy_(iDim,rMP,1,xrMP,1)
call daxpy_(iDim,One,xnrMP,1,xrMP,1)
call dCopy_(iDim,xrMP,1,xxrMP,1)

iPrint = 1
do l=lMax-1,0,-1
  call CutOff_Error(l,lMax,xrMP,xxrMP,nij,EC,C_o_C,nElem,Work(iScratch_New),Work(iScratch_Org),nAtoms,iPrint,rms)
end do

call Free_Work(iScratch_New)
call Free_Work(iScratch_Org)
call Sphere_Free()

!-- Print warning if non-ANO basis sets have been used

if (.not. get_BasisType('ANO')) then
  write(6,*)
  write(6,*) 'WARNING: The calculation were performed with at '
  write(6,*) '         least one non-ANO basis! The results'
  write(6,*) '         might therefore be erroneous.'
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return
! Avoid unused argument warnings
if (.false.) call Unused_real_array(ChPol)

end subroutine Print_Local
