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

subroutine Prop(Short,qplab,cen1,cen2,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,lpole,labs,tmat,temp,ifallorb)
!***********************************************************************
!                                                                      *
!     purpose: preprocessing of tables for different tensor            *
!              properties                                              *
!                                                                      *
!     Short           logical option for either Short (Short           *
!                     =.true., the total electronic contribution)      *
!                     long (Short=.false., orbital contributions)      *
!                     output                                           *
!     oplab           operator label as defined in SEWARD              *
!     cen1(1:3)       coordinates of centre no.1 (operator)            *
!     cen2(1:3)       coordinates of centre no.2 (gauge)               *
!     nIrrep          the number of irreducible representations        *
!     nBas            the number of functions in each representa-      *
!     (0:nIrrep-1)    tion                                             *
!     nTot            the total number of elements supplied for        *
!                     each component of the property tensor; equal     *
!                     either to 1 (total electronic and total nuc-     *
!                     lear contributions) or to the dimension of       *
!                     the basis set.                                   *
!     Occ(1:nTot)     occupation numbers for all eigenvectors,         *
!                     a dummy for Short outputs                        *
!     ThrSV           threshold for occupation numbers; If             *
!                     Occ(i).le.ThrSV the contribution will not        *
!                     be printed                                       *
!     PrEl(1:nTot,    matrix elements for all components 1,2,...,      *
!          1:maxlab)  maxlab, nTot entries for each component          *
!                     maxlab=(lpole+1)*(lpole+2)/2                     *
!     PrNu(1:maxlab)  nuclear contributions for each component         *
!     labs(1:maxlab)  labels for each component                        *
!     temp(1:maxlab)  auxiliary storage area                           *
!     ifallorb        logical option for whether the property of       *
!                     all orbitals are printed (and not weighted by    *
!                     occupation number)in property calculation when   *
!                     short=.false. (S.S.Dong, 2018)                   *
!                                                                      *
!     lpole is the value for l in l-pole moments                       *
!     allocate integer storage area of the size appropriate            *
!     for the actual lpole value: size as below                        *
!                                                                      *
! 2000 Dept. of Chem. Phys., Univ. of Lund, Sweden                     *
! Modified by S.S.Dong, 2018, Univ. of Minnesota                       *
! - Enable properties to be printed for all orbitals                   *
! (including virtuals) and not weighted by occupation numbers          *
!***********************************************************************

use Constants, only: Zero, One, Angstrom, Debye
use Definitions, only: wp, iwp, u6

implicit none
logical(kind=iwp), intent(in) :: Short, ifallorb
character(len=8), intent(in) :: qplab
integer(kind=iwp), intent(in) :: nIrrep, nBas(0:nIrrep-1), nTot, lpole
real(kind=wp), intent(in) :: cen1(3), cen2(3), Occ(nTot), ThrSV, PrEl(nTot,(lpole+1)*(lpole+2)/2), PrNu((lpole+1)*(lpole+2)/2)
character(len=16), intent(out) :: labs((lpole+1)*(lpole+2)/2)
real(kind=wp), intent(out) :: tmat((lpole+1)*(lpole+2)/2,(lpole+1)*(lpole+2)/2), temp((lpole+1)*(lpole+2)/2)
integer(kind=iwp) :: i, icen, icen1, icen2, ilab, inp, iOcc, ip_, iPL, ipPrTot, iSt, iTol, iTol_E0, iTol_E1, ix, ixx, iy, iyy, iz, &
                     izz, j, jMax, maxlab
real(kind=wp) :: Fact, Molecular_Charge = Zero, PrElAug(nTot,(lpole+1)*(lpole+2)/2+1), PrNuAug((lpole+1)*(lpole+2)/2+1), sig, tmp, &
                 X_Coor, Y_Coor, Z_Coor
logical(kind=iwp) :: StoreInfo
integer(kind=iwp), parameter :: lmax = 16
character(len=lmax) :: lab, labsAug((lpole+1)*(lpole+2)/2+1)
character(len=80) :: Line
character(len=8) :: oplab
character(len=5) :: lab5
character(len=4) :: lab4
character(len=3) :: lab3
integer(kind=iwp), external :: Cho_X_GetTol, iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
#include "hfc_logical.fh"
#include "WrkSpc.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
!                                                                      *
!***********************************************************************
!                                                                      *
! Some properties can be calculated at many, many centers,
! and writing all values through Add_Info affects performance,
! those will be deactivated and the summation is printed by the
! calling function.
! This only affects the call to Add_Info, otherwise the properties
! are calculated and printed normally
!
StoreInfo = .true.
!                                                                      *
!***********************************************************************
!                                                                      *

! decipher the operator label:oplab
oplab = qplab
call UpCase(oplab)
lab4 = oplab(1:4)
if (lab4 == 'MLTP') then
  if (lpole < 0) then
    return
  else
    Go To 200
  end if
end if
if (lab4(1:2) == 'EF') Go To 300
if (lab4 == 'DMS ') Go To 500
if (lab4(1:3) == 'PAM') Go To 600
if (lab4(1:3) == 'CNT') Go To 700
if (lab4 == 'VELO') return
if (lab4 == 'ANGM') return

! invalid label supplied

write(u6,'(//1x,a,a8,a//)') ' The label: *',oplab,'* is not a valid label ... stop'
call Quit_OnUserError()
!                                                                      *
!***********************************************************************
!                                                                      *
! Multipole moment section ... generate labels for printing

200 continue
if (lPole > lMax) then
  write(u6,*) 'Prop: lPole.gt.lMax'
  write(u6,*) 'lPole=',lPole
  write(u6,*) 'Increase lMax and recompile!'
  call Abend()
end if
ilab = 0
do ix=lpole,0,-1
  do iy=lpole-ix,0,-1
    iz = lpole-ix-iy
    ilab = ilab+1
    lab = '                '
    do izz=lmax,lmax-iz+1,-1
      lab(izz:izz) = 'Z'
    end do
    do iyy=lmax-iz,lmax-iz-iy+1,-1
      lab(iyy:iyy) = 'Y'
    end do
    do ixx=lmax-iz-iy,lmax-lpole+1,-1
      lab(ixx:ixx) = 'X'
    end do
    labs(ilab) = lab
  end do
end do

maxlab = ilab
call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)

! Print cartesian moments

!----------------------------------------------------------------------*
if ((iPL == 2) .and. Short) then
!----------------------------------------------------------------------*

  do i=1,MaxLab
    Work(ipPrTot+i-1) = PrNu(i)-PrEl(1,i)
  end do

  ! New style output

  Line = ' '
  if (lPole == 0) then
    Line = 'Charge (e):'
    Fact = One
  else if (lPole == 1) then
    Line = 'Dipole Moment (Debye):'
    Fact = Debye
  else if (lPole == 2) then
    Line = 'Quadrupole Moment (Debye*Ang):'
    Fact = Debye*Angstrom
  else if (lPole == 3) then
    Line = 'Octupole Moment (Debye*Ang**2):'
    Fact = Debye*Angstrom**2
  else if (lPole == 4) then
    Line = 'Hexadecapole Moment (Debye*Ang**3):'
    Fact = Debye*Angstrom**3
  else
    Line = ''
    if (lpole <= 9) then
      write(Line,'(I1)') lPole
      iSt = 2
    else
      write(Line,'(I2)') lPole
      iSt = 3
    end if
    Line(iSt:iSt+26) = 'th-pole Moment (Debye*Ang**'
    iSt = iSt+27
    if (lpole <= 10) then
      write(Line(iSt:iSt),'(I1)') lpole-1
      iSt = iSt+1
    else
      write(Line(iSt:iSt+1),'(I2)') lpole-1
      iSt = iSt+2
    end if
    Line(iSt:iSt+2) = '):'
    Fact = Debye*Angstrom**(lPole-1)
  end if
  write(u6,'(6X,A)') trim(Line)
  if (lpole > 0) then
    write(u6,'(6X,A,3F10.4)') 'Origin of the operator (Ang)=',(cen1(i)*Angstrom,i=1,3)
  end if
  if (lPole == 0) then
    tmp = Work(ipPrTot)
    write(u6,'(6X,A,A,F10.4)') labs(1),'=',Work(ipPrTot)*Fact
    Molecular_Charge = Work(ipPrTot)*Fact
  else if (lPole == 1) then
    tmp = sqrt(Work(ipPrTot)**2+Work(ipPrTot+1)**2+Work(ipPrTot+2)**2)
    write(u6,'(4X,4(A,A,ES12.4))') labs(1),'=',Work(ipPrTot)*Fact,labs(2),'=',Work(ipPrTot+1)*Fact,labs(3),'=', &
                                   Work(ipPrTot+2)*Fact,'           Total','=',tmp*Fact
    if (abs(Molecular_Charge) > 0.9_wp) then
      write(u6,'(6X,A)') 'Center of Charge (Ang)'
      X_Coor = Angstrom*(Work(ipPrTot)/Molecular_Charge)
      Y_Coor = Angstrom*(Work(ipPrTot+1)/Molecular_Charge)
      Z_Coor = Angstrom*(Work(ipPrTot+2)/Molecular_Charge)
      write(u6,'(6X,3(A,A,F14.8))') labs(1),'=',X_Coor,labs(2),'=',Y_Coor,labs(3),'=',Z_Coor
      Molecular_Charge = Zero
    end if
    call Put_DArray('Dipole moment',Work(ipPrTot),3)
    !call peek_iScalar('xml opened',isopen)
    !if (isopen == 1) then
    call xml_dDump('dipole','Dipole moment','Debye',1,Work(ipPrTot),3,1)
    !end if
  else if (lPole >= 2) then
    ip_ = ipPrTot
    tmp = Zero
    do i=0,Maxlab-1
      tmp = max(tmp,abs(Work(ipPrTot+i)))
    end do
    ip_ = ipPrTot
    do i=1,maxlab,4
      jMax = min(maxlab-i,3)
      write(u6,'(4X,4(A,A,ES12.4))') (labs(i+j),'=',Work(ip_+j)*Fact,j=0,jMax)
      ip_ = ip_+4
    end do
  end if
  if ((lpole >= 2) .and. (lpole <= 4)) then

    ! Transform cartesian moments to multipole moments

    inp = 0
    call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
    inp = 1
    call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
    do i=1,MaxLab
      Work(ipPrTot+i-1) = PrNu(i)-PrEl(1,i)
    end do

    if (lPole >= 3) then
      write(u6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang**',lPole-1,')'
    else
      write(u6,'(6X,A,I1,A)') 'In traceless form (Debye*Ang)'
    end if
    ip_ = ipPrTot
    do i=1,maxlab,4
      jMax = min(maxlab-i,3)
      write(u6,'(4X,4(A,A,ES12.4))') (labs(i+j),'=',Work(ip_+j)*Fact,j=0,jMax)
      ip_ = ip_+4
    end do

  end if

!----------------------------------------------------------------------*
else if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
!----------------------------------------------------------------------*

  ! Old style output

  if (lpole == 0) lab5 = ' 0-th'
  if (lpole == 1) lab5 = ' 1-st'
  if (lpole == 2) lab5 = ' 2-nd'
  if (lpole == 3) lab5 = ' 3-rd'
  if (lpole > 3) then
    lab3 = '-th'
    write(lab5,'(i2,a3)') lpole,lab3
  end if
  write(u6,'(//6x,a5,a,3(f12.8,a))') lab5,' cartesian moments: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
  write(u6,'(6x,76(''-''))')
  sig = -One
  call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,0,ifallorb)
  if (lpole == 1) then
    write(u6,'(6x,76(''-''))')
    write(u6,'(6x,a,3f16.8,3x,a)') 'Total             ',(Work(ipPrTot+j)*Debye,j=0,2),'Debye'
    call Put_DArray('Dipole moment',Work(ipPrTot),3)
  end if

  if ((lpole >= 2) .and. (lpole <= 4)) then
    write(u6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,'-pole moment: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,76(''-''))')

    ! Transform cartesian moments to multipole moments

    inp = 0
    call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
    inp = 1
    call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)

    ! Print l-pole cartesian moments

    sig = -One
    call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,0,ifallorb)
  end if

!----------------------------------------------------------------------*
else
!----------------------------------------------------------------------*
  do i=1,MaxLab
    Work(ipPrTot+i-1) = PrNu(i)-PrEl(1,i)
  end do
  if (lpole == 1) call Put_DArray('Dipole moment',Work(ipPrTot),3)
  if ((lpole >= 2) .and. (lpole <= 4)) then

    ! Transform cartesian moments to multipole moments

    inp = 0
    call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
    inp = 1
    call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)
    do i=1,MaxLab
      Work(ipPrTot+i-1) = PrNu(i)-PrEl(1,i)
    end do
  end if
!----------------------------------------------------------------------*
end if ! iPL
!----------------------------------------------------------------------*
! Prop is also called in other programs where MAG_X2C could
! be uninitialized, if a test is required in such case
! please initialize MAG_X2C to false in related programs
if (MAG_X2C) StoreInfo = .false.
Go To 999
!                                                                      *
!***********************************************************************
!                                                                      *
! electric field section

300 continue
read(oplab,'(a3,i5)') lab3,icen

if (lab3 == 'EF0') then
  MaxLab = 1
  ! set labels
  labs(1) = '                '
  ! Print cartesian components of the electric field gradient
  ! tensor at the given centre

  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
    write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric potential:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,72(''-''))')
  else
    if (icen == 1) then
      write(u6,'(//6X,A)') 'Electric potential:'
    end if
  end if
  sig = -One
end if
if (lab3 == 'EF1') then
  MaxLab = 3
  ! set labels
  labs(1) = '               X'
  labs(2) = '               Y'
  labs(3) = '               Z'
  ! Print cartesian components of the electric field vector
  ! at the given centre

  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
    write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric field:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,72(''-''))')
  else
    if (icen == 1) then
      write(u6,'(//6X,A)') 'Electric field:'
      write(u6,'(5X,6A16)') (labs(i),i=1,MaxLab)
    end if
  end if
  sig = +One
end if
if (lab3 == 'EF2') then
  ! Seven components because we add r*r to the "normal" 6
  MaxLab = 7
  ! set labels
  labs(1) = '  (2*XX-YY-ZZ)/2'
  labs(2) = '          1.5*XY'
  labs(3) = '          1.5*XZ'
  labs(4) = '  (2*YY-ZZ-XX)/2'
  labs(5) = '          1.5*YZ'
  labs(6) = '  (2*ZZ-XX-YY)/2'
  ! Actually, r*r is already stored as the 6th component
  ! we now move it to the 7th using the "augmented" arrays
  do j=1,5
    do iOcc=1,nTot
      PrElAug(iOcc,j) = PrEl(iOcc,j)
    end do
    PrNuAug(j) = PrNu(j)
    labsAug(j) = labs(j)
  end do
  do iOcc=1,nTot
    ! Generate the actual value of the 6th component
    PrElAug(iOcc,6) = -PrEl(iOcc,1)-PrEl(iOcc,4)
    PrElAug(iOcc,7) = PrEl(iOcc,6)
  end do
  PrNuAug(6) = PrNu(6)
  PrNuAug(7) = Zero
  labsAug(6) = labs(6)
  labsAug(7) = '     RR=XX+YY+ZZ'

  ! Print cartesian components of the electric field gradient
  ! tensor at the given centre

  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
    write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Electric field gradient:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
    write(u6,'(6x,78(''-''))')
  else
    if (icen == 1) then
      write(u6,'(//6X,A)') 'Electric field gradient:'
      write(u6,'(5X,6A16)') (labsAug(i),i=1,MaxLab)
    end if
  end if
  sig = +One
end if
call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)

! Print the values using "augmented" arrays if needed

if (Maxlab == 7) then
  call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrElAug,PrNuAug,maxlab,labsAug,Work(ipPrTot),iPL,icen,ifallorb)
  MaxLab = 6 ! Reset so call to Add_Info is correct!

else
  call Prout(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,icen,ifallorb)
end if

if (lab3 == 'EF2') then
  Tmp = Work(ipPrTot+2)
  Work(ipPrTot+2) = Work(ipPrTot+3)
  Work(ipPrTot+3) = Tmp
  if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) call Print_EigenValues(Work(ipPrTot),3)
  Tmp = Work(ipPrTot+2)
  Work(ipPrTot+2) = Work(ipPrTot+3)
  Work(ipPrTot+3) = Tmp
end if

! do not write the different electric field components through Add_Info
StoreInfo = .false.
Go To 999
!                                                                      *
!***********************************************************************
!                                                                      *
! diamagnetic shielding section

500 continue
read(oplab,'(a4,i2,i2)') lab4,icen2,icen1

maxlab = 9
call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)

! set labels
labs(1) = '   (YO*YG+ZO*ZG)'
labs(2) = '         (XO*YG)'
labs(3) = '         (XO*ZG)'
labs(4) = '         (YO*XG)'
labs(5) = '   (XO*XG+ZO*ZG)'
labs(6) = '         (YO*ZG)'
labs(7) = '         (ZO*XG)'
labs(8) = '         (ZO*YG)'
labs(9) = '   (XO*XG+YO*YG)'

! Print cartesian components of the diamagnetic shielding
! tensor at the given centre (O=cen1) with the gauge origin
! at centre G=cen2

write(u6,'(//6x,a,i3,1x,a,3(f12.8,a)/1x,a,3(f12.8,a))') ' Diamagnetic shielding:   centre no.',icen1,'(',cen1(1),',',cen1(2),',', &
                                                        cen1(3),')','                       gauge origin at (',cen2(1),',', &
                                                        cen2(2),',',cen2(3),')'
write(u6,'(6x,78(''-''))')
sig = One
call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,0,ifallorb)
!                                                                      *
!***********************************************************************
!                                                                      *
! PAM section ... generate labels for printing

600 continue
ilab = 0
do ix=lpole,0,-1
  do iy=lpole-ix,0,-1
    iz = lpole-ix-iy
    ilab = ilab+1
    lab = '                '
    do izz=lmax,lmax-iz+1,-1
      lab(izz:izz) = 'Z'
    end do
    do iyy=lmax-iz,lmax-iz-iy+1,-1
      lab(iyy:iyy) = 'Y'
    end do
    do ixx=lmax-iz-iy,lmax-lpole+1,-1
      lab(ixx:ixx) = 'X'
    end do
    labs(ilab) = lab
  end do
end do

maxlab = ilab
call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)

! Print cartesian moments

if (lpole == 0) lab5 = ' 0-th'
if (lpole == 1) lab5 = ' 1-st'
if (lpole == 2) lab5 = ' 2-nd'
if (lpole == 3) lab5 = ' 3-rd'
if (lpole > 3) then
  lab3 = '-th'
  write(lab5,'(i2,a3)') lpole,lab3
end if
write(u6,'(//6x,a5,a,3(f12.8,a))') lab5,' cartesian moments: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
write(u6,'(6x,76(''-''))')
sig = One
call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,0,ifallorb)

if ((lpole >= 2) .and. (lpole <= 4)) then
  write(u6,'(//6x,a,i2,a,3(f12.8,a))') 'Cartesian ',lpole,'-pole moment: origin at (',cen1(1),',',cen1(2),',',cen1(3),')'
  write(u6,'(6x,76(''-''))')

  inp = 0
  call Tmltpl(inp,lpole,maxlab,labs,1,PrNu,tmat,temp)
  inp = 1
  call Tmltpl(inp,lpole,maxlab,labs,nTot,PrEl,tmat,temp)

  ! Print 0-pole PAM integrals (sig=1 in opposite multipole moments)

  sig = One
  call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,0,ifallorb)
end if
!go to 999
!                                                                      *
!***********************************************************************
!                                                                      *
! Contact term section

700 continue
read(oplab,'(a3,i5)') lab3,icen

maxlab = 1
call GetMem('PrTot','Allo','Real',ipPrTot,maxlab)

! set labels
labs(1) = '      Delta(R-C)'

if ((iPL >= 3) .or. ((.not. Short) .and. (iPL == 2))) then
  write(u6,'(//6x,a,i5,1x,a,3(f12.8,a))') ' Contact term:  centre no.',icen,'(',cen1(1),',',cen1(2),',',cen1(3),')'
  write(u6,'(6x,78(''-''))')
else
  if (icen == 1) then
    write(u6,'(//6X,A)') 'Contact term:'
    write(u6,'(5X,6A16)') (labs(i),i=1,MaxLab)
  end if
end if
sig = +One
call PrOut(Short,sig,nIrrep,nBas,nTot,Occ,ThrSV,PrEl,PrNu,maxlab,labs,Work(ipPrTot),iPL,icen,ifallorb)

! do not write the contact term through Add_Info
StoreInfo = .false.
!                                                                      *
!***********************************************************************
!                                                                      *
999 continue
iTol = 5
iTol_E0 = 8
iTol_E1 = Cho_X_GetTol(iTol_E0)
iTol = int(real(iTol*iTol_E1,kind=wp)/real(iTol_E0,kind=wp))
if (StoreInfo) call Add_Info(OpLab,Work(ipPrTot),maxlab,iTol)
call GetMem('PrTot','Free','Real',ipPrTot,maxlab)

return
#ifdef _WARNING_WORKAROUND_
if (.false.) call Unused_integer(icen2)
#endif

end subroutine Prop
