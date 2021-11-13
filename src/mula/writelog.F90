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
! Copyright (C) 1995, Niclas Forsberg                                  *
!***********************************************************************

!module InOutMod

!  Contains:
!    ReadInp
!    WriteLog
!    WriteFC
!    WriteHeader
!    WrMold
!
!  Uses:
!    Constants
!    InputData
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

!contains

subroutine WriteLog(PotCoef,AtomLbl,AtCoord,Mass,InterVec,stand_dev,max_err,energy,Hess,G,V,harmfreq,qMat,Bond,nBond,xvec,D3,D4, &
                    PED,x_anharm,anharmfreq,max_term,nState,ForceField,NumOfAt,nOsc)
!  Purpose:
!    Write results to log file.
!
!  Input:
!    AtomLbl    : Array of character - contains the labels for the
!                 atoms.
!    AtCoord    : Two dimensional Real*8 array - contains
!                 the cartesian coordinates of the atoms.
!    Mass       : Real*8 array - contains the mass of the
!                 atoms.
!    InterVec   : Integer array - contains the atoms that are used
!                 in the calculations of each internal coordinate.
!    ipow       : Two dimesional integer array - terms of polynomial.
!    PotCoef    : Real*8 two dimensional array - coefficients
!                 of terms given in ipow.
!    stand_dev  : Real*8 variable - standard deviation of fitted
!                 values of polynomial.
!    max_err    : Real*8 variable - maximum error of fitted
!                 values of polynomial.
!    energy     : Real*8 variable - energy in minimum.
!    Hess       : Real*8 array -  contains the force
!                 constants expressed in internal coordinates.
!    G          : Real*8 array.
!    V          : Real*8 array - contains the eigenvectors
!                 of G*F as columns.
!    harmfreq   : Real*8 array - harmonic frequencies.
!    qMat       : Real*8 array - contains the cartesian
!                 displacements of the atoms.
!    Bond       : Integer array - contains atom pairs that are to be
!                 bonded together in a plot.
!    nBond      : Integer - dim of Bond, i.e. 2*(number of bonds).
!    xvec       : Real*8 array - contains the internal coordinates
!                 at equilibrium.
!    D3         : Real*8 array - contains the third derivatives.
!    D4         : Real*8 array - contains the fourth derivatives.
!    PED        : Real*8 three dimensional array - potential
!                 energy distribution.
!    x_anharm   : Real*8 two dimensional array - anharmonicity
!                 constants.
!    anharmfreq : Real*8 array - anharmonic frequencies.
!    max_term   : Integer - highest power of term in fitted polynomial.
!    nState     : Integer - 1 or 2, depending upon which state it is.
!
!  Output:
!    Log file
!
!  Written by:
!    Niclas Forsberg,
!    Dept. of Theoretical Chemistry, Lund University, 1995.

use mula_global, only: Huge_Print, ipow, MaxNumAt, nPolyTerm, VibModPlot
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: auTocm
use Definitions, only: u6

implicit real*8(a-h,o-z)
integer VibPlotUnit
real*8 energy
real*8 Hess(nOsc,nOsc), G(nOsc,nOsc), V(nOsc,nOsc)
real*8 D3(nOsc,nOsc,nOsc)
real*8 D4(nOsc,nOsc,nOsc,nOsc)
real*8 qMat(3*NumOfAt,nOsc)
real*8 harmfreq(nOsc)
real*8 xvec(nOsc)
integer InterVec(MaxNumAt*15)
real*8 Mass(MaxNumAt)
real*8 AtCoord(3,NumOfAt)
character*4 AtomLbl(MaxNumAt)
real*8 PED(nOsc,nOsc,nOsc)
real*8 x_anharm(nOsc,nOsc)
real*8 anharmfreq(nOsc)
integer Bond(2*MaxNumAt)
!logical Plot
logical ForceField, cont
real*8 stand_dev, max_err
character*11 VibPlotFile
character*8 BondString
real*8 PotCoef(nPolyTerm,1)
integer, allocatable :: aNormModes(:)
integer, external :: isfreeunit

NumInt = nOsc
call mma_allocate(aNormModes,NumInt,label='aNormModes')
do i=1,NumInt
  aNormModes(i) = i
end do

if (nState == 1) then
  write(u6,*)
  write(u6,*)
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ','|           Results for the first state           |'
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,*)
end if
if (nState == 2) then
  write(u6,*)
  write(u6,*)
  write(u6,*)
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ','|           Results for the second state          |'
  write(u6,'(a27,a)') ' ','|                                                 |'
  write(u6,'(a27,a)') ' ',' ================================================='
  write(u6,*)
end if

! Write coefficients of polynomial fit, standard deviation and
! maximum error to log.
if (.not. ForceField) then
  write(u6,*)
  write(u6,'(A)') '  Coefficients of different terms in polynomial:'
  write(u6,'(A)') '  ----------------------------------------------'
  do i=1,nPolyTerm
    write(u6,'(1x,f18.12,10(3x,i2))') PotCoef(i,1),ipow(i,:)
  end do
  write(u6,*)
  write(u6,*)
  write(u6,'(A,es14.6)') '  Standard deviation of fitted values:',stand_dev
  write(u6,'(A)') '  ------------------------------------'
  write(u6,*)
  write(u6,'(A,es14.6)') '  Maximum error of fitted values:     ',max_err
  write(u6,'(A)') '  -------------------------------'

  ! Write energy in minimum.
  write(u6,*)
  write(u6,*)
  write(u6,'(A,F15.8,A)') '  Energy in minimum:',energy
  write(u6,'(A)') '  ------------------'
end if

! Write force constant matrix to log file.
if (Huge_Print) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Force Constant Matrix:'
  write(u6,*) ' ','======================'
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (nInt,nInt=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (Hess(mInt,nInt),nInt=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do

  ! Write matrix G to log file.
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Inverse Mass Tensor :'
  write(u6,*) ' ','====================='
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (nInt,nInt=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (G(mInt,nInt),nInt=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do
  write(u6,*)

  ! Eigenvectors.
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Eigenvectors :'
  write(u6,*) ' ','=============='
  write(u6,*)
  maxCol = 10
  nRow = NumInt/maxCol+1
  nCol = mod(NumInt,maxCol)
  m1 = 1
  m2 = 1
  cont = .true.
  do i=1,nRow
    if ((i == nRow) .or. ((i == 1) .and. (NumInt < maxCol))) then
      if (i == 1) then
        m2 = m2+nCol-1
      else
        m2 = m2+nCol
      end if
    else
      if (i == 1) then
        m2 = m2+maxCol-1
      else
        m2 = m2+maxCol
      end if
    end if
    if ((nCol == 0) .and. (i == nRow)) cont = .false.
    if (cont) then
      write(u6,'(10I12)') (nInt,nInt=m1,m2)
      do mInt=1,NumInt
        write(u6,'(10F12.6)') (V(mInt,nInt),nInt=m1,m2)
      end do
      write(u6,*)
      write(u6,*)
      m1 = m2+1
    end if
  end do
  write(u6,*)
end if

! Third derivatives.
if (max_term > 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Cubic force constants:'
  write(u6,'(A)') '  ----------------------'
  do i=1,NumInt
    do j=i,NumInt
      do k=j,NumInt
        write(u6,'(a,i2,i2,i2,f15.8)') '  ',i,j,k,D3(i,j,k)
      end do
    end do
  end do
end if

! Fourth derivatives.
if (max_term > 3) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Quartic force constants:'
  write(u6,'(A)') '  ------------------------'
  do i=1,NumInt
    do j=i,NumInt
      do k=j,NumInt
        do l=k,NumInt
          write(u6,'(a,i2,i2,i2,i2,f15.8)') '  ',i,j,k,l,D4(i,j,k,l)
        end do
      end do
    end do
  end do
end if

! Anharmonicity constants.
if (max_term > 2) then
  write(u6,*)
  write(u6,*)
  write(u6,'(A)') '  Anharmonicity constants:'
  write(u6,'(A)') '  ------------------------'
  do i=1,NumInt
    write(u6,'(a,20f15.8)') ' ',(auTocm*x_anharm(i,j),j=1,i)
  end do
end if

! Write labels, coordinates and masses to log file.
call WriteCartCoord(AtomLbl,AtCoord,Mass,NumOfAt)

! Write internal coordinates to log file.
call WriteIntCoord(InterVec,AtomLbl,xvec,NumInt)

! Harmonic frequencies in reciprocal cm, GHz and hartrees.
call WriteFreq(harmfreq,aNormModes,NumInt,'Harmonic frequencies')

! Fundamental frequencies in reciprocal cm, GHz and hartrees.
if (max_term > 2) call WriteFreq(anharmfreq,aNormModes,NumInt,'Anharmonic frequencies')

! Write Potential Energy Distribution for each mode to log
! if the molecule is not too large.
if (Huge_Print .and. (NumInt <= 10)) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Potential Energy Distribution :'
  write(u6,*) ' ','==============================='
  do imode=1,NumInt
    write(u6,*)
    write(u6,*) ' ','Mode',imode
    write(u6,*) ' ','-------'
    do k=1,NumInt
      if (NumInt < 21) then
        write(u6,'(60F6.2)') (PED(k,l,imode),l=1,NumInt)
      else
        write(u6,'(60F5.2)') (PED(k,l,imode),l=1,NumInt)
      end if
    end do
  end do
end if

! Find length of longest displacement vector.
vMax = Dnrm2_(3,qMat(1,1),1)
do nAtom=1,NumOfAt
  do j=1,NumInt
    vLength = Dnrm2_(3,qMat((nAtom-1)*3+1,j),1)
    if (vLength > vMax) then
      vMax = vLength
    end if
  end do
end do

! Print displacement vectors scaled with length of longest vector
! if the molecule is not too large.
if (Huge_Print .and. (NumInt <= 10)) then
  write(u6,*)
  write(u6,*)
  write(u6,*) ' ','Scaled displacement vectors :'
  write(u6,*) ' ','============================='
  write(u6,*) ' ','(length of longest vector set equal to one)'
  write(u6,*)
  maxCol = 10
  do nAtom=1,NumOfAt
    write(u6,*) ' ',' Atom: ',AtomLbl(nAtom)
    write(u6,*) ' ',' ---------'
    nRow = NumInt/maxCol
    nCol = mod(NumInt,maxCol)
    m1 = 1
    m2 = 0
    do i=0,nRow
      if (i == nRow) then
        m2 = m2+nCol
      else
        m2 = m2+maxCol
      end if
      if (m2 /= 0) then
        write(u6,'(10I12)') (j,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  x ',(qMat((nAtom-1)*3+1,j)/vMax,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  y ',(qMat((nAtom-1)*3+2,j)/vMax,j=m1,m2)
        write(u6,'(A4,10F12.8)') '  z ',(qMat((nAtom-1)*3+3,j)/vMax,j=m1,m2)
        write(u6,*)
        m1 = m2+1
      end if
    end do
  end do
  write(u6,*)
  write(u6,*)
end if

! Generate MOLDEN input:

if (nState == 1) then
  call WrMold('MOLDEN-1',NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
else
  call WrMold('MOLDEN-2',NumOfAt,AtomLbl,AtCoord,NumInt,HarmFreq,QMat)
end if

! Print geometry of molecule, together with cartesian displacements,
! in a file.
if (VibModPlot) then
  VibPlotUnit = isfreeunit(21)
  if (nState == 1) then
    VibPlotFile = 'plot.modes1'
  end if
  if (nState == 2) then
    VibPlotFile = 'plot.modes2'
  end if
  call molcas_Open(vibplotunit,vibplotfile)
  !open(VibPlotUnit,VibPlotFile)
  write(VibPlotUnit,*) 2*NumOfAt
  write(VibPlotUnit,*) nBond/2
  write(VibPlotUnit,*) '1'
  write(VibPlotUnit,*) NumInt
  iBond = 1
  do while (iBond+1 <= nBond)
    write(BondString,'(i4,i4)') Bond(iBond),Bond(iBond+1)
    write(VibPlotUnit,'(a8)') BondString
    iBond = iBond+2
  end do
  do nAtom=1,NumOfAt
    x1 = AtCoord(1,nAtom)
    x2 = AtCoord(2,nAtom)
    x3 = AtCoord(3,nAtom)
    write(VibPlotUnit,'(a,3f12.3)') '>',x2,x3,x1
  end do
  do iInt=1,NumInt
    do nAtom=1,NumOfAt
      x1 = AtCoord(1,nAtom)
      x2 = AtCoord(2,nAtom)
      x3 = AtCoord(3,nAtom)
      write(VibPlotUnit,'(a,3f12.3)') AtomLbl(nAtom),x2+qMat((nAtom-1)*3+2,iInt)/vMax,x3+qMat((nAtom-1)*3+3,iInt)/vMax, &
                                      x1+qMat((nAtom-1)*3+1,iInt)/vMax
    end do
  end do
  close(VibPlotUnit)
end if

call mma_deallocate(aNormModes)

end subroutine WriteLog
