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

subroutine PRIMO(Header,PrOcc,PrEne,ThrOcc,ThrEne,nSym,nBas,nOrb,BName,Ene,Occ,CMO,iPrForm)
!***********************************************************************
!                                                                      *
! Purpose: print a set of orbitals                                     *
!                                                                      *
! Calling parameters:                                                  *
! Header : Header line which is printed together with the              *
!          section title                                               *
! nSym   : number of irreducible representations                       *
! nOrb   : Total number of orbitals per irred. rep.                    *
! nBas   : Number of basis functions per irred. rep.                   *
! BName  : Center and function type label per basis function           *
! Ene    : Orbital Energies                                            *
! CMO    : Orbital coefficients                                        *
! Occ    : Orbital Occupations                                         *
! PrEne  : Switch to turn printing of orbital Energies ON/OFF          *
! PrOcc  : Switch to turn printing of orbital Occupatios ON/OFF        *
! ThrOcc : Threshold for Occupation number to be printed               *
!          Orbitals with Occupation number less than ThrOcc are        *
!          not printed                                                 *
! ThrEne : Threshold for orbitals Energies to be printed               *
!          Orbitals with Energies larger than ThrEne are               *
!          not printed                                                 *
! iPrForm: print level: -1 - Old way, 0- None, 1-list, 2-Short,        *
!          3-long                                                      *
!                                                                      *
!***********************************************************************

use Constants, only: Zero, One, Ten, Half
use Definitions, only: wp, iwp, u6

implicit none
#include "Molcas.fh"
character(len=*), intent(in) :: Header
logical(kind=iwp), intent(in) :: PrOcc, PrEne
real(kind=wp), intent(in) :: ThrOcc, ThrEne, Ene(*), Occ(*), CMO(*)
integer(kind=iwp), intent(in) :: nSym, nBas(*), nOrb(*), iPrForm
character(len=LenIn8), intent(in) :: BName(*)
integer(kind=iwp) :: i, iB, IEO, IEO2, Inc, IO, iPL, iPos, ISB, ISCMO, ISO, ISX, iSym, itmp, jB, jSB, nB, nCol, NCOLS, nO, nTot
real(kind=wp) :: Cff_Mx, Energy, Occupation, tmp, Shift(10)
logical(kind=iwp) :: Large, Go, Header_Done
integer(kind=iwp), parameter :: Magic = 5+1+LenIn8+1+6+3
character(len=Magic) :: ChCMO
character(len=180) :: Line
character(len=24) :: FMT0, FMT1, FMT2, LABEL0, LABEL1, LABEL2
character(len=20) :: Fmt_s, Fmt_l, Fmt_f, Frmt
character(len=3) :: IrrepName(8)
character(len=4) :: Star(10)
integer(kind=iwp), external :: iPrintLevel
logical(kind=iwp), external :: Reduce_Prt
character(len=LenIn8), external :: Clean_BName

!                                                                      *
!----------------------------------------------------------------------*
!                                                                      *
do i=1,10
  Star(i) = ' '
end do
iPL = iPrintLevel(-1)
if (Reduce_Prt() .and. (iPL < 3)) iPL = 0
if (iPL <= 1) return

if (iPrForm == 0) return

!----------------------------------------------------------------------*
! Print Header                                                         *
!----------------------------------------------------------------------*

write(u6,*)
call CollapseOutput(1,'   Molecular orbitals:')
write(u6,'(6X,A)') '-------------------'
write(u6,*)
write(u6,'(6X,2A)') 'Title: ',trim(Header)
!write(u6,*)
!write(u6,*) 'test print out'
!----------------------------------------------------------------------*
! Legend (general)                                                     *
!----------------------------------------------------------------------*
if (iPrForm >= 4) then
  write(u6,'(6X,A)') 'LEGEND'
  write(u6,'(6X,A)') '============================================================================================================='
  write(u6,'(6X,A)') 'A basis set label has the structure "nl+m" or "nl-m" where n is the principal quantum number for atomic basis'
  write(u6,'(6X,A)') 'functions (a star, "*", denotes a primitive, diffuse or a polarization function), l denotes the angular,'
  write(u6,'(6X,A)') 'and m the magnetic quantum number. The basis functions are normally real Spherical harmonics. For Cartesian'
  write(u6,'(6X,A)') 'type basis functions we use the notation "ijk" to denote the power of the x, y, and z component in the radial'
  write(u6,'(6X,A)') 'part. For p-functions we always use the notation px, py, and pz.'
end if
!----------------------------------------------------------------------*
! Define standard output format                                        *
!----------------------------------------------------------------------*
NCOLS = 10
LABEL0 = 'Orbital '
LABEL1 = 'Energy  '
LABEL2 = 'Occ. No.'
FMT0 = '(10X,A12,3X,10(I5,A,1X))'
FMT1 = '(10X,A12,2X,10F10.4)'
FMT2 = '(5X,I4,1X,A,10F10.4)'

!----------------------------------------------------------------------*
! Set up list of pointers                                              *
!----------------------------------------------------------------------*

ISX = 1
ISB = 1
ISCMO = 1

!----------------------------------------------------------------------*
! Start loop over all irreducible representations                      *
!----------------------------------------------------------------------*

nTot = 0
do iSym=1,nSym
  nTot = nTot+nBas(iSym)
end do
Large = .false.
if (iPrForm == -1) Large = nTot > 256
if (iPrForm == 1) Large = .false.
if (iPrForm == 2) Large = .true.
if (iPrForm == 3) Large = .false.
if (iPrForm == 4) Large = .false.
call Get_cArray('Irreps',IrrepName,24)

do iSym=1,nSym
  nB = nBas(iSym)
  nO = nOrb(iSym)
  if (nO == 0) cycle

  Header_Done = .false.

  !--------------------------------------------------------------------*
  ! Start loop over columns                                            *
  !--------------------------------------------------------------------*

  if (Large) then
    FMT_s = '(I5,1X,A,A,F6.3,A)'
    FMT_l = '(I5,1X,A,A,F6.2,A)'
    FMT_f = '(I5,1X,A,A,F6.1,A)'

    do iSO=0,nO-1
      if (PrEne .and. PrOcc) then
        Go = (Ene(ISX+IsO) <= ThrEne) .and. (Occ(ISX+IsO) >= ThrOcc)
      else if (PrEne) then
        Go = Ene(ISX+IsO) <= ThrEne
      else if (PrOcc) then
        Go = Occ(ISX+IsO) >= ThrOcc
      else
        Go = .false.
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Go=',Go
#     endif
      if (Go) then
        if (.not. Header_Done) then
          write(u6,'(/6X,A,I2,A,A)') 'Molecular orbitals for symmetry species',iSym,': ',IrrepName(iSym)
          write(u6,*)

          !------------------------------------------------------------*
          ! Start by listing the basis of this irrep                   *
          !------------------------------------------------------------*

          write(u6,*)
          write(u6,'(6X,A)') 'Basis set list:'
          nCol = min((nB+9)/10,5)
          Inc = (nB+nCol-1)/nCol
          jSB = iSB-1
          do iB=1,Inc
            write(u6,'(4X,5(I5,1X,A,9X))') (jB,Clean_BName(BName(jSB+jB),LenIn),jB=iB,nB,Inc)
          end do
          write(u6,*)

          write(u6,*)
          write(u6,*) ' Index Energy  Occupation Coefficients ...'
          write(u6,*)
          Header_Done = .true.
        end if
        ! This will occupy the first 25 positions
        if (PrEne) then
          Energy = Ene(ISX+iSO)
        else
          Energy = Zero
        end if
        if (PrOcc) then
          Occupation = Occ(ISX+iSO)
        else
          Occupation = Zero
        end if
        write(Line,'(I5,2F10.4)') iSO+1,Energy,Occupation
#       ifdef _DEBUGPRINT_
        write(u6,'(A,A)') 'Line=',Line
#       endif
        iPos = 1
        Cff_Mx = Zero
        do iB=0,NB-1
          if (abs(CMO(isCMO+iSO*nB+iB)) > Cff_Mx) Cff_Mx = abs(CMO(isCMO+iSO*nB+iB))
        end do
#       ifdef _DEBUGPRINT_
        write(u6,*) 'Cff_Mx=',Cff_Mx
#       endif
        do iB=0,NB-1
#         ifdef _DEBUGPRINT_
          write(u6,*) ' iB=',iB
#         endif
          if (abs(CMO(isCMO+iSO*nB+iB)) >= Cff_Mx*Half) then
#           ifdef _DEBUGPRINT_
            write(u6,*) ' Process iB=',iB
#           endif
            if (abs(CMO(isCMO+iSO*nB+iB)) >= 100.0_wp) then
              Frmt = Fmt_f
            else if (abs(CMO(isCMO+iSO*nB+iB)) >= Ten) then
              Frmt = Fmt_l
            else
              Frmt = Fmt_s
            end if
#           ifdef _DEBUGPRINT_
            write(u6,*) ' Frmt=',Frmt
#           endif
            write(ChCMO,Frmt) iB+1,Clean_BName(BName(ISB+IB),LenIn),'(',CMO(isCMO+iSO*nB+iB),'), '
#           ifdef _DEBUGPRINT_
            write(u6,'(A)') ChCMO
#           endif
            if (iPos == 1) then
              Line(2+Magic:2+Magic*2-1) = ChCMO
              iPos = 2
            else if (iPos == 2) then
              Line(2+Magic*2:2+Magic*3-1) = ChCMO
              iPos = 3
            else if (iPos == 3) then
              Line(2+Magic*3:2+Magic*4-1) = ChCMO
              iPos = 4
            else if (iPos == 4) then
              Line(2+Magic*4:2+Magic*5-1) = ChCMO
              iPos = 1
              write(u6,'(A)') Line
              Line = ' '
            end if
          end if
        end do
        if (iPos /= 1) write(u6,'(A)') Line
        write(u6,*)
      end if
    end do

  else

    do ISO=0,NO-1,NCOLS
      IEO = ISO+NCOLS-1
      IEO = min((NO-1),IEO)
      IEO2 = IEO

      !vv NAG
      !do IO=ISO,IEO2
      !  if (PrEne .and. PrOcc) then
      !    if ((Ene(ISX+IO) > ThrEne) .or. (Occ(ISX+IO) < ThrOcc)) IEO = IEO-1
      !  else if (PrEne) Then
      !    if (Ene(ISX+IO) > ThrEne ) IEO = IEO-1
      !  else if (PrOcc) Then
      !    if (Occ(ISX+IO) < ThrOcc ) IEO = IEO-1
      !  end if
      !end do

      if (PrEne .and. PrOcc) then
        do IO=ISO,IEO2
          if ((Ene(ISX+IO) > ThrEne) .or. (Occ(ISX+IO) < ThrOcc)) IEO = IEO-1
        end do
      else if (PrEne) then
        do IO=ISO,IEO2
          if (Ene(ISX+IO) > ThrEne) IEO = IEO-1
        end do
      else if (PrOcc) then
        do IO=ISO,IEO2
          if (Occ(ISX+IO) < ThrOcc) IEO = IEO-1
        end do
      end if
      if (IEO >= ISO) then
        if (.not. Header_Done) then
          write(u6,'(/6X,A,I2,A,A)') 'Molecular orbitals for symmetry species',iSym,': ',IrrepName(iSym)
          Header_Done = .true.
        end if
        if (PrEne) then
          tmp = Zero
          do IO=ISO,IEO
            Shift(IO-ISO+1) = One
            tmp = max(abs(Ene(ISX+IO)),tmp)
          end do
          if (tmp > 1.0e3_wp) then
            write(u6,*)
            write(u6,'(10X,A)') 'Some orbital energies have been scaled by powers of 10, the power is written right after the &
                               &orbital index'
            write(u6,*)
            do IO=ISO,IEO
              tmp = abs(Ene(ISX+IO))
              if (tmp > 1.0e3_wp) then
                tmp = tmp*1.0e-2_wp
                itmp = int(log10(tmp))
                Shift(IO-ISO+1) = Ten**itmp
                write(Star(IO-ISO+1),'(A,I1,A)') ' (',itmp,')'
              else
                Star(IO-ISO+1) = '    '
              end if
            end do
          else
            do IO=ISO,IEO
              Star(IO-ISO+1) = '    '
            end do
          end if
        end if
        write(u6,*)
        write(u6,FMT0) LABEL0,(1+IO,Star(IO-ISO+1),IO=ISO,IEO)
        if (PrEne) write(u6,FMT1) LABEL1,(Ene(ISX+IO)/Shift(IO-ISO+1),IO=ISO,IEO)
        if (PrOcc) write(u6,FMT1) LABEL2,(Occ(ISX+IO),IO=ISO,IEO)
        write(u6,*)
        if (iPrForm /= 1) then
          do IB=0,NB-1
            write(u6,FMT2) IB+1,Clean_BName(BName(ISB+IB),LenIn),(CMO(ISCMO+IO*NB+IB),IO=ISO,IEO)
          end do
        end if
      end if

      !----------------------------------------------------------------*
      ! End of loop over columns                                       *
      !----------------------------------------------------------------*

    end do
  end if

  !--------------------------------------------------------------------*
  ! Update list of pointers                                            *
  !--------------------------------------------------------------------*

  ISX = ISX+NO
  ISB = ISB+NB
  ISCMO = ISCMO+NO*NB

  !--------------------------------------------------------------------*
  ! End of loop over all irreducible representations                   *
  !--------------------------------------------------------------------*

end do

call CollapseOutput(0,'   Molecular orbitals:')
write(u6,*)

return

end subroutine PRIMO
