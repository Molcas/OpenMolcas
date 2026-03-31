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

subroutine TSHop(CI1,CI2)

use rassi_aux, only: ipglob
use rassi_global_arrays, only: JBNUM, LROOT
use Cntrl, only: ChkHop, ISTATE1, ISTATE2, iTOC15, JBNAME, LuIph, nCI1, nCI2
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, Quart
use Definitions, only: wp, iwp, u6

implicit none
real(kind=wp) :: CI1(NCI1), CI2(NCI2)
integer(kind=iwp) :: file1, file2, I, IAD, IAD3, IADR3(3), IDISK, JOB1, JOB2, LROOT1, maxHop, nHop
real(kind=wp) :: prdct(2,2)
logical(kind=iwp) :: fexist, lAllowHop, lHop, lHopped, lMaxHop
character(len=80) :: filnam, filother
real(kind=wp), allocatable :: CI1pr(:), CI2pr(:)

! Skip the test if a hop has occurred
!  (this should happen when testing for state n+1
!  right after a hop to state n-1)
call Get_iScalar('Relax CASSCF root',I)
if (I /= ISTATE1) then
  if (IPGLOB >= 2) write(u6,'(6X,A)') 'A hop has just been detected, skipping this state.'
  return
end if

call mma_allocate(CI1pr,NCI1,Label='CI1pr')
call mma_allocate(CI2pr,NCI2,Label='CI2pr')

! Initialization
CI1pr(:) = Zero
CI2pr(:) = Zero
file1 = 83
file2 = 83
lHopped = .false.

! Get the CI coefficients for current state

! Open JOBIPH file:
JOB1 = JBNUM(ISTATE1)
call DANAME(LUIPH,JBNAME(JOB1))
! Read table of contents on this JOBIPH file:
IAD = 0
call IDAFILE(LUIPH,2,ITOC15,15,IAD)
! Read CI coefficients from interface.
IDISK = ITOC15(4)
LROOT1 = LROOT(ISTATE1)
do I=1,LROOT1-1
  call DDAFILE(LUIPH,0,CI1,NCI1,IDISK)
end do
call DDAFILE(LUIPH,2,CI1,NCI1,IDISK)
call DACLOS(LUIPH)

! Get the CI coefficients for state2

! Open JOBIPH file:
JOB2 = JBNUM(ISTATE2)
call DANAME(LUIPH,JBNAME(JOB2))
! Read table of contents on this JOBIPH file:
IAD = 0
call IDAFILE(LUIPH,2,ITOC15,15,IAD)
! Read CI coefficients from interface.
IDISK = ITOC15(4)
LROOT1 = LROOT(ISTATE2)
do I=1,LROOT1-1
  call DDAFILE(LUIPH,0,CI2,NCI2,IDISK)
end do
call DDAFILE(LUIPH,2,CI2,NCI2,IDISK)
call DACLOS(LUIPH)

! Check if it is a hop up or hop down

I = ISTATE1-ISTATE2
if (I > 0) then
  if (IPGLOB >= 2) write(u6,*) 'Checking for a hop to a root lower in energy.'
  filnam = 'CIVECTOR'
  filother = 'CIVECTUP'
else if (I < 0) then
  if (IPGLOB >= 2) write(u6,*) 'Checking for a hop to a root higher in energy.'
  filnam = 'CIVECTUP'
  filother = 'CIVECTOR'
else
  write(u6,*) 'Unknown problem'
  write(u6,*) 'ISTATE1 = ',ISTATE1
  write(u6,*) 'ISTATE2 = ',ISTATE2
end if

! Open the file with the previous CI-vectors

call f_inquire(filnam,fexist)
if (fexist) then
  call DANAME(file1,filnam)
  if (IPGLOB >= 3) write(u6,*) trim(filnam)//' file exists.'
else
  ! If the file does not exist, create a new one with the current vectors
  call DANAME(file1,filnam)
  ! Dummy table of contents is written
  do I=1,3
    IADR3(I) = 0
  end do
  IAD3 = 0
  call IDAFILE(file1,1,IADR3,3,IAD3)
  IADR3(1) = IAD3
  ! Current CI coefficients are written
  call DDAFILE(file1,1,CI1,NCI1,IAD3)
  IADR3(2) = IAD3
  call DDAFILE(file1,1,CI2,NCI2,IAD3)
  IADR3(3) = IAD3
  ! Write the real table of contents
  IAD3 = 0
  call IDAFILE(file1,1,IADR3,3,IAD3)
  if (IPGLOB >= 3) write(u6,*) trim(filnam)//' file created.'
end if

! Check for surface hop if the energy difference is smaller than
! the threshold.

if (ChkHop) then

  ! Read table of contents on this file
  IAD3 = 0
  call IDAFILE(file1,2,IADR3,3,IAD3)
  ! Read the CI coefficients of ISTATE1 from previuos step
  IAD3 = IADR3(1)
  call DDAFILE(file1,2,CI1pr,NCI1,IAD3)
  ! Read the CI coefficients of ISTATE2 from previuos step
  IAD3 = IADR3(2)
  call DDAFILE(file1,2,CI2pr,NCI2,IAD3)

  ! Calculate the scalar product of the CI coefficient vectors.
  prdct(1,1) = Zero
  prdct(1,2) = Zero
  prdct(2,1) = Zero
  prdct(2,2) = Zero
  if (IPGLOB >= 3) write(u6,'(4(A16))') 'CI1','CI1pr','CI2','CI2pr'
  do i=1,NCI1
    if (IPGLOB >= 3) write(u6,'(4(3X,ES13.6))') CI1(i),CI1pr(i),CI2(i),CI2pr(i)
    prdct(1,1) = prdct(1,1)+CI1pr(i)*CI1(i)
    prdct(1,2) = prdct(1,2)+CI1pr(i)*CI2(i)
    prdct(2,1) = prdct(2,1)+CI2pr(i)*CI1(i)
    prdct(2,2) = prdct(2,2)+CI2pr(i)*CI2(i)
  end do
  if (IPGLOB >= 2) then
    write(u6,'(6X,A)') 'The scalar products of the CI-vectors:'
    write(u6,3000) 'CIpr(state1) * CI(state1) =',prdct(1,1)
    write(u6,3000) 'CIpr(state1) * CI(state2) =',prdct(1,2)
    write(u6,3000) 'CIpr(state2) * CI(state1) =',prdct(2,1)
    write(u6,3000) 'CIpr(state2) * CI(state2) =',prdct(2,2)
    write(u6,*)
  end if
  ! Check the conditions for a surface hop
  if ((abs(prdct(1,2)) >= Quart) .and. (abs(prdct(2,1)) >= Quart)) then
    write(u6,'(6X,3A)') '+',repeat('-',78),'+'
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,A1,T35,A,T86,A1)') '|','A HOP event is detected!','|'
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,A1,T32,2(A,I3,4X),T86,A1)') '|','From state:',ISTATE1,'To state:',ISTATE2,'|'
    ! Check if the number of Hops is limited:
    call qpg_iScalar('MaxHops',lMaxHop)
    if (lMaxHop) then
      call Get_iScalar('MaxHops',maxHop)
      if (maxHop < 1) lMaxHop = .false.
    end if
    if (lMaxHop) then
      call qpg_iScalar('Number of Hops',lHop)
      if (lHop) then
        call Get_iScalar('Number of Hops',nHop)
      else
        nHop = 0
      end if
      if (maxHop <= nHop) then
        lAllowHop = .false.
        write(u6,'(6X,A1,T40,A,T86,A1)') '|','maxHop > nHop','|'
        write(u6,'(6X,A1,T31,A,T86,A1)') '|','This surface HOP is not allowed','|'
        write(u6,'(6X,A1,T24,A,T86,A1)') '|','because the number of allowed Hops is exceeded','|'
      else
        lAllowHop = .true.
      end if
    else
      lAllowHop = .true.
    end if
    ! Bottom of the printed box
    write(u6,'(6X,A1,T86,A1)') '|','|'
    write(u6,'(6X,3A,//)') '+',repeat('-',78),'+'
    ! Set the numbers of Hops
    if (lAllowHop) then
      call Put_iScalar('Relax CASSCF root',ISTATE2)
      call Put_iScalar('NumGradRoot',ISTATE2)
      nHop = nHop+1
      call Put_iScalar('Number of Hops',nHop)
      lHopped = .true.
    end if
  end if
end if
if (IPGLOB >= 3) then
  write(u6,'(2(6X,A8,I3))') 'ISTATE1=',ISTATE1,'ISTATE2=',ISTATE2
  do i=1,NCI1
    write(u6,'(6X,ES12.5,8X,ES12.5)') CI1(i),CI2(i)
  end do
end if

call mma_deallocate(CI1pr)
call mma_deallocate(CI2pr)

! Save the CI coefficients to the right file,
! dependending on whether or not a hop has occurred

if (lHopped) then

  ! Delete the file with the CI-vectors
  call DaEras(file1)

  ! Note that, if a hop occurred, the vectors are written
  ! in the *other* file, and in reversed order
  call f_inquire(filother,fexist)
  if (fexist) then
    call DANAME(file2,filother)
    IAD3 = 0
    call IDAFILE(file2,2,IADR3,3,IAD3)
    IAD3 = IADR3(1)
    call DDAFILE(file2,1,CI2,NCI2,IAD3)
    IAD3 = IADR3(2)
    call DDAFILE(file2,1,CI1,NCI1,IAD3)
    call DACLOS(file2)
  else
    ! The file does not exist, create a new one
    call DANAME(file2,filother)
    IAD3 = 0
    call IDAFILE(file2,1,IADR3,3,IAD3)
    IADR3(1) = IAD3
    call DDAFILE(file2,1,CI2,NCI2,IAD3)
    IADR3(2) = IAD3
    call DDAFILE(file2,1,CI1,NCI1,IAD3)
    IADR3(3) = IAD3
    IAD3 = 0
    call IDAFILE(file2,1,IADR3,3,IAD3)
    call DACLOS(file2)
  end if
else

  ! Write the CI-vectors normally if no hop occurred
  IAD3 = 0
  call IDAFILE(file1,2,IADR3,3,IAD3)
  IAD3 = IADR3(1)
  call DDAFILE(file1,1,CI1,NCI1,IAD3)
  IAD3 = IADR3(2)
  call DDAFILE(file1,1,CI2,NCI2,IAD3)
  call DACLOS(file1)
end if

return

3000 format(6X,A,F7.2)

end subroutine TSHop
