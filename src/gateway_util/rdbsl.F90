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
! Copyright (C) 1991, Bjorn O. Roos                                    *
!               1991, Roland Lindh                                     *
!               Valera Veryazov                                        *
!***********************************************************************

subroutine Rdbsl(BasDir,BSLbl,bType,nCGTO,mCGTO,lAng,lCGTO,lUnit,iAtmNr,BasisTypes,ExtBasDir)
!***********************************************************************
! Object: Decode the basis set label and read the basis set            *
!         from the library                                             *
!                                                                      *
! Called from: GetBS                                                   *
! Subroutines called: Decode                                           *
!                                                                      *
! Author: Bjoern Roos, Theoretical Chemistry, Chemical Centre          *
!         University of Lund, Lund, Sweden                             *
!         February 1991                                                *
! Patched: Valera Veryazov                                             *
!***********************************************************************

use Definitions, only: wp, iwp, u6

implicit none
character(len=*), intent(inout) :: BasDir
character(len=80), intent(inout) :: BSLbl
character(len=80), intent(out) :: bType
integer(kind=iwp), intent(in) :: lCGTO, lUnit
integer(kind=iwp), intent(out) :: nCGTO(0:lCGTO), mCGTO(0:lCGTO), lAng, iAtmNr, BasisTypes(4)
character(len=*), intent(in) :: ExtBasDir
#include "getlnqoe.fh"
#include "angtp.fh"
integer(kind=iwp) :: i, i1, iLast1, iLast2, iLast3, iLast4, iLast_JR, irecl, istatus, j, k, lAngm, n
real(kind=wp) :: cg
logical(kind=iwp) :: Do_Cycle, Exists, Hit, is_error, lStop
character(len=256) :: BasLoc, FileName, TmpString
character(len=180) :: BSLB
character(len=80) :: atom, atomb, author, aux, aux2, basis, CGTO, CGTOm, string
character(len=16) :: Str16
character :: kAng(0:iTabMx)
#ifdef _DEBUGPRINT_
#define _TEST_ .true.
#else
#define _TEST_ .true.
#endif
logical(kind=iwp), parameter :: IfTest = .false.
integer(kind=iwp), external :: Lbl2Nr
character(len=180), external :: Get_Ln_Quit

if (iTabMx < lCGTO) then
  call WarningMessage(2,'RdBsL: iTabMx < lCGTO;Update the code!')
  call Abend()
else
  do i=0,iTabMx
    kAng(i) = Angtp(i)
    call UpCase(kAng(i))
  end do
end if
if (IfTest) then
  write(u6,*) ' Enter Rdbsl'
  write(u6,'(2a)') ' BsLbl=',BsLbl
  write(u6,'(2a)') ' BasDir=',BasDir
  write(u6,'(2a)') ' ExtDir=',ExtBasDir
end if
call ResetErrorLine()
call UpCase(BSLBl)

iLast1 = len_trim(BasDir)
BasLoc = BasDir
if (ExtBasDir /= ' ') then
  Hit = .true.
  call Decode(BSLBl(1:80),bType,2,Hit)
end if
call find_basis_set(BasLoc,ExtBasDir,bType)
iLast_JR = len_trim(BasLoc)
if (index(BasDir,'c_Basis') /= 0) then
  BasLoc(iLast_JR+1:iLast_JR+8) = '/c_Basis'
  iLast_JR = iLast_JR+8
else if (index(BasDir,'j_Basis') /= 0) then
  BasLoc(iLast_JR+1:iLast_JR+8) = '/j_Basis'
  iLast_JR = iLast_JR+8
else if (index(BasDir,'jk_Basis') /= 0) then
  BasLoc(iLast_JR+1:iLast_JR+9) = '/jk_Basis'
  iLast_JR = iLast_JR+9
end if
call BasisTbl(BSLBl,BasLoc(1:iLast_JR))

Hit = .true.
call Decode(BSLBl(1:80),atom,1,Hit)
if (IfTest) write(u6,'(1x,a,a)') 'Atom=',atom
iAtmNr = Lbl2Nr(atom)

Hit = .true.
call Decode(BSLBl(1:80),bType,2,Hit)
if (IfTest) write(u6,'(1x,a,a)') 'Type=',bType

Hit = .true.
call Decode(BSLBl(1:80),author,3,Hit)
if (IfTest) write(u6,'(1x,a,a)') 'Author=',author

Hit = .true.
call Decode(BSLBl(1:80),basis,4,Hit)
if (IfTest) write(u6,'(1x,a,a)') 'Basis=',basis

Hit = .true.
call Decode(BSLBl(1:80),CGTO,5,Hit)
if (IfTest) write(u6,'(1x,a,a)') 'CGTO=',CGTO

Hit = .false.
call Decode(BSLBl(1:80),Aux,6,Hit)
if (.not. Hit) Aux = ' '
if (IfTest) write(u6,'(1x,a,a)') 'Aux=',Aux

Hit = .false.
call Decode(BSLBl(1:80),Aux2,7,Hit)
if (.not. Hit) Aux2 = ' '
if (IfTest) write(u6,'(1x,a,a)') 'Aux2=',Aux2

if ((CGTO == '') .and. ((bType == 'ECP') .or. (bType == 'PSD') .or. (bType == 'ANO'))) then
  if (IfTest) write(u6,*) ' Early exit'
  call WarningMessage(2,'Abend in RdBsl:No CGTO basis set provided in basis label')
  call Quit_OnUserError()
end if

! We do not need more than 10 dots

j = 0
do i=1,10
  j = j+index(BsLbl(j+1:80),'.')
  if (j == 0) exit
end do
if (j > 0) BsLbl = BsLbl(1:j)

! Open basis library

TmpString = bType
call Upcase(TmpString)
BasDir = adjustl(BasDir)
TmpString = adjustl(TmpString)
iLast1 = len_trim(BasDir)
iLast2 = len_trim(TmpString)
if (IfTest) then
  write(u6,'(I3,A)') iLast1,BasDir
  write(u6,'(I3,A)') iLast2,TmpString
end if
Filename = BasLoc(1:iLast_JR)//'/'//TmpString(1:iLast2)
TmpString = Author
call Upcase(TmpString)
iLast4 = len_trim(Filename)
iLast2 = len_trim(TmpString)
TmpString = Filename(1:iLast4)//'.'//TmpString(1:iLast2)
call f_Inquire(TmpString,Exists)
if (Exists) Filename = TmpString
iLast3 = len_trim(Filename)
if (IfTest) write(u6,'(A,A)') 'Filename=',Filename
call f_Inquire(Filename,Exists)
if (.not. Exists) then
  ! Try to find name in trans.tbl file
  call TransTbl(Filename)
  call f_Inquire(Filename,Exists)
  if (.not. Exists) then
    iLast3 = len_trim(Filename)
    call WarningMessage(2,'Basis set file '//Filename(1:iLast3)//' does not exist!;;(1) For a valence basis set: check the '// &
                        'spelling of the basis set label and that the basis set file is present in the basis set library '// &
                        'directory.;;(2) For an external auxiliary basis set: check that the basis set file is present in the '// &
                        'appropiate basis set library subdirectory.')
    call Quit_OnUserError()
  end if
end if
! check basistype
call BasisType(Filename,0,BasisTypes)
call molcas_open_ext2(lUnit,FileName,'sequential','formatted',istatus,.false.,irecl,'unknown',is_error)
!open(unit=lUnit,file=Filename,form='FORMATTED',iostat=istatus)
if (istatus /= 0) then
  iLast3 = len_trim(Filename)
  call WarningMessage(2,' Problems opening basis set file '//Filename(1:iLast3))
  call Quit_OnUserError()
end if
rewind(lUnit)

! loop over the basis set library to find the correct label

if (IfTest) write(u6,*) ' Locate basis set label in library'

Do_Cycle = .true.
do while (Do_Cycle)

  BSLB = Get_Ln_Quit(lUnit,0)
  if (Quit_On_Error) then
    iLast3 = len_trim(BsLbl)
    call WarningMessage(2,'The requested basis set label: '//BsLbl(:iLast3)//';'//'was not found in basis library: '//Filename)
    call Abend()
  end if
  call UpCase(BSLB)
  if (BSLB(1:1) /= '/') cycle
  n = index(BSLB,' ')
  do i=n,80
    BSLB(i:i) = '.'
  end do
  Hit = .true.
  call Decode(BSLB(2:80),atomb,1,Hit)
  if (atomb /= atom) cycle

  if (bType /= '') then
    Hit = .true.
    call Decode(BSLB(2:80),string,2,Hit)
    if (string /= bType) cycle
  end if

  if (author /= '') then
    Hit = .true.
    call Decode(BSLB(2:80),string,3,Hit)
    if (string /= author) cycle
  end if

  if (basis /= '') then
    Hit = .true.
    call Decode(BSLB(2:80),string,4,Hit)
    if (string /= basis) cycle
  end if

  ! If a contraction sequence has been specified it must be identical
  ! to what is in the library file if the basis set type does not
  ! not allow any other contraction sequence.

  if ((CGTO /= '') .and. (bType(1:3) /= 'ANO') .and. (bType /= 'ECP') .and. (bType /= 'PSD') .and. (bType /= 'RYDBERG') .and. &
      (Aux /= 'ECP')) then
    Hit = .true.
    call Decode(BSLB(2:80),string,5,Hit)
    if (string /= CGTO) cycle
  end if

  Hit = .true.
  call Decode(BSLB(2:80),string,6,Hit)
  if (Aux /= '') then
    if (string /= aux) cycle
  else if (string == 'ECP') then
    ! If the library basis (BSLB) has ECP, add it to the BsLbl
    j = 0
    do i=1,5
      j = j+index(BsLbl(j+1:80),'.')
    end do
    BsLbl = BsLbl(1:j)//'ECP'//BsLbl(j+1:)
    Aux = 'ECP'

  end if

  if (Aux2 /= '') then
    Hit = .true.
    call Decode(BSLB(2:80),string,7,Hit)
    if (string /= aux2) cycle
  end if

  Do_Cycle = .false.
end do

if (IfTest) write(u6,*) ' Process library label'
Hit = .true.
call Decode(BSLB(2:80),CGTOm,5,Hit)
if (CGTO == '') CGTO = CGTOm

! Here when the basis set label has been identified on the file
! Now decode the CGTO label

lAng = -1
i1 = 1
do i=1,80
  if (CGTO(i:i) == '') exit
  do k=0,lCGTO
    if (CGTO(i:i) == kAng(k)) then
      !read(CGTO(i1:i-1),*) nCGTO(k)
      Str16 = CGTO(i1:i-1)
      read(Str16,'(F16.0)') cg
      nCGTO(k) = nint(cg)
      i1 = i+1
      if (k /= lAng+1) then
        call WarningMessage(2,'RdBsl: Error in contraction label')
        write(u6,*) 'Conflict for ',kAng(k),' shell'
        write(u6,*) 'Erroneous label:',CGTO
        call Abend()
      end if
      lAng = max(lAng,k)
      exit
    end if
  end do
end do
#ifdef _DEMO_
if (lAng > 1) then
  call WarningMessage(2,'Demo version can handle only s and p basis functions')
  call Quit_OnUserError()
end if
#endif
!write(u6,*) ' lAng=',lAng

! Check for size of contracted basis set

lAngm = 0
i1 = 1
do i=1,80
  if (CGTOm(i:i) == '') exit
  do k=0,lCGTO
    if (CGTOm(i:i) == kAng(k)) then
      !read(CGTOm(i1:i-1),*) mCGTO(k)
      Str16 = CGTOm(i1:i-1)
      read(Str16,'(F16.0)') cg
      mCGTO(k) = nint(cg)
      i1 = i+1
      lAngM = max(lAngM,k)
      exit
    end if
  end do
end do
if (IfTest) then
  write(u6,'(2a)') 'Type=',bType
  write(u6,*) 'nCGTO=',(nCGTO(k),k=0,lCGTO)
  write(u6,*) 'mCGTO=',(mCGTO(k),k=0,lCGTO)
end if
!write(u6,*) ' lAngm=',lAngm
if (lAngM < lAng) then
  call WarningMessage(2,'Abend in RdBsl:Too high angular momentum in basis set input')
  call Quit_OnUserError()
end if
lStop = .false.
if (bType == 'ANO') then
  ! Gen.cont.: never more contracted functions than available
  do i=0,lang
    if (IfTest) write(u6,*) 'i,nCGTO(i),mCGTO(i)=',i,nCGTO(i),mCGTO(i)
    if (nCGTO(i) > mCGTO(i)) then
      write(u6,4000) kAng(i),nCGTO(i),mCGTO(i)
      lStop = .true.
    end if
  end do
else if (bType == 'ECP') then
  ! Segmented basis set: never more functions than primitives
  ! (to be checked later (getbs.f)) and never less functions than available
  do i=0,lang
    if (IfTest) write(u6,*) 'i,nCGTO(i),mCGTO(i)=',i,nCGTO(i),mCGTO(i)
    if (nCGTO(i) < mCGTO(i)) then
      write(u6,4100) kAng(i),nCGTO(i),mCGTO(i)
      lStop = .true.
    end if
  end do
end if
if (lStop) then
  call WarningMessage(2,'Abend in RdBsl:Requested basis inconsistent with library')
  call Quit_OnUserError()
end if

return

4000 format(/1x,'Too many CGTOs of ',a,'-type ',I3,' Max=',I3)
4100 format(/1x,'Too few segmented CGTOs of ',a,'-type ',I3,' Min=',I3)

end subroutine Rdbsl
