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

subroutine Rdbsl(BasDir,BSLbl,type,nCGTO,mCGTO,lAng,lCGTO,lUnit,iAtmNr,BasisTypes,ExtBasDir)
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

implicit real*8(A-H,O-Z)
#include "getlnqoe.fh"
dimension nCGTO(0:lCGTO), mCGTO(0:lCGTO)
character*80 BSLBl, BSLB*180, string, atom, type, author, basis, blank, CGTO, CGTOm, atomb, aux, aux2, Get_Ln_Quit*180
character*(*) BasDir, ExtBasDir
! In case anyone would wonder: basis set path fixed to 100 char max,
! (but no problem in increasing this later on) JR 2005
character*256 BasLoc
#include "angtp.fh"
character*1 kAng(0:iTabMx), Str16*16, FileName*256, TmpString*256
integer StrnLn
external StrnLn
logical lStop, Hit, IfTest, Exist, is_error
integer irecl
integer iLast_JR, iLast1
integer BasisTypes(4)
data IfTest/.false./

#ifdef _DEBUGPRINT_
IfTest = .true.
#endif
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
  write(6,*) ' Enter Rdbsl'
  write(6,'(2a)') ' BsLbl=',BsLbl
  write(6,'(2a)') ' BasDir=',BasDir
  write(6,'(2a)') ' ExtDir=',ExtBasDir
end if
call ResetErrorLine()
call UpCase(BSLBl)

iLast1 = StrnLn(BasDir)
BasLoc = BasDir
if (ExtBasDir /= ' ') then
  Hit = .true.
  call Decode(BSLBl(1:80),type,2,Hit)
end if
call find_basis_set(BasLoc,ExtBasDir,type)
iLast_JR = StrnLn(BasLoc)
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
if (IfTest) write(6,'(1x,a,a)') 'Atom=',atom
iAtmNr = Lbl2Nr(atom)

Hit = .true.
call Decode(BSLBl(1:80),type,2,Hit)
if (IfTest) write(6,'(1x,a,a)') 'Type=',type

Hit = .true.
call Decode(BSLBl(1:80),author,3,Hit)
if (IfTest) write(6,'(1x,a,a)') 'Author=',author

Hit = .true.
call Decode(BSLBl(1:80),basis,4,Hit)
if (IfTest) write(6,'(1x,a,a)') 'Basis=',basis

Hit = .true.
call Decode(BSLBl(1:80),CGTO,5,Hit)
if (IfTest) write(6,'(1x,a,a)') 'CGTO=',CGTO

Hit = .false.
call Decode(BSLBl(1:80),Aux,6,Hit)
if (.not. Hit) Aux = ' '
if (IfTest) write(6,'(1x,a,a)') 'Aux=',Aux

Hit = .false.
call Decode(BSLBl(1:80),Aux2,7,Hit)
if (.not. Hit) Aux2 = ' '
if (IfTest) write(6,'(1x,a,a)') 'Aux2=',Aux2

blank = ' '
if ((CGTO == blank) .and. ((type == 'ECP') .or. (type == 'PSD') .or. (type == 'ANO'))) then
  if (IfTest) write(6,*) ' Early exit'
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

TmpString = type
call Upcase(TmpString)
call LeftAd(BasDir)
call LeftAd(TmpString)
iLast1 = StrnLn(BasDir)
iLast2 = StrnLn(TmpString)
if (IfTest) then
  write(6,'(I3,A)') iLast1,BasDir
  write(6,'(I3,A)') iLast2,TmpString
end if
Filename = BasLoc(1:iLast_JR)//'/'//TmpString(1:iLast2)
TmpString = Author
call Upcase(TmpString)
iLast4 = StrnLn(Filename)
iLast2 = StrnLn(TmpString)
TmpString = Filename(1:iLast4)//'.'//TmpString(1:iLast2)
call f_Inquire(TmpString,Exist)
if (Exist) Filename = TmpString
iLast3 = StrnLn(Filename)
if (IfTest) write(6,'(A,A)') 'Filename=',Filename
call f_Inquire(Filename,Exist)
if (.not. Exist) then
  ! Try to find name in trans.tbl file
  call TransTbl(Filename)
  call f_Inquire(Filename,Exist)
  if (.not. Exist) then
    iLast3 = StrnLn(Filename)
    call WarningMessage(2,'Basis set file '//Filename(1:iLast3)//' does not exist!;;(1) For a valence basis set: check the '// &
                        'spelling of the basis set label and that the basis set file is present in the basis set library '// &
                        'directory.;;(2) For an external auxiliary basis set: check that the basis set file is present in the '// &
                        'appropiate basis set library subdirectory.')
    call Quit_OnUserError()
  end if
end if
! check basistype
call BasisType(Filename,0,BasisTypes)
call molcas_open_ext2(lUnit,FileName,'sequential','formatted',iostat,.false.,irecl,'unknown',is_error)
!open(unit=lUnit,file=Filename,form='FORMATTED',iostat=IOStat)
if (IOStat /= 0) then
  iLast3 = StrnLn(Filename)
  call WarningMessage(2,' Problems opening basis set file '//Filename(1:iLast3))
  call Quit_OnUserError()
end if
rewind(lUnit)

! loop over the basis set library to find the correct label

if (IfTest) write(6,*) ' Locate basis set label in library'
10 BSLB = Get_Ln_Quit(lUnit,0)
if (Quit_On_Error) then
  iLast3 = StrnLn(BsLbl)
  call WarningMessage(2,'The requested basis set label: '//BsLbl(:iLast3)//';'//'was not found in basis library: '//Filename)
  call Abend()
end if
call UpCase(BSLB)
if (BSLB(1:1) /= '/') Go To 10
n = index(BSLB,' ')
do i=n,80
  BSLB(i:i) = '.'
end do
Hit = .true.
call Decode(BSLB(2:80),atomb,1,Hit)
if (atomb /= atom) Go To 10

if (type /= blank) then
  Hit = .true.
  call Decode(BSLB(2:80),string,2,Hit)
  if (string /= type) Go To 10
end if

if (author /= blank) then
  Hit = .true.
  call Decode(BSLB(2:80),string,3,Hit)
  if (string /= author) Go To 10
end if

if (basis /= blank) then
  Hit = .true.
  call Decode(BSLB(2:80),string,4,Hit)
  if (string /= basis) Go To 10
end if

! If a contraction sequence has been specified it must be identical
! to what is in the library file if the basis set type does not
! not allow any other contraction sequence.

if ((CGTO /= blank) .and. (type(1:3) /= 'ANO') .and. (type /= 'ECP') .and. (type /= 'PSD') .and. (type /= 'RYDBERG') .and. &
    (Aux /= 'ECP')) then
  Hit = .true.
  call Decode(BSLB(2:80),string,5,Hit)
  if (string /= CGTO) Go To 10
end if

Hit = .true.
call Decode(BSLB(2:80),string,6,Hit)
if (Aux /= blank) then
  if (string /= aux) Go To 10
else if (string == 'ECP') then
  ! If the library basis (BSLB) has ECP, add it to the BsLbl
  j = 0
  do i=1,5
    j = j+index(BsLbl(j+1:80),'.')
  end do
  BsLbl = BsLbl(1:j)//'ECP'//BsLbl(j+1:)
  Aux = 'ECP'

end if

if (Aux2 /= blank) then
  Hit = .true.
  call Decode(BSLB(2:80),string,7,Hit)
  if (string /= aux2) Go To 10
end if

if (IfTest) write(6,*) ' Process library label'
Hit = .true.
call Decode(BSLB(2:80),CGTOm,5,Hit)
if (CGTO == blank) CGTO = CGTOm

! Here when the basis set label has been identified on the file
! Now decode the CGTO label

lAng = -1
i1 = 1
do i=1,80
  if (CGTO(i:i) == blank) go to 21
  do k=0,lCGTO
    if (CGTO(i:i) == kAng(k)) then
      !read(CGTO(i1:i-1),*) nCGTO(k)
      Str16 = CGTO(i1:i-1)
      read(Str16,'(F16.0)') cg
      nCGTO(k) = nint(cg)
      i1 = i+1
      if (k /= lAng+1) then
        call WarningMessage(2,'RdBsl: Error in contraction label')
        write(6,*) 'Conflict for ',kAng(k),' shell'
        write(6,*) 'Erroneous label:',CGTO
        call Abend()
      end if
      lAng = max(lAng,k)
      Go To 20
    end if
  end do
20 continue
end do
21 continue
#ifdef _DEMO_
if (lAng > 1) then
  call WarningMessage(2,'Demo version can handle only s and p basis functions')
  call Quit_OnUserError()
end if
#endif
!write(6,*) ' lAng=',lAng

! Check for size of contracted basis set

lAngm = 0
i1 = 1
do i=1,80
  if (CGTOm(i:i) == blank) go to 31
  do k=0,lCGTO
    if (CGTOm(i:i) == kAng(k)) then
      !read(CGTOm(i1:i-1),*) mCGTO(k)
      Str16 = CGTOm(i1:i-1)
      read(Str16,'(F16.0)') cg
      mCGTO(k) = nint(cg)
      i1 = i+1
      lAngM = max(lAngM,k)
      go to 30
    end if
  end do
30 continue
end do
31 continue
if (IfTest) then
  write(6,'(2a)') 'Type=',type
  write(6,*) 'nCGTO=',(nCGTO(k),k=0,lCGTO)
  write(6,*) 'mCGTO=',(mCGTO(k),k=0,lCGTO)
end if
!write(6,*) ' lAngm=',lAngm
if (lAngM < lAng) then
  call WarningMessage(2,'Abend in RdBsl:Too high angular momentum in basis set input')
  call Quit_OnUserError()
end if
lStop = .false.
if (type == 'ANO') then
  ! Gen.cont.: never more contracted functions than available
  do i=0,lang
    if (IfTest) write(6,*) 'i,nCGTO(i),mCGTO(i)=',i,nCGTO(i),mCGTO(i)
    if (nCGTO(i) > mCGTO(i)) then
      write(6,4000) kAng(i),nCGTO(i),mCGTO(i)
4000  format(/1x,'Too many CGTOs of ',a,'-type ',I3,' Max=',I3)
      lStop = .true.
    end if
  end do
else if (type == 'ECP') then
  ! Segmented basis set: never more functions than primitives
  ! (to be checked later (getbs.f)) and never less functions than available
  do i=0,lang
    if (IfTest) write(6,*) 'i,nCGTO(i),mCGTO(i)=',i,nCGTO(i),mCGTO(i)
    if (nCGTO(i) < mCGTO(i)) then
      write(6,4100) kAng(i),nCGTO(i),mCGTO(i)
4100  format(/1x,'Too few segmented CGTOs of ',a,'-type ',I3,' Min=',I3)
      lStop = .true.
    end if
  end do
end if
if (lStop) then
  call WarningMessage(2,'Abend in RdBsl:Requested basis inconsistent with library')
  call Quit_OnUserError()
end if

return

end subroutine Rdbsl
