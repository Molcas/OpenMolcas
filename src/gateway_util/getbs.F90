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
! Copyright (C) 1990,2020,  Roland Lindh                               *
!               1990, IBM                                              *
!***********************************************************************
!***********************************************************************
!                                                                      *
!    Object: to read basis set Exponents and Contraction Coefficients  *
!            from a library file.                                      *
!            The contraction coefficients will at this point be radial *
!            normalized.                                               *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
!                                                                      *
!            DDName is the path to the library directory               *
!                                                                      *
!                                                                      *
! Calling    : RecPrt, Rdbsl                                           *
!                                                                      *
!     Author: Roland Lindh, IBM Almaden Research Center, San Jose, CA  *
!***********************************************************************

subroutine GetBS(DDname,BSLbl,iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,iSTDINP,L_STDINP,Expert,ExtBasDir)

use Basis_Info, only: dbsc, Extend_Shells, nCnttp, Shells
use DKH_Info, only: iRELMP
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u5, u6

implicit none
#include "Molcas.fh"
character(len=*), intent(in) :: DDname, ExtBasDir
character(len=80), intent(inout) :: BSLbl
integer(kind=iwp), intent(inout) :: iShll, iSTDINP
character(len=180), intent(out) :: Ref(2)
logical(kind=iwp), intent(in) :: UnNorm, L_STDINP, Expert
integer(kind=iwp), intent(in) :: LuRd
integer(kind=iwp), intent(out) :: BasisTypes(4)
character(len=180), intent(in) :: STDINP(MxAtom*2)
#include "angtp.fh"
integer(kind=iwp) :: i, iAdded, iAIMP, iAng, iDominantSet, iEnd, iErr, iFlgOne, iFrst, iMPShll, iNow, iPrevNow, iPrim, iPrint, &
                     iPrSh, iValSh, j, j1, j2, jNow, jPrSh, jValSh, lAng, lUnit, LUQRP, mCGTO(0:iTabMx), mDel, mSOC, mVal, nAdded, &
                     nAIMP, nCGTO(0:iTabMx), nCntrc, nEorb, nPrim, nProj, Nwords
real(kind=wp) :: Coeff, RatioThres
character(len=263) :: Filename
character(len=256) :: Basis_Lib, DirName
character(len=180) :: Line
character(len=80) :: Atom, Filenm, bType
character(len=24) :: Words(2)  ! CGGn
character(len=20) :: MPLbl
logical(kind=iwp) :: Cart(0:iTabMx), Found, Hit, IfTest, inLn1, inLn2, inLn3, isEorb, isFock, UnContracted
real(kind=wp), allocatable :: ExpMerged(:), Temp(:,:)
integer(kind=iwp), external :: Lbl2Nr, StrnLn
character(len=180), external :: Get_Ln
character(len=*), parameter :: DefNm = 'basis_library'
character(len=*), parameter :: KWord(18) = ['END ','VALE','CORE','EXTE','EXCH','1STO','NOPA','NOP1','NOP2','NOP3','NOPF','RESC', &
                                            'RA0H','RA0F','RAIH','MIXE','SOC ','DKSO']

!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
IfTest = .true.
iPrint = 99
#else
IfTest = .false.
iPrint = 5
#endif
if (IfTest) iPrint = 99
dbsc(nCnttp)%FOp = .true.
lAng = 0

Cart(:) = .false.
if (IfTest) write(u6,'(A,A)') 'DDName=',DDName
Line = DDName
call UpCase(Line)
if ((index(Line,'INLINE') /= 0) .and. (index(Line,'EXPERT') == 0)) then
  if (.not. Expert) then
    write(u6,*) 'INLINE keyword is detected'
    write(u6,*) 'EXPERT keyword is not set...'
    write(u6,*) '       We better abort right now!'
    call abend()
  end if
end if
! Locate the first word
iFrst = 1
call NxtWrd(DDName,iFrst,iEnd)
Filenm = DDName(iFrst:iEnd)
DirName = DDName(iFrst:iEnd)
call UpCase(Filenm)
if (index(Filenm,'INLINE') /= 0) then
  inLn1 = .true.
  inLn2 = .true.
  inLn3 = .true.
else if (index(Filenm,'MM') /= 0) then
  dbsc(nCnttp)%IsMM = 1
  inLn1 = .true.
  inLn2 = .false.
  inLn3 = .false.
else
  inLn1 = .false.
  inLn2 = .false.
  inLn3 = .false.
end if
iFrst = iEnd+1
iEnd = len(DDName)
! Check the 2nd and the 3rd field
if (index(DDName(iFrst:iEnd),'INLINE') /= 0) then
  ! Here if T T or F T
  inLn3 = .true.
  ! Locate the 2nd field
  call NxtWrd(DDName,iFrst,iEnd)
  ! Check the 2nd field
  if (index(DDName(iFrst:iEnd),'INLINE') /= 0) then
    inLn2 = .true.
  else
    inLn2 = .false.
  end if
end if
if (IfTest) write(u6,*) inLn1,inLn2,inLn3

if (.not. inLn1) then
  lUnit = 11

  ! Find and decode basis set label

  call Rdbsl(DirName,BSLbl,bType,nCGTO,mCGTO,lAng,Itabmx,lUnit,dbsc(nCnttp)%AtmNr,BasisTypes,ExtBasDir)
  Line = Get_Ln(lUnit)
  Ref(1) = Line
  Line = Get_Ln(lUnit)
  Ref(2) = Line
else
  Hit = .true.
  call Decode(BSLbl,atom,1,Hit)
  dbsc(nCnttp)%AtmNr = Lbl2Nr(atom)
  lUnit = LuRd
  Ref(1) = ''
  Ref(2) = ''
  Hit = .true.
  call Decode(BSLBl(1:80),bType,2,Hit)
  Basis_Lib = ' '
  i = 1
  Basis_Lib = DefNm
  call Find_Basis_Set(Basis_Lib,' ',' ')
  i = StrnLn(Basis_Lib)
  call BasisType(Basis_Lib(1:i)//'/'//bType,0,BasisTypes)
end if
Line(1:3) = bType(1:3)
call UpCase(Line(1:3))
if (Line(1:3) == 'AUX') dbsc(nCnttp)%Aux = .true.
if (IfTest) then
  write(u6,'(A,A)') 'Ref(1):',trim(Ref(1))
  write(u6,'(A,A)') 'Ref(2):',trim(Ref(2))
end if
Uncontracted = BasisTypes(1) == 6
if (dbsc(nCnttp)%IsMM == 1) then
  lAng = 0
  dbsc(nCnttp)%Charge = Zero
  return
end if

! begin parsing options
isEorb = .false.
isFock = .false.
if (L_STDINP .and. inLn1) then  ! CGGn
  iSTDINP = iSTDINP+1           ! CGGn
  Line = STDINP(iSTDINP)        ! CGGn
else                            ! CGGn
  Line = Get_Ln(lUnit)
end if                          ! CGGn
if (Line == 'Options') then
  do
    Line = Get_Ln(lUnit)
    if (Line /= 'EndOptions') then
      if (Line == 'OrbitalEnergies') then
        if (IfTest) write(u6,*) 'Orbital energies are included'
        isEorb = .true.
      else if (Line == 'FockOperator') then
        if (IfTest) write(u6,*) 'Fock operator is included'
        isEorb = .true.
        isFock = .true.
      else if (Line(1:9) == 'Cartesian') then
        Line = Line(10:)
        call LoCase(Line)
        ! Only d or higher shells are tested
        if (index(Line,'all') /= 0) then
          Cart(2:) = .true.
        else
          do i=2,iTabMx
            if (index(Line,AngTp(i)) /= 0) Cart(i) = .true.
          end do
        end if
      else
        write(u6,*) 'Illegal option: ',Line
        call Abend()
      end if
    else
      exit
    end if
  end do
  Line = Get_Ln(lUnit)
end if
! end parsing options
if (IfTest) write(u6,'(A,A)') 'Line=',Line
if (L_STDINP .and. inLn1) then                         ! CGGn
  call Pick_Words(Line,2,Nwords,Words)                 ! CGGn
  if (Nwords /= 2) call Abend()                        ! CGGn
  call Get_dNumber(Words(1),dbsc(nCnttp)%Charge,iErr)  ! CGGn
  if (iErr /= 0) call Abend()                          ! CGGn
  call Get_iNumber(Words(2),lAng,iErr)                 ! CGGn
  if (iErr /= 0) call Abend()                          ! CGGn
else                                                   ! CGGn
  call get_f1(1,dbsc(nCnttp)%Charge)
  if (inLn1) call get_i1(2,lAng)
end if                                                 ! CGGn
if (iPrint >= 99) then
  write(u6,*) 'lAng, Charge=',lAng,dbsc(nCnttp)%Charge
  write(u6,*) ' Start reading valence basis'
end if
if (lAng > iTabMx) then
  write(u6,*) 'GetBS: lAng > iTabMx'
  write(u6,*) 'lAng,iTabMx=',lAng,iTabMx
  call Abend()
end if
! Loop over each shell type (s,p,d,etc....)
iValSh = iShll
dbsc(nCnttp)%nVal = lAng+1
mVal = 0
dbsc(nCnttp)%iVal = iShll+1
do iAng=0,lAng
  if (IfTest) write(u6,*) 'iAng=',iAng
  iShll = iShll+1
  if (iShll > size(Shells)) call Extend_Shells()
  if (IfTest) then
    write(u6,'(A,A)') 'Line=',Line
    write(u6,*) L_STDINP,inLn1
  end if
  if (L_STDINP .and. inLn1) then            ! CGGn
    iSTDINP = iSTDINP+1                     ! CGGn
    Line = STDINP(iSTDINP)                  ! CGGn
    call Pick_Words(Line,2,Nwords,Words)    ! CGGn
    if (Nwords /= 2) call Abend()           ! CGGn
    call Get_iNumber(Words(1),nPrim,iErr)   ! CGGn
    if (iErr /= 0) call Abend()             ! CGGn
    call Get_iNumber(Words(2),nCntrc,iErr)  ! CGGn
    if (iErr /= 0) call Abend()             ! CGGn
    nCGTO(iAng) = 0                         ! CGGn
    mCGTO(iAng) = 0                         ! CGGn
  else                                      ! CGGn
    Line = Get_Ln(lUnit)
    if (IfTest) write(u6,'(A,A)') 'Line=',Line
    call Get_i1(1,nPrim)
    if (inLn1) then
      call Get_i1(2,nCntrc)
      nCGTO(iAng) = 0
      mCGTO(iAng) = 0
    else
      nCntrc = nCGTO(iAng)
    end if
  end if                                    ! CGGn
  if (IfTest) write(u6,*) ' nPrim, nCntrc=',nPrim,nCntrc

  Shells(iShll)%nExp = nPrim
  Shells(iShll)%nBasis_c = nCntrc
  if (Cart(iAng)) then
    Shells(iShll)%Transf = .false.
    Shells(iShll)%Prjct = .false.
  end if
  call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
  ! Read gaussian exponents
  if (nPrim > 0) then
    if (IfTest) write(u6,*) 'Read gaussian exponents'
    call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,Ierr)
    if (Ierr /= 0) then
      call WarningMessage(2,'GetBS: Error while reading the exponents')
      call Quit_OnUserError()
    end if
    if (IfTest) write(u6,*) 'Done with exponents'
  end if
  if (iPrint >= 99) call RecPrt(' Exponents',' ',Shells(iShll)%Exp,nPrim,1)

  ! Storage of coefficients for both contracted and uncontracted case.

  call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
  Shells(iShll)%Cff_c(:,:,1) = Zero
  call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
  Shells(iShll)%nBasis = nCntrc
  call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
  ! Read contraction coefficients
  ! Observe that the matrix will have nPrim rows and nCntrc columns
  if (IfTest) then
    write(u6,'(2A)') ' Type=',bType
    write(u6,*) 'mCGTO(iAng)=',mCGTO(iAng)
    write(u6,*) 'nCGTO(iAng)=',nCGTO(iAng)
    write(u6,*) 'nCntrc=',nCntrc
  end if
  if (IfTest) write(u6,*) ' Read/Process coefficients'

  if ((inLn1 .or. (mCGTO(iAng) == nCntrc)) .or. (nCntrc == 0)) then
    ! Read in coeffs. in GC format, as the standard case
    if (IfTest) write(u6,*) ' Standard case'
    Shells(iShll)%Cff_c(:,:,:) = Zero
    if (UnContracted) then
      do i=1,nPrim
        Shells(iShll)%Cff_c(i,i,1) = One
      end do
    else
      do iPrim=1,nPrim
        call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),iPrim,nCntrc*nPrim,nPrim,Ierr)
        if (Ierr /= 0) then
          call WarningMessage(2,'GetBS: Error reading coeffs in GC format')
          call Quit_OnUserError()
        end if
      end do
    end if

  else

    ! Here, we want nCntrc generally contracted functions (gcf)
    ! of nPrim primitives resulting from the mCGTO(iAng) gcfs
    ! of the library plus the nCntrc-mCGTO(iAng) outermost
    ! primitive functions of the unextended basis set, that is,
    ! of the basis set which excludes any added polarization or
    ! diffuse functions.
    ! (Note that everything is in generally contraction format.)
    ! Example:  library:.6s.3s.  + input:.6s.5s.   -->   result
    !                  x x 0                          x x 0 0 0
    !                  x x 0                          x x 0 0 0
    !                  x x 0                          x x 0 0 0
    !                  x x 0                          x x 1 0 0
    !                  x x 0                          x x 0 1 0
    !                  0 0 1                          0 0 0 0 1

    if (IfTest) write(u6,*) ' Initial GC + outermost primitives'
    if (nCntrc > nPrim) then
      call WarningMessage(2,'Number of contracted more than the number of primitive: correct the basis set label!')
      call Quit_OnUserError()
    end if
    call mma_allocate(Temp,nPrim,max(nCntrc,mCGTO(iAng)),Label='Temp')
    Temp(:,:) = Zero
    ! read the block in the library as it is
    do iPrim=1,nPrim
      call Read_v(lUnit,Temp,iPrim,nPrim*mCGTO(iAng),nPrim,Ierr)
      if (Ierr /= 0) then
        call WarningMessage(2,'GetBS: Error reading the block')
        call Quit_OnUserError()
      end if
    end do

    ! Order the exponents

    call OrdExp1(nPrim,Shells(iShll)%Exp,mCGTO(iAng),Temp)

    ! identify the presence of added polarization and diffusse functions;
    iAdded = 0
    ! Examine all the contracted functions starting with the last
    Outer: do jNow=mCGTO(iAng),1,-1
      iFlgOne = 0
      ! Examine the primitives
      do iNow=1,nPrim
        Coeff = Temp(iNow,jNow)
        ! Stop if it is obvious that this is not a diffuse
        ! or polarization function.
        if ((Coeff /= Zero) .and. (Coeff /= One)) exit Outer
        if (Coeff == One) then
          if (iFlgOne == 1) cycle Outer
          iFlgOne = 1
        end if
      end do
      iAdded = iAdded+1
      if (IfTest) write(u6,*) 'function',jNow,' is an added one'
    end do Outer

    nAdded = iAdded
    if (nAdded == mCGTO(iAng)) nAdded = 0
    if (IfTest) write(u6,*) ' nAdded=',nAdded
    if (nAdded > 0) then
      ! shift the added polarization and diffuse functions to the right
      do jNow=1,nAdded
        j1 = nCntrc-jNow+1
        j2 = mCGTO(iAng)-jNow+1
        do iNow=1,nPrim
          Temp(iNow,j1) = Temp(iNow,j2)
        end do
      end do
    end if
    ! insert/append the outermost primitives (in GC format)
    do jNow=mCGTO(iAng)+1-nAdded,nCntrc-nAdded
      if (IfTest) write(u6,*) 'jNow=',jNow
      Temp(:,jNow) = Zero
      j = jNow-(mCGTO(iAng)-nAdded)
      iPrevNow = nPrim-nAdded-(nCntrc-mCGTO(iAng))
      iNow = iPrevNow+j
      Temp(iNow,jNow) = One
    end do
    Shells(iShll)%Cff_c(:,:,1) = Temp(:,1:nCntrc)
    call mma_deallocate(Temp)
  end if

  if (IfTest) write(u6,*) ' Done! Now Process.'

  ! Order the exponents

  call OrdExp(nPrim,Shells(iShll)%Exp,nCntrc,Shells(iShll)%Cff_c(1,1,1))
  if (nPrim*nCntrc /= 0) mVal = mVal+1

  ! Decontract if integrals required in the primitive basis

  if (nPrim /= 0) then
    Shells(iShll)%Cff_p(:,:,1) = Zero
    do i=1,nPrim
      Shells(iShll)%Cff_p(i,i,1) = One
    end do

    ! Save the contraction coefficients once more after the coefficients.
    ! The second set will not be normalized!

    Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
    Shells(iShll)%Cff_c(:,:,2) = Shells(iShll)%Cff_c(:,:,1)
    Shells(iShll)%Cff_p(:,:,2) = Shells(iShll)%Cff_p(:,:,1)

    ! The normalization coefficients are assumed to be for
    ! normalized Gaussians. In Nrmlz the contraction coefficients are
    ! multiplied with the normalization coefficient of each primitive
    ! Gaussian. The contracted Gaussian are then normalized with respect
    ! the radial overlap.

    if (.not. UnNorm) then
      call Nrmlz(Shells(iShll)%Exp,nPrim,Shells(iShll)%Cff_c(1,1,1),nCntrc,iAng)
      call Nrmlz(Shells(iShll)%Exp,nPrim,Shells(iShll)%Cff_p(1,1,1),nPrim,iAng)
    end if

    if (iPrint >= 99) then
      nPrim = Shells(iShll)%nExp
      nCntrc = Shells(iShll)%nBasis_C
      call RecPrt(' Coefficients (normalized)',' ',Shells(iShll)%Cff_c(1,1,1),nPrim,nCntrc)
      call RecPrt(' Coefficients (unnormalized)',' ',Shells(iShll)%Cff_c(1,1,2),nPrim,nCntrc)
    end if
  end if
  if (nPrim == 0) cycle
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  ! Begin read orbital energies

  call mma_allocate(Shells(iShll)%FockOp,nCntrc,nCntrc,Label='FockOp')
  Shells(iShll)%nFockOp = nCntrc
  Shells(iShll)%FockOp(:,:) = Zero

  if (isFock) then
    !dbsc(nCnttp)%FOp = dbsc(nCnttp)%FOp .and. .true.
    Line = Get_Ln(lUnit)
    call Get_i1(1,nEorb)
    call mma_allocate(Temp,nEorb,nEorb,Label='Temp')
    do i=1,nEorb
      Line = Get_Ln(lUnit)
      call Get_F(1,Temp(1,i),nEOrb)
    end do
    Shells(iShll)%FockOp(1:min(nEorb,nCntrc),1:min(nEorb,nCntrc)) = Temp(1:min(nEorb,nCntrc),1:min(nEorb,nCntrc))
    call mma_deallocate(Temp)
#   ifdef _DEBUGPRINT_
    call RecPrt('Fock',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#   endif
  else if (isEorb) then
    !dbsc(nCnttp)%FOp = dbsc(nCnttp)%FOp .and. .true.
    Line = Get_Ln(lUnit)
    call Get_i1(1,nEorb)
    call mma_allocate(Temp,nEorb,1,Label='Temp')
    if (nEorb > 0) then
      Line = Get_Ln(lUnit)
      call Get_F(1,Temp(1,1),nEorb)
    end if
    do i=1,min(nEOrb,nCntrc)
      Shells(iShll)%FockOp(i,i) = Temp(i,1)
    end do
    call mma_deallocate(Temp)
#   ifdef _DEBUGPRINT_
    call RecPrt('Eorb',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#   endif
  else
    dbsc(nCnttp)%FOp = .false.
#   ifdef _DEBUGPRINT_
    call RecPrt('Empty',' ',Shells(iShll)%FockOp,nCntrc,nCntrc)
#   endif
  end if

  !  End read orbital energies
  !                                                                    *
  !*********************************************************************
  !                                                                    *

end do
if (mVal == 0) dbsc(nCnttp)%nVal = 0
!***********************************************************************
!
! If PAM basis set read the potentials and coefficient!

if (inLn2 .and. (.not. inLn1)) then
  if (IfTest) write(u6,*) ' Close library and start to read from standard input'
  close(lUnit)
  lUnit = u5
end if
if (index(BSLBl,'.PAM.') /= 0) then
  if (IfTest) write(u6,*) ' Process PAM'
  dbsc(nCnttp)%lPAM2 = .true.
  if (iPrint >= 99) write(u6,*) ' Start reading PAMs'
  call GetPAM(lUnit,nCnttp)

  if (inLn3 .and. (.not. inLn2)) then
    close(lUnit)
    lUnit = u5
  end if
end if
!***********************************************************************
!
! If FRAGMENT basis set read the fragment's coordinates, basis sets,
! orbital energies and density matrix

if (inLn2 .and. (.not. inLn1)) then
  if (IfTest) write(u6,*) ' Close library and start to read from standard input'
  close(lUnit)
  lUnit = u5
end if
if (index(BSLBl,'.FRAGMENT.') /= 0) then
  if (IfTest) write(u6,*) ' Process FRAGMENT'
  if (iPrint >= 99) write(u6,*) ' Start reading fragment data'
  call GetFragment(lUnit,nCnttp)

  if (inLn3 .and. (.not. inLn2)) then
    close(lUnit)
    lUnit = u5
  end if
end if
!***********************************************************************
!
! If ECP basis set read the rest!

if (inLn2 .and. (.not. inLn1)) then
  if (IfTest) write(u6,*) ' Close library and start to read from standard input'
  close(lUnit)
  lUnit = u5
end if

nProj = -1
nAIMP = -1
mSOC = -1
if ((index(BSLBl,'.ECP.') /= 0) .or. (index(BSLBl,'.REL.') /= 0)) then
  if (IfTest) write(u6,*) ' Process ECPs/RELs'
  dbsc(nCnttp)%ECP = .true.
  iPrSh = iShll
  if (iPrint >= 99) write(u6,*) ' Start reading ECPs/RELs'
  dbsc(nCnttp)%iPrj = iShll+1
  call GetECP(lUnit,iShll,nProj,UnNorm)
  dbsc(nCnttp)%nPrj = nProj+1

  if (inLn3 .and. (.not. inLn2)) then
    close(lUnit)
    lUnit = u5
  end if

  ! Now read the spectral resolvent basis set
  !Line = GetLn(lUnit)
  Line = Get_Ln(lUnit)
  call UpCase(Line)
  if (index(Line,'SPEC') == 0) then
    call WarningMessage(2,'ERROR: Keyword SPECTRAL expected, offending line : '//Line)
    call Quit_OnUserError()
  end if

  iMPShll = iShll
  do
    Line = Get_Ln(lUnit)
    call UpCase(Line)
    select case (Line(1:4))

      case default
        call WarningMessage(2,' Invalid keyword in GetBS;'//Line)
        call Quit_OnUserError()

      case (KWord(1)) ! END
        exit

      case (KWord(2)) ! VALE
        ! Valence basis set

        if (nAIMP /= -1) then
          call WarningMessage(2,' SR basis set is already defined!')
          call Quit_OnUserError()
        end if
        nAIMP = lAng
        dbsc(nCnttp)%nSRO = nAIMP+1
        dbsc(nCnttp)%iSRO = iShll+1
        jValSh = iValSh
        do iAIMP=0,nAIMP
          iShll = iShll+1
          if (iShll > size(Shells)) call Extend_Shells()
          jValSh = jValSh+1
          call mma_allocate(Shells(iShll)%Exp,Shells(jValSh)%nExp,Label='Exp')
          Shells(iShll)%Exp(:) = Shells(jValSh)%Exp(:)
          Shells(iShll)%nExp = Shells(jValSh)%nExp
          Shells(iShll)%nBasis = 0
        end do

      case (KWord(3)) ! CORE
        ! Core basis set

        if (nAIMP /= -1) then
          call WarningMessage(2,' SR basis set is already defined!')
          call Quit_OnUserError()
        end if
        nAIMP = nProj
        dbsc(nCnttp)%nSRO = nAIMP+1
        dbsc(nCnttp)%iSRO = iShll+1
        jPrSh = iPrSh
        do iAIMP=0,nAIMP
          iShll = iShll+1
          if (iShll > size(Shells)) call Extend_Shells()
          jPrSh = jPrSh+1
          call mma_allocate(Shells(iShll)%Exp,Shells(jPrSh)%nExp,Label='Exp')
          Shells(iShll)%Exp(:) = Shells(jPrSh)%Exp(:)
          Shells(iShll)%nExp = Shells(jPrSh)%nExp
          Shells(iShll)%nBasis = 0
        end do

      case (KWord(4)) ! EXTE
        ! External basis set

        if (nAIMP /= -1) then
          call WarningMessage(2,' SR basis set is already defined!')
          call Quit_OnUserError()
        end if
        !Line = GetLn(lUnit)
        Line = Get_Ln(lUnit)
        call Get_i1(1,nAIMP)
        dbsc(nCnttp)%nSRO = nAIMP+1
        dbsc(nCnttp)%iSRO = iShll+1
        do iAIMP=0,nAIMP
          iShll = iShll+1
          if (iShll > size(Shells)) call Extend_Shells()
          Line = Get_Ln(lUnit)
          call Get_i1(1,nPrim)
          call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
          Shells(iShll)%nExp = nPrim
          Shells(iShll)%nBasis = 0

          if (nPrim > 0) then
            call read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,Ierr)
            if (Ierr /= 0) then
              call WarningMessage(2,'GetBS: Error reading SRO exponents')
              call Quit_OnUserError()
            end if
          end if

        end do

      case (KWord(5)) ! NOPA
        ! Exchange operator

        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,0)

      case (KWord(6)) ! 1STO
        ! 1st order relativistic correction

        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,1)
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,2)
        Line = Get_Ln(lUnit)
        MPLbl = Line(1:20)

      case (KWord(7)) ! NOPA
        ! one-centre no-pair operators

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 0)) call Error()
        IRELMP = 0
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(8)) ! NOP1
        ! one-centre no-pair operators (DK1)

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 1)) call Error()
        IRELMP = 1
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(9)) ! NOP2
        ! one-centre no-pair operators (DK2)

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 2)) call Error()
        IRELMP = 2
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(10)) ! NOP3
        ! one-centre no-pair operators (DK3)

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 3)) call Error()
        IRELMP = 3
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(11)) ! NOPF
        ! one-centre no-pair operators (DK3)

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 4)) call Error()
        IRELMP = 4
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(12)) ! RESC
        ! one-centre RESC operators

        dbsc(nCnttp)%NoPair = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 11)) call Error()
        IRELMP = 11
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(13)) ! RA0H
        ! one-centre ZORA operators

        dbsc(nCnttp)%NoPair = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 21)) call Error()
        IRELMP = 21
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(14)) ! RA0F
        ! one-centre ZORA-FP operators

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 22)) call Error()
        IRELMP = 22
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (Kword(15)) ! RAIH
        ! one-centre IORA operators

        dbsc(nCnttp)%NoPair = .true.
        dbsc(nCnttp)%SODK = .true.
        if ((IRELMP >= 0) .and. (IRELMP /= 23)) call Error()
        IRELMP = 23
        dbsc(nCnttp)%nOpt = ibset(dbsc(nCnttp)%nOpt,3)

      case (KWord(16)) ! MIXE
        ! Mixed basis set (valence + core), with dominance of the
        ! valence basis set, i.e. the valence basis set plus those
        ! exponents in the core primitives not too close to the
        ! valence ones.

        if (nAIMP /= -1) then
          call WarningMessage(2,' SR basis set is already defined!')
          call Quit_OnUserError()
        end if

        ! Threshold for the ratio between consecutive exponents
        ! (if the ratio between two consecutive exponents is lower
        ! than the threshold, the exponent of the core basis set
        ! will be removed)

        Line = Get_Ln(lUnit)
        call Get_F1(1,RatioThres)

        nAIMP = lAng
        dbsc(nCnttp)%iSRO = iShll+1
        dbsc(nCnttp)%nSRO = nAIMP+1
        jValSh = iValSh
        jPrSh = iPrSh

        do iAIMP=0,nAIMP
          iShll = iShll+1
          if (iShll > size(Shells)) call Extend_Shells()

          jValSh = jValSh+1
          nCntrc = Shells(jValSh)%nBasis

          if (iAIMP <= nProj) then
            jPrSh = jPrSh+1

            iDominantSet = 2
            call mma_allocate(ExpMerged,Shells(jPrSh)%nExp+Shells(jValSh)%nExp,Label='ExpMerged')
            call MergeBS(Shells(jPrSh)%Exp,Shells(jPrSh)%nExp,Shells(jValSh)%Exp,Shells(jValSh)%nExp,ExpMerged,Shells(iShll)%nExp, &
                         RatioThres,iDominantSet)
            call mma_allocate(Shells(iShll)%Exp,Shells(iShll)%nExp,Label='Exp')
            Shells(iShll)%Exp(:) = ExpMerged(1:Shells(iShll)%nExp)
            call mma_deallocate(ExpMerged)

          else

            Shells(iShll)%nExp = Shells(jValSh)%nExp
            call mma_allocate(Shells(iShll)%Exp,Shells(iShll)%nExp,Label='Exp')
            Shells(iShll)%Exp(:) = Shells(jValSh)%Exp(:)

          end if

          Shells(iShll)%nBasis = 0

        end do

      case (KWord(17)) ! SOC
        ! SOC basis set

        if (mSOC /= -1) then
          call WarningMessage(2,' SOC basis set is already defined!')
          call Quit_OnUserError()
        end if

        dbsc(nCnttp)%iSOC = iShll+1
        Line = Get_Ln(lUnit)
        call Get_I1(1,mSOC)
        dbsc(nCnttp)%nSOC = mSOC+1
        if (IfTest) write(u6,'(A,I4)') 'dbsc(nCnttp)%nSOC =',dbsc(nCnttp)%nSOC
        if (mSOC >= 0) then
          do iAng=0,mSOC
            if (IfTest) write(u6,'(A,I4)') ' iAng=',iAng
            iShll = iShll+1
            if (iShll > size(Shells)) call Extend_Shells()
            Line = Get_Ln(lUnit)
            call Get_I1(1,nPrim)
            call Get_I1(2,nCntrc)
            call Get_I1(3,mDel)
            dbsc(nCnttp)%kDel(iAng) = mDel
            if (IfTest) write(u6,*) 'nPrim = ',nPrim,' nCntrc = ',nCntrc
            if (IfTest) write(u6,*) 'nDeleted = ',mDel
            call mma_allocate(Shells(iShll)%Exp,nPrim,Label='Exp')
            Shells(iShll)%nExp = nPrim
            Shells(iShll)%nBasis = nCntrc
            if (IfTest) write(u6,*) 'getBS: ishll,nCntrc',ishll,nCntrc
            if (IfTest) write(u6,'(A)') ' Reading Exponents'
            if (nPrim > 0) call Read_v(lUnit,Shells(iShll)%Exp,1,nPrim,1,ierr)
            if (IfTest) call RecPrt('Exponents',' ',Shells(iShll)%Exp,1,nPrim)
            call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
            call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
            Shells(iShll)%nBasis = nCntrc
            call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
            Shells(iShll)%Cff_p(:,:,:) = Zero
            if (IfTest) write(u6,'(A)') ' Reading coefficients'
            do iPrim=1,nPrim
              call Read_v(lUnit,Shells(iShll)%Cff_c(1,1,1),iPrim,nPrim*nCntrc,nPrim,ierr)
            end do

            Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
            Shells(iShll)%Cff_c(:,:,2) = Shells(iShll)%Cff_c(:,:,1)
            Shells(iShll)%Cff_p(:,:,2) = Shells(iShll)%Cff_p(:,:,1)

            if (IfTest) call RecPrt('Coefficients',' ',Shells(iShll)%Cff_c(1,1,1),nPrim,nCntrc)

          end do
        end if

      case (KWord(18)) ! DKSO
        ! Use DKSO on request

        dbsc(nCnttp)%SODK = .true.

    end select
  end do
  if (btest(dbsc(nCnttp)%nOpt,1) .and. btest(dbsc(nCnttp)%nOpt,3)) then
    call WarningMessage(2,' 1st order relativistic correction and no-pair approximation cannot be used simultaneously!')
    call Quit_OnUserError()
  end if
  if (nAIMP >= 0) then
    Basis_Lib = DefNm
    call Find_Basis_Set(Basis_Lib,' ',' ')
    Filename = trim(Basis_Lib)//'/QRPLIB'
    call f_Inquire(Filename,Found)
    if (.not. Found) then
      write(u6,*) 'File '//trim(Filename)//' not found'
      call abend()
    end if
    LUQRP = 33
    call molcas_open(LUQRP,Filename)
    !open(LUQRP,file='QRPLIB',form='formatted')
    call CalcAMt(dbsc(nCnttp)%nOpt,LUQRP,MPLbl,nAIMP,iMPShll+1,nProj,iPrSh+1,real(dbsc(nCnttp)%AtmNr,kind=wp))
    close(LUQRP)
  end if
end if

lAng = max(lAng,nProj,nAIMP)
if (.not. inLn3) close(lUnit)

return

contains

subroutine Error()

  call WarningMessage(2,'Mixing of iRELMP values')
  call Untested('GetBS')

end subroutine Error

end subroutine GetBS
