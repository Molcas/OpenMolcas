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
! Copyright (C) Roland Lindh                                           *
!***********************************************************************

subroutine Mk_RI_Shells(LuRd)
!***********************************************************************
!                                                                      *
!    Objective: To expand the data for the auxiliary functions         *
!                                                                      *
! Called from: RdCtl                                                   *
!                                                                      *
! Calling    : GetBS                                                   *
!                                                                      *
!     Author: Roland Lindh                                             *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, Max_Shells, nCnttp, Shells
use Sizes_of_Seward, only: S
use RICD_Info, only: iRI_Type
use Gateway_Info, only: UnNorm
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuRd
#include "Molcas.fh"
#include "print.fh"
#include "getlnqoe.fh"
integer(kind=iwp) :: BasisTypes(4), i, iAng, ib, iCnttp, iEnd, iEnds, Ierr, iLast3, Indx, iPrint, iRout, iSh, iShll, iSph, iStrt, &
                     j, jShll, lAng, lSTDINP, Lu_lib, mCnttp, mdc, n, nCnt, nCntrc, nn, nPrim, nSet
logical(kind=iwp) :: Hit, IfTest
character(len=256) :: Basis_lib, Fname
character(len=180) :: BSLB, Line, Ref(2)
character(len=80) :: atom, atomb, author, Aux, basis, BSLbl, btype, CGTO
character(len=180), allocatable :: STDINP(:) !CGGn
integer(kind=iwp), external :: IsFreeUnit, StrnLn
character(len=180), external :: Get_Ln, Get_Ln_Quit
character(len=*), parameter :: DefNm = 'basis_library'

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)

! Temporary setup of symmetry information

call mma_allocate(STDINP,mxAtom*2,label='STDINP')
IfTest = .false.
!IfTest = .true.

! Add the auxiliary basis set

BasisTypes(:) = 0
iShll = S%Mx_Shll-1
lSTDINP = 0
mCnttp = nCnttp

! Branch to special loop for reading external RICD basis sets
! since these have a different infrastructure.

if (iRI_Type == 5) then
  !                                                                    *
  !*********************************************************************
  !*********************************************************************
  !                                                                    *
  ! Specially designed loop to read a RICD auxiliary basis set from
  ! an external library.

  Lu_lib = 17
  Lu_lib = IsFreeUnit(Lu_lib)
  call molcas_open(Lu_lib,'RICDLIB')

  do iCnttp=1,mCnttp
    if (dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nVal == 0)) cycle
    mdc = dbsc(iCnttp)%mdci

    Hit = .true.
    call Decode(dbsc(iCnttp)%Bsl_old,atom,1,Hit)
    btype = ' '
    Author = ' '
    basis = ' '
    CGTO = ' '
    Aux = ' '
    if (IfTest) then
      write(u6,*) 'Bsl_Old=',dbsc(iCnttp)%Bsl_old
      write(u6,*) 'Atom=',Atom
    end if

    Indx = index(dbsc(iCnttp)%Bsl_old,' ')
    BSLbl = ' '
    BSLbl(1:Indx-1) = dbsc(iCnttp)%Bsl_old(1:Indx-1)

    ! Find the basis set

    rewind(Lu_lib)

    ! Loop over the basis set library to find the correct label

    if (IfTest) write(u6,*) ' Locate basis set label in library'
    do
      BSLB = Get_Ln_Quit(Lu_lib,0)
      if (Quit_On_Error) then
        iLast3 = StrnLn(BsLbl)
        call WarningMessage(2,'The requested basis set label: '//BsLbl(:iLast3)//';was not found in basis library: RICDLIB')
        call Abend()
      end if

      call UpCase(BSLB)
      if (BSLB(1:1) /= '/') cycle
      if (IfTest) write(u6,*) 'BSLB=',BSLB
      n = index(BSLB,' ')
      do i=n,80
        BSLB(i:i) = '.'
      end do
      Hit = .true.
      call Decode(BSLB(2:80),atomb,1,Hit)
      if (atomb == atom) exit
    end do
    if (IfTest) write(u6,*) 'atomb=',atomb

    ! Now we should have found the correct basis set label!

    nSet = -1

    do while (nSet /= 0)
      Line = Get_Ln(Lu_lib)
      if (IfTest) then
        write(u6,*) 'nSet=',nSet
        write(u6,*) 'Line=',Line
      end if
      call Get_I1(2,lAng)
      if (nSet == -1) call Get_I1(3,nSet)
      if (IfTest) write(u6,*) 'lAng,nSet=',lAng,nSet

      Line = Get_Ln(Lu_lib)
      Line = Get_Ln(Lu_lib)

      nCnttp = nCnttp+1
      if (nCnttp > Mxdbsc) then
        call WarningMessage(2,'Error in Mk_RI_Shells')
        write(u6,*) 'Mk_RI_Shells: Increase Mxdbsc'
        call Abend()
      end if
      if (Show .and. (iPrint >= 6)) then
        write(u6,*)
        write(u6,*)
        write(u6,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLb
        write(u6,'(1X,A)') 'Basis set is read from the workdir.'
      end if

      dbsc(nCnttp)%Bsl = BSLB(2:80)
      dbsc(nCnttp)%Bsl_old = dbsc(nCnttp)%Bsl

      ! Loop over the angular shells

      jShll = iShll
      do iAng=0,lAng
        iShll = iShll+1
        Line = Get_Ln(Lu_lib)
        call Get_I1(1,nPrim)
        call Get_I1(2,nCntrc)
        call Get_I1(3,iSph)
        if (IfTest) then
          write(u6,*) 'iAng=',iAng
          write(u6,*) 'nPrim=',nPrim
          write(u6,*) 'nCntrc=',nCntrc
          write(u6,*) 'iSph=',iSph
        end if

        ! Read Gaussian exponents

        call mma_Allocate(Shells(iShll)%Exp,nPrim,Label='ExpRI')
        Shells(iShll)%nExp = nPrim
        Shells(iShll)%nBasis_C = nCntrc
        iEnd = iStrt-1
        if (nPrim > 0) then
          if (IfTest) write(u6,*) ' Read gaussian exponents'
          call Read_v(Lu_lib,Shells(iShll)%Exp,1,nPrim,1,Ierr)
          if (Ierr /= 0) then
            call WarningMessage(2,'GetBS: Error while reading the exponents')
            call Quit_OnUserError()
          end if
          if (IfTest) write(u6,*) ' Done with exponents'
          if ((iPrint >= 99) .or. IfTest) call RecPrt(' Exponents',' ',Shells(iShll)%Exp,nPrim,1)
        end if
        iStrt = iEnd+1

        ! Read contraction coefficients. Storage of coefficients
        ! for both contracted and uncontracted case.

        call mma_allocate(Shells(iShll)%Cff_c,nPrim,nCntrc,2,Label='Cff_c')
        call mma_allocate(Shells(iShll)%pCff,nPrim,nCntrc,Label='pCff')
        Shells(iShll)%nBasis = nCntrc
        call mma_allocate(Shells(iShll)%Cff_p,nPrim,nPrim,2,Label='Cff_p')
        iEnds = iEnd
        iEnd = iStrt-1
        ! Read contraction coefficients
        ! Observe that the matrix will have nPrim rows and
        ! nCntrc columns
        if (IfTest) write(u6,*) ' Read coefficients'

        ! Read in coeffs. in GC format, as the standard case

        if (IfTest) write(u6,*) ' Standard case'
        if (nPrim*nCntrc > 0) then
          Shells(iShll)%Cff_c(:,:,:) = Zero

          ! Note that we now change the order!!!
          read(Lu_lib,*) ((Shells(iShll)%Cff_c(i,j,2),j=1,nCntrc),i=1,nPrim)
          if (Ierr /= 0) then
            call WarningMessage(2,'GetBS: Error reading coeffs in GC format')
            call Quit_OnUserError()
          end if
          if (IfTest) call RecPrt(' Coeffs',' ',Shells(iShll)%Cff_c(1,1,2),nPrim,nCntrc)
          Shells(iShll)%Cff_c(:,:,1) = Shells(iShll)%Cff_c(:,:,2)
          if (IfTest) call RecPrt(' Coeffs',' ',Shells(iShll)%Cff_c(1,1,1),nPrim,nCntrc)

          ! Put in unit matrix of uncontracted set

          Shells(iShll)%Cff_p(:,:,1) = Zero
          do i=1,nPrim
            Shells(iShll)%Cff_p(i,i,1) = One
          end do

          Shells(iShll)%Cff_p(:,:,2) = Shells(iShll)%Cff_p(:,:,1)
          call Nrmlz(Shells(iShll)%Exp,nPrim,Shells(iShll)%Cff_p(1,1,1),nPrim,iAng)

          Shells(iShll)%Cff_c(:,:,2) = Shells(iShll)%Cff_c(:,:,1)
          call Fix_Coeff(nPrim,nCntrc,Shells(iShll)%Cff_c(1,1,2),Shells(iShll)%Cff_p(1,1,1),'F')
          Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
        end if

        iEnd = iEnds
        if (iSph == 0) then
          Shells(iShll)%Transf = .false.
          Shells(iShll)%Prjct = .false.
        else if (iSph == 1) then
          Shells(iShll)%Transf = .false.
          Shells(iShll)%Prjct = .true.
        else if (iSph == 2) then
          Shells(iShll)%Transf = .true.
          Shells(iShll)%Prjct = .false.
        else
          Shells(iShll)%Transf = .true.
          Shells(iShll)%Prjct = .true.
        end if

        Shells(iShll)%nBasis = Shells(iShll)%nBasis_C
        Shells(iShll)%Aux = .true.

      end do ! iAng

      dbsc(nCnttp)%Aux = .true.
      S%iAngMx = max(S%iAngMx,lAng)

      dbsc(nCnttp)%iVal = jShll+1
      dbsc(nCnttp)%nVal = lAng+1
      dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal

      nCnt = dbsc(iCnttp)%nCntr
      dbsc(nCnttp)%nCntr = nCnt
      dbsc(nCnttp)%mdci = mdc
      dbsc(nCnttp)%Parent_iCnttp = iCnttp
      ! Create a pointer to the actual coordinates.
      dbsc(nCnttp)%Coor => dbsc(iCnttp)%Coor(1:3,1:nCnt)

      S%Mx_Shll = iShll+1
      Max_Shells = S%Mx_Shll
      S%Mx_mdc = mdc

      nSet = nSet-1
      if (nSet /= 0) Line = Get_Ln(Lu_lib)
    end do ! do while (nSet /= 0)

  end do ! iCnttp

  close(Lu_lib)

else

  do iCnttp=1,mCnttp
    if (dbsc(iCnttp)%Frag .or. (dbsc(iCnttp)%nVal == 0)) cycle
    mdc = dbsc(iCnttp)%mdci
    nCnttp = nCnttp+1

    if (nCnttp > Mxdbsc) then
      call WarningMessage(2,'Error in Mk_RI_Shells')
      write(u6,*) 'Mk_RI_Shells: Increase Mxdbsc'
      call Abend()
    end if

    ! Resolve the name of the valence basis and find the name of
    ! the appropriate auxiliary basis set.

    dbsc(nCnttp)%Bsl = dbsc(iCnttp)%Bsl_old

    Hit = .true.
    call Decode(dbsc(nCnttp)%Bsl,atom,1,Hit)
    Hit = .true.
    call Decode(dbsc(nCnttp)%Bsl,btype,2,Hit)
    Hit = .true.
    call Decode(dbsc(nCnttp)%Bsl,author,3,Hit)
    Hit = .true.
    call Decode(dbsc(nCnttp)%Bsl,basis,4,Hit)
    Hit = .true.
    call Decode(dbsc(nCnttp)%Bsl,CGTO,5,Hit)
    Hit = .false.
    call Decode(dbsc(nCnttp)%Bsl,Aux,6,Hit)
    if (.not. Hit) Aux = ' '

    n = index(Atom,' ')-1
    dbsc(nCnttp)%Bsl(1:n+1) = atom(1:n)//'.'
    nn = n+1

    n = index(btype,' ')-1
    dbsc(nCnttp)%Bsl(nn+1:nn+n+5) = btype(1:n)//'.....'

    ! Modify basis set library correctly

    Indx = index(dbsc(nCnttp)%Bsl,' ')
    BSLbl = ' '
    BSLbl(1:Indx-1) = dbsc(nCnttp)%Bsl(1:Indx-1)
    call WhichMolcas(Basis_lib)
    if (Basis_lib(1:1) /= ' ') then
      ib = index(Basis_lib,' ')-1
      if (ib < 1) call SysAbendMsg('rdCtl','Too long PATH to MOLCAS',' ')
      if (iRI_Type == 1) then
        Fname = Basis_lib(1:ib)//'/basis_library/j_Basis'
      else if (iRI_Type == 2) then
        Fname = Basis_lib(1:ib)//'/basis_library/jk_Basis'
      else if (iRI_Type == 3) then
        Fname = Basis_lib(1:ib)//'/basis_library/c_Basis'
      else
        call WarningMessage(2,'Error in Mk_RI_Shells')
        write(u6,*) 'Wrong iRI_Type!'
        write(u6,*) 'iRI_Type=',iRI_Type
        call Abend()
      end if
    else
      if (iRI_Type == 1) then
        Fname = DefNm//'/j_Basis'
      else if (iRI_Type == 2) then
        Fname = DefNm//'/jk_Basis'
      else if (iRI_Type == 3) then
        Fname = DefNm//'/c_Basis'
      else
        call WarningMessage(2,'Error in Mk_RI_Shells')
        write(u6,*) 'Wrong iRI_Type!'
        write(u6,*) 'iRI_Type=',iRI_Type
        call Abend()
      end if
    end if

    if (Show .and. (iPrint >= 6)) then
      write(u6,*)
      write(u6,*)
      write(u6,'(1X,A,I5,A,A)') 'Basis Set ',nCnttp,' Label: ',BSLbl(1:Indx-1)
      write(u6,'(1X,A,A)') 'Basis set is read from library:',Fname(1:index(Fname,' '))
    end if

    jShll = iShll
    dbsc(nCnttp)%Bsl_old = dbsc(nCnttp)%Bsl
    call GetBS(Fname,dbsc(nCnttp)%Bsl,iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,lSTDINP,.false.,.true.,' ')

    dbsc(nCnttp)%Aux = .true.
    dbsc(nCnttp)%Charge = Zero

    if (Show .and. (iPrint >= 6) .and. (Ref(1) /= '') .and. (Ref(2) /= '')) then
      write(u6,'(1x,a)') 'Basis Set Reference(s):'
      if (Ref(1) /= '') write(u6,'(5x,a)') trim(Ref(1))
      if (Ref(2) /= '') write(u6,'(5x,a)') trim(Ref(2))
      write(u6,*)
      write(u6,*)
    end if
    dbsc(nCnttp)%ECP = dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nPP+dbsc(nCnttp)%nM1+dbsc(nCnttp)%nM2 /= 0

    lAng = max(dbsc(nCnttp)%nVal,dbsc(nCnttp)%nSRO,dbsc(nCnttp)%nPrj)-1
    S%iAngMx = max(S%iAngMx,lAng)
    ! No transformation needed for s and p shells
    Shells(jShll+1)%Transf = .false.
    Shells(jShll+1)%Prjct = .false.
    Shells(jShll+2)%Transf = .false.
    Shells(jShll+2)%Prjct = .false.
    dbsc(nCnttp)%pChrg = dbsc(iCnttp)%pChrg
    dbsc(nCnttp)%Fixed = dbsc(iCnttp)%Fixed
    dbsc(nCnttp)%Parent_iCnttp = iCnttp
    dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal+dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nPP

    do iSh=jShll+1,iShll
      Shells(iSh)%nBasis = Shells(iSh)%nBasis_c
      call mma_deallocate(Shells(iShll)%pCff)
      call mma_allocate(Shells(iShll)%pCff,Shells(iSh)%nExp,Shells(iSh)%nBasis,Label='pCff')
      Shells(iShll)%pCff(:,:) = Shells(iShll)%Cff_c(:,:,1)
      Shells(iSh)%Aux = .true.
    end do

    nCnt = dbsc(iCnttp)%nCntr
    dbsc(nCnttp)%nCntr = nCnt
    dbsc(nCnttp)%mdci = mdc
    ! Create a pointer to the actual coordinates of the parent dbsc
    dbsc(nCnttp)%Coor => dbsc(iCnttp)%Coor(1:3,1:nCnt)

    S%Mx_Shll = iShll+1
    Max_Shells = S%Mx_Shll
    S%Mx_mdc = mdc

  end do

end if
!                                                                      *
!***********************************************************************
!                                                                      *
! Add the final DUMMY SHELL!

call Mk_Dummy_Shell()
call mma_deallocate(STDINP)
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Mk_RI_Shells
