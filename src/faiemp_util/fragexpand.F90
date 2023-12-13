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
! Copyright (C) Ben Swerts                                             *
!               2020, Roland Lindh                                     *
!***********************************************************************

subroutine FragExpand(LuRd)
!***********************************************************************
!                                                                      *
!    Objective: To expand the data for the fragments and append them   *
!               to regular arrays in Basis_Info, Center_Info           *
!                                                                      *
! Called from: Input                                                   *
!                                                                      *
! Calling    : GetBS                                                   *
!                                                                      *
!     Author: Ben Swerts                                               *
!                                                                      *
!***********************************************************************

use Basis_Info, only: dbsc, Max_Shells, nCnttp, Shells
use Center_Info, only: dc, n_dc
use Sizes_of_Seward, only: S
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: LuRd
#include "Molcas.fh"
integer(kind=iwp), parameter :: storageSize = 200, LineWords = storageSize/8
real(kind=wp) :: eqBasis(LineWords), x1, y1, z1
integer(kind=iwp) :: BasisTypes(4), LenLbl, iAtom, ib, iBas, iCnttp, iCntr, ii, Indx, iSh, iShll, jShll, lAng, Last, LenBSL, &
                     lSTDINP, mCnttp, mdc, ndc
logical(kind=iwp) :: UnNorm
character(len=4) :: label
character(len=13) :: DefNm = 'basis_library'
character(len=180) :: Ref(2)
character(len=storageSize) :: sBasis
character(len=256) :: Basis_lib, Fname
character(len=180), allocatable :: STDINP(:)
! external functions and procedures
integer(kind=iwp), external :: iMostAbundantIsotope, iCLast
real(kind=wp), external :: NucExp, rMass
#ifdef _DEBUGPRINT_
integer(kind=iwp) :: i
#endif

!                                                                      *
!***********************************************************************
!                                                                      *
call mma_allocate(STDINP,mxAtom*2,label='STDINP')
UnNorm = .false.
LenLbl = 0
mdc = dbsc(nCnttp)%mdci+dbsc(nCnttp)%nCntr
BasisTypes(:) = 0
iShll = S%Mx_Shll-1
lSTDINP = 0
mCnttp = nCnttp
#ifdef _DEBUGPRINT_
write(u6,*) 'nCnttp, iShll, mdc = ',nCnttp,iShll,mdc
#endif

#ifdef _DEBUGPRINT_
write(u6,'(A,i6)') 'FragExpand: just before the ''Do iCnttp'''
write(u6,'(A,i6)') 'FragExpand:       mdc          = ',mdc
write(u6,'(A,i6)') 'FragExpand:    mCnttp          = ',mCnttp
write(u6,'(A)') ' dbsc(nCnttp)%mdci   dbsc(nCnttp)%nCntr  nFragType(nCnttp)   nFragCoor(nCnttp)  '
do i=1,mCnttp
  write(u6,'(4(3X,I6,11X))') dbsc(i)%mdci,dbsc(i)%nCntr,dbsc(i)%nFragType,dbsc(i)%nFragCoor
end do
#endif

! Loop over distrinct basis set centers (dbsc)

ndc = 0 ! distinct center index
do iCnttp=1,mCnttp

  ! Skip if this is not a fragment dbsc to expand.

  if (dbsc(iCnttp)%nFragType <= 0) then
    ndc = ndc+dbsc(iCnttp)%nCntr
    cycle
  end if

  ! Loop over the centers associated with this dbsc.

  do iCntr=1,dbsc(iCnttp)%nCntr
    ndc = ndc+1

    ! Loop over the fragments associated with this dbsc.

    do iAtom=1,dbsc(iCnttp)%nFragCoor

      ! Create a new basis set center

      nCnttp = nCnttp+1
      if (nCnttp > Mxdbsc) then
        write(u6,*) ' Increase Mxdbsc'
        call Quit_OnUserError()
      end if

      ! Read the associated basis set in sBasis

      iBas = int(dbsc(iCnttp)%FragCoor(1,iAtom))
      call dcopy_(LineWords,dbsc(iCnttp)%FragType(1,iBas),1,eqBasis,1)
      sBasis = transfer(eqBasis,sBasis) ! ???

      ! Get the basis set directory

      LenBSL = len(sBasis)
      Last = iCLast(sBasis,LenBSL)
      Indx = index(sBasis,'/')
      if (Indx == 0) then
        call WhichMolcas(Basis_lib)
        if (Basis_lib(1:1) /= ' ') then
          ib = index(Basis_lib,' ')-1
          if (ib < 1) call SysAbendMsg('fragexpand','Too long PATH to MOLCAS',' ')
          Fname = Basis_lib(1:ib)//'/basis_library'
        else
          Fname = DefNm
        end if
        Indx = Last+1
        dbsc(nCnttp)%Bsl = trim(sBasis)
      else
        Fname = sBasis(Indx+2:Last)
        if (Fname == ' ') then
          write(u6,*) ' No basis set library specified for'
          write(u6,'(A,A)') 'Fname=',Fname
          call Quit_OnUserError()
        end if
        Fname = adjustl(Fname)
        dbsc(nCnttp)%Bsl = sBasis(1:Indx-1)
      end if
#     ifdef _DEBUGPRINT_
      write(u6,*) 'Setting Bsl(',nCnttp,') to ',dbsc(nCnttp)%Bsl
      write(u6,*) 'Fname = ',Fname
#     endif
      ! Now Fname contains the basis set directory and dbsc(.)%Bsl
      ! contains the basis set label

      jShll = iShll
      dbsc(nCnttp)%mdci = mdc
      call GetBS(Fname,sBasis(1:Indx-1),iShll,Ref,UnNorm,LuRd,BasisTypes,STDINP,lSTDINP,.false.,.true.,' ')
      lAng = max(dbsc(nCnttp)%nVal,dbsc(nCnttp)%nSRO,dbsc(nCnttp)%nPrj)-1
      S%iAngMx = max(S%iAngMx,lAng)
      Shells(jShll+1)%Transf = .false.
      Shells(jShll+1)%Prjct = .false.
      Shells(jShll+2)%Transf = .false.
      Shells(jShll+2)%Prjct = .false.
      dbsc(nCnttp)%Fixed = .true.
      dbsc(nCnttp)%nShells = dbsc(nCnttp)%nVal+dbsc(nCnttp)%nPrj+dbsc(nCnttp)%nSRO+dbsc(nCnttp)%nSOC+dbsc(nCnttp)%nPP
      dbsc(nCnttp)%CntMass = rMass(dbsc(nCnttp)%AtmNr)
      do iSh=jShll+1,iShll
        Shells(iSh)%Frag = .true.
      end do
      if (index(sBasis(1:Indx-1),'6-31G') /= 0) then
        iSh = jShll+3
        if (iSh <= iShll) then
          Shells(iSh)%Prjct = .false.
          Shells(iSh)%Transf = .false.
        end if
      end if
      dbsc(nCnttp)%Frag = .true.

      ! add the coordinates (1 atom / basis set center)
      dbsc(nCnttp)%nCntr = 1

      mdc = mdc+1
      n_dc = max(mdc,n_dc)
      if (mdc > MxAtom) then
        write(u6,*) ' FragExpand: Increase MxAtom'
        write(u6,*) '        MxAtom=',MxAtom
        call Quit_OnUserError()
      end if
      ! get the relative coordinates
      x1 = dbsc(iCnttp)%FragCoor(2,iAtom)
      y1 = dbsc(iCnttp)%FragCoor(3,iAtom)
      z1 = dbsc(iCnttp)%FragCoor(4,iAtom)
      ! make them absolute
      x1 = x1+dbsc(iCnttp)%Coor(1,iCntr)
      y1 = y1+dbsc(iCnttp)%Coor(2,iCntr)
      z1 = z1+dbsc(iCnttp)%Coor(3,iCntr)
#     ifdef _DEBUGPRINT_
      write(u6,'(a,i3,3(a,F12.7))') 'FragExpand: Center ',nCnttp,' Coordinates:  x =',x1,' y=',y1,' z=',z1
#     endif
      ! store them
      call mma_allocate(dbsc(nCnttp)%Coor_Hidden,3,1,Label='dbsc:C')
      dbsc(nCnttp)%Coor => dbsc(nCnttp)%Coor_Hidden(:,:)
      dbsc(nCnttp)%Coor(1,1) = x1
      dbsc(nCnttp)%Coor(2,1) = y1
      dbsc(nCnttp)%Coor(3,1) = z1
      ! store the Mulliken charge
      dbsc(nCnttp)%FragCharge = dbsc(iCnttp)%FragCoor(5,iAtom)
      ! create custom (hopefully) unique label
      LenLbl = index(sBasis,'.')-1
      label = sBasis(1:LenLbl)
      do ii=LenLbl+1,4
        label(ii:ii) = '_'
      end do
#     ifdef _DEBUGPRINT_
      write(u6,'(2A)') 'Label=',label
#     endif
      ! LenIn possible BUG
      dc(mdc)%LblCnt = label
      if (mdc < 10) then
        write(label,'(a3,i1)') '___',mdc
      else if (mdc < 100) then
        write(label,'(a2,i2)') '__',mdc
      else if (mdc < 1000) then
        write(label,'(a1,i3)') '_',mdc
      else
        write(label,'(i4)') mdc
      end if
      dc(mdc)%LblCnt(5:LenIn2) = label
#     ifdef _DEBUGPRINT_
      write(u6,'(2A)') 'Label=',label
      write(u6,'(2A)') 'LblCnt(mdc)=',dc(mdc)%LblCnt
#     endif
      call Chk_LblCnt(dc(mdc)%LblCnt,mdc-1)
      ! store a reference to the originating fragment placeholder
      ! misuse nFragCoor for this purpose: it will not overwrite anything, but
      ! beware of redimensioning this array to anything other than Mxdbsc
      ! To signify this we store the negative value such that we can identify
      ! the that the actual number of centers is 0 and that the corresponding
      ! size of dbsc()%FragCoor is 0 and nothing else.
      dbsc(nCnttp)%nFragCoor = -ndc  ! DO NOT CHANGE THIS!!!!

      if (dbsc(nCnttp)%ExpNuc < Zero) dbsc(nCnttp)%ExpNuc = NucExp(iMostAbundantIsotope(dbsc(nCnttp)%AtmNr))
    end do  ! iAtom
  end do    ! iCntr
end do      ! iCnttp
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _DEBUGPRINT_
write(u6,'(A,i6)') 'FragExpand: After the ''Do iCnttp'''
write(u6,'(A,i6)') 'FragExpand:       mdc          = ',mdc
write(u6,'(A,i6)') 'FragExpand:    nCnttp          = ',nCnttp
write(u6,'(A)') ' dbsc(nCnttp)%mdci   dbsc(nCnttp)%nCntr  nFragType(nCnttp)   nFragCoor(nCnttp)  '
do i=1,nCnttp
  write(u6,'(4(3X,I6,11X))') dbsc(i)%mdci,dbsc(i)%nCntr,dbsc(i)%nFragType,dbsc(i)%nFragCoor
end do
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
S%Mx_Shll = iShll+1
Max_Shells = S%Mx_Shll
!                                                                      *
!***********************************************************************
!                                                                      *
call mma_deallocate(STDINP)

return

end subroutine FragExpand
