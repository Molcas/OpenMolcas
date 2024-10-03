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
! Copyright (C) 1995, Martin Schuetz                                   *
!***********************************************************************

subroutine ChkOrt(iD,OffMx)
!***********************************************************************
!                                                                      *
!     purpose: Check orthogonality of CMOs                             *
!                                                                      *
!     output:                                                          *
!       OffMx   : maximal off diagonal element                         *
!                                                                      *
!***********************************************************************

use Index_Functions, only: nTri_Elem
use InfSCF, only: CMO, MaxBas, nBas, nOrb, nSym, Ovrlp
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero, One
use Definitions, only: wp, iwp, u6

implicit none
integer(kind=iwp), intent(in) :: iD
real(kind=wp), intent(out) :: OffMx
integer(kind=iwp) :: i, iCMO, iDgNo1, ij, iOff, iOffMx, iSym, j, jOffMx, nBs, nOr
real(kind=wp) :: DgNo1
logical(kind=iwp) :: termin
real(kind=wp), allocatable :: Aux(:), OvlS(:)
real(kind=wp), parameter :: OrtThr = 1.0e-9_wp

call mma_allocate(OvlS,MaxBas**2,Label='OvlS')
call mma_allocate(Aux,MaxBas**2,Label='Aux')

ij = 1
iCMO = 1
termin = .false.
do iSym=1,nSym
  OffMx = Zero
  DgNo1 = Zero
  nBs = nBas(iSym)
  nOr = nOrb(iSym)
  if (nOr > 0) then
    call Square(Ovrlp(ij),OvlS,1,nBs,nBs)
    call DGEMM_('N','N',nBs,nOr,nBs, &
                One,OvlS,nBs, &
                CMO(iCMO,iD),nBs, &
                Zero,Aux,nBs)
    call DGEMM_('T','N',nOr,nOr,nBs, &
                One,CMO(iCMO,iD),nBs, &
                Aux,nBs, &
                Zero,OvlS,nOr)
    ! get largest non zero off diag element
    do i=1,nOr
      do j=1,i-1
        iOff = (j-1)*nOr+i
        OffMx = max(OffMx,abs(OvlS(iOff)))
      end do
    end do
    ! get diag element most different from one
    do i=1,nOr
      iOff = (i-1)*nOr+i
      DgNo1 = max(DgNo1,abs(OvlS(iOff)-One))
    end do
    ! check, if orthogonality violated
    if ((OffMx > OrtThr) .or. (DgNo1 > OrtThr)) then
      ! Ooooops...
      !write(u6,*) 'WARNING: reorthonormalizing MOs...',OffMx,DgNo1

      ! try to re-orthonormalize...
      call Orthox(OvlS,CMO(iCMO,iD),nOr,nBs)

      ! and test again...
      OffMx = Zero
      DgNo1 = Zero
      call Square(Ovrlp(ij),OvlS,1,nBs,nBs)
      call DGEMM_('N','N',nBs,nOr,nBs, &
                  One,OvlS,nBs, &
                  CMO(iCMO,iD),nBs, &
                  Zero,Aux,nBs)
      call DGEMM_('T','N',nOr,nOr,nBs, &
                  One,CMO(iCMO,iD),nBs, &
                  Aux,nBs, &
                  Zero,OvlS,nOr)
      ! get largest non zero off diag element
      do i=1,nOr
        do j=1,i-1
          iOff = (j-1)*nOr+i
          OffMx = max(OffMx,abs(OvlS(iOff)))
        end do
      end do
      ! get diag element most different from one
      do i=1,nOr
        iOff = (i-1)*nOr+i
        DgNo1 = max(DgNo1,abs(OvlS(iOff)-One))
      end do
      ! check, if orthogonality violated
      if ((OffMx > OrtThr) .or. (DgNo1 > OrtThr)) then
        termin = .true.
        if (OffMx > OrtThr) then
          ! off diag element too large. Now we have time, since
          ! program will terminate anyway -> Go through matrix
          ! again, this time we want element indices...
          outer: do i=1,nOr
            do j=1,i-1
              iOff = (j-1)*nOr+i
              if (abs(OvlS(iOff)) >= OffMx) then
                ! found
                iOffMx = i
                jOffMx = j
                OffMx = OvlS(iOff)
                ! exit loop
                exit outer
              end if
            end do
          end do outer
          !call WarningMessage(0,'Orthogonality violated')
          write(u6,*) ' iSym =',iSym
          write(u6,*) ' largest off diag element:',' [',iOffMx,',',jOffMx,']',' = ',OffMx
        end if
        if (DgNo1 > OrtThr) then
          ! diag element too different from One
          do i=1,nOr
            iOff = (i-1)*nOr+i
            if (abs(OvlS(iOff)-One) >= DgNo1) then
              ! found
              iDgNo1 = i
              DgNo1 = OvlS(iOff)
              ! exit loop
              exit
            end if
          end do
          !call WarningMessage(0,'Orthogonality violated')
          write(u6,*)
          write(u6,*) ' ***** Orthogonality violated *****'
          write(u6,*) ' iSym =',iSym
          write(u6,*) ' diag element most different from 1.0:',' [',iDgNo1,',',iDgNo1,']',' = ',DgNo1
        end if
      end if
    end if
  end if
  ij = ij+nTri_Elem(nBs)
  iCMO = iCMO+nBs*nOr
end do

call mma_deallocate(Aux)
call mma_deallocate(OvlS)
if (termin) then
  call WarningMessage(0,'Orthogonality cannot be recovered\n Basis set problem???')
  call Abend()
end if

return

end subroutine ChkOrt
