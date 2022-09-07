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
! Copyright (C) Francesco Aquilante                                    *
!***********************************************************************

subroutine SORT_mat(irc,nDim,nVec,iD_A,nSym,lu_A0,mode,lScr,Scr,Diag)
!***********************************************************************
!
!     Author:  F. Aquilante
!
!***********************************************************************

implicit real*8(a-h,o-z)
integer irc, nSym, lScr
integer iD_A(*), nDim(nSym), nVec(nSym), lu_A0(nSym)
real*8 Scr(lScr)
character*7 mode
character Name_A*6
real*8, optional :: Diag(*)

!write(6,*) 'Mode=',Mode
irc = 0
if (mode == 'GePivot') then  ! returns iD_A
  if (.not. present(Diag)) call Abend()
  is = 1
  ! 19112013VVP: The threshold changed from 1.d-14 to 1.d-12
  Thr = 1.0D-12
  ! The original threshold:
  !Thr = 1.0D-14
  do iSym=1,nSym
    if (nDim(iSym) == 0) cycle
    lu_A = 7
    write(Name_A,'(A4,I2.2)') 'ZMAT',iSym-1
    call DaName_MF_WA(lu_A,Name_A)
    !call RecPrt('Diag',' ',Diag(iS),1,nDim(iSym))
    call get_pivot_idx(Diag(is),nDim(iSym),nVec(iSym),lu_A0(iSym),lu_A,iD_A(is),Scr,lScr,Thr)
    call DaEras(lu_A) ! we do not need it
    !call RecPrt('Diag',' ',Diag(iS),1,nDim(iSym))
    is = is+nDim(iSym)
  end do
else if (mode == 'DoPivot') then ! store full-pivoted UT A-matrix
  is = 1
  do iSym=1,nSym
    if (nVec(iSym) /= 0) then
      lu_A = 7
      write(Name_A,'(A4,I2.2)') 'AMAT',iSym-1
      call DaName_MF_WA(lu_A,Name_A)
      call Pivot_mat(nDim(iSym),nVec(iSym),lu_A0(iSym),lu_A,iD_A(is),Scr,lScr)
      call DaEras(lu_A0(iSym))
      lu_A0(iSym) = lu_A
    end if
    is = is+nDim(iSym)
  end do

else if (mode == 'Restore') then !store squared Q-mat (col. piv.)
  is = 1
  do iSym=1,nSym
    if (nVec(iSym) /= 0) then
      lu_A = 7
      write(Name_A,'(A4,I2.2)') 'QVEC',iSym-1
      call DaName_MF_WA(lu_A,Name_A)
      call Restore_mat(nDim(iSym),nVec(iSym),lu_A0(iSym),lu_A,iD_A(is),Scr,lScr,.false.)
      call DaEras(lu_A0(iSym))
      lu_A0(iSym) = lu_A
    end if
    is = is+nDim(iSym)
  end do

else
  write(6,*) ' SORT_mat: invalid mode! '
  irc = 66
end if

return

end subroutine SORT_mat
