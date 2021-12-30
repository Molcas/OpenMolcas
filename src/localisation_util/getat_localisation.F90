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

subroutine GetAt_Localisation(X,nBas,m,XAt,nAtoms,iOpt,nBas_per_Atom,nBas_Start,Norm)

implicit real*8(a-h,o-z)
real*8 X(nBas,m), XAt(nAtoms,*)
integer nBas_per_Atom(nAtoms), nBas_Start(nAtoms)
character*3 Norm  ! 'MAX' or 'FRO'
character*3 myNorm

if ((nBas < 1) .or. (nAtoms < 1)) return

myNorm = Norm
call UpCase(myNorm)

if (iOpt == 1) then
  call dCopy_(nAtoms*m,[0.0d0],0,XAt,1)
  if (myNorm == 'MAX') then
    do j=1,m
      do iAt=1,nAtoms
        i1 = nBas_Start(iAt)
        i2 = i1+nBas_per_Atom(iAt)-1
        do i=i1,i2
          XAt(iAt,j) = max(abs(X(i,j)),XAt(iAt,j))
        end do
      end do
    end do
  else if (myNorm == 'FRO') then
    do j=1,m
      do iAt=1,nAtoms
        i1 = nBas_Start(iAt)
        i2 = i1+nBas_per_Atom(iAt)-1
        do i=i1,i2
          XAt(iAt,j) = XAt(iAt,j)+X(i,j)**2
        end do
        XAt(iAt,j) = sqrt(XAt(iAt,j))
      end do
    end do
  end if
else
  if (m /= nBas) then
    call SysAbendMsg('GetAt_Localisation','Fatal error','m != nBas')
  end if
  call dCopy_(nAtoms*nAtoms,[0.0d0],0,XAt,1)
  if (myNorm == 'MAX') then
    do jAt=1,nAtoms
      j1 = nBas_Start(jAt)
      j2 = j1+nBas_per_Atom(jAt)-1
      do j=j1,j2
        do iAt=1,nAtoms
          i1 = nBas_Start(iAt)
          i2 = i1+nBas_per_Atom(iAt)-1
          do i=i1,i2
            XAt(iAt,jAt) = max(abs(X(i,j)),XAt(iAt,jAt))
          end do
        end do
      end do
    end do
  else if (myNorm == 'FRO') then
    do jAt=1,nAtoms
      j1 = nBas_Start(jAt)
      j2 = j1+nBas_per_Atom(jAt)-1
      do j=j1,j2
        do iAt=1,nAtoms
          i1 = nBas_Start(iAt)
          i2 = i1+nBas_per_Atom(iAt)-1
          do i=i1,i2
            XAt(iAt,jAt) = XAt(iAt,jAt)+X(i,j)**2
          end do
        end do
      end do
    end do
    do jAt=1,nAtoms
      do iAt=1,nAtoms
        XAt(iAt,jAt) = sqrt(XAt(iAt,jAt))
      end do
    end do
  end if
end if

end subroutine GetAt_Localisation
