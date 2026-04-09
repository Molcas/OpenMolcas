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

subroutine RdAB()

use Index_Functions, only: nTri_Elem
use MCLR_Data, only: ChDisp, CMO, lDisp
use input_mclr, only: ESCF, iMethod, McKinley, nBas, nDel, nDisp, nIsh, nIsh, nOrb, nSym, ntBSqr, ntIsh, ntISqr, ntITri, ntPert, &
                      Perturbation, PT2
use stdalloc, only: mma_allocate
use Definitions, only: iwp, u6

implicit none
integer(kind=iwp) :: iDum, iOpt, iRC, iSym, Length
character(len=8) :: Label

Perturbation = 'NONE'

if (MCKINLEY) then
  LABEL = 'TDISP'
  iRc = -1
  iOpt = 0
  call RdMck(iRC,iOpt,Label,idum,ntpert,idum)
  if (iRC /= 0) then
    write(u6,*) 'RdAB: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  LABEL = 'PERT'
  iRc = -1
  iOpt = 0
  call cRdMck(iRC,iOpt,Label,idum,Perturbation,idum)
  if (iRC /= 0) then
    write(u6,*) 'RdAB: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
end if

if (iMethod == 1) then
  call Get_dScalar('Last energy',ESCF)
  call Get_iArray('nIsh',nIsh,nSym)
  call Get_iArray('nDel',nDel,nSym)
  !--------------------------------------------------------------------*
  !   Precompute the total sum of variables and size of matrices       *
  !--------------------------------------------------------------------*
  ntIsh = 0
  ntItri = 0
  ntIsqr = 0
  ntBsqr = 0
  Length = 0
  do iSym=1,nSym
    ntIsh = ntIsh+nIsh(iSym)
    ntItri = ntItri+nTri_Elem(nIsh(iSym))
    ntIsqr = ntIsqr+nIsh(iSym)*nIsh(iSym)
    ntBsqr = ntBsqr+nbas(isym)**2
    norb(isym) = nbas(isym)-ndel(isym)
    Length = Length+nBas(iSym)*nOrb(iSym)
  end do
  call mma_allocate(CMO,Length,Label='CMO')
  call Get_dArray_chk('Last orbitals',CMO,Length)
end if

if (MCKINLEY) then
  Label = 'ldisp'
  iRc = -1
  iOpt = 0
  call RdMck(iRc,iOpt,label,idum,ldisp,idum)
  if (iRC /= 0) then
    write(u6,*) 'RdAB: Error reading MCKINT'
    write(u6,'(A,A)') 'Label=',Label
    call Abend()
  end if
  nDisp = sum(lDisp(1:nSym))
  if (ndisp /= 0) then
    Label = 'Chdisp'
    iRc = -1
    iOpt = 0
    call cRdmck(iRc,iOpt,label,idum,ChDisp(1),idum)
    if (iRC /= 0) then
      write(u6,*) 'RdAB: Error reading MCKINT'
      write(u6,'(A,A)') 'Label=',Label
      call Abend()
    end if
  end if
end if

if (PT2) then
  ldisp(1) = 1
  ldisp(2:nsym) = 0
end if

return

end subroutine RdAB
