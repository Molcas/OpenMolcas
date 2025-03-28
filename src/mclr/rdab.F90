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

use MCLR_Data, only: CMO
use MCLR_Data, only: ChDisp, lDisp
use input_mclr, only: nSym, nIsh, nBas, nOrb, Perturbation, McKinley, iMethod, nIsh, ntITri, ntISqr, ntBSqr, nDisp, PT2, ESCF, &
                      nDel, ntPert, ntIsh
use stdalloc, only: mma_allocate
use Definitions, only: u6

implicit none
character(len=8) Label
integer iRC, iOpt, Length, iSym, iS, iDum

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
    ntItri = ntItri+nIsh(iSym)*(nIsh(iSym)+1)/2
    ntIsqr = ntIsqr+nIsh(iSym)*nIsh(iSym)
    ntbSQR = ntbsqr+nbas(isym)**2
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
  nDisp = 0
  do iS=1,nSym
    nDisp = nDisp+lDisp(iS)
  end do
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
  call icopy(nsym,[0],0,ldisp,1)
  ldisp(1) = 1
end if

return

end subroutine RdAB
