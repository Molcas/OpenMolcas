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

subroutine Build_Mp2Dens_Old(ip_TriDens,ip_Density,CMO,mSym,nOrbAll,nOccAll,Diagonalize)

implicit real*8(a-h,o-z)
integer ip_AOTriBlock
real*8 CMO(*)
integer ipSymRec(8)
integer ipSymTri(8)
integer ipSymLin(8)
integer nOrbAll(8), nOccAll(8), ip_Density(8)
logical Diagonalize
character*30 Note
#include "WrkSpc.fh"
#include "corbinf.fh"
iTri(i,j) = max(i,j)*(max(i,j)-3)/2+i+j

nOrbAllTot = nOrbAll(1)
nOrbAllMax = nOrbAll(1)
lRecTot = nOrbAll(1)*nOrbAll(1)
do iSym=2,mSym
  nOrbAllTot = nOrbAllTot+nOrbAll(iSym)
  nOrbAllMax = max(nOrbAllMax,nOrbAll(iSym))
  lRecTot = lRecTot+nOrbAll(iSym)*nOrbAll(iSym)
end do

! A blockmatrix of the size Orb(iSym) X Orb(iSym) is
! allocated and set to zero
call GetMem('AORecBlock','Allo','Real',ip_AORecBlock,nOrbAllMax**2)
call GetMem('TmpRecBlock','Allo','Real',ip_TmpRecBlock,nOrbAllMax**2)
call GetMem('AOTriBlock','Allo','Real',ip_AOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2)
if (Diagonalize) then
  call GetMem('MOTriBlock','Allo','Real',ip_MOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2)
  call GetMem('EigenVecBlock','Allo','Real',ip_EigenVecBlock,nOrbAllMax*nOrbAllMax)
  call GetMem('EigenValBlock','Allo','Real',ip_EigenValBlock,nOrbAllMax)
  call GetMem('EigenVectors','Allo','Real',ip_EigenVecTot,lRecTot)
  call GetMem('EigenValues','Allo','Real',ip_EigenValTot,nOrbAllTot)
  call GetMem('Energies','Allo','Real',ip_Energies,nOrbAllTot)
  call GetMem('IndT','Allo','Inte',ip_IndT,7*mSym)
  call FZero(Work(ip_Energies),nOrbAllTot)
end if

call FZero(Work(ip_AORecBlock),nOrbAllMax**2)
call FZero(Work(ip_TmpRecBlock),nOrbAllMax**2)
call FZero(Work(ip_AOTriBlock),nOrbAllMax*(nOrbAllMax+1)/2)

! Setup a pointer to a symmetryblock in rect or tri representation.
ipSymRec(1) = 0
ipSymTri(1) = 0
ipSymLin(1) = 0
do iSym=2,8
  ipSymTri(iSym) = ipSymTri(iSym-1)+(nOrbAll(iSym-1))*(nOrbAll(iSym-1)+1)/2
  ipSymRec(iSym) = ipSymRec(iSym-1)+(nOrbAll(iSym-1))**2
  ipSymLin(iSym) = ipSymLin(iSym-1)+(nOrbAll(iSym-1))
end do

do iSym=1,mSym

  if (nOrbAll(iSym) /= 0) then

    ! Set initial values for the eigenvectors to the CMO-matrix.
    ! If we transform the CMO-matrix in the same way as the MO-matrix
    ! needs to be treated to be diagonal we get a transformation from
    ! AO to canonical basis. Also sets the energies to zero since they
    ! have no physical relevance in this basis.
    if (Diagonalize) then
      do i=1,nOrbAll(iSym)**2
        Work(ip_EigenVecBlock+i-1) = CMO(ipSymRec(iSym)+i)
      end do
    end if
    ! Transform the symmetryblock to AO-density
    call DGEMM_('N','N',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),1.0d0,CMO(ipSymRec(iSym)+1),nOrbAll(iSym), &
                Work(ip_Density(iSym)),nOrbAll(iSym),0.0d0,Work(ip_TmpRecBlock),nOrbAll(iSym))
    call DGEMM_('N','T',nOrbAll(iSym),nOrbAll(iSym),nOrbAll(iSym),1.0d0,Work(ip_TmpRecBlock),nOrbAll(iSym),CMO(ipSymRec(iSym)+1), &
                nOrbAll(iSym),0.0d0,Work(ip_AORecBlock),nOrbAll(iSym))
    !call RecPrt('AODens:','(20F8.5)',Work(ip_AORecBlock),nOrb(iSym),nOrb(iSym))
    !call RecPrt('MODens:','(20F8.5)',Work(ip_MORecBlock),nOrb(iSym), nOrb(iSym))
    call Fold_Mat(1,nOrbAll(iSym),Work(ip_AORecBlock),Work(ip_AOTriBlock))
    call dcopy_(nOrbAll(iSym)*(nOrbAll(iSym)+1)/2,Work(ip_AOTriBlock),1,Work(ip_TriDens+ipSymTri(iSym)),1)

    if (Diagonalize) then
      ! Make a normal folded matrix

      index = 0
      do i=1,nOrbAll(iSym)
        do j=1,i
          Work(ip_MOTriBlock+index) = Work(ip_Density(iSym)+j-1+(i-1)*(nOrbAll(iSym)))
          index = index+1
        end do
      end do

      call NIDiag(Work(ip_MOTriBlock),Work(ip_EigenVecBlock),nOrbAll(iSym),nOrbAll(iSym))

      do i=1,nOrbAll(iSym)
        Work(ip_EigenValBlock+i-1) = Work(ip_MOTriBlock+iTri(i,i)-1)
      end do

      call JacOrd3(Work(ip_EigenValBlock),Work(ip_EigenVecBlock),nOrbAll(iSym),nOrbAll(iSym))

#     ifdef _DEBUGPRINT_
      write(6,*) 'The sorted eigenvalues are'
      do i=1,nOrbAll(iSym)
        write(6,*) Work(ip_EigenValBlock+i-1)
      end do
      write(6,*) 'Eigenvectors sorted'
      do i=1,nOrbAll(iSym)**2
        write(6,*) Work(ip_EigenVecBlock+i-1)
      end do
#     endif

      call dcopy_(nOrbAll(iSym)**2,Work(ip_EigenVecBlock),1,Work(ip_EigenVecTot+ipSymRec(iSym)),1)
      call dcopy_(nOrbAll(iSym),Work(ip_EigenValBlock),1,Work(ip_EigenValTot+ipSymLin(iSym)),1)

    end if

  end if
end do

! Put the MP2 natural orbitals we just produced on disk to the file MP2ORB.

if (Diagonalize) then
  LuMP2 = 50
  LuMP2 = IsFreeUnit(LuMP2)
  ! Build the TypeIndex array
  iOff = ip_IndT
  do iSym=1,mSym
    iWork(iOff+0) = nFro(iSym)
    iWork(iOff+1) = nOcc(iSym)
    iWork(iOff+2) = 0
    iWork(iOff+3) = 0
    iWork(iOff+4) = 0
    iWork(iOff+5) = nOrb(iSym)-nFro(iSym)-nOcc(iSym)-nDel(iSym)
    iWork(iOff+6) = nDel(iSym)
    iOff = iOff+7
  end do
  Note = '*  Natural MP2 orbitals'
  call WrVec('MP2ORB',LuMP2,'COEI',mSym,nOrbAll,nOrbAll,Work(ip_EigenVecTot),Work(ip_EigenValTot),Work(ip_Energies), &
             iWork(ip_IndT),Note)
  ! Create a molden-file
  iUHF = 0
  call Molden_Interface(iUHF,'MP2ORB','MD_MP2')

end if

call GetMem('AORecBlock','Free','Real',ip_AORecBlock,nOrbAllMax**2)
call GetMem('TmpRecBlock','Free','Real',ip_TmpRecBlock,nOrbAllMax**2)
call GetMem('AOTriBlock','Free','Real',ip_AOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2)

if (Diagonalize) then
  call GetMem('MOTriBlock','Free','Real',ip_MOTriBlock,nOrbAllMax*(nOrbAllMax+1)/2)
  call GetMem('EigenVecBlock','Free','Real',ip_EigenVecBlock,nOrbAllMax*nOrbAllMax)
  call GetMem('EigenValBlock','Free','Real',ip_EigenValBlock,nOrbAllMax)
  call GetMem('EigenVectors','Free','Real',ip_EigenVecTot,lRecTot)
  call GetMem('EigenValues','Free','Real',ip_EigenValTot,nOrbAllTot)
  call GetMem('Energies','Free','Real',ip_Energies,nOrbAllTot)
  call GetMem('IndT','Free','Inte',ip_IndT,7*mSym)
end if

return
! Avoid unused argument warnings
if (.false.) call Unused_integer_array(nOccAll)

end subroutine Build_Mp2Dens_Old
