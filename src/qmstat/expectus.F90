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

subroutine Expectus(QMMethod,HmatOld,Vmat,VpolMat,Smat,MxDim,iVEC,nDim,lEig,iEig,ip_ExpVal)

implicit real*8(a-h,o-z)
#include "numbers.fh"
#include "WrkSpc.fh"
#include "warnings.h"
dimension HmatOld(MxDim), Vmat(MxDim), VpolMat(MxDim), Smat(MxDim)
character QMMethod*5
logical lEig

! Take different path for different QM-method.

if (QMMethod(1:5) == 'RASSI') then

  ! For how many roots are the eigenvalues to be computed.

  if (lEig) then
    nRoots = iEig
  else
    nRoots = nDim
  end if

  ! Loop over roots and compute expectation values according to
  ! well-known formulas.

  nDTri = nDim*(nDim+1)/2
  call GetMem('DenTemp','Allo','Real',iDTmp,nDTri)
  call GetMem('ExpVals','Allo','Real',ip_ExpVal,4*nRoots)
  do iRoot=1,nRoots

    ! Generate density matrix for relevant root.

    call DensiSt(Work(iDTmp),Work(iVEC),iRoot,nDim,nDim)

    ! Expectation values.

    Work(ip_ExpVal+4*(iRoot-1)+0) = Ddot_(nDTri,Work(iDTmp),iOne,HmatOld,iOne)
    Work(ip_ExpVal+4*(iRoot-1)+1) = Ddot_(nDTri,Work(iDTmp),iOne,Vmat,iOne)
    Work(ip_ExpVal+4*(iRoot-1)+2) = Ddot_(nDTri,Work(iDTmp),iOne,Vpolmat,iOne)
    Work(ip_ExpVal+4*(iRoot-1)+3) = Ddot_(nDTri,Work(iDTmp),iOne,Smat,iOne)
  end do
  call GetMem('DenTemp','Free','Real',iDTmp,nDTri)

else if (QMMethod(1:5) == 'SCF  ') then
  ! If it's SCF we are running.

  nDTri = nDim*(nDim+1)/2
  call GetMem('DenTemp','Allo','Real',iDTmp,nDTri)
  call GetMem('ExpVals','Allo','Real',ip_ExpVal,4)
  call Densi_MO(Work(iDTmp),Work(iVEC),1,iEig,nDim,nDim)

  ! Expectation values.

  Work(ip_ExpVal+0) = Ddot_(nDTri,Work(iDTmp),iOne,HmatOld,iOne)
  Work(ip_ExpVal+1) = Ddot_(nDTri,Work(iDTmp),iOne,Vmat,iOne)
  Work(ip_ExpVal+2) = Ddot_(nDTri,Work(iDTmp),iOne,Vpolmat,iOne)
  Work(ip_ExpVal+3) = Ddot_(nDTri,Work(iDTmp),iOne,Smat,iOne)
  call GetMem('DenTemp','Free','Real',iDTmp,nDTri)

else
  ! Shit happens.

  write(6,*)
  write(6,*) ' Now how did this happen, says Expectus!'
  call Quit(_RC_INTERNAL_ERROR_)
end if

! What's you major malfunction, numb nuts!

return

end subroutine Expectus
