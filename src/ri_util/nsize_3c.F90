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

integer function nSize_3C(kS,lS,nShBf,nShell,nIrrep,iOff,nBas_Aux)
!***********************************************************************
!                                                                      *
!     Compute the size of ({nu,mu}|K) and the offsets to the           *
!     different symmetry blocks.                                       *
!                                                                      *
!***********************************************************************

integer nShBf(0:nIrrep-1,nShell), iOff(3,0:nIrrep-1), nBas_Aux(0:nIrrep-1)

nSize_3C = 0

if (nIrrep == 1) then
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call IZero(iOff,3)
  if (kS /= lS) then
    nK = nShBf(0,kS)
    nL = nShBf(0,lS)
    nKL = nK*nL
  else
    nK = nShBf(0,kS)
    nKL = nK*(nK+1)/2
  end if
  iOff(1,0) = nKL

  nJ = nBas_Aux(0)-1
  iOff(2,0) = nJ
  nSize_3C = nJ*nkl
  !                                                                    *
  !*********************************************************************
  !                                                                    *
else
  !                                                                    *
  !*********************************************************************
  !                                                                    *
  call IZero(iOff,3*nIrrep)
  do klIrrep=0,nIrrep-1
    iOff(3,klIrrep) = nSize_3C

    nKL = 0
    if (kS /= lS) then
      do kIrrep=0,nIrrep-1
        nK = nShBf(kIrrep,kS)
        lIrrep = ieor(klIrrep,kIrrep)
        nL = nShBf(lIrrep,lS)
        nKL = nKL+nK*nL
      end do
    else
      do kIrrep=0,nIrrep-1
        nK = nShBf(kIrrep,kS)
        lIrrep = ieor(klIrrep,kIrrep)
        nL = nShBf(lIrrep,lS)

        if (kIrrep > lIrrep) then
          nKL = nKL+nK*nL
        else if (kIrrep == lIrrep) then
          nKL = nKL+nK*(nK+1)/2
        else
          nKL = nKL+0
        end if

      end do
    end if
    iOff(1,klIrrep) = nKL

    nJ = nBas_Aux(klIrrep)
    if (klIrrep == 0) nJ = nJ-1
    iOff(2,klIrrep) = nJ
    nSize_3C = nSize_3C+nJ*nKL

  end do
  !                                                                    *
  !*********************************************************************
  !                                                                    *
end if
!                                                                      *
!***********************************************************************
!                                                                      *
return

end function nSize_3C
