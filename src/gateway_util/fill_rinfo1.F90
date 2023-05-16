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

subroutine Fill_rInfo1()

use Basis_Info, only: dbsc, nCnttp, Shells
use Definitions, only: iwp, u6

implicit none
#include "Molcas.fh"
integer(kind=iwp) :: iAng, icnt, iCnttp, jSh, kCof, kExp, krBas, krCnt, krCof, krExp, nExpj
#include "rinfo.fh"

!                                                                      *
!***********************************************************************
!                                                                      *
! Generate stuff for rinfo.fh

krExp = 0
krCof = 0
krBas = 0
krCnt = 0

! Loop over basis sets
do iCnttp=1,nCnttp
  ! Loop over distinct centers
  do icnt=1,dbsc(iCnttp)%nCntr
    krCnt = krCnt+1
    nAngr(krCnt) = dbsc(iCnttp)%nVal-1

    ! Start with s type shells
    jSh = dbsc(iCnttp)%iVal
    do iAng=0,dbsc(iCnttp)%nVal-1

      krBas = krBas+1
      if (krBas > MxAO) then
        call WarningMessage(2,'Too many shells')
        write(u6,*) 'MORE THAN ',MxAO,' SHELLS'
        write(u6,*) 'Increase MxAO in Molcas.fh and recompile the code!'
        call Abend()
      end if
      nExpj = Shells(jSh)%nExp
      nPrimr(krBas) = nExpj
      nBasisr(krBas) = Shells(jSh)%nBasis_C

      if (krExp+nExpj > MxPrim) then
        call WarningMessage(2,'Too many primitives')
        write(u6,*) 'MORE THAN ',MxPrim,' PRIMITIVES'
        write(u6,*) 'Increase MxPrim in rinfo.fh and recompile the code!'
        call Abend()
      end if
      do kExp=1,nExpj
        krExp = krExp+1
        rExp(krExp) = Shells(jSh)%Exp(kExp)
      end do

      ! Pointer to the untouched contraction matrix as after input.

      if (krCof+nExpj*Shells(jSh)%nBasis > MxrCof) then
        call WarningMessage(2,'Too many contraction coefficients')
        write(u6,*) 'MORE THAN ',MxrCof,' CONTRACTION COEFFICIENTS'
        write(u6,*) 'Increase MxrCof in rinfo.fh and recompile the code!'
        call Abend()
      end if
      do kCof=1,Shells(jSh)%nBasis_C
        do kExp=1,nExpj
          krCof = krCof+1
          rCof(krCof) = Shells(jSh)%Cff_c(kExp,kCof,2)
        end do
      end do

      jSh = jSh+1
    end do
  end do
end do
!                                                                      *
!***********************************************************************
!                                                                      *
return

end subroutine Fill_rInfo1
