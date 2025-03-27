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
! Copyright (C) 1997, Anders Bernhardsson                              *
!***********************************************************************

subroutine Niclas(H,coor,LUT)

use Basis_Info
use Center_Info
use Symmetry_Info, only: nIrrep, iChTbl
use stdalloc, only: mma_allocate, mma_deallocate
use Constants, only: Zero
use Definitions, only: wp

implicit real*8(a-h,o-z)
#include "SysDef.fh"
real*8 H(*)
character*40 Label
integer nDeg(200), ldisp(0:7)
integer inddsp(100,0:7)
logical, external :: TF
real*8 Coor(*)
real*8 Dummy(1)
real*8, allocatable :: Htmp(:), Tmp(:)
! Statement functions
itri(i,j) = max(i,j)*(max(i,j)-1)/2+min(i,j)
irec(i,j) = nd*(j-1)+i-1

idsp = 0
call iCOPY(nirrep,[0],0,ldisp,1)
do iIrrep=0,nIrrep-1
  mdc = 0
  do iCnttp=1,nCnttp
    nCnti = dbsc(iCnttp)%nCntr
    do iCnt=1,nCnti
      mdc = mdc+1
      IndDsp(mdc,iIrrep) = idsp
      do iCar=0,2
        iComp = 2**iCar
        if (TF(mdc,iIrrep,iComp)) then
          idsp = idsp+1
          ldisp(iirrep) = ldisp(iirrep)+1
          ndeg(idsp) = nIrrep/dc(mdc)%nStab
        end if
      end do
    end do
  end do
end do

!***********************************************************************
!
! Steady
!
! Make the symmetrized Hessian correct for degenerated geometries
!
!***********************************************************************

nD = 0
do i=0,nIrrep-1
  nD = ldisp(i)+nd
end do
call mma_allocate(TMP,nd**2,Label='Tmp')
call mma_allocate(HTMP,nd**2,Label='Htmp')
Htmp(:) = Zero
ii = 0
iii = 0
do iS=1,Nirrep
  do i=1,ldisp(iS-1)
    do j=1,i
      Tmp(itri(iii+i,iii+j)) = sqrt(real(nDeg(i+iii)*nDeg(j+iii),kind=wp))*H(ii+itri(i,j))
      !write(u6,*) H(ii+itri(i,j)),Tmo(itri(iii+i,iii+j))
    end do
  end do
  ii = ii+ldisp(is-1)*(ldisp(is-1)+1)/2
  iii = iii+ldisp(is-1)
end do

!***********************************************************************
!
! Go
!
!***********************************************************************

call FCOOR(LUT,Coor)
mdc = 0
iPERT = 0
do iCnttp=1,nCnttp
  nCnti = dbsc(iCnttp)%nCntr
  do iCnt=1,nCnti
    mdc = mdc+1

    nCenti = nIrrep/dc(mdc)%nStab

    ndc = 0
    jPERT = 0
    do jCnttp=1,nCnttp
      nCntj = dbsc(jCnttp)%nCntr
      do jCnt=1,nCntj
        ndc = ndc+1

        nCentj = nIrrep/dc(ndc)%nStab
        do iIrrep=0,nIrrep-1
          iDsp = IndDsp(mdc,iIrrep)
          do iCar=0,2
            iComp = 2**iCar
            if (TF(mdc,iIrrep,iComp)) then
              idsp = idsp+1
              jDsp = IndDsp(ndc,iIrrep)
              do jCar=0,2
                jComp = 2**jCar
                if (TF(ndc,iIrrep,jComp)) then
                  jdsp = jdsp+1
                  HE = Tmp(itri(idsp,jdsp))
                  do iCo=0,Ncenti-1
                    do jCo=0,Ncentj-1
                      i = iPert+ico*3+icar+1
                      j = jPert+jco*3+jcar+1
                      kop_m = dc(mdc)%iCoSet(iCo,0)
                      nop_m = nropr(kop_m)
                      kop_n = dc(ndc)%iCoSet(jCo,0)
                      nop_n = nropr(kop_n)
                      riPh = real(iPrmt(nop_m,icomp)*iChTbl(iIrrep,nop_m),kind=wp)/sqrt(real(nCENTI,kind=wp))
                      rjPh = real(iPrmt(nop_n,jcomp)*ichtbl(iirrep,nop_n),kind=wp)/sqrt(real(nCENTJ,kind=wp))
                      Htmp(1+irec(i,j)) = Htmp(1+irec(i,j))+riph*rjph*HE
                    end do ! jco
                  end do ! ico
                end if
              end do ! jcar
            end if
          end do ! icar
        end do ! irrep
        jPert = jpert+ncentj*3
      end do ! jcnt
    end do ! jcnttp
    iPert = ipert+ncenti*3
  end do ! icnt
end do ! icnttp

Label = 'Unsymmetrized Hessian'
write(LUT,'(A)') Label
write(LUT,'(A)') '*BEGIN HESSIAN'
write(LUT,'(A,I5)') '*Number of pert. ',nd
call WRH(LUT,1,[nd],[nd],Htmp,Dummy,0,Label)
write(LUT,'(A)') '*END HESSIAN'

call Put_dArray('FC-Matrix',Htmp,nd**2)

call mma_deallocate(HTMP)
call mma_deallocate(TMP)

return

end subroutine Niclas
