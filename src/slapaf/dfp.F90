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

subroutine DFP(B,nDim,Bd,Delta,Gamma)

implicit real*8(a-h,o-z)
real*8 B(nDim,nDim), Bd(nDim), gamma(nDim), Delta(nDim)
real*8, parameter :: Thr = 1.0D-8

!                                                                      *
!***********************************************************************
!                                                                      *
!#define _DEBUGPRINT_
!                                                                      *
!***********************************************************************
!                                                                      *
!
#ifdef _DEBUGPRINT_
call RecPrt('DFP: B',' ',B,nDim,nDim)
!call RecPrt('DFP: Bd',' ',Bd,1,nDim)
call RecPrt('DFP: Gamma',' ',Gamma,1,nDim)
call RecPrt('DFP: Delta',' ',Delta,1,nDim)
#endif
call DGEMM_('N','N',nDim,1,nDim,1.0d0,B,nDim,Delta,nDim,0.0d0,Bd,nDim)
gd = DDot_(nDim,Gamma,1,Delta,1)
dBd = DDot_(nDim,Delta,1,Bd,1)
#ifdef _DEBUGPRINT_
call RecPrt('DFP: Bd',' ',Bd,1,nDim)
write(6,*) 'gd=',gd
write(6,*) 'dBd=',dBd
write(6,*) 'Thr=',Thr
#endif
if (gd < 0.0d0) then
  call MSP(B,Gamma,Delta,nDim)
else

  do i=1,nDim
    do j=1,nDim
      B(i,j) = B(i,j)+(gamma(i)*gamma(j))/gd-(Bd(i)*Bd(j))/dBd
    end do
  end do
end if

#ifdef _DEBUGPRINT_
call RecPrt('DFP: B',' ',B,nDim,nDim)
#endif

return

end subroutine DFP
