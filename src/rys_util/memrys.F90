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
! Copyright (C) 1990, Roland Lindh                                     *
!               1990, IBM                                              *
!***********************************************************************

subroutine MemRys(iAnga,MemPrm)

implicit real*8(a-h,o-z)
! This routine will compute the memory requirement of RYS
! Memory requirement is per primitive!
#include "itmax.fh"
#include "print.fh"
#include "FMM.fh"
!gh - stuff for short range integrals
#include "srint.fh"
integer iAnga(4)
! Statement function for canonical index, etc.
nabSz(ixyz) = (ixyz+1)*(ixyz+2)*(ixyz+3)/6-1

iRout = 13
iPrint = nPrint(iRout)
la = iAnga(1)
lb = iAnga(2)
lc = iAnga(3)
ld = iAnga(4)
nRys = (la+lb+lc+ld+2)/2
labMin = nabSz(max(la,lb)-1)+1
labMax = nabSz(la+lb)
lcdMin = nabSz(max(lc,ld)-1)+1
lcdMax = nabSz(lc+ld)
labcd = (labMax-labMin+1)*(lcdMax-lcdMin+1)
if (iPrint >= 99) then
  write(6,*) ' labMin=',labMin
  write(6,*) ' labMax=',labMax
  write(6,*) ' lcdMin=',lcdMin
  write(6,*) ' lcdMax=',lcdMax
end if
MemPrm = 0
! [a0|c0]
MemPrm = MemPrm+labcd
!                                                                      *
!***********************************************************************
!                                                                      *
! For FMM, we only want short-range integrals, using twice the memory
! to store full and long-range components (which are subtracted)
! -same for MOLPRO shortrange

if (FMM_shortrange .or. shortrange) MemPrm = MemPrm+labcd
!                                                                      *
!***********************************************************************
!                                                                      *
nabMax = la+lb
!nabMin = max(la,lb)
ncdMax = lc+ld
!ncdMin = max(lc,ld)
nabcd = (nabMax+1)*(ncdMax+1)
lB10 = max(min(nabMax-1,1),0)
lB01 = max(min(ncdMax-1,1),0)
lB00 = max(min(min(nabMax,ncdMax),1),0)
! Normalization
MemPrm = MemPrm+1
! 2D-Integrals
MemPrm = MemPrm+nabcd*3*nRys
! Coefficients for recurrence relations
MemPrm = MemPrm+3*nRys+3*nRys+3*nRys*(lB10+lB01+lB00)
! Roots
MemPrm = MemPrm+nRys
! The inverse of the arguments
MemPrm = MemPrm+1
! Arguments
MemPrm = MemPrm+1
! Expanded versions of Zeta, ZetInv, Eta, EtaInv, rKapab, rKapcd, P and Q
MemPrm = MemPrm+12
if (iPrint >= 99) then
  write(6,*) ' [e0|f0] integrals   :',labcd
  write(6,*) ' Normalization factor:',1
  write(6,*) ' 2D-integrals        :',nabcd*3*nRys
  write(6,*) ' PAQP vector         :',3*nRys
  write(6,*) ' QCPQ vector         :',3*nRys
  write(6,*) ' B10 coefficients    :',nRys*3*lB10
  write(6,*) ' B00 coefficients    :',nRys*3*lB00
  write(6,*) ' B01 coefficients    :',nRys*3*lB01
  write(6,*) ' Roots               :',nRys
  write(6,*) ' Inverse arguments   :',1
  write(6,*) ' Arguments           :',1
end if

return

end subroutine MemRys
