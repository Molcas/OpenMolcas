************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
************************************************************************
      subroutine kindiag(TKIN,TKINTRIA,ndim,evec,eval,breit)
      implicit real*8 (a-h,o-z)
cbs   determines eigenvectors and -values of TKIN
      dimension tkin(ndim,ndim),
     *TKINTRIA((ndim*ndim+ndim)/2),eval(ndim),evec(ndim,ndim)
      logical breit
cbs   move symmetric matrix to triangular matrix
      itria=1
      do irun2=1,ndim
      do irun1=1,irun2
      TKINTRIA(itria)=TKIN(irun1,irun2)
      itria=itria+1
      enddo
      enddo
      do irun2=1,ndim
      do irun1=1,ndim
      evec(irun1,irun2)=0d0
      enddo
      enddo
      do irun1=1,ndim
      evec(irun1,irun1)=1d0
      enddo
cbs   now diagonalize
            CALL Jacob(TKINTRIA,evec,ndim,ndim)
cbs   get the eigenvalues
      do irun=1,ndim
      eval(irun)=TKINTRIA((irun*irun+irun)/2)
      enddo
      if (breit) then
      do irun=1,ndim
      eval(irun)=0d0
      enddo
      endif
cbs   ensure normalization of the vectors.
      do IRUN=1,ndim
      fact=0d0
      do JRUN=1,ndim
      fact=fact+evec(JRUN,IRUN)*evec(JRUN,IRUN)
      enddo
      fact=1d0/sqrt(fact)
      do JRUN=1,ndim
      evec(JRUN,IRUN)=fact*evec(JRUN,IRUN)
      enddo
      enddo
      return
      end
