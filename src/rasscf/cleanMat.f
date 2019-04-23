************************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
*                                                                      *
* Copyright (C) 2016, Giovanni Li Manni                                *
************************************************************************
      Subroutine CleanMat(MAT)
************* by G. Li Manni Stuttgart April 2016 *************
*
* MAT: One-body density matrix in MO basis as passed by QMC calculation.

* It could well be an average matrix in SA calculation.
*
* It has following shape:
*        11
*        12 22
*        ** ** 33
*        ** ** ** 44
*        ** ** ** 45 55
*        ** ** ** 46 56 66
*        ** ** ** 47 57 67 77
*        ** ** ** ** ** ** ** 88
*        ** ** ** ** ** ** ** 89  99
*        ** ** ** ** ** ** ** 810 910 1010
*        """""""""""""""""""""""""""""""""""
* mimicking a system with (2 0 0 1 4 3 0 0)  actice orbitals (blocked by Irreps)

*           DMAT will be destroyed and replaced with a positive semi-definite one.
*           N-representability will be preserved.

      implicit none
#include "WrkSpc.fh"
#include "rasdim.fh"
#include "rasscf.fh"
#include "general.fh"
* NACPAR = NAC*(NAC+1)/2 with NAC total number of active orbitals
      integer rc,LEVC,j,i,iTmp,iTmp2
      real*8 MAT(NACPAR), trace
      Character*12 routine
      Parameter (routine = 'CleanMat')

      Call qEnter(routine)

      rc = 0
      If (nacpar .lt. 1) then
       rc= -1
       write(6,*) 'matrix size < 1.'
       Go To 10
      end if

* Allocate memory for eigenvectors and new DMAT
      Call GetMem('EVC','Allo','Real',LEVC,NAC**2)
* Initialize eigenvectors
      Call dCopy_(NAC**2,[0.0d0],0,Work(LEVC),1)
* set eigenvector array to identity for this version of JACOB
      Call dCopy_(NAC,[1.0d0],0,Work(LEVC),NAC+1)

* Step 1: Diagonalize MAT. Eigenvalues are stored in diagonal of MAT
      trace = 0.0d0
      DO I=1,NAC
         trace = trace+ MAT(I*(I+1)/2)
      END DO
      CALL JACOB(MAT,Work(LEVC),NAC,NAC)

#ifdef _DEBUG_
      write(6,*) 'eigenvalues: '
      DO I=1,NAC
         write(6,*) MAT(I*(I+1)/2)
      END DO
      write(6,*) 'eigenvectors: '
      do i = 0, nac -1
       write(6,*) (Work(LEVC+i*NAC+j), j=0,NAC - 1)
      end do
#endif
* Set to zero negative eigenvalue and to TWO values larger than 2.0d0.
      do j = 1, nac
       if(MAT(j*(j+1)/2).gt.2.0d0)     MAT(j*(j+1)/2) = 2.0d0
       if(MAT(j*(j+1)/2).lt.1.0d-12)   MAT(j*(j+1)/2) = 0.0d0
      end do

      trace = 0.0d0
      DO I=1,NAC
         trace = trace + MAT(I*(I+1)/2)
      END DO
      write(6,*) 'trace after removing negative eigenvalues =', trace

* Combine pieced to form the output MAT
* blas routine for square*triangular operation
      Call GetMem('Scr','Allo','Real',iTmp,nac*nac)
      Call GetMem('Scr2','Allo','Real',iTmp2,nac*nac)
      Call dCopy_(nac*nac,[0.0d0],0,Work(iTmp),1)
      Call dCopy_(nac*nac,[0.0d0],0,Work(iTmp2),1)
c     call DTRMM('R','L','N','n',nac,nac,1.0d0,MAT,nac,Work(iTmp))
      do i = 0, nac-1
          do j = 0, nac-1
           work(iTmp+i*nac+j) = Work(LEVC+i*NAC+j)*MAT((I+1)*(I+2)/2)
          end do
      end do
      Call DGEMM_('N','T',nac,nac,nac,
     &            1.0d0,Work(iTmp),nac,Work(LEVC),nac,
     &            0.0d0,Work(iTmp2),nac)
* Copy back to MAT
      do i = 1, nac
        do j = 1, i
         MAT((i-1)*i/2+j) = Work(iTmp2+(i-1)*nac + j-1)
        end do
      end do
#ifdef _DEBUG_
      write(6,*) 'trace after recombination:'
      trace = 0.0d0
      DO I=1,NAC
         trace = trace+ MAT(I*(I+1)/2)
      END DO
#endif
* Release memory
      Call GetMem('Scr2','Free','Real',iTmp2,nac*nac)
      Call GetMem('Scr','Free','Real',iTmp,nac*nac)
      Call GetMem('EVC','Free','Real',LEVC,nac*nac)
****************** Exit ****************
10    Continue
      Call qExit(routine)
      return
      end subroutine
