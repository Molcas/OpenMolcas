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
      SUBROUTINE TRDTMP(DPT2,NDPT2)
      use definitions, only: iwp, wp
      use constants, only: Zero
      USE Para_Info, ONLY: King
      use stdalloc, only: mma_allocate, mma_deallocate
      use EQSOLV, only: iVecc
      use caspt2_module, only: nSym, nIsh, nAsh, nOrb, nAES, nAshT
      IMPLICIT None


      integer(kind=iwp), intent(in)::  NDPT2
      real(kind=wp), intent(inout):: DPT2(NDPT2)

      real(kind=wp), Allocatable:: DTemp(:,:)
      integer(kind=iwp) idpt, iofdpt, isym, it, itabs, itq, iu, iuabs,
     &                  iuq, na, ndtemp, ni, no
      real(kind=wp) Value

      if(nasht==0) return


      ndtemp=nasht**2
      call mma_allocate(dtemp,nAshT,nAshT,Label='DTemp')
      DTemp(:,:)=Zero
CSVC: trdact is still serial and expects to work on the LUSOLV file,
C which is only on the master node. As long as the MKWW subroutines are
C not functioning in parallel, this part should be done only on the
C master node:
      IF (KING()) call trdact(IVECC,IVECC,dtemp)

      call GADSUM(DTEMP,NDTEMP)

      iofdpt=0
      do isym=1,nsym
        ni=nish(isym)
        na=nash(isym)
        no=norb(isym)
        do it=1,na
          itq=ni+it
          itabs=naes(isym)+it
          do iu=1,na
            iuq=ni+iu
            iuabs=naes(isym)+iu
            value=dtemp(itabs,iuabs)
            idpt=iofdpt+itq+no*(iuq-1)
            DPT2(idpt)=DPT2(idpt)+value
          end do
        end do
        iofdpt=iofdpt+no**2
      end do
      call mma_deallocate(dtemp)

      end SUBROUTINE TRDTMP
