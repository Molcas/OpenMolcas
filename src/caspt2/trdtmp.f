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
      SUBROUTINE TRDTMP(DPT2)
      USE Para_Info, ONLY: King
      IMPLICIT REAL*8 (A-H,O-Z)


#include "rasdim.fh"
#include "caspt2.fh"
#include "eqsolv.fh"
#include "WrkSpc.fh"

      DIMENSION DPT2(*)
      if(nasht.eq.0) return

C Test print:
CTEST      WRITE(*,*)' At beginning of TRDTMP. The DPT2 array:'
CTEST      iofdpt=0
CTEST      do isym=1,nsym
CTEST        WRITE(*,*)' DPT2, symmetry block ISYM=',ISYM
CTEST        no=norb(isym)
CTEST        do i=1,no
CTEST          WRITE(*,'(1x,6f12.6)')(dpt2(iofdpt+i+no*(j-1)),j=1,no)
CTEST        end do
CTEST        iofdpt=iofdpt+no**2
CTEST      end do

      ndtemp=nasht**2
      call getmem('dtemp','allo','real',ldtemp,ndtemp)
      call dcopy_(ndtemp,[0.0d0],0,work(ldtemp),1)
CTEST      WRITE(*,*)' Memory list:'
CTEST      call getmem('Point A','list','real',ldum,ndum)
CTEST      WRITE(*,*)' Calling TRDACT.'
CSVC: trdact is still serial and expects to work on the LUSOLV file,
C which is only on the master node. As long as the MKWW subroutines are
C not functioning in parallel, this part should be done only on the
C master node:
      IF (KING()) call trdact(IVECC,IVECC,work(ldtemp))
      call GADSUM(WORK(LDTEMP),NDTEMP)
CTEST      WRITE(*,*)' Back from TRDACT.'
CTEST      WRITE(*,*)' After TRDACT, the DTEMP array:'
CTEST      do i=1,nasht
CTEST        WRITE(*,'(1x,6f12.6)')
CTEST     &       (work(ldtemp-1+i+(j-1)*nasht),j=1,nasht)
CTEST      end do
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
            value=work(ldtemp-1+itabs+nasht*(iuabs-1))
            idpt=iofdpt+itq+no*(iuq-1)
            dpt2(idpt)=dpt2(idpt)+value
          end do
        end do
        iofdpt=iofdpt+no**2
      end do
C Test print:
CTEST      WRITE(*,*)' At end of TRDTMP. The DPT2 array:'
CTEST      iofdpt=0
CTEST      do isym=1,nsym
CTEST        WRITE(*,*)' DPT2, symmetry block ISYM=',ISYM
CTEST        no=norb(isym)
CTEST        do i=1,no
CTEST          WRITE(*,'(1x,6f12.6)')(dpt2(iofdpt+i+no*(j-1)),j=1,no)
CTEST        end do
CTEST        iofdpt=iofdpt+no**2
CTEST      end do
CTEST      WRITE(*,*)' At end of TRDTMP. Memory list:'
CTEST      call getmem('Point B','list','real',ldum,ndum)
      call getmem('dtemp','free','real',ldtemp,ndtemp)
CTEST      WRITE(*,*)' Leaving TRDTMP.'
      return
      end
