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
      Subroutine RBufF_tra2(LUHLFX,W,LL,LBuf,NOTU,KKTU,IST,IADXS,MEMX)
      Implicit real*8(a-h,o-z)
      Integer  MEMX,BLKSZ,NBLCK,NPASS,BPASS,NRST
      Integer  i,j

#include "SysDef.fh"
#include "stdalloc.fh"
      Real*8 W(*)

      Real*8, Allocatable:: BUF(:,:)

      BLKSZ=(NOTU-1)*IADXS+LBuf
      NBLCK=MEMX/BLKSZ

      Call mma_allocate(BUF,BLKSZ,NBLCK,LABEL='BUF')

      BPASS=(LL+LBuf-1)/Lbuf
      NPASS=(BPASS+NBLCK-1)/NBLCK

C      WRITE(6,*) "LL=",LL
C      WRITE(6,*) "LBUF=",LBUF
C      WRITE(6,*) "MEMX=",MEMX
C      WRITE(6,*) "BLKSZ=",BLKSZ
C      WRITE(6,*) "NBLCK=",NBLCK
C      WRITE(6,*) "BPASS=",BPASS
C      WRITE(6,*) "NPASS=",NPASS

      IADX=(KKTU-1)*IADXS
      IST=1

      DO I=1,NPASS-1
         CALL dDAFILE(LUHLFX,2,BUF,SIZE(BUF),IADX)

         DO J=1,NBLCK
            call dcopy_(LBuf,BUF(:,J),1,W(IST),1)
            IST=IST+LBuf
         END DO

         BPASS=BPASS-NBLCK
      END DO

      NRST=(BPASS-1)*BLKSZ+mod(LL,LBuf)
      Call dDAFILE(LUHLFX,2,BUF,NRST,IADX)

      DO J=1,BPASS-1
         call dcopy_(LBuf,BUF(:,J),1,W(IST),1)
         IST=IST+LBuf
      END DO

      NRST=mod(LL,LBuf)

      call dcopy_(NRST,BUF(:,BPASS),1,W(IST),1)

      IST=1

      Call mma_deallocate(BUF)

      Return
      End
