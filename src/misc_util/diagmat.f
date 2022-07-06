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
* Copyright (C) 2022, Jie J. Bao                                       *
************************************************************************
* ****************************************************************
* history:                                                       *
* Jie J. Bao, on Jul 01, 2022, created this file.               *
* ****************************************************************

      Subroutine DiagMat(Mat,EigVal,nDim,nElem)
      use stdalloc, only : mma_allocate, mma_deallocate
      INTEGER nDim,nElem
      Real*8  Mat(nElem),EigVal(nDim)

      Real*8,DIMENSION(:),Allocatable::Scr
      INTEGER nScr

      CALL GetDiagScr(nScr,Mat,EigVal,nDim)

      CALL mma_allocate(Scr,nScr)

      CALL DiagMat_Inner(Mat,EigVal,nDim,nElem,Scr,nScr)

      CALL mma_deallocate(Scr)

      End Subroutine



      Subroutine GetDiagScr(nScr,Mat,EigVal,nDim)
      INTEGER nScr,nDim,INFO
      Real*8 Mat(nDim**2)
      Real*8 EigVal(nDim)
      Real*8 Scr(2)

      CALL DSYEV_('V','U',nDim,Mat,nDim,EigVal,Scr,-1,INFO)
      NScr=INT(Scr(1))
      RETURN
      End Subroutine


      Subroutine DiagMat_Inner(Mat,EigVal,nDim,nElem,Scr,nScr)
      INTEGER nDim,nElem,nScr,INFO
      Real*8  Mat(nElem),EigVal(nDim),Scr(nScr)

      CALL DSYEV_('V','U',nDim,Mat,nDim,EigVal,Scr,nScr,INFO)

      IF(INFO.gt.0) write(6,*)
     &  'DSYEV failed to converge in DiagMat_Inner'
      IF(INFO.lt.0) THEN
       write(6,'(6X,A46)')
     &   'The following element in the matrix is illegal'
       write(6,'(A13,2X,I6,2X,A5,F9.6)')
     &  'Element Index',-INFO,'Value',Mat(-INFO)
      END IF

      RETURN

      End Subroutine

