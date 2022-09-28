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

************************************************************************
* History:                                                             *
* Jie J. Bao on May 09, 2022, created this file                        *
************************************************************************



       Subroutine ShiftDiag(Mat,RShift,lShift,nDim,Digit)

       INTEGER nDim,Digit
       Real*8,DIMENSION(nDim**2)::Mat
       Real*8 RShift
       Logical lShift


       Real*8,DIMENSION(nDim)::RDiag
       Real*8 MaxElem
       INTEGER I,iShift

       DO I=1,nDim
        RDiag(I)=Mat((I-1)*nDim+I)
       END DO

       MaxElem=maxval(RDiag)

*       write(6,*) 'maximum of diagonal elements',MaxElem
*       CALL RecPrt(' ',' ',RDiag,1,nDim)

       IF(abs(MaxElem).lt.Real(Digit,8)) THEN
        lShift=.false.
        RETURN
       ELSE
        lShift=.true.
        IShift=Int(MaxElem,8)/Digit*Digit
        RShift=Real(IShift,8)
*        write(6,*) lShift,IShift,RShift
       END IF

       RETURN
       End Subroutine
