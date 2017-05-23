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
*--------------------------------------------------------------------------*
* Compute rotation matrix for this specific solvent.                       *
*--------------------------------------------------------------------------*
      Subroutine TransRot(Cordst,i,Rot,xt,yt,zt,Ax,Ay,Az)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "warnings.fh"

      Dimension Rot(3,3),Cordst(MxPut*MxCen,3)

      XO=CORDST(I,1)-AX
      YO=CORDST(I,2)-AY
      ZO=CORDST(I,3)-AZ
      XH=CORDST(I+1,1)-AX
      YH=CORDST(I+1,2)-AY
      ZH=CORDST(I+1,3)-AZ
      XA=CORDST(I+2,1)-AX
      YA=CORDST(I+2,2)-AY
      ZA=CORDST(I+2,3)-AZ
      DELX=(XH+XA)/2.-XO
      DELY=(YH+YA)/2.-YO
      DELZ=(ZH+ZA)/2.-ZO
      DELR=DELX*DELX+DELY*DELY+DELZ*DELZ
      DELR=DELR-1.225449d0
      IF(ABS(DELR).GT..0001)THEN !This is a check of the water geometry.
        Write(6,*)'Molecule',((i-1)/5)+1  !If we enter here, something
        WRITE(6,*)' WARNING IN TRANSROT ', 'delr',delr  !is wrong.
        WRITE(6,*)' O',XO,YO,ZO
        WRITE(6,*)' H',XH,YH,ZH
        WRITE(6,*)' A',XA,YA,ZA
        Call Quit(_RC_GENERAL_ERROR_)
      ENDIF
      XT=XO+.3d0/1.107d0*DELX
      YT=YO+.3d0/1.107d0*DELY
      ZT=ZO+.3d0/1.107d0*DELZ
      ROT(1,3)=(XO-XT)/.3d0
      ROT(2,3)=(YO-YT)/.3d0
      ROT(3,3)=(ZO-ZT)/.3d0
      ROT(1,2)=(XH-XA)/2.86d0
      ROT(2,2)=(YH-YA)/2.86d0
      ROT(3,2)=(ZH-ZA)/2.86d0
      ANORM=1./SQRT(ROT(1,3)**2+ROT(2,3)**2+ROT(3,3)**2)
      ROT(1,3)=ROT(1,3)*ANORM
      ROT(2,3)=ROT(2,3)*ANORM
      ROT(3,3)=ROT(3,3)*ANORM
      ANORM=1./SQRT(ROT(1,2)**2+ROT(2,2)**2+ROT(3,2)**2)
      ROT(1,2)=ROT(1,2)*ANORM
      ROT(2,2)=ROT(2,2)*ANORM
      ROT(3,2)=ROT(3,2)*ANORM
      ROT(1,1)=1.-ROT(1,3)**2-ROT(1,2)**2
      IF(ROT(1,1).LT.0.)ROT(1,1)=0d0
      ROT(1,1)=SQRT(ROT(1,1))
      ROT(2,1)=1.-ROT(2,2)**2-ROT(2,3)**2
      IF(ROT(2,1).LT.0.)ROT(2,1)=0d0
      ROT(2,1)=SQRT(ROT(2,1))
      ROT(3,1)=1.-ROT(3,3)**2-ROT(3,2)**2
      IF(ROT(3,1).LT.0.)ROT(3,1)=0d0
      ROT(3,1)=SQRT(ROT(3,1))

      IFLAG=0
732   Continue
        TAL=ROT(1,1)*ROT(1,2)+ROT(2,1)*ROT(2,2)+ROT(3,1)*ROT(3,2)
        ROT(1,1)=ROT(1,1)-TAL*ROT(1,2)
        ROT(2,1)=ROT(2,1)-TAL*ROT(2,2)
        ROT(3,1)=ROT(3,1)-TAL*ROT(3,2)
        TAL=ROT(1,1)*ROT(1,3)+ROT(2,1)*ROT(2,3)+ROT(3,1)*ROT(3,3)
        ROT(1,1)=ROT(1,1)-TAL*ROT(1,3)
        ROT(2,1)=ROT(2,1)-TAL*ROT(2,3)
        ROT(3,1)=ROT(3,1)-TAL*ROT(3,3)
        A=1d0/SQRT(ROT(1,1)**2+ROT(2,1)**2+ROT(3,1)**2)
        ROT(1,1)=ROT(1,1)*A
        ROT(2,1)=ROT(2,1)*A
        ROT(3,1)=ROT(3,1)*A
        IFLAG=IFLAG+1
        IF(IFLAG.GT.3) THEN
          WRITE(6,*) ' STOP IN TRANSROT'
          Call Quit(_RC_GENERAL_ERROR_)
        ENDIF
      IF(A.GT.10.) GO TO 732

      Return
      End
