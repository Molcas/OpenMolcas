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
* Copyright (C) 2015, Luis Manuel Frutos                               *
*               2015, Alessio Valentini                                *
************************************************************************

C This module calculates an external force applied to the system.

      Subroutine extf(ireturn)

#include "real.fh"
#include "constants.fh"
      integer aN, atomNumberx3, nsAtom, efatom1, efatom2, ext
      real*8  gradient(1000), ExtGrad(1000), modgrad(1000)
      real*8  coord(1000), posvect12(3)
      real*8  efmodul, efmodulAU, norm
      real*8  nnewt
      parameter (nnewt = CONV_AU_TO_KJ_/CONST_BOHR_RADIUS_IN_SI_*1D12)
      logical linear
      CHARACTER*180  Key, Line
      CHARACTER*180  Get_Ln
      EXTERNAL       Get_Ln

C get initial values

C      call Mem_Info('EXTF')

      call Get_iScalar('Unique atoms',nsAtom)

      atomNumberx3 = 3 * nsAtom

      call get_darray('GRAD',gradient,atomNumberx3)
      call Get_dArray('Unique Coordinates',coord,atomNumberx3)
* read molcas input

      LuSpool=isfreeunit(21)
      call SpoolInp(LuSpool)

      REWIND(LuSpool)
      call RdNLst(LuSpool,'extf')
  999 CONTINUE
      Key = Get_Ln(LuSpool)
      Line = Key
      call UpCase(Line)
      IF (Line(1:4).EQ.'LINE') GOTO 1100
      IF (Line(1:3).EQ.'END')  GOTO 9000

*>>>>>>>>>>>>>>>>>>>> FORCe *<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 1100 CONTINUE
      linear=.true.
      write(6,*) 'Linear forces between two atoms selected'
      Line = Get_Ln(LuSpool)
      call Get_I1(1, efatom1)
      Line = Get_Ln(LuSpool)
      call Get_I1(1, efatom2)
      Line = Get_Ln(LuSpool)
      call Get_F1(1, efmodul)
      Line = Get_Ln(LuSpool)
      call Get_I1(1, ext)
      write(6,*) 'atom1:', efatom1
      write(6,*) 'atom2:', efatom2
      write(6,*) 'Force:', efmodul, ' nN'
      if (ext.eq.1) then
         write(6,*) 'Compression force'
      else
         write(6,*) 'Extension force'
      end if
      GOTO 999
*>>>>>>>>>>>>>>>>>>>> END <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 9000 CONTINUE
C>>>>>>>>>>>>>>>>>>>>> LINEAR CODE <<<<<<<<<<<<<<<<<<<<<<<<<<<
      IF (linear) then
      aN=1
      write(6,*) 'Gradient Found:'
      do i=1, atomNumberx3, 3
      write(6,*) aN, gradient(i), gradient(i+1), gradient(i+2)
      aN=aN+1
      end do
      write(6,*) '  '

C from nN to atomic units
      efmodulAU=efmodul/nnewt

C initialize the external gradient vector
      do i=1, atomNumberx3
       ExtGrad(i)=Zero
      end do

C posvect is the position vector from atom2 respect atom1
      posvect12(1)=coord((efatom2-1)*3+1)-coord((efatom1-1)*3+1)
      posvect12(2)=coord((efatom2-1)*3+2)-coord((efatom1-1)*3+2)
      posvect12(3)=coord((efatom2-1)*3+3)-coord((efatom1-1)*3+3)

C creating "norm": the norm of the posvect12 vector
      norm=sqrt((posvect12(1))**2+(posvect12(2))**2+(posvect12(3))**2)
      posvect12(1)=posvect12(1)/norm
      posvect12(2)=posvect12(2)/norm
      posvect12(3)=posvect12(3)/norm

      ExtGrad((efatom1-1)*3+1)=posvect12(1)*efmodulAU
      ExtGrad((efatom1-1)*3+2)=posvect12(2)*efmodulAU
      ExtGrad((efatom1-1)*3+3)=posvect12(3)*efmodulAU
      ExtGrad((efatom2-1)*3+1)=-posvect12(1)*efmodulAU
      ExtGrad((efatom2-1)*3+2)=-posvect12(2)*efmodulAU
      ExtGrad((efatom2-1)*3+3)=-posvect12(3)*efmodulAU
C Extension gradient vector created.

C Checking if it is a compression or extension force. In the former case
C the direction of the ExtGrad vector will be changed
      if (ext.eq.1) then
        ExtGrad((efatom1-1)*3+1)=-ExtGrad((efatom1-1)*3+1)
        ExtGrad((efatom1-1)*3+2)=-ExtGrad((efatom1-1)*3+2)
        ExtGrad((efatom1-1)*3+3)=-ExtGrad((efatom1-1)*3+3)
        ExtGrad((efatom2-1)*3+1)=-ExtGrad((efatom2-1)*3+1)
        ExtGrad((efatom2-1)*3+2)=-ExtGrad((efatom2-1)*3+2)
        ExtGrad((efatom2-1)*3+3)=-ExtGrad((efatom2-1)*3+3)
      end if

      write(6,*) '  '
      aN=1
      write(6,*) 'External Force'
      do i=1, atomNumberx3, 3
      write(6,*) aN, -ExtGrad(i), -ExtGrad(i+1), -ExtGrad(i+2)
      aN=aN+1
      end do
      write(6,*) '  '


C Creating the final modified external gradient vector (modgrad)
      do i=1,atomNumberx3
       modgrad(i)=gradient(i)+ExtGrad(i)
      end do

      end if
C>>>>>>>>>>>>>>>>>>>>>>> end of linear code <<<<<<<<<<<

      write(6,*) '  '
      aN=1
      write(6,*) 'Gradient after force application:'
      do i=1, atomNumberx3, 3
      write(6,*) aN, modgrad(i), modgrad(i+1), modgrad(i+2)
      aN=aN+1
      end do
      write(6,*) '  '

      call put_darray('GRAD',modgrad,atomNumberx3)


      CLOSE(LuSpool)


      ireturn=0
      Return
      End

