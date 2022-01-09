!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!***********************************************************************

subroutine Print_OpInfo()

#ifdef _EFP_
use EFP_Module
use EFP
#endif
use External_Centers
use Symmetry_Info, only: nIrrep

implicit real*8(A-H,O-Z)
#include "print.fh"
#include "rmat.fh"
character*72 tempStr
character*14 Format_XF
logical PrintOperators
real*8 A(3)
integer iStb(0:7), jCoSet(8,8)

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
LuWr = 6
!                                                                      *
!***********************************************************************
!                                                                      *
PrintOperators = .false.
PrintOperators = PrintOperators .or. (nEF /= 0)
PrintOperators = PrintOperators .or. (nDMS /= 0)
PrintOperators = PrintOperators .or. (nWel /= 0)
PrintOperators = PrintOperators .or. allocated(XF)
PrintOperators = PrintOperators .or. RMat_On
if (PrintOperators) then
  write(LuWr,*)
  call CollapseOutput(1,'   Operator info:')
  write(LuWr,'(3X,A)') '   --------------'
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nEF /= 0) then
  if (nOrdEF == 0) then
    write(LuWr,'(2X,A,1X,I8)') 'Centers for electric potential option:',nEF
  else if (nOrdEF == 1) then
    write(LuWr,'(2X,A,1X,I8)') 'Centers for electric field option:',nEF
  else if (nOrdEF == 2) then
    write(LuWr,'(2X,A,1X,I8)') 'Centers for electric field gradient and contact option:',nEF
  end if
  do i=1,nEF
    write(LuWr,'(4X,I8,3(1X,F14.8))') i,(EF_Centers(j,i),j=1,3)
  end do
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nDMS /= 0) then
  call RecPrt(' Gauge Origin for diamagnetic shielding',' ',Dxyz,1,3)
  call RecPrt(' Centers for diamagnetic shielding',' ',DMS_Centers,3,nDMS)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nWel /= 0) then
  write(LuWr,*)
  write(LuWr,*) ' Spherical well specification in au'
  write(LuWr,*) ' =================================='
  write(LuWr,*) '   Coeff.      Exp.        R0      '
  do iWel=1,nWel
    write(LuWr,'(3(F10.6,2x))') Wel_Info(3,iWel),Wel_Info(2,iWel),Wel_Info(1,iWel)
  end do
  write(LuWr,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(XF)) then

  if (nPrint(2) < 6) Go To 666
  if (iXPolType > 0) then
    tempStr = '       a(xx)       a(xy)       a(xz)       a(yy)       a(yz)       a(zz)'
  else
    tempStr = ' '
  end if
  write(LuWr,*)
  write(LuWr,*) ' External field specification in au'
  write(LuWr,*) ' =================================='
  if (nOrd_XF == 0) then
    write(LuWr,*) '     x           y           z           Z'//tempStr
  elseif (nOrd_XF == 1) then
    write(LuWr,*) '     x           y           z           Z         my(x)       my(y)       my(z)'//tempStr
  elseif (nOrd_XF == 2) then
    write(LuWr,*) '     x           y           z           Z         my(x)       my(y)       my(z)'// &
                  '       Q(xx)       Q(xy)       Q(xz)       Q(yy)       Q(yz)       Q(zz)'//tempStr
  elseif (nOrd_XF == -1) then
    write(LuWr,*) '     x           y           z '//tempstr
  else
    call WarningMessage(2,'Option not implemented yet!')
    call Abend()
  end if

666 continue

  write(Format_XF,'(A,I2.2,A)') '(',nData_XF,'(F10.6,2x))'
  XnetCharg = 0.0
  do iXF=1,nXF
    A(1:3) = XF(1:3,iXF)
    Charge_iXF = XF(4,iXF)
    iChxyz = iChAtm(A)
    iDum = 0
    call Stblz(iChxyz,nStab_iXF,iStb,iDum,jCoSet)
    if (nPrint(2) >= 6) write(LuWr,Format_XF) (XF(i,iXF),i=1,nData_XF)
    XnetCharg = XnetCharg+dble(nIrrep/nStab_iXF)*Charge_iXF
  end do
  write(LuWr,*)
  write(LuWr,*) ' Net charge from external field: ',XnetCharg
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_On) then
  write(LuWr,*) ' Parameters for radial integration (R-matrix option)'
  write(LuWr,*) ' ==================================================='
  write(LuWr,'(A,G12.5)') '   rmatr     :',RmatR
  write(LuWr,'(A,G12.5)') '   epsabs    :',Epsabs
  write(LuWr,'(A,G12.5)') '   epsrel    :',Epsrel
  write(LuWr,'(A,G12.5)') '   qcoul     :',qCoul
  write(LuWr,'(A,G12.5)') '   dipol(1)  :',dipol(1)
  write(LuWr,'(A,G12.5)') '   dipol(2)  :',dipol(2)
  write(LuWr,'(A,G12.5)') '   dipol(3)  :',dipol(3)
  write(LuWr,'(A,G12.5)') '   epsq      :',epsq
  write(LuWr,'(A,G12.5)') '   bparm     :',bParm
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _EFP_
if (nEFP_fragments /= 0) then
  call EFP_PRINT_BANNER()
  write(LuWr,*)
  write(LuWr,*) ' Specification of Effective Fragment Potentials'
  write(LuWr,*)
  if (Coor_Type == XYZABC_type) then
    write(LuWr,*) 'In XYZABC format'
  elseif (Coor_Type == POINTS_type) then
    write(LuWr,*) 'In Points format'
  elseif (Coor_Type == ROTMAT_type) then
    write(LuWr,*) 'In RotMat format'
  else
    write(LuWr,*) 'Illegal Coor_type:',Coor_Type
    call Abend()
  end if
  do i=1,nEFP_Fragments
    write(LuWr,*)
    write(LuWr,*) 'Fragment:',FRAG_TYPE(i)
    if (Coor_Type == XYZABC_type) then
    elseif (Coor_Type == POINTS_type) then
      do j=1,3
        write(LuWr,'(A10,3F20.10)') ABC(j,i)(1:10),(EFP_Coors((j-1)*3+k,i),k=1,3)
      end do
    elseif (Coor_Type == ROTMAT_type) then
    end if
  end do
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
if (PrintOperators) then
  call CollapseOutput(0,'   Operator info:')
  write(LuWr,*)
end if

return

end subroutine Print_OpInfo
