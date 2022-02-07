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
use EFP_Module, only: ABC, Coor_Type, EFP_Coors, FRAG_TYPE, nEFP_fragments, POINTS_type, ROTMAT_type, XYZABC_type
use EFP, only: EFP_PRINT_BANNER
#endif
use External_Centers, only: DMS_Centers, Dxyz, EF_Centers, iXPolType, nData_XF, nDMS, nEF, nOrd_XF, nOrdEF, nWel, nXF, Wel_Info, XF
use Symmetry_Info, only: nIrrep
use Constants, only: Zero
use Definitions, only: wp, iwp, u6

implicit none
#include "print.fh"
#include "rmat.fh"
integer(kind=iwp) :: i, iChxyz, iDum, iPrint, iRout, iStb(0:7), iWel, iXF, j, jCoSet(8,8), nSTab_iXF
#ifdef _EFP_
integer(kind=iwp) :: k
#endif
real(kind=wp) :: A(3), Charge_iXF, XnetCharg
logical(kind=iwp) :: PrintOperators
character(len=72) :: tempStr
character(len=14) :: Format_XF
integer(kind=iwp), external :: iChAtm

!                                                                      *
!***********************************************************************
!                                                                      *
iRout = 2
iPrint = nPrint(iRout)
if (iPrint == 0) return
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
  write(u6,*)
  call CollapseOutput(1,'   Operator info:')
  write(u6,'(3X,A)') '   --------------'
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (nEF /= 0) then
  if (nOrdEF == 0) then
    write(u6,'(2X,A,1X,I8)') 'Centers for electric potential option:',nEF
  else if (nOrdEF == 1) then
    write(u6,'(2X,A,1X,I8)') 'Centers for electric field option:',nEF
  else if (nOrdEF == 2) then
    write(u6,'(2X,A,1X,I8)') 'Centers for electric field gradient and contact option:',nEF
  end if
  do i=1,nEF
    write(u6,'(4X,I8,3(1X,F14.8))') i,(EF_Centers(j,i),j=1,3)
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
  write(u6,*)
  write(u6,*) ' Spherical well specification in au'
  write(u6,*) ' =================================='
  write(u6,*) '   Coeff.      Exp.        R0      '
  do iWel=1,nWel
    write(u6,'(3(F10.6,2x))') Wel_Info(3,iWel),Wel_Info(2,iWel),Wel_Info(1,iWel)
  end do
  write(u6,*)
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (allocated(XF)) then

  if (nPrint(2) >= 6) then
    if (iXPolType > 0) then
      tempStr = '       a(xx)       a(xy)       a(xz)       a(yy)       a(yz)       a(zz)'
    else
      tempStr = ' '
    end if
    write(u6,*)
    write(u6,*) ' External field specification in au'
    write(u6,*) ' =================================='
    if (nOrd_XF == 0) then
      write(u6,*) '     x           y           z           Z'//tempStr
    else if (nOrd_XF == 1) then
      write(u6,*) '     x           y           z           Z         my(x)       my(y)       my(z)'//tempStr
    else if (nOrd_XF == 2) then
      write(u6,*) '     x           y           z           Z         my(x)       my(y)       my(z)'// &
                  '       Q(xx)       Q(xy)       Q(xz)       Q(yy)       Q(yz)       Q(zz)'//tempStr
    else if (nOrd_XF == -1) then
      write(u6,*) '     x           y           z '//tempstr
    else
      call WarningMessage(2,'Option not implemented yet!')
      call Abend()
    end if
  end if

  write(Format_XF,'(A,I2.2,A)') '(',nData_XF,'(F10.6,2x))'
  XnetCharg = Zero
  do iXF=1,nXF
    A(1:3) = XF(1:3,iXF)
    Charge_iXF = XF(4,iXF)
    iChxyz = iChAtm(A)
    iDum = 0
    call Stblz(iChxyz,nStab_iXF,iStb,iDum,jCoSet)
    if (nPrint(2) >= 6) write(u6,Format_XF) (XF(i,iXF),i=1,nData_XF)
    XnetCharg = XnetCharg+real(nIrrep/nStab_iXF,kind=wp)*Charge_iXF
  end do
  write(u6,*)
  write(u6,*) ' Net charge from external field: ',XnetCharg
end if
!                                                                      *
!***********************************************************************
!                                                                      *
if (RMat_On) then
  write(u6,*) ' Parameters for radial integration (R-matrix option)'
  write(u6,*) ' ==================================================='
  write(u6,'(A,G12.5)') '   rmatr     :',RmatR
  write(u6,'(A,G12.5)') '   epsabs    :',Epsabs
  write(u6,'(A,G12.5)') '   epsrel    :',Epsrel
  write(u6,'(A,G12.5)') '   qcoul     :',qCoul
  write(u6,'(A,G12.5)') '   dipol(1)  :',dipol(1)
  write(u6,'(A,G12.5)') '   dipol(2)  :',dipol(2)
  write(u6,'(A,G12.5)') '   dipol(3)  :',dipol(3)
  write(u6,'(A,G12.5)') '   epsq      :',epsq
  write(u6,'(A,G12.5)') '   bparm     :',bParm
end if
!                                                                      *
!***********************************************************************
!                                                                      *
#ifdef _EFP_
if (nEFP_fragments /= 0) then
  call EFP_PRINT_BANNER()
  write(u6,*)
  write(u6,*) ' Specification of Effective Fragment Potentials'
  write(u6,*)
  if (Coor_Type == XYZABC_type) then
    write(u6,*) 'In XYZABC format'
  else if (Coor_Type == POINTS_type) then
    write(u6,*) 'In Points format'
  else if (Coor_Type == ROTMAT_type) then
    write(u6,*) 'In RotMat format'
  else
    write(u6,*) 'Illegal Coor_type:',Coor_Type
    call Abend()
  end if
  do i=1,nEFP_Fragments
    write(u6,*)
    write(u6,*) 'Fragment:',FRAG_TYPE(i)
    if (Coor_Type == XYZABC_type) then
    else if (Coor_Type == POINTS_type) then
      do j=1,3
        write(u6,'(A10,3F20.10)') ABC(j,i)(1:10),(EFP_Coors((j-1)*3+k,i),k=1,3)
      end do
    else if (Coor_Type == ROTMAT_type) then
    end if
  end do
end if
#endif
!                                                                      *
!***********************************************************************
!                                                                      *
if (PrintOperators) then
  call CollapseOutput(0,'   Operator info:')
  write(u6,*)
end if

return

end subroutine Print_OpInfo
