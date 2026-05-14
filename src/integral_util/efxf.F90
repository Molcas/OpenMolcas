!***********************************************************************
! This file is part of OpenMolcas.                                     *
!                                                                      *
! OpenMolcas is free software; you can redistribute it and/or modify   *
! it under the terms of the GNU Lesser General Public License, v. 2.1. *
! OpenMolcas is distributed in the hope that it will be useful, but it *
! is provided "as is" and without any express or implied warranties.   *
! For more details see the full text of the license in the file        *
! LICENSE or in <http://www.gnu.org/licenses/>.                        *
!                                                                      *
! Copyright (C) 2004, Par Soderhjelm                                   *
!***********************************************************************

subroutine EFXF(coord,XF,nXF,nOrd_XF,iXPolType,dEF,XMolnr,nXMolnr,iGrid,scal14)
!***********************************************************************
!                                                                      *
!     Object:  Calculate electric field in one point                   *
!              from XFIELD multipoles                                  *
!              Note: Ignores symmetry totally!                         *
!                                                                      *
!     Authors: P. Soderhjelm                                           *
!              Dept. of Theor. Chem., Univ. of Lund, Sweden.           *
!                                                                      *
!              November 2004                                           *
!***********************************************************************

use Index_Functions, only: nTri_Elem1
use Constants, only: Zero, One, Two, Three, Half
use Definitions, only: wp, iwp

implicit none
real(kind=wp), intent(in) :: coord(3), XF(*), Scal14
integer(kind=iwp), intent(in) :: nXF, nOrd_XF, iXPolType, nXMolnr, XMolnr(nXMolnr,nXF), iGrid
real(kind=wp), intent(inout) :: dEF(3)
integer(kind=iwp) :: i, iFD, Inc, iOrdOp
real(kind=wp) :: Dax, Day, Daz, QaSum, Qaxx, Qaxy, Qaxz, Qayy, Qayz, Qazz, R12, Scal, x, y, z, ZA
logical(kind=iwp) :: LExcl

if (nOrd_XF < 0) return
! Calculate number of entries per XFIELD point
Inc = 3
do iOrdOp=0,nOrd_XF
  Inc = Inc+nTri_Elem1(iOrdOp)
end do
if (iXPolType > 0) Inc = Inc+6

! Loop over XF points
do iFd=1,nXF
  scal = One
  if ((iXPolType > 0) .and. (iGrid <= nXF)) then
    LExcl = .false.
    if (iFd == iGrid) LExcl = .true.
    do i=1,nXMolnr
      if (XMolnr(1,iGrid) == XMolnr(i,iFd)) LExcl = .true.
      if (XMolnr(1,iGrid) == -XMolnr(i,iFd)) scal = scal14
    end do
    if (LExcl) cycle
    !if (LExcl) then
    !  write(u6,*) 'EXCLUDE ',iFd,' from field at ',iGrid
    !  cycle
    !else if (scal < One) then
    !  write(u6,*) 'SCALE ',iFd,' from field at ',iGrid,' with',scal
    !end if
  end if
  ZA = Zero
  DAx = Zero
  DAy = Zero
  DAz = Zero
  QAxx = Zero
  QAxy = Zero
  QAxz = Zero
  QAyy = Zero
  QAyz = Zero
  QAzz = Zero
  if (nOrd_XF == 0) then
    ZA = XF((iFd-1)*Inc+4)*scal
  else if (nOrd_XF == 1) then
    ZA = XF((iFd-1)*Inc+4)*scal
    DAx = XF((iFd-1)*Inc+5)*scal
    DAy = XF((iFd-1)*Inc+6)*scal
    DAz = XF((iFd-1)*Inc+7)*scal
  else if (nOrd_XF == 2) then
    ZA = XF((iFd-1)*Inc+4)*scal
    DAx = XF((iFd-1)*Inc+5)*scal
    DAy = XF((iFd-1)*Inc+6)*scal
    DAz = XF((iFd-1)*Inc+7)*scal
    QAxx = XF((iFd-1)*Inc+8)*scal
    QAxy = XF((iFd-1)*Inc+9)*scal
    QAxz = XF((iFd-1)*Inc+10)*scal
    QAyy = XF((iFd-1)*Inc+11)*scal
    QAyz = XF((iFd-1)*Inc+12)*scal
    QAzz = XF((iFd-1)*Inc+13)*scal
  else
    call WarningMessage(2,'Efxf: Option not implemented yet!')
    call Abend()
  end if
  x = XF((iFd-1)*Inc+1)-coord(1)
  y = XF((iFd-1)*Inc+2)-coord(2)
  z = XF((iFd-1)*Inc+3)-coord(3)
  r12 = sqrt(x**2+y**2+z**2)

  ! Z field
  dEF(1) = dEF(1)-ZA*x/r12**3
  dEF(2) = dEF(2)-ZA*y/r12**3
  dEF(3) = dEF(3)-ZA*z/r12**3

  if (nOrd_XF < 1) cycle

  ! D field
  dEF(1) = dEF(1)+Three*(DAx*x+DAy*y+DAz*z)*x/r12**5-DAx/r12**3
  dEF(2) = dEF(2)+Three*(DAx*x+DAy*y+DAz*z)*y/r12**5-DAy/r12**3
  dEF(3) = dEF(3)+Three*(DAx*x+DAy*y+DAz*z)*z/r12**5-DAz/r12**3

  if (nOrd_XF < 2) cycle

  ! Q field
  QAsum = (QAxx*x*x+QAyy*y*y+QAzz*z*z+Two*(QAxy*x*y+QAxz*x*z+QAyz*y*z))

  dEF(1) = dEF(1)+Half*(-15.0_wp/r12**7*x*QAsum+Three/r12**5*(Three*QAxx*x+Two*QAxy*y+Two*QAxz*z+QAyy*x+QAzz*x))

  dEF(2) = dEF(2)+Half*(-15.0_wp/r12**7*y*QAsum+Three/r12**5*(QAxx*y+Two*QAxy*x+Three*QAyy*y+Two*QAyz*z+QAzz*y))

  dEF(3) = dEF(3)+Half*(-15.0_wp/r12**7*z*QAsum+Three/r12**5*(QAxx*z+Two*QAxz*x+QAyy*z+Two*QAyz*y+Three*QAzz*z))

  ! These formulas gives the corresponding energy terms in drvn0
  !eZD = ZA*(DRBx*x+DRBy*y+DRBz*z)/r12**3
  !eDD = (DAx*DRBx+DAy*DRBy+DAz*DRBz)/r12**3-Three*(DAx* x+DAy* y+DAz *z)*(DRBx*x+DRBy*y+DRBz*z)/r12**5c
  !eQD = -Half*(-15.0_wp/r12**7*(DRBx*x+DRBy*y+DRBz*z)*QAsum+Three/r12**5* &
  !             (Three*DRBx*QAxx*x+DRBy*QAxx*y+DRBz*QAxx*z+Two*DRBx*QAxy*y+Two*DRBy*QAxy*x+Two*DRBx*QAxz*z+ &
  !              Two*DRBz*QAxz*x+DRBx*QAyy*x+Three*DRBy*QAyy*y+DRBz*QAyy*z+Two*DRBy*QAyz*z+Two*DRBz*QAyz*y+DRBx*QAzz*x+ &
  !              DRBy*QAzz*y+Three*DRBz*QAzz*z))

end do   !iFd

return
end subroutine EFXF
