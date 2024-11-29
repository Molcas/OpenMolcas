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
Module MCLR_Data

! Stuff from Pointers.fh
Integer ipMat(8,8),ipMatLT(8,8),ipCM(8),ipMC(8,8),ipmatba(8,8),ipMO(8,8,8),iADMO(8,8,8)
Integer nDens,nDensLT,nCMO,nscrtch,ipCI,nDens2,ipf0,ndensc,na(8),nb(8)
Integer nmba,nna,n1dens,n2dens,nconf1,nacpar,nacpr2

! Stuff from spin_mclr.fh
real*8     rms,rbetaa,rbetas,ralpha

! Stuff from genop.fh
!. Type of operator in action
!        I12        :         2:Both on and two electron integrals
!        IST        :
!        Square     :         Integrals square/triangular
Logical square
Integer I12,IST

!Stuff from machine.fh
integer  nrec

!Stuff from csfsd.fh
Integer i1,iAnders,lconf,lldet,iAllo

!Stuff from incdia.fh
integer ipdia

! Print flags from MCLR_Data.fh
Integer IPRSTR,IPRCIX,IPRORB,IPRDIA,IPRXT,IPRANA

!Stuff from sa.fh
Logical Fancy_Preconditioner,SA,esterr,isNAC,override
Logical isMECIMSPD
Integer istate,irlxroot,NACstates(2),NSSA(2)

End Module MCLR_Data
