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
      Subroutine Rd1Int_MCLR
************************************************************************
*                                                                      *
*     Read header and matrices from the one-electron integral file     *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*                                                                      *
*----------------------------------------------------------------------*
*                                                                      *
*     history: none                                                    *
*                                                                      *
************************************************************************
      Implicit Real*8 (A-H,O-Z)
#include "Input.fh"
#include "Pointers.fh"
#include "WrkSpc.fh"
      Character*8 Method
*
*----------------------------------------------------------------------*
*     Start                                                            *
*----------------------------------------------------------------------*
*     debug=.true.
*
      Call Get_cArray('Relax Method',Method,8)
      If ( Method.eq.'RHF-SCF ' ) then
       iMethod=1
      Else if ( Method.eq.'RASSCF  ' .or.
     &          Method.eq.'CASSCF  ') then
       iMethod=2
      Else if ( Method.eq.'RASSCFSA' .or.
     &          Method.eq.'CASSCFSA' ) then
       iMethod=2
      Else if ( Method.eq.'CASPT2  ' ) then
       iMethod=3
      Else If ( Method.eq.'MBPT2   ' ) then
       iMethod=4
      Else If ( Method.eq.'MCPDFT  ' ) then
       iMethod=2
      End If
*---  read file header  -----------------------------------------------*
      Call Get_cArray('Seward Title',Header1I,144)
*---  read number of symm. species ------------------------------------*
      Call get_iScalar('nSym',nSym)
*---  read number of basis functions per symmetry species -------------*
      Call Get_iArray('nBas',nBas,nSym)
*---  read nuclear potential ------------------------------------------*
*     Call Get_PotNuc(PotNuc)
      Call Get_dScalar('PotNuc',Potnuc)
*---  read number of atoms --------------------------------------------*
      Call Get_iScalar('Unique atoms',nAtoms)
*---  read atom labels ------------------------------------------------*
      Call Get_cArray('Unique Atom Names',AtLbl,LENIN*nAtoms)
*---  read atom coordinates -------------------------------------------*
      Call Get_dArray('Unique Coordinates',Coor,3*nAtoms)
*---  read labels of the irreps ---------------------------------------*
      Call Get_cArray('Irreps',ChIrr,24)
*----------------------------------------------------------------------*
*     Precompute the total sum of variables and size of matrices       *
*----------------------------------------------------------------------*
      ntBas=0
      ntBtri=0
      ntBsqr=0
      Do iSym=1,nSym
         ntBas=ntBas+nBas(iSym)
         ntBtri=ntBtri+nBas(iSym)*(nBas(iSym)+1)/2
         ntBsqr=ntBsqr+nBas(iSym)*nBas(iSym)
      End Do
*----------------------------------------------------------------------*
*     Exit                                                             *
*----------------------------------------------------------------------*
      Return
      End
