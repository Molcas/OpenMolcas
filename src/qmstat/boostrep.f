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
      Subroutine BoostRep(AddRep,SmatPure,iVecs,nSize,InCutOff)
      Implicit Real*8 (a-h,o-z)

#include "maxi.fh"
#include "numbers.fh"
#include "qminp.fh"
#include "WrkSpc.fh"

      Dimension SmatPure(*)
      Logical InCutOff

*
*-- Enter.
*
*
*-- Common section.
*
      Scalar=0
*
*-- Take different route for different QM-method.
*
      If(QmType(1:3).eq.'SCF') then
*
*-- Repulsion term added to the Energy
*-- Calculated as S_i*S_i*CPsi_m* CPis_n
*-- S_i is the overlap integral (with the solvent molecule)
*-- for the occupied orbitals of the quantum system
*-- CPsi_m and CPsi_n are the transformation coeficients
*-- obtained from the diagonalization procedure of the Fock matrix
*-- to go from the original wavefunction to the final wavefunction
*-- after the SCF procedure. These coeficientes run over all basis set
*
        Do 801, iO1=1,nSize
          Do 802, iO2=1,nSize
            Do 803, i=1,iOcc1
                kaunter=i*(i+1)/2
                ind1=nSize*(iO1-1)+i-1
                ind2=nSize*(iO2-1)+i-1
                Scalar=Scalar+(Work(iVecs+ind1)*Work(iVecs+ind2)
     &                       *SmatPure(kaunter))
803         Continue
802       Continue
801     Continue
        AddRep=exrep4*abs(Scalar)**2+exrep6*abs(Scalar)**3
     &                       +exrep10*abs(Scalar)**5
      Elseif(QmType(1:4).eq.'RASS') then
        Do 813, i=1,nSize
          Do 814, j=1,nSize
            If(i.ge.j) then
              kaunter=i*(i+1)/2-i+j
            Else
              kaunter=j*(j+1)/2-j+i
            Endif
            ind1=nSize*(nEqState-1)+i-1
            ind2=nSize*(nEqState-1)+j-1
            Scalar=Scalar+Work(iVecs+ind1)*Work(iVecs+ind2)
     &                   *SmatPure(kaunter)
814       Continue
813     Continue
        AddRep=exrep4*abs(Scalar)**2+exrep6*abs(Scalar)**3
     &                       +exrep10*abs(Scalar)**5
      Endif

*
*-- Crazy energy added if inner cut-off has been passed. Ensure reject.
*
      If(InCutOff) AddRep=1D+20


      Return
      End
