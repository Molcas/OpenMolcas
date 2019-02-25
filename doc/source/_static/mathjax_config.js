/***********************************************************************
* This file is part of OpenMolcas.                                     *
*                                                                      *
* OpenMolcas is free software; you can redistribute it and/or modify   *
* it under the terms of the GNU Lesser General Public License, v. 2.1. *
* OpenMolcas is distributed in the hope that it will be useful, but it *
* is provided "as is" and without any express or implied warranties.   *
* For more details see the full text of the license in the file        *
* LICENSE or in <http://www.gnu.org/licenses/>.                        *
***********************************************************************/

//MathJax.Hub.Config({
window.MathJax = {
  "HTML-CSS": {
    scale: 90,
    preferredFont: "STIX",
    webFont: "STIX-Web",
  },
  TeX: {
    extensions: ["AMSmath.js","AMSsymbols.js","mhchem.js"],
    Macros: {
      mat: ["\\boldsymbol{\#1}", 1],
      sign: ["\\operatorname{sign}", 0],
      Tr: ["\\operatorname{Tr}", 0],
      abs: ["\\operatorname{abs}", 0],
      bra: ["\\left<\#1\\right|", 1],
      ket: ["\\left|\#1\\right>", 1],
      braket: ["\\left<\#1\\middle|\#2\\right>", 2],
      braopket: ["\\left<\#1\\middle|\#2\\middle|\#3\\right>", 3],
    },
  },
  SVG: {
    scale: 90,
    font: "STIX-Web",
  },
//});
};
