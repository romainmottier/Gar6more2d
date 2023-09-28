!! Copyright INRIA. Contributors : Julien DIAZ and Abdelaaziz EZZIANI
!! 
!! Julien.Diaz@inria.fr and Abdelaaziz.Ezziani@univ-pau.fr
!! 
!! This software is a computer program whose purpose is to
!! compute the analytical solution of problems of waves propagation in two 
!! layered media such as
!! - acoustic/acoustic
!! - acoustic/elastodynamic
!! - acoustic/porous
!! - porous/porous,
!! based on the Cagniard-de Hoop method.
!! 
!! This software is governed by the CeCILL license under French law and
!! abiding by the rules of distribution of free software.  You can  use, 
!! modify and/ or redistribute the software under the terms of the CeCILL
!! license as circulated by CEA, CNRS and INRIA at the following URL
!! "http://www.cecill.info". 
!! 
!! As a counterpart to the access to the source code and  rights to copy,
!! modify and redistribute granted by the license, users are provided only
!! with a limited warranty  and the software's author,  the holder of the
!! economic rights,  and the successive licensors  have only  limited
!! liability. 
!! 
!! In this respect, the user's attention is drawn to the risks associated
!! with loading,  using,  modifying and/or developing or reproducing the
!! software by the user in light of its specific status of free software,
!! that may mean  that it is complicated to manipulate,  and  that  also
!! therefore means  that it is reserved for developers  and  experienced
!! professionals having in-depth computer knowledge. Users are therefore
!! encouraged to load and test the software's suitability as regards their
!! requirements in conditions enabling the security of their systems and/or 
!! data to be ensured and,  more generally, to use and operate it in the 
!! same conditions as regards security. 
!! 
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL license and that you accept its terms.
!! ========================================================================

FUNCTION derivee(x,y,h,V1,V2,pp)
  USE m_const
  IMPLICIT NONE
  REAL*8,INTENT(in) ::h,x,y,V1,V2
  COMPLEX*16,INTENT(in)  :: pp
  COMPLEX*16 :: derivee,z1,z2
  IF (ABS(AIMAG(pp)).EQ.1.D0/V1**2) THEN
     derivee=0;
  ELSE
     IF ((ABS(imag(pp))>1/V1).AND.(REAL(pp).EQ.0.d0)) THEN
        z1=ima*SQRT(-1.D0/V1**2-pp**2)
     ELSE
        z1=zsqrt(1.D0/V1**2+pp**2)
     ENDIF
     IF ((ABS(imag(pp))>1/V2).AND.(REAL(pp).EQ.0.d0)) THEN
        z2=ima*SQRT(-1.D0/V2**2-pp**2)
     ELSE
        z2=zsqrt(1.D0/V2**2+pp**2)
     ENDIF
     IF ((ABS(imag(pp))>1.D0/V1).AND.(REAL(pp).EQ.0.d0)) THEN
        derivee=-y*pp/z2+h*pp/z1+ima*x;
     ELSE
        derivee=-y*pp/z2+h*pp/z1+ima*x;
     END IF
     IF (derivee.NE.0.d0) THEN
        derivee=1.D0/derivee;
     ELSE
        derivee=0;
     END IF
  END IF
END FUNCTION derivee
