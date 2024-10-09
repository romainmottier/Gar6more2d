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

FUNCTION calccoeffP_ee(pp)
  USE m_const
  USE m_phys
  IMPLICIT NONE
  COMPLEX*16,INTENT(in) :: pp
  INTEGER ::Info,IPIV(4)
  COMPLEX*16 :: kp1,ks1,ks2,kp2
  COMPLEX*16 :: Mat(4,4),coeff(4),calccoeffP_ee(4),X(2)
  IF(REAL(pp).NE.0.d0) THEN
     kp1=SQRT(pp**2+1/Vp1**2)
     ks1=SQRT(pp**2+1/Vs1**2)
     kp2=SQRT(pp**2+1/Vp2**2)
     ks2=SQRT(pp**2+1/Vs2**2)
  ELSE
     IF (AIMAG(pp)**2>1/Vp1**2) THEN
        kp1=ima*SQRT(-pp**2-1/Vp1**2)
     ELSE
        kp1=SQRT(pp**2+1/Vp1**2)
     END IF
     IF (AIMAG(pp)**2>1/Vs1**2) THEN
        ks1=ima*SQRT(-pp**2-1/Vs1**2)
     ELSE
        ks1=SQRT(pp**2+1/Vs1**2)
     END IF
     IF (AIMAG(pp)**2>1/Vp2**2) THEN
        kp2=ima*SQRT(-pp**2-1/Vp2**2)
     ELSE
        kp2=SQRT(pp**2+1/Vp2**2)
     END IF
     IF (AIMAG(pp)**2>1/Vs2**2) THEN
        ks2=ima*SQRT(-pp**2-1/Vs2**2)
     ELSE
        ks2=SQRT(pp**2+1/Vs2**2)
     END IF
  END IF
  MAT(1,1)=-ima*pp
  MAT(2,1)=-kp1
  MAT(3,1)=2*mu1*kp1*ima*pp
  MAT(4,1)=lambda1/Vp1**2+2*mu1*kp1**2


  MAT(1,2)=-ks1
  MAT(2,2)=ima*pp
  MAT(3,2)=mu1*(ks1**2+pp**2)
  MAT(4,2)=-2*mu1*ks1*ima*pp

  MAT(1,3)=ima*pp
  MAT(2,3)=-kp2
  MAT(3,3)=2*mu2*kp2*ima*pp
  MAT(4,3)=-lambda2/Vp2**2-2*mu2*kp2**2

  MAT(1,4)=-ks2
  MAT(2,4)=-ima*pp
  MAT(3,4)=-mu2*(ks2**2+pp**2)
  MAT(4,4)=-2*mu2*ks2*ima*pp

  Coeff(1)=ima*pp/(2*kp1)
  Coeff(2)=-1/2D0
  Coeff(3)=mu1*ima*pp
  Coeff(4)=-(lambda1/Vp1**2+2*mu1*kp1**2)/(2*kp1)
  Coeff=Coeff/Vp1**2
  CALL ZGESV( 4, 1, MAT, 4, IPIV, Coeff, 4, INFO )
  calccoeffP_ee=Coeff
END FUNCTION calccoeffP_ee
