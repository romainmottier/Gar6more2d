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

FUNCTION source(t)
  USE m_source
  USE m_const
  IMPLICIT NONE
  REAL*8,INTENT(in) :: t
  REAL*8 :: alpha,t1,source
  INTEGER :: sign

  alpha=-pi**2*f0**2
  sign=(-1)**(deriv+1)

  IF ((t.GT.0).AND.(t.LT.2*tdelay)) THEN
    t1=t-tdelay
    ! Rickert:
    SELECT CASE (deriv)
      CASE(0)
        source=sign*dexp(alpha*t1**2)
      CASE(1)
        source=sign*(2*alpha*t1)*dexp(alpha*t1**2)
      CASE(2)
        source=sign*(2*alpha+4*alpha**2*t1**2)*dexp(alpha*t1**2)
      CASE(3)
        source=sign*(12*alpha**2*t1+8*alpha**3*t1**3)*dexp(alpha*t1**2)
      CASE(4)
        source=sign*(12*alpha**2+48*alpha**3*t1**2+16*alpha**4*t1**4)*dexp(alpha*t1**2)
      CASE DEFAULT
        source=0.D0
    END SELECT
  ELSE
     source=0.D0
  END IF
END FUNCTION source
