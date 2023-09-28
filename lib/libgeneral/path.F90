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

FUNCTION path(x,y,h,V1,V2,t)
  USE m_const
  IMPLICIT NONE
#include "m_inter.h90"
  REAL*8,INTENT(in) ::h,x,y,V1,V2,t
  REAL*8 :: delta,test(4),scale
  REAL*8:: VS(4,4),Rwork(8),WR(4),WI(4),fv1(4),fv2(4),fv3(4)
  REAL*16:: alpha,coe,dh,dx,dy,dV1,dV2,dt
  COMPLEX*16 ::Eig(4),path,WC(4),Mat(4,4),Work(200),VL(4,4),VR(4,4)
  INTEGER ::Info,ii,Sdim,ILO,IHI,I
  LOGICAL Bwork(4)
  Mat=0.D0
  dh=h
  dx=x
  dy=y
  dv1=V1
  dv2=V2
  dt=t
  IF (y.NE.0) THEN
     IF (x.NE.0d0) THEN
!!$      alpha=(x**2+y**2+h**2)**2-4*h**2*y**2;
!!$      Mat(1,1)=-4*t*x*(x**2+y**2+h**2)/alpha
!!$      Mat(1,2)=2*(y**2/V2**2*(x**2+y**2-h**2)+h**2/V1**2*(x**2+h**2-y**2)-t**2*(3*x**2+y**2+h**2))/alpha
!!$      Mat(1,3)=4*t*x*(y**2/V2**2+h**2/V1**2-t**2)/alpha
!!$      Mat(1,4)=-((y**2/V2**2-h**2/V1**2)**2+t**2*(t**2-2*y**2/V2**2-2*h**2/V1**2))/alpha
!!$      Mat(2,1)=1.D0
!!$      Mat(3,2)=1.D0
!!$      Mat(4,3)=1.D0
!!$write(6,*) Mat
!!$sdim=400
!!$WR=0.D0
!!$WI=0.D0
!!$info=0
!!$Bwork=1
!!$VS=0.D0
!!$  CALL DGEES( 'N', 'N', 0, 4, Mat, 4, SDIM, WR, WI,&
!!$             &VS, 4, WORK, 136, BWORK, INFO )
!!$  eig=wr+ima*wi
!!$write(6,*) eig,sdim,info,bwork,VS
        
        alpha=(dx**2+dy**2+dh**2)**2-4.Q0*dh**2*dy**2;
        coe=4.Q0*dt*dx*(dx**2+dy**2+dh**2)/(2.Q0*(dy**2/V2**2*(dx**2+dy**2-dh**2)+&
             &dh**2/dV1**2*(dx**2+dh**2-dy**2)-dt**2*(3.Q0*dx**2+dy**2+dh**2))) 
        Mat=0.D0
        Mat(1,1)=-coe*4.Q0*dt*dx*(dx**2+dy**2+dh**2)/alpha
        Mat(1,2)=2.Q0*coe**2*(dy**2/dV2**2*(dx**2+dy**2-dh**2)+dh**2/dV1**2*(dx&
             &**2+dh**2-dy**2)-dt**2*(3.Q0*dx**2+dy**2+dh**2))/alpha
        Mat(1,3)=4.Q0*coe**3*dt*dx*(dy**2/dV2**2+dh**2/dV1**2-dt**2)/alpha
        Mat(1,4)=-coe**4*((dy**2/dV2**2-dh**2/dV1**2)**2+dt**2*(dt**2-2.Q0*dy**2/dV2**2-2.Q0*dh**2/dV1**2))/alpha 
        !      Mat(1,1)=-1e3*4*t*x*(x**2+y**2+h**2)/alpha
        !      Mat(1,2)=1.e6*2*(y**2/V2**2*(x**2+y**2-h**2)+&
        !           &h**2/V1**2*(x**2+h**2-y**2)-t**2*(3*x**2+y**2+h**2))/alpha
        !      Mat(1,3)=4*1.e9*t*x*(y**2/V2**2+h**2/V1**2-t**2)/alpha
        !      Mat(1,4)=-1.e12*((y**2/V2**2-h**2/V1**2)**2+t**2*(t**2-2*y**2/V2**2-2*h**2/V1**2))/alpha
        Mat(2,1)=1.D0
        Mat(3,2)=1.D0
        Mat(4,3)=1.D0
        !write(6,*) Mat
        sdim=1200
        WR=0.D0
        WI=0.D0
        info=0
        Bwork=.TRUE. 
        VS=0.D0
!        WRITE(6,*) MAT
!       CALL  ZGEBAL( 'B', 4, Mat, 4, ILO, IHI, SCALE, INFO )
        CALL ZGEEV( 'V', 'V',  4, Mat, 4,  WC, &
             &VR,4,VL,4, WORK, 200,RWORK,  INFO )
!      call  cg(4,4,real(mat),aimag(mat),wr,wi,mat,vr,vl,fv1,fv2,fv3,info)
!        DO I=1,4
!           WRITE(6,*) SUM(Mat(I,:)*VR(:,2))/VR(I,2)
!           END DO
!        WRITE(6,*) WC
!        WRITE(6,*) work(1)
        !  eig=wr+ima*wi
        !write(6,*) eig,sdim,info,bwork,VS
        !
        !      write(6,*) 'toto',eig(1)**4,Mat(1,1)*eig(1)**3,Mat(1,2)*eig(1)**2,Mat(1,3)*eig(1),Mat(1,4)
        EIG=-AIMAG(WC)+ima*REAL(WC)!-WI+ima*WR
        !      test=eig**4+Mat(1,1)*eig**3+Mat(1,2)*eig**2+Mat(1,3)*eig+Mat(1,4)
        eig=eig/coe
        test=ABS(-y*SQRT(1/V2**2+eig**2)+h*SQRT(1/V1**2+eig**2)+ima*x*eig-t)
        !      write(6,*) 'toto',-y*sqrt(1/V2**2+eig(4)**2),h*sqrt(1/V1**2+eig(4)**2),ima*x*eig(4)-t,test
        ii=MINLOC(test,dim=1)
        path=eig(ii)
!WRITE(6,*) path
        !write(6,*) eig
        !write(6,*) aimag(derivee(-abs(x),y,h,V1,V2,eig(1)))
        !write(6,*) aimag(derivee(-abs(x),y,h,V1,V2,eig(2)))
        !write(6,*) aimag(derivee(-abs(x),y,h,V1,V2,eig(3)))
        !write(6,*) aimag(derivee(-abs(x),y,h,V1,V2,eig(4)))
        IF(ABS(REAL(path)/AIMAG(path)).LT.1e-8) THEN
           IF (AIMAG(derivee(-ABS(x),y,h,V1,V2,path)).LT.0.D0) THEN
              test(ii)=100
              ii=MINLOC(test,dim=1)
              path=eig(ii)
           END IF
        END IF
        IF (REAL(path).LT.0.D0) THEN
           path=-REAL(path)+ima*AIMAG(path)
        END IF
     ELSE IF (ABS(y).NE.h) THEN
        alpha=(x**2+y**2+h**2)**2-4*h**2*y**2; 
        Mat(1,1)=2*(y**2/V2**2*(x**2+y**2-h**2)+h**2/V1**2*(x**2+h**2-y**2)-t**2*\
        (3*x**2+y**2+h**2))/alpha 
        Mat(1,2)=((y**2/V2**2-h**2/V1**2)**2+t**2*(t**2-2*y**2/V2**2-2*h**2/V1**\
2       ))/alpha
        delta=Mat(1,1)**2-4*Mat(1,2)
        IF (delta<0) THEN
           Eig(1)=SQRT((-Mat(1,1)+ima*SQRT(-delta))/2.)
           Eig(2)=-SQRT((-Mat(1,1)+ima*SQRT(-delta))/2.)
           Eig(3)=SQRT((-Mat(1,1)-ima*SQRT(-delta))/2.)
           Eig(4)=-SQRT((-Mat(1,1)-ima*SQRT(-delta))/2.)
        ELSE
           Eig(1)=SQRT((-Mat(1,1)+SQRT(delta))/2.+0*ima)
           Eig(2)=-SQRT((-Mat(1,1)+SQRT(delta))/2.+0*ima)
           Eig(3)=SQRT((-Mat(1,1)-SQRT(delta))/2.+0*ima)
           Eig(4)=-SQRT((-Mat(1,1)-SQRT(delta))/2.+0*ima)
        END IF
        test=ABS(-y*SQRT(1/V2**2+eig**2)+h*SQRT(1/V1**2+eig**2)+ima*x*eig-t) 
        ii=MINLOC(test,dim=1) 
        path=eig(ii) 
        IF(ABS(REAL(path)/AIMAG(path)).LT.1e-8) THEN 
           IF (AIMAG(derivee(-ABS(x),y,h,V1,V2,path)).LT.0.D0) THEN 
              test(ii)=100 
              ii=MINLOC(test,dim=1) 
              path=eig(ii) 
           END IF
        END IF
        IF (REAL(path).LT.0.D0) THEN 
           path=-REAL(path)+ima*AIMAG(path) 
        END IF

     ELSE
        path=SQRT(t**2*(t**2-2*h**2/V2**2-2*h**2/V1**2)/(4*t**2*h**2)+0*ima)
     END IF
  ELSE
     delta=h**2*(t**2-(h**2+x**2)/V1**2)
     IF (delta<0) THEN
        IF(x>0) THEN
           path=-ima*x*t/(h**2+x**2)+ima*SQRT(-delta)/(h**2+x**2);
        ELSE
           path=-ima*x*t/(h**2+x**2)-ima*SQRT(-delta)/(h**2+x**2);
        END IF
     ELSE
        path=-ima*x*t/(h**2+x**2)+SQRT(delta)/(h**2+x**2);
     END IF
  END IF
END FUNCTION path

