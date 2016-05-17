!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TURBTWO_CONVECT_Upwind
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1

    INTEGER:: Lxyz
    INTEGER:: NII,NJJ,NKK
    INTEGER:: NdI,NdJ,NdK
    INTEGER:: In2,Jn2,Kn2
    INTEGER:: In1,Jn1,Kn1
    INTEGER:: In0,Jn0,Kn0
    INTEGER:: Ip1,Jp1,Kp1
    INTEGER:: Ip2,Jp2,Kp2

    REAL:: UL,UR,VLl,VRr,WL,WR,PL,PR,TL,TR
    REAL:: rcpcv,DENL,DENR
    REAL:: sav1,sav2,sav3,sav,sav1n,sav2n,sav3n
    REAL:: ulnormal,urnormal,rulnormal,rurnormal
    REAL:: small
    REAL:: eigTKO,TKL,TKR,dTK,TOL,TOR,dTO


    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1



    DO Lxyz=1,3
        NII=NI1;  NJJ=NJ1;  NKK=NK1
        NdI=0;    NdJ=0;    NdK=0
        IF (Lxyz==1) THEN
             NII=NI;  NdI=1
        ELSEIF (Lxyz==2) THEN
             NJJ=NJ;  NdJ=1
        ELSE
             NKK=NK;  NdK=1
        ENDIF
        
        DO K=1,NKK
        DO J=1,NJJ
        DO I=1,NII
	     !In0=I  ;		Jn0=J;		    Kn0=K   
	     In1=I-NdI;		Jn1=J-NdJ;      Kn1=K-NdK
	     !In2=In1-NdI;	Jn2=Jn1-NdJ;    Kn2=Kn1-NdK
	     !Ip1=I+NdI;		Jp1=J+NdJ;      Kp1=K+NdK
	     !Ip2=Ip1+NdI;	Jp2=Jp1+NdJ;    Kp2=Kp1+NdK



       UL=VL(Lxyz,1,i,j,k)
	   UR=VR(Lxyz,1,i,j,k)
         !!!!!!!!!!!!!!!!!
         VLl=VL(Lxyz,2,i,j,k)
	   VRr=VR(Lxyz,2,i,j,k)
         !!!!!!!!!!!!!!!!!
         WL=VL(Lxyz,3,i,j,k)
	   WR=VR(Lxyz,3,i,j,k)

         PL=VL(Lxyz,4,i,j,k)
	   PR=VR(Lxyz,4,i,j,k)

	   TL=VL(Lxyz,5,i,j,k)
	   TR=VR(Lxyz,5,i,j,k)
         !!!!!!!!!!!!!!!!!

	   rcpcv=ppf !??
         DENL=PL/(rcpcv*TL)
         DENR=PR/(rcpcv*TR)

!         surface area vectors
!
          sav1=SD(Lxyz,1,I,J,K)  !aix(i,j,k)
          sav2=SD(Lxyz,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(Lxyz,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(Lxyz,I,J,K)   !aim(i,j,k)+small
          sav1n=sav1/sav
          sav2n=sav2/sav
          sav3n=sav3/sav


!          inviscid flux
   
           ulnormal=(sav1*UL+sav2*VLl+sav3*WL)                      ! ulnormal=(sav1*ul+sav2*vl+sav3*wl+vib(i,j,k))
           urnormal=(sav1*UR+sav2*VRr+sav3*WR)                      ! urnormal=(sav1*ur+sav2*vr+sav3*wr+vib(i,j,k))
           rulnormal=DENL*ulnormal                                   ! rulnormal=rl*ulnormal
           rurnormal=DENR*urnormal                                   ! rurnormal=rr*urnormal
 
		 small = 1.e-10

           eigTKO=0.5*max(abs(ulnormal),abs(urnormal))  !max(eig12,eig13)+max(eig22,eig23)
           TKL=max(VL(Lxyz,6,i,j,k),small) !tel=max(teil(i,j,k),small)
           TKR=max(VR(Lxyz,6,i,j,k),small) !ter=max(teir(i,j,k),small)
           dTK=-0.5*(rulnormal*TKL+rurnormal*TKR &
     &                 -eigtko*(DENR*TKR-DENL*TKL))   !drtei= -0.5*(rulnormal*tel+rurnormal*ter-(eigto)*(rr*ter-rl*tel))
	
           TOL=max(VL(Lxyz,7,i,j,k),small)  !oql=max(oqil(i,j,k),small)
           TOR=max(VR(Lxyz,7,i,j,k),small)  !oqr=max(oqir(i,j,k),small)
           dTO=-0.5*(rulnormal*TOL+rurnormal*TOR    &
     &	 	        -eigtko*(DENR*TOR-DENL*TOL))    !droqi= -0.5*(rulnormal*oql+rurnormal*oqr-(eigto)*(rr*oqr-rl*oql))

  
	    D(6,I,J,K)=D(6,I,J,K)-dTK
          D(6,In1,Jn1,Kn1)=D(6,In1,Jn1,Kn1)+dTK

	    D(7,I,J,K)=D(7,I,J,K)-dTO
          D(7,In1,Jn1,Kn1)=D(7,In1,Jn1,Kn1)+dTO


    	  ENDDO
    	  ENDDO
	      ENDDO !kji

    ENDDO !Lxyz


    RETURN
END SUBROUTINE TURBTWO_CONVECT_Upwind