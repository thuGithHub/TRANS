!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TIME_LUSGS
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


!	INTEGER IL,JL,KL
	REAL::	cvs,  RDSt  ,stept

    REAL:: AinRo,AinRx,AinRy,AinRz,AinRe
    REAL:: pRo,roinv,pVx,pVy,pVz,pE,pPP,pTT,pGam
    REAL:: AAc,UVW2,ALL,ubetac2

    REAL:: DRoIm1,DRvxIm1,DRvyIm1,DRvzIm1,DReIm1
    REAL:: DRoJm1,DRvxJm1,DRvyJm1,DRvzJm1,DReJm1
    REAL:: DRoKm1,DRvxKm1,DRvyKm1,DRvzKm1,DReKm1
    REAL:: AinX,AinY,AinZ,AinP
    REAL:: aSDx,aSDy,aSDz
    REAL:: DifRo,DifRvx,DifRvy,DifRvz,DifRe
    real::vibc,vjbc,vkbc


    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1




!        if(ThisBlock%ID_Present_Blk==3)write(*,*)"LUSGS U00=",KI,D(1,1,1,1), F(1,1,1,1)

!prblems now:
!is dtm right in subiterations? (maybe different definitions in Ryx ver)
!smooth is necssary?
! residuals what to do?  ->  in smooth
! why when i=0 DRoIm1.... = 0?


! update Fluxes, as time series

! niterc should do?   iteration in LU

	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1
		!stept = Dtm(I,J,K) 
!	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
!          step(i,j,k)=2.*ste/(3.*ste/tstep+2.)


!           open(unit=120,file='rho record-0.txt')
 !            write(120,*) I, J, K,  D(1,I,J,K),D(5,I,J,K)
!if(ThisBlock%ID_Present_Blk==1) write(*,*)"Uopbef=",Ki,I,j,k,D(1:5,i,j,k)
  !           continue
         
		 do L=1,5
!		   Q(L,I,J,K) = -(Q(L,I,J,K) - D(L,I,J,K))   ! should be + or -?
	       D(L,I,J,K) = D(L,I,J,K) / Vol(I,J,K)         &
     &  -(3.0*V(L,I,J,K)-4.0*Wt0(L,I,J,K)+Wt1(L,I,J,K))*0.5/      &
     &	DT_LUSGS_Main   !stept
	!   should save for Wt0 and Wt1,  Vol, DTM right?  Dtm=step?
         enddo
         
        
!if(ThisBlock%ID_Present_Blk==1) write(*,*)"Uopaft=",Ki,I,j,k,D(1:5,i,j,k)
!if(ThisBlock%ID_Present_Blk==3.and.KI>11.and.k>130)write(*,*)"LUDIni=",KI,i,j,k,D(1:5,i,j,k)
!endif
                 
	enddo
	enddo
	enddo

!     implicit operators
!	temporarily use F as diag
	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1
		!stept = Dtm(I,J,K) 
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)

		F(1,I,J,K) = 1.0+(Rds(1,I,J,K)+Rds(2,I,J,K)+Rds(3,I,J,K))       &
     &               *(stept/Vol(I,J,K))
!		cvs = F(1,I,J,K)
!		do L=1,5      !should do this??
!			Q(L,I,J,K)=Q(L,I,J,K)*cvs
!		enddo


	enddo
	enddo
	enddo




	IF (IF_PRECONDITION) THEN
!	IF (0) THEN

!		write(*,*) 'Do ', (D(L,1,9,1),L=1,5)


	    do K=1,NK1
	    do J=1,NJ1
	    do I=1,NI1

			    AinRo=D(1,I,J,K)
			    AinRx=D(2,I,J,K)
			    AinRy=D(3,I,J,K)
			    AinRz=D(4,I,J,K)
			    AinRe=D(5,I,J,K)

			    pRo = V(1,I,J,K)
			    roinv = 1.0/pRo
			    pVx = V(2,I,J,K)*roinv
			    pVy = V(3,I,J,K)*roinv
			    pVz = V(4,I,J,K)*roinv
			    pE  = V(5,I,J,K)*roinv
			    pPP = PP(I,J,K)
			    pTT = T(I,J,K)

		    ! with Inverse of transform matrix
                !AAc =SQRT(gam(I,J,K)*pPP*roinv)
                AAc =SQRT(1.4*pPP*roinv)    !by ydd
			    UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

		     ALL = min(Alagm(1,I,J,K),Alagm(2,I,J,K),Alagm(3,I,J,K))
	         ubetac2 = 0.0
		     CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I,J,K), &
         &				ALL,ubetac2)

		     CALL PrecondDW (AinRo,AinRx,AinRy,AinRz,AinRe, &
         &			pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
    !     >			AinP,AinX,AinY,AinZ,AinT)
			    D(1,I,J,K)=AinRo
			    D(2,I,J,K)=AinRx
			    D(3,I,J,K)=AinRy
			    D(4,I,J,K)=AinRz
			    D(5,I,J,K)=AinRe

	    enddo
	    enddo
	    enddo
	ENDIF




!L OP



!	III=IF_PRECONDITION
!	IF_PRECONDITION=0

	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1

		!stept = Dtm(I,J,K) 
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		cvs=stept/Vol(I,J,K) !beta

!if(ThisBlock%ID_Present_Blk==4)write(*,*)"iniIDR",KI,i,j,k,D(1:5,i-1,j,k),D(1:5,i,j,k)
		!I dir
		IF (I==1) THEN !??? why
			DRoIm1 =0.0
			DRvxIm1=0.0
			DRvyIm1=0.0
			DRvzIm1=0.0
			DReIm1 =0.0
		ELSE
		 
			AinRo=D(1,I-1,J,K)
			AinRx=D(2,I-1,J,K)
			AinRy=D(3,I-1,J,K)
			AinRz=D(4,I-1,J,K)
			AinRe=D(5,I-1,J,K)
            
          
          

			pRo = V(1,I-1,J,K)
			roinv = 1.0/pRo
			pVx = V(2,I-1,J,K)*roinv
			pVy = V(3,I-1,J,K)*roinv
			pVz = V(4,I-1,J,K)*roinv
			pE  = V(5,I-1,J,K)*roinv

			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
            !Real Gas
            
            !AinP=(gam(I-1,J,K)-1.0)*(AinRe-0.5*AinRo*        &
			AinP=(1.4-1.0)*(AinRe-0.5*AinRo*        &       !by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))     


		    aSDx = 0.5*(SD(1,1,I,J,K)+SD(1,1,I-1,J,K))
		    aSDy = 0.5*(SD(1,2,I,J,K)+SD(1,2,I-1,J,K))
		    aSDz = 0.5*(SD(1,3,I,J,K)+SD(1,3,I-1,J,K))
!added by ydd
            vibc=0.5*(vibn(i,j,k)+vibn(i-1,j,k))

			pPP = PP(I-1,J,K)
			pTT = T(I-1,J,K)
            !pGam=gam(I-1,J,K)
            pGam=1.4    !by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
			!--CALL LUSGS_ADQ
            CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vibc)

		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
			AAc =SQRT(1.4*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I-1,J,K),    Alagm(2,I-1,J,K),Alagm(3,I-1,J,K))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I-1,J,K),ALL,ubetac2)

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)

		 ENDIF


			RDSt=Rds(1,I-1,J,K)
			DRoIm1 =0.5*(DifRo + RDSt*D(1,I-1,J,K))
			DRVxIm1=0.5*(DifRvx+ RDSt*D(2,I-1,J,K))
			DRVyIm1=0.5*(DifRvy+ RDSt*D(3,I-1,J,K))
			DRVzIm1=0.5*(DifRvz+ RDSt*D(4,I-1,J,K))
			DReIm1 =0.5*(DifRe + RDSt*D(5,I-1,J,K))
            
             
        
        ENDIF
        
!if(ThisBlock%ID_Present_Blk==4)write(*,*)"IDR",KI,i,j,k,DifRo,D(1:5,i-1,j,k),DRoIm1,DRVxIm1,DRVyIm1,DRVzIm1
		!J dir
		IF (J==1) THEN !??? why
			DRoJm1 =0.0
			DRvxJm1=0.0
			DRvyJm1=0.0
			DRvzJm1=0.0
			DReJm1 =0.0
		ELSE
			AinRo=D(1,I,J-1,K)
			AinRx=D(2,I,J-1,K)
			AinRy=D(3,I,J-1,K)
			AinRz=D(4,I,J-1,K)
			AinRe=D(5,I,J-1,K)

!        if(ThisBlock%ID_Present_Blk==2)write(*,*)"DJ-1=",KI,i,j,k,D(1:5,i,j-1,k)
			pRo = V(1,I,J-1,K)
			roinv = 1.0/pRo
			pVx = V(2,I,J-1,K)*roinv
			pVy = V(3,I,J-1,K)*roinv
			pVz = V(4,I,J-1,K)*roinv
			pE  = V(5,I,J-1,K)*roinv
			
			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
            !AinP=(gam(I,J-1,K)-1.)*(AinRe-0.5*AinRo*     &
			AinP=(1.4-1.)*(AinRe-0.5*AinRo*     &       !by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))

		    aSDx = 0.5*(SD(2,1,I,J,K)+SD(2,1,I,J-1,K))
		    aSDy = 0.5*(SD(2,2,I,J,K)+SD(2,2,I,J-1,K))
		    aSDz = 0.5*(SD(2,3,I,J,K)+SD(2,3,I,J-1,K))
!added by ydd
            vjbc=0.5*(vjbn(i,j,k)+vjbn(j-1,j,k))

			pPP = PP(I,J-1,K)
			pTT = T(I,J-1,K)
      !      pGam=gam(I,J-1,K)
            pGam=1.4    !by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
            CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vjbc)


		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
			AAc =SQRT(1.4*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I,J-1,K),     &
     &		 Alagm(2,I,J-1,K),Alagm(3,I,J-1,K))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I,J-1,K),		ALL,ubetac2)

!	IF(I==1.and.J==9.and.K==1)THEN
!			CALL PrecondDW2 (DifRo,DifRvx,DifRvy,DifRvz,DifRe,
!     >			pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)
!	ELSE

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,		pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)
!	ENDIF
		 ENDIF

			RDSt=RDS(2,I,J-1,K)
			DRoJm1 =0.5*(DifRo + RDSt*D(1,I,J-1,K))
			DRVxJm1=0.5*(DifRvx+ RDSt*D(2,I,J-1,K))
			DRVyJm1=0.5*(DifRvy+ RDSt*D(3,I,J-1,K))
			DRVzJm1=0.5*(DifRvz+ RDSt*D(4,I,J-1,K))
			DReJm1 =0.5*(DifRe + RDSt*D(5,I,J-1,K))
		ENDIF

!if(ThisBlock%ID_Present_Blk==3.and.KI==13)write(*,*)"JDR",KI,i,j,k,DifRo,D(1:5,i,j-1,k),DRoJm1,DRVxJm1,DRVyJm1,DRVzJm1


		!K dir
		IF (K==1) THEN !??? why
			DRoKm1 =0.0
			DRvxKm1=0.0
			DRvyKm1=0.0
			DRvzKm1=0.0
			DReKm1 =0.0
		ELSE
			AinRo=D(1,I,J,K-1)
			AinRx=D(2,I,J,K-1)
			AinRy=D(3,I,J,K-1)
			AinRz=D(4,I,J,K-1)
			AinRe=D(5,I,J,K-1)

			pRo = V(1,I,J,K-1)
			roinv = 1.0/pRo
			pVx = V(2,I,J,K-1)*roinv
			pVy = V(3,I,J,K-1)*roinv
			pVz = V(4,I,J,K-1)*roinv
			pE  = V(5,I,J,K-1)*roinv
			
			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
!           AinP=(gam(I,J,K-1)-1.0)*(AinRe-0.5*AinRo*    &
            AinP=(1.4-1.0)*(AinRe-0.5*AinRo*    &   ! by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))
		
		    aSDx = 0.5*(SD(3,1,I,J,K)+SD(3,1,I,J,K-1))
		    aSDy = 0.5*(SD(3,2,I,J,K)+SD(3,2,I,J,K-1))
		    aSDz = 0.5*(SD(3,3,I,J,K)+SD(3,3,I,J,K-1))
!added by ydd
            vkbc=0.5*(vkbn(i,j,k)+vkbn(i,j,k-1))

			pPP = PP(I,J,K-1)
			pTT = T(I,J,K-1)
            !pGam=gam(I,J,K-1)
            pGam=1.4    !by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
           CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vkbc)

		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
            AAc =SQRT(1.4*pPP*roinv)    ! by ydd
      !      AAc =SQRT(gam(I,J,K)*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I,J,K-1), Alagm(2,I,J,K-1),Alagm(3,I,J,K-1))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I,J,K-1),ALL,ubetac2)

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)

		 ENDIF

			RDSt=RDS(3,I,J,K-1)
			DRoKm1 =0.5*(DifRo + RDSt*D(1,I,J,K-1))
			DRVxKm1=0.5*(DifRvx+ RDSt*D(2,I,J,K-1))
			DRVyKm1=0.5*(DifRvy+ RDSt*D(3,I,J,K-1))
			DRVzKm1=0.5*(DifRvz+ RDSt*D(4,I,J,K-1))
			DReKm1 =0.5*(DifRe + RDSt*D(5,I,J,K-1))
		ENDIF
!if(ThisBlock%ID_Present_Blk==3.and.KI==13)write(*,*)"KDR",KI,i,j,k,DifRo,D(1:5,i,j,k-1),DRoKm1,DRVxKm1,DRVyKm1,DRVzKm1



!if(ThisBlock%ID_Present_Blk==4)write(*,*)"LUDKbef=",KI,i,j,k,D(1:5,i,j,k)
!        if(ThisBlock%ID_Present_Blk==2)write(*,*)"DJbef=",KI,i,j,k,D(1,i,j,k)
		! update -> for U Op
		!stept = Dtm(I,J,K) 
  	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		D(1,I,J,K)=(D(1,I,J,K)*stept +      (DRoIm1+DRoJm1+DRoKm1)*cvs) / F(1,I,J,K)
		D(2,I,J,K)=(D(2,I,J,K)*stept +      (DRvxIm1+DRvxJm1+DRvxKm1)*cvs) / F(1,I,J,K)
		D(3,I,J,K)=(D(3,I,J,K)*stept +      (DRvyIm1+DRvyJm1+DRvyKm1)*cvs) / F(1,I,J,K)
		D(4,I,J,K)=(D(4,I,J,K)*stept +      (DRvzIm1+DRvzJm1+DRvzKm1)*cvs) / F(1,I,J,K)
		D(5,I,J,K)=(D(5,I,J,K)*stept +      (DReIm1+DReJm1+DReKm1)*cvs) / F(1,I,J,K)
        
           
          
        
          
!if(ThisBlock%ID_Present_Blk==4)write(*,*)"LUDK=",KI,i,j,k,D(1:5,i,j,k),DRoIm1,DRoJm1,DRoKm1

        enddo
	enddo
	enddo




!	U Op      
	do K=NK1,1,-1
	do J=NJ1,1,-1
	do I=NI1,1,-1

		!stept = Dtm(I,J,K) 
	   stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		cvs=stept/Vol(I,J,K) !beta

		!I dir
		IF (I==NI1) THEN !??? why
			DRoIm1 =0.0
			DRvxIm1=0.0
			DRvyIm1=0.0
			DRvzIm1=0.0
			DReIm1 =0.0
		ELSE
			AinRo=D(1,I+1,J,K)
			AinRx=D(2,I+1,J,K)
			AinRy=D(3,I+1,J,K)
			AinRz=D(4,I+1,J,K)
			AinRe=D(5,I+1,J,K)
            
            

			pRo = V(1,I+1,J,K)
			roinv = 1.0/pRo
			pVx = V(2,I+1,J,K)*roinv
			pVy = V(3,I+1,J,K)*roinv
			pVz = V(4,I+1,J,K)*roinv
			pE  = V(5,I+1,J,K)*roinv
			
			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
!           AinP=(gam(I+1,J,K)-1.0)*(AinRe-0.5*AinRo*    &
			AinP=(1.4-1.0)*(AinRe-0.5*AinRo*    &   !by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))

		    aSDx = 0.5*(SD(1,1,I+2,J,K)+SD(1,1,I+1,J,K))
		    aSDy = 0.5*(SD(1,2,I+2,J,K)+SD(1,2,I+1,J,K))
		    aSDz = 0.5*(SD(1,3,I+2,J,K)+SD(1,3,I+1,J,K))
!added by ydd
            vibc=0.5*(vibn(i+2,j,k)+vibn(i+1,j,k))
			pPP = PP(I+1,J,K)
			pTT = T(I+1,J,K)
            !pGam=gam(I+1,J,K)
            pGam=1.4    ! by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
			!--CALL LUSGS_ADQ
            CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vibc)

		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
			AAc =SQRT(1.4*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I+1,J,K), Alagm(2,I+1,J,K),Alagm(3,I+1,J,K))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I+1,J,K),ALL,ubetac2)

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)

		 ENDIF

!if(ThisBlock%ID_Present_blk==2) write(*,*)i,j,k,KI,RDS(1,I+1,J,K),DifRo,D(1,I+1,J,K)
			RDSt=RDS(1,I+1,J,K)
			DRoIm1 =0.5*(DifRo - RDSt*D(1,I+1,J,K))
			DRVxIm1=0.5*(DifRvx- RDSt*D(2,I+1,J,K))
			DRVyIm1=0.5*(DifRvy- RDSt*D(3,I+1,J,K))
			DRVzIm1=0.5*(DifRvz- RDSt*D(4,I+1,J,K))
			DReIm1 =0.5*(DifRe - RDSt*D(5,I+1,J,K))
		ENDIF


		!J dir
		IF (J==NJ1) THEN !??? why
			DRoJm1 =0.0
			DRvxJm1=0.0
			DRvyJm1=0.0
			DRvzJm1=0.0
			DReJm1 =0.0
		ELSE
			AinRo=D(1,I,J+1,K)
			AinRx=D(2,I,J+1,K)
			AinRy=D(3,I,J+1,K)
			AinRz=D(4,I,J+1,K)
			AinRe=D(5,I,J+1,K)

			pRo = V(1,I,J+1,K)
			roinv = 1.0/pRo
			pVx = V(2,I,J+1,K)*roinv
			pVy = V(3,I,J+1,K)*roinv
			pVz = V(4,I,J+1,K)*roinv
			pE  = V(5,I,J+1,K)*roinv
			
			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
            !AinP=(gam(I,J+1,K)-1.0)*(AinRe-0.5*AinRo*        &
            AinP=(1.4-1.0)*(AinRe-0.5*AinRo*        &       !by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))

		    aSDx = 0.5*(SD(2,1,I,J+2,K)+SD(2,1,I,J+1,K))
		    aSDy = 0.5*(SD(2,2,I,J+2,K)+SD(2,2,I,J+1,K))
		    aSDz = 0.5*(SD(2,3,I,J+2,K)+SD(2,3,I,J+1,K))
!added by ydd
            vjbc=0.5*(vjbn(i,j+2,k)+vjbn(i,j+1,k))

			pPP = PP(I,J+1,K)
			pTT = T(I,J+1,K)
            !pGam= gam(I,J+1,K)
            pGam= 1.4       ! by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
			!--CALL LUSGS_ADQ
            CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vjbc)

		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
			AAc =SQRT(1.4*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I,J+1,K), Alagm(2,I,J+1,K),Alagm(3,I,J+1,K))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I,J+1,K),ALL,ubetac2)

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)

		 ENDIF
           
			RDSt=RDS(2,I,J+1,K)
			DRoJm1 =0.5*(DifRo - RDSt*D(1,I,J+1,K))
			DRVxJm1=0.5*(DifRvx- RDSt*D(2,I,J+1,K))
			DRVyJm1=0.5*(DifRvy- RDSt*D(3,I,J+1,K))
			DRVzJm1=0.5*(DifRvz- RDSt*D(4,I,J+1,K))
			DReJm1 =0.5*(DifRe - RDSt*D(5,I,J+1,K))
		ENDIF


		!K dir
		IF (K==NK1) THEN !??? why
			DRoKm1 =0.0
			DRvxKm1=0.0
			DRvyKm1=0.0
			DRvzKm1=0.0
			DReKm1 =0.0
		ELSE
			AinRo=D(1,I,J,K+1)
			AinRx=D(2,I,J,K+1)
			AinRy=D(3,I,J,K+1)
			AinRz=D(4,I,J,K+1)
			AinRe=D(5,I,J,K+1)

			pRo = V(1,I,J,K+1)
			roinv = 1.0/pRo
			pVx = V(2,I,J,K+1)*roinv
			pVy = V(3,I,J,K+1)*roinv
			pVz = V(4,I,J,K+1)*roinv
			pE  = V(5,I,J,K+1)*roinv
			
			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
            !AinP=(gam(I,J,K+1)-1.0)*(AinRe-0.5*AinRo*        &
			AinP=(1.4-1.0)*(AinRe-0.5*AinRo*        &   !by ydd
     &             (pVx**2+pVy**2+pVz**2) - pRo*    &
     &			  (AinX*pVx+AinY*pVy+AinZ*pVz))

		    aSDx = 0.5*(SD(3,1,I,J,K+2)+SD(3,1,I,J,K+1))
		    aSDy = 0.5*(SD(3,2,I,J,K+2)+SD(3,2,I,J,K+1))
		    aSDz = 0.5*(SD(3,3,I,J,K+2)+SD(3,3,I,J,K+1))
!added by ydd
            vkbc=0.5*(vkbn(i,j,k+2)+vkbn(i,j,k+1))

			pPP = PP(I,J,K+1)
			pTT = T(I,J,K+1)
            !pGam= gam(I,J,K+1)
            pGam= 1.4   !by ydd

			!include "..\common\ns_LUSGS_ADQ.for"
			!--CALL LUSGS_ADQ
             CALL LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
                &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
                &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vkbc)

		 IF(IF_PRECONDITION) THEN

		! with Inverse of transform matrix
            AAc =SQRT(1.4*pPP*roinv)    !by ydd
!            AAc =SQRT(gam(I,J,K)*pPP*roinv)
			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz

			ALL = min(Alagm(1,I,J,K+1), Alagm(2,I,J,K+1),Alagm(3,I,J,K+1))
			ubetac2 = 0.0
			CALL getPrecondPara(UVW2,AAc*AAc,pRo,pTT,Rmiu(I,J,K+1),ALL,ubetac2)

			CALL PrecondDW (DifRo,DifRvx,DifRvy,DifRvz,DifRe,pRo,pVx,pVy,pVz,pE,   AAc,ubetac2,pPP,pTT) !,
!     >			AinP,AinX,AinY,AinZ,AinT)

		 ENDIF

			RDSt=RDS(3,I,J,K+1)           
			DRoKm1 =0.5*(DifRo - RDSt*D(1,I,J,K+1))
			DRVxKm1=0.5*(DifRvx- RDSt*D(2,I,J,K+1))
			DRVyKm1=0.5*(DifRvy- RDSt*D(3,I,J,K+1))
			DRVzKm1=0.5*(DifRvz- RDSt*D(4,I,J,K+1))
			DReKm1 =0.5*(DifRe - RDSt*D(5,I,J,K+1))
		ENDIF

		! update
		D(1,I,J,K)=D(1,I,J,K) -     (DRoIm1+DRoJm1+DRoKm1)*cvs/F(1,I,J,K)
		D(2,I,J,K)=D(2,I,J,K) -     (DRvxIm1+DRvxJm1+DRvxKm1)*cvs/F(1,I,J,K)
		D(3,I,J,K)=D(3,I,J,K) -     (DRvyIm1+DRvyJm1+DRvyKm1)*cvs/F(1,I,J,K)
		D(4,I,J,K)=D(4,I,J,K) -     (DRvzIm1+DRvzJm1+DRvzKm1)*cvs/F(1,I,J,K)
		D(5,I,J,K)=D(5,I,J,K) -     (DReIm1+DReJm1+DReKm1)*cvs/F(1,I,J,K)
        
           
	enddo
	enddo
	enddo





!	IF_PRECONDITION=III






END SUBROUTINE TIME_LUSGS







SUBROUTINE LUSGS_ADQ(pGam,pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv &
    &   ,AinRo,AinRx,AinRy,AinRz,AinRe  &
    &   ,DifRo,DifRvx,DifRvy,DifRvz,DifRe,vbcc)
    IMPLICIT NONE
    
    REAL:: pVx,pVy,pVz,aSDx,aSDy,aSDz,pPP,roinv,vbcc
    REAL:: AinRo,AinRx,AinRy,AinRz,AinRe
    REAL:: DifRo,DifRvx,DifRvy,DifRvz,DifRe

    REAL:: VP,Ut2
    REAL:: aSDt
    REAL:: pGam,Gm1,Gm2
    REAL:: AA2,gUt2
    REAL:: GEt1, GEt2


			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz+vbcc
			Ut2  = pVx*pVx+pVy*pVy+pVz*pVz
!			VP          Vc   = Vxi*Dx  + Vyi*Dy  + Vzi*Dz
!			dVP= AinX*aSDx +AinY*aSDy +AinZ*aSDz
!			F(3,I,J,K) = dVP

!			AinRo=D(1,I-1,J,K)
!			AinRx=D(2,I-1,J,K)
!			AinRy=D(3,I-1,J,K)
!			AinRz=D(4,I-1,J,K)
!			AinRe=D(5,I-1,J,K)
			aSDt = 0.0
!			pRo, pVx, pVy, pVz, pE
!			needed? pPP,pTT
!           Real Gas
			Gm1 = pGam - 1.0
			Gm2 = pGam - 2.0
              AA2  = pGam*pPP*roinv
			gUt2 = 0.5*Gm1*Ut2


!
			DifRo =  AinRo * aSDt                   +       &
     &                 AinRx * aSDx                               +     &
     &                 AinRy * aSDy                               +     &
     &                 AinRz * aSDz                               +     &
     &                 AinRe * 0.0

			DifRvx = AinRo * (aSDx*gUt2  - pVx*VP)		        +     &
     &                 AinRx * (aSDt+VP  - aSDx*Gm2*pVx)          +     &
     &                 AinRy * (aSDy*pVx - aSDx*Gm1*pVy)          +     &
     &                 AinRz * (aSDz*pVx - aSDx*Gm1*pVz)          +     &
     &                 AinRe * (           aSDx*Gm1    )

			DifRvy = AinRo * (aSDy*gUt2  - pVy*VP)		        +     &
     &                 AinRx * (aSDx*pVy - aSDy*Gm1*    pVx)          +     &
     &                 AinRy * (aSDt+VP  - aSDy*Gm2*pVy)          +     &
     &                 AinRz * (aSDz*pVy - aSDy*Gm1*pVz)          +     &
     &                 AinRe * (           aSDy*Gm1    )

			DifRvz = AinRo * (aSDz*gUt2  - pVz*VP)		        +     &
     &                 AinRx * (aSDx*pVz - aSDz*Gm1*pVx)          +     &
     &                 AinRy * (aSDy*pVz - aSDz*Gm1*pVy)          +     &
     &                 AinRz * (aSDt+VP  - aSDz*Gm2*pVz)          +     &
     &                 AinRe * (           aSDz*Gm1    )

			GEt1 = 0.5*Gm2*Ut2 - AA2/Gm1
			GEt2 = 0.5*Ut2 + AA2/Gm1

			DifRe  = AinRo * VP* GEt1            		        +     &
     &                 AinRx * (aSDx*GEt2 - Gm1*pVx*VP)           +     &
     &                 AinRy * (aSDy*GEt2 - Gm1*pVy*VP)           +     &
     &                 AinRz * (aSDz*GEt2 - Gm1*pVz*VP)           +     &
     &                 AinRe * (aSDt      + pGam*VP    )


END SUBROUTINE LUSGS_ADQ
