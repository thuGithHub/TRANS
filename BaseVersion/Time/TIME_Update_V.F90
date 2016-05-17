SUBROUTINE TIME_Update_V
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1

	INTEGER:: I1,J1,K1
	INTEGER:: I2,J2,K2
    INTEGER:: NM
    
    REAL:: Aee,AE1
    REAL:: DRSa1,DRSa2,DRSa3,DRSa4,DRSa5
    REAL:: DUa1,DUa2,DUa3,DUa4,DUa5
    REAL:: DR1,DR2,DR3,DR4,DR5
    REAL:: Ima1,Ima2,Ima3,Ima4,Ima5
    REAL:: Jma1,Jma2,Jma3,Jma4,Jma5
    REAL:: Kma1,Kma2,Kma3,Kma4,Kma5
    
    REAL:: DRSk,DRSo
    REAL:: DUk,DUo
    REAL:: DRk,DRo
    INTEGER:: Imk,Jmk,Kmk,Imo,Jmo,Kmo
    
    REAL:: DRSg
    REAL:: DUg
    REAL:: DRg
    INTEGER:: Img,Jmg,Kmg
    
    REAL:: DRS,DUa,DR,DRSa
    INTEGER:: Ima,Jma,Kma

	REAL:: PPP_Tiny, Den_Tiny, Eng_Tiny
	REAL:: VV
	
	REAL, ALLOCATABLE:: AA(:),BB(:),CC(:,:)

	REAL:: roinv  !1./ro
	REAL:: pRo, pVx, pVy, pVz, pE  !primitive vars
	REAL:: pRx, pRy, pRz, pRe, pPP,pTT
	REAL:: AinRo, AinRx, AinRy, AinRz, AinRe !Flux in rho vars
	REAL:: AinX, AinY, AinZ, AinP    !Flux in primitive vars
    
	REAL:: aSDx, aSDy, aSDz  ! SD in center point
	
    !REAL:: RXM2	
    REAL:: DV
    
    REAL:: AinRTk, AinRTO
    REAL:: ARTk,ARTO
    REAL:: DistN

	REAL:: Wij,Wjk,Wik,Sii,Sjj,Skk,Sij,Sik,Sjk,Vort,Skl
	REAL:: Den,vK,vOmega,TT,RmiuL,Dist,FunC2
	REAL:: Rmut1,Rmut2,Rmut3,aLt,Delta,criterion,Fhybrid,Rmiu_org,Rmiu_oth
	
	REAL:: Dudx,Dudy,Dudz
	REAL:: Dvdx,Dvdy,Dvdz
	REAL:: Dwdx,Dwdy,Dwdz
	REAL:: Div,Votg,prod
	
	REAL:: Cmua,arg2
	
	REAL:: GMXX
	
	REAL:: Rg_LMT_min, Rg_LMT_max


    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1

    NM=max0(NI,NJ,NK)

	ALLOCATE(AA(NM))
	ALLOCATE(BB(NM))
	ALLOCATE(CC(ML,NM))






!	 goto 120


       !smooth

	  Aee=0.666666666666667
!
        NM=max0(NI,NJ,NK)
        BB(1)=1.0+2.0*Aee
        AA(1)=Aee/BB(1)
        DO I=2,NM
           BB(I)=BB(1)-Aee*Aee/BB(I-1)
           AA(I)=Aee/BB(I)
        ENDDO
!!!!!!
        AE1=Aee+1.
!!!!!!
       ! DRS=0.0
        DRSa=0.0
        DUa=0.0

        DO K=1,NK1
           DO J=1,NJ1
              DO I=1,NI1
!                 DO L=1,ML-2
!	              Q(L,I,J,K)=DU_Total(L,I,J,K)
!	           ENDDO
!
                 DR=ABS(D(1,I,J,K)/V(1,I,J,K))
                 
!                open(unit=111,file='rho record.txt')
!              write(111,*) I, J, K, V(5,I,J,K)
 !            continue
         
                 
 !                write(*,*) I,J,K
                 IF(ISNAN(D(1,I,J,K)))then
                     write(*,*) I,J,K,"BlockID=",ThisBlock%ID_Present_blk
                     write(*,*) KI
                     write(*,*)"conservative Vars=",V(1:5,I,J,K)
                     write(*,*)"FLux D=",D(1,I,J,K),i,j,k
                        stop
                     continue
                 endif
                 
                 

                 IF(DUa<=DR)THEN
                    DUa=DR
                    Ima=I
                    Jma=J
                    Kma=K
	           ENDIF
                 DRSa=DRSa+DR*DR
	        ENDDO 
	     ENDDO
        ENDDO  

        DRSa=SQRT(DRSa/NI1/NJ1/NK1)

        if(ISNAN(DRSa)) then
        write(*,*)"Residual of blk",ThisBlock%id_present_blk,"is NAN"
        stop

        endif
!!!!!!
        DO K=1,NK1
           DO J=1,NJ1
              DO L=1,5
                 CC(L,1)=AE1*D(L,1,J,K)
!                 D(L,NI1,J,K)=AE1*D(L,NI1,J,K)
	        ENDDO
              DO I=2,NI1
                 I1=I-1
                 DO L=1,5
                    CC(L,I)=D(L,I,J,K)+AA(I1)*CC(L,I1)
	           ENDDO
	        ENDDO
              DO L=1,5
!                 D(L,NI1,J,K)=CC(L,NI1)/BB(NI1)
	        ENDDO
              DO I1=1,NI1-1
                 I=NI1-I1
                 I2=I+1
                 DO L=1,5
!                    D(L,I,J,K)=CC(L,I)/BB(I)+AA(I)*D(L,I2,J,K)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO
!!!!!!
        DO K=1,NK1
           DO I=1,NI1
              DO L=1,5
                 CC(L,1)=AE1*D(L,I,1,K)
!                 D(L,I,NJ1,K)=AE1*D(L,I,NJ1,K)
	        ENDDO
              DO J=2,NJ1
                 J1=J-1
                 DO L=1,5
                    CC(L,J)=D(L,I,J,K)+AA(J1)*CC(L,J1)
	           ENDDO
	        ENDDO
              DO L=1,5
!                 D(L,I,NJ1,K)=CC(L,NJ1)/BB(NJ1)
	        ENDDO
              DO J1=1,NJ1-1
                 J=NJ1-J1
                 J2=J+1
                 DO L=1,5
!                    D(L,I,J,K)=CC(L,J)/BB(J)+AA(J)*D(L,I,J2,K)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO
!!!!!!
        DO J=1,NJ1
           DO I=1,NI1
              DO L=1,5
                 CC(L,1)=AE1*D(L,I,J,1)
!                 D(L,I,J,NK1)=AE1*D(L,I,J,NK1)
	        ENDDO
              DO K=2,NK1
                 K1=K-1
                 DO L=1,5
                    CC(L,K)=D(L,I,J,K)+AA(K1)*CC(L,K1)
	           ENDDO
	        ENDDO
              DO L=1,5
!                 D(L,I,J,NK1)=CC(L,NK1)/BB(NK1)
	        ENDDO
              DO K1=1,NK1-1
                 K=NK1-K1
                 K2=K+1
                 DO L=1,5
!                    D(L,I,J,K)=CC(L,K)/BB(K)+AA(K)*D(L,I,J,K2)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO


!	 goto 120
       !smooth

	  Aee=0.66666666666667
!
        NM=max0(NI,NJ,NK)
        BB(1)=1.0+2.0*Aee
        AA(1)=Aee/BB(1)
        DO I=2,NM
           BB(I)=BB(1)-Aee*Aee/BB(I-1)
           AA(I)=Aee/BB(I)
        ENDDO
!!!!!!
        AE1=Aee+1.0
!!!!!!
        DRSk=0.0
        DUk=0.0

        DO K=1,NK1
           DO J=1,NJ1
              DO I=1,NI1
!                 DO L=1,ML-2
!	              Q(L,I,J,K)=DU_Total(L,I,J,K)
!	           ENDDO
!
                 DRk=ABS(D(6,I,J,K)/V(1,I,J,K))

                 IF(DUk<=DRk)THEN
                    DUk=DRk
                    Imk=I
                    Jmk=J
                    Kmk=K
	           ENDIF
                 DRSk=DRSk+DRk*DRk
	        ENDDO 
	     ENDDO
        ENDDO  

        DRSk=SQRT(DRSk/NI1/NJ1/NK1)

!!!!!!
        DRSo=0.0
        DUo=0.0

        DO K=1,NK1
           DO J=1,NJ1
              DO I=1,NI1
!                 DO L=1,ML-2
!	              Q(L,I,J,K)=DU_Total(L,I,J,K)
!	           ENDDO
!
                 DRo=ABS(D(7,I,J,K)/V(1,I,J,K))

                 IF(DUo<=DRo)THEN
                    DUo=DRo
                    Imo=I
                    Jmo=J
                    Kmo=K
	           ENDIF
                 DRSo=DRSo+DRo*DRo
	        ENDDO 
	     ENDDO
        ENDDO  

        DRSo=SQRT(DRSo/NI1/NJ1/NK1)

!!!!!!
        DO K=1,NK1
           DO J=1,NJ1
              DO L=6,7
                 CC(L,1)=AE1*D(L,1,J,K)
!                 D(L,NI1,J,K)=AE1*D(L,NI1,J,K)
	        ENDDO
              DO I=2,NI1
                 I1=I-1
                 DO L=6,7
                    CC(L,I)=D(L,I,J,K)+AA(I1)*CC(L,I1)
	           ENDDO
	        ENDDO
              DO L=6,7
!                 D(L,NI1,J,K)=CC(L,NI1)/BB(NI1)
	        ENDDO
              DO I1=1,NI1-1
                 I=NI1-I1
                 I2=I+1
                 DO L=6,7
!                    D(L,I,J,K)=CC(L,I)/BB(I)+AA(I)*D(L,I2,J,K)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO
!!!!!!
        DO K=1,NK1
           DO I=1,NI1
              DO L=6,7
                 CC(L,1)=AE1*D(L,I,1,K)
!                 D(L,I,NJ1,K)=AE1*D(L,I,NJ1,K)
	        ENDDO
              DO J=2,NJ1
                 J1=J-1
                 DO L=6,7
                    CC(L,J)=D(L,I,J,K)+AA(J1)*CC(L,J1)
	           ENDDO
	        ENDDO
              DO L=6,7
!                 D(L,I,NJ1,K)=CC(L,NJ1)/BB(NJ1)
	        ENDDO
              DO J1=1,NJ1-1
                 J=NJ1-J1
                 J2=J+1
                 DO L=6,7
!                    D(L,I,J,K)=CC(L,J)/BB(J)+AA(J)*D(L,I,J2,K)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO
!!!!!!
        DO J=1,NJ1
           DO I=1,NI1
              DO L=6,7
                 CC(L,1)=AE1*D(L,I,J,1)
!                 D(L,I,J,NK1)=AE1*D(L,I,J,NK1)
	        ENDDO
              DO K=2,NK1
                 K1=K-1
                 DO L=6,7
                    CC(L,K)=D(L,I,J,K)+AA(K1)*CC(L,K1)
	           ENDDO
	        ENDDO
              DO L=6,7
!                 D(L,I,J,NK1)=CC(L,NK1)/BB(NK1)
	        ENDDO
              DO K1=1,NK1-1
                 K=NK1-K1
                 K2=K+1
                 DO L=6,7
!                    D(L,I,J,K)=CC(L,K)/BB(K)+AA(K)*D(L,I,J,K2)
	           ENDDO
	        ENDDO
	     ENDDO
	  ENDDO

!	 goto 120


       !smooth

	  Aee=0.6666666666666667
!
        NM=max0(NI,NJ,NK)
        BB(1)=1.0+2.0*Aee
        AA(1)=Aee/BB(1)
        DO I=2,NM
           BB(I)=BB(1)-Aee*Aee/BB(I-1)
           AA(I)=Aee/BB(I)
        ENDDO
!!!!!!
        AE1=Aee+1.
!!!!!!
        DRSg=0.0
        DUg=0.0

!        DO K=1,NK1 ! by ydd
!           DO J=1,NJ1
!              DO I=1,NI1
!                 DO L=1,ML-2
!	              Q(L,I,J,K)=DU_Total(L,I,J,K)
!	           ENDDO
!
!                 DRg=ABS(D(8,I,J,K)/V(1,I,J,K))

!                 IF(DUg<=DRg)THEN
!                    DUg=DRg
!                    Img=I
!                    Jmg=J
!                    Kmg=K
!	           ENDIF
!                 DRSg=DRSg+DRg*DRg
!	        ENDDO 
!	     ENDDO
!        ENDDO  

!        DRSg=SQRT(DRSg/NI1/NJ1/NK1)


!!!!!!
!        DO K=1,NK1
!           DO J=1,NJ1
!              DO L=8,8
!                 CC(L,1)=AE1*D(L,1,J,K)
!                 D(L,NI1,J,K)=AE1*D(L,NI1,J,K)
!	        ENDDO
!              DO I=2,NI1
!                 I1=I-1
!                 DO L=8,8
!                    CC(L,I)=D(L,I,J,K)+AA(I1)*CC(L,I1)
!	           ENDDO
!	        ENDDO
!              DO L=8,8
!!                 D(L,NI1,J,K)=CC(L,NI1)/BB(NI1)
!	        ENDDO
!              DO I1=1,NI1-1
!                 I=NI1-I1
!                 I2=I+1
!                 DO L=8,8
!!                    D(L,I,J,K)=CC(L,I)/BB(I)+AA(I)*D(L,I2,J,K)
!	           ENDDO
!	        ENDDO
!	     ENDDO
!	  ENDDO
!!!!!!
!        DO K=1,NK1
!           DO I=1,NI1
!              DO L=8,8
!                 CC(L,1)=AE1*D(L,I,1,K)
!                 D(L,I,NJ1,K)=AE1*D(L,I,NJ1,K)
!	        ENDDO
!              DO J=2,NJ1
!                 J1=J-1
!                 DO L=8,8
!                    CC(L,J)=D(L,I,J,K)+AA(J1)*CC(L,J1)
!	           ENDDO
!	        ENDDO
!              DO L=8,8
!!                 D(L,I,NJ1,K)=CC(L,NJ1)/BB(NJ1)
!	        ENDDO
!              DO J1=1,NJ1-1
!                 J=NJ1-J1
!                 J2=J+1
!                 DO L=8,8
!!                    D(L,I,J,K)=CC(L,J)/BB(J)+AA(J)*D(L,I,J2,K)
!	           ENDDO
!	        ENDDO
!	     ENDDO
!	  ENDDO
!!!!!!
!        DO J=1,NJ1
!           DO I=1,NI1
!              DO L=8,8
!                 CC(L,1)=AE1*D(L,I,J,1)
!!                 D(L,I,J,NK1)=AE1*D(L,I,J,NK1)
!	        ENDDO
!              DO K=2,NK1
!                 K1=K-1
!                 DO L=8,8
!                    CC(L,K)=D(L,I,J,K)+AA(K1)*CC(L,K1)
!	           ENDDO
!	        ENDDO
!              DO L=8,8
!!                 D(L,I,J,NK1)=CC(L,NK1)/BB(NK1)
!	        ENDDO
!              DO K1=1,NK1-1
!                 K=NK1-K1
!                 K2=K+1
!                 DO L=8,8
!!                    D(L,I,J,K)=CC(L,K)/BB(K)+AA(K)*D(L,I,J,K2)
!	           ENDDO
!	        ENDDO
!	     ENDDO
!	  ENDDO



120    CONTINUE


!!!!!!        
	  Rg_LMT_min=0.0
	  Rg_LMT_max=0.999999999
!!!!!!

      ! update VARs
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
        
       
        
        
		IF (IF_PRECONDITION) THEN !IF_PRECONDITION) THEN !Low Ma
!
!		! with Inverse of transform matrix
!			AAc =SQRT(1.4*PP(I,J,K)*roinv)
!			UVW2 = pVx*pVx+pVy*pVy+pVz*pVz
!
!		 TTL=T(I,J,K)
!		 ALL = min(Alagm(1,I,J,K),Alagm(2,I,J,K),Alagm(3,I,J,K))
!	     ubetac2 = 0.0
!		 CALL getPrecondPara(UVW2,AAc*AAc,pRo,TTL,Rmiu(I,J,K),
!     >				ALL,ubetac2)


!		!K transform
!			AinRe = AinRe - (pVx*AinRx+pVy*AinRy+pVz*AinRz)
!     >			-(pE/1.4-1.2*(pVx*pVx+pVy*pVy+pVz*pVz))*AinRo
!			AinRx = AinRx - pVx*AinRo
!			AinRy = AinRy - pVy*AinRo
!			AinRz = AinRz - pVz*AinRo
!			AinRo = AinRo

!		!beta transform
!			rCp = 1.4/0.4*PPF * pRo
!			!a_det = rCp/(AAc*AAc*ubetac2)   !=rCp/a11

!			a11=AAc*AAc*ubetac2
!			a15=-0.4*ubetac2/rCp  !-0.4/AAc/AAc/a_det
!			b11=1.0*a11/rCp !1.0/a_det
!			b15=(1.0+0.4*ubetac2)/rCp

!			AinP=a11*AinRo+a15*AinRe
!			AinX=AinRx*roinv
!			AinY=AinRy*roinv
!			AinZ=AinRz*roinv
!			AinT=b11*AinRo+b15*AinRe

!			AinRo=AinP/AinT/PPF
!			AinRo=(PP(I,J,K)+AinP)/(T(I,J,K)+AinT)/PPF-pRo


!		if (AinP.ge.0.) then
!			PP(I,J,K)=PP(I,J,K)+AinP
!		else
!			PP(I,J,K)=PP(I,J,K)*exp(AinP/(PP(I,J,K)+1.d-10))
!		endif
!		PP(I,J,K)=ABS(PP(I,J,K)) ! shouldn't be calced using other vars?


!		if (AinT.ge.0.) then
!			T(I,J,K)=T(I,J,K)+AinT
!		else
!			T(I,J,K)=T(I,J,K)*exp(AinT/(T(I,J,K)+1.d-10))
!		endif
!		T(I,J,K)=ABS(T(I,J,K)) ! shouldn't be calced using other vars?

!		T(I,J,K)=PP(I,J,K)/(pRo*rcpcv)  !???
!		RXM2=Xm*Xm*1.4
!          pRo=RXM2*PP(I,J,K)/T(I,J,K)


		pRo = V(1,I,J,K)
		roinv = 1.0/pRo
		pRx = V(2,I,J,K)
		pRy = V(3,I,J,K)
		pRz = V(4,I,J,K)
		pRe = V(5,I,J,K)
		pRo = pRo+AinRo
		pRo = ABS(pRo) ! pRo is new updated value

		pRx = pRx+AinRx
		pRy = pRy+AinRy
		pRz = pRz+AinRz

		pRe = pRe+AinRe
		pRe = ABS(pRe)
        

		DV=0.5*(pRx**2+pRy**2+pRz**2)*roinv
        !real Gas
!        pPP=(pRe-DV)*(gam(I,J,K)-1.0)
        pPP = (pRe-DV)*(1.4-1.0)    !by ydd

		if (pPP .le. 0.0) then
			pPP=ABS(pPP)
            !pRe=pPP/(gam(I,J,K)-1.0)+DV
            pRe=pPP/(1.4-1.0)+DV    !by ydd
		endif

		!RXM2=Xm*Xm*1.4
          pTT=RXM2*PP(I,J,K)/pRo

		V(1,I,J,K)=pRo
		V(2,I,J,K)=pRx
		V(3,I,J,K)=pRy
		V(4,I,J,K)=pRz
		V(5,I,J,K)=pRe
		PP(I,J,K)=pPP
		T(I,J,K)=pTT

   ELSE        !no preconditioned

			AinX=(AinRx-pVx*AinRo)*roinv
			AinY=(AinRy-pVy*AinRo)*roinv
			AinZ=(AinRz-pVz*AinRo)*roinv
			!gamma=1.4
            !Real Gas
            !AinP=(gam(I,J,K)-1.0)*(AinRe-0.5*AinRo*        &   !by ydd
            AinP=(1.4-1.0)*(AinRe-0.5*AinRo*        &
     &          (pVx**2+pVy**2+pVz**2) - pRo*       &
     &		  (AinX*pVx+AinY*pVy+AinZ*pVz))
            
 !       open(unit=112,file='rho record.txt')
  !       write(112,*) I, J, K, AinP
  !       continue

	
		if (AinRo.ge.0.) then
			pRo = pRo+AinRo
		else
			pRo = pRo * exp(AinRo/(pRo + 1.e-10))
		endif
		pRo = ABS(pRo) ! pRo is new updated value

		if (AinP.ge.0.) then
			PP(I,J,K)=PP(I,J,K)+AinP
		else
			PP(I,J,K)=PP(I,J,K)*exp(AinP/(PP(I,J,K)+1.e-10))
		endif
		PP(I,J,K)=ABS(PP(I,J,K)) ! shouldn't be calced using other vars?

!		T(I,J,K)=PP(I,J,K)/(pRo*rcpcv)  !???
		!RXM2=Xm*Xm*1.4
          T(I,J,K)=RXM2*PP(I,J,K)/pRo
          
           !!!!!Real Gas!!!!!
!        IF(IF_RealGas == 1) then   !by ydd
!                TTD=T(I,J,K)*Tinf
!                gam(I,J,K)=(B1+2*B2*TTD+3*B3*TTD**2+4*B4*TTD**3+5*B5*TTD**4)    &
!     &           /(B1-Rconst+2*B2*TTD+3*B3*TTD**2+4*B4*TTD**3+5*B5*TTD**4)
!                if(gam(I,J,K)>1.4)then
!                    gam(I,J,K)=1.4
!                elseif(gam(I,J,K)<1.25)then
!                    gam(I,J,K)=1.25
!                endif
!        ELSEIF(IF_RealGas == 0) then
!            gam(I,J,K)=1.4
!        ENDIF




		pVx = pVx+AinX
		pVy = pVy+AinY
		pVz = pVz+AinZ

		V(1,I,J,K)=pRo

		V(2,I,J,K)=pRo*pVx
		V(3,I,J,K)=pRo*pVy
		V(4,I,J,K)=pRo*pVz


		DV=0.5*(pVx**2+pVy**2+pVz**2)
		
        !V(5,I,J,K)=PP(I,J,K)/(1.4-1.0)+DV*pRo
        !V(5,I,J,K)=PP(I,J,K)/(gam(I,J,K)-1.0)+DV*pRo
        V(5,I,J,K)=PP(I,J,K)/(1.4-1.0)+DV*pRo-0.5*(omega(1)*Rad(i,j,k))**2.0*pRo    !by ydd

		ENDIF


!	  PPP_Tiny = PPf / 25.
!	  Den_Tiny = 1.E-6
!	  Eng_Tiny = 1.E-6
 !	       if(V(1,I,J,K)<=Den_Tiny) V(1,I,J,K)=Den_Tiny
 !	       if(V(5,I,J,K)<=Eng_Tiny) V(5,I,J,K)=Eng_Tiny
!           VV=U(2,I,J,K)*U(2,I,J,K)+U(3,I,J,K)*U(3,I,J,K)+
!     >        U(4,I,J,K)*U(4,I,J,K)
!           PP(I,J,K)=0.4*U(5,I,J,K)-0.2*VV/U(1,I,J,K)
!	     if(PP(I,J,K)<=PPP_Tiny) PP(I,J,K)=PPP_Tiny
!            T(I,J,K)=RXM2*PP(I,J,K)/U(1,I,J,K)



         !Laminar update,  is it necessary?  or is it suitable to update at here?
         !Rmiu(I,J,K) =(1.+Csthlnd)/(T(I,J,K)+Csthlnd)*T(I,J,K)**1.5

	if (IF_euler ==0 .and. IF_turb == 1 .and. KI >= KItm) then  
	else
		GOTO 220
	endif                 



		AinRTk = D(6,I,J,K)
		AinRTO = D(7,I,J,K)

		ARTk = V(6,I,J,K)
		ARTO = V(7,I,J,K)
	
		if (AinRTk.ge.0.0) then
			ARTk = ARTk+AinRTk
		else
			ARTk = ARTk * exp(AinRTk/(ARTk + 1.e-10))
		endif
		ARTk = max(ARTk, 1.d-10) 

		if (AinRTO.ge.0.0) then
			ARTO = ARTO+AinRTO
		else
			ARTO = ARTO * exp(AinRTO/(ARTO + 1.e-10))
		endif
		ARTO = max(ARTO, 1.d-10) 




		pRo = V(1,I,J,K)
		roinv = 1.0/pRo

         !Laminar update,  is it necessary?  or is it suitable to update at here?

!     new  viscosity

			VV=V(2,I,J,K)**2+V(3,I,J,K)**2+V(4,I,J,K)**2
!			PPP=0.4*V(5,I,J,K)-0.2*VV/pRo
            !PPP=(gam(I,J,K)-1.0)*V(5,I,J,K)-0.5*(gam(I,J,K)-1.0)*(VV/pRo- &
            PPP=(1.4-1.0)*V(5,I,J,K)-0.5*(1.4-1.0)*(VV/pRo- &     !by ydd
                &       (omega(1)*rad(i,j,k))**2.0*pRo)
			!RXM2=Xm*Xm*1.4
			TT=RXM2*PPP/pRo
			RmiuL =(1.0+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
			DistN = Dst(I,J,K)





			Dudx = DqDxyz(1,I,J,K)
			Dudy = DqDxyz(2,I,J,K)
			Dudz = DqDxyz(3,I,J,K)
			Dvdx = DqDxyz(4,I,J,K)
			Dvdy = DqDxyz(5,I,J,K)
			Dvdz = DqDxyz(6,I,J,K)
			Dwdx = DqDxyz(7,I,J,K)
			Dwdy = DqDxyz(8,I,J,K)
			Dwdz = DqDxyz(9,I,J,K)

		Div  = DuDx + DvDy + DwDz
		Votg = (dudy-dvdx)**2+(dvdz-dwdy)**2+(dwdx-dudz)**2
	    prod = (2.0*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+    &
     &	 (dudy+dvdx)**2+(dvdz+dwdy)**2+(dwdx+dudz)**2-  &
     &	  2.0/3.0*Div*Div)



   
		
        
        if (Kind_model == 1 ) then !right???


			arg2 = 2.0 * sqrt(ARTk/pRo)/Beta_Star/  &
     &				(ARTo/pRo)/DistN/Ref
			arg2 = max(arg2,    &
     &			500.0*RmiuL/DistN/DistN/ARTO/Ref/Ref)
			arg2 = tanh(arg2**2)

			
			Sii = DuDx
			Sjj = DvDy
			Skk = DwDz
			Sij = (DuDy+DvDx)/2.0
			Sik = (DuDz+DwDx)/2.0
			Sjk = (DvDz+DwDy)/2.0

	!           Vort=SQRT(2.*(2.*Wij*Wij+2.*Wik*Wik+2.*Wjk*Wjk))
		    Skl =SQRT(2.0*(1.0*Sii*Sii+1.0*Sjj*Sjj+1.0*Skk*Skk      &
     &                  +2.0*Sij*Sij+2.0*Sik*Sik+2.0*Sjk*Sjk))


			Cmua = min(1.0, a1*ARTO/pRo/((Skl)*arg2+1.e-10)*Ref) !?
             
             
             
             ARTO=max(ARTO,      &
     &      (pRo*cmua*abs(prod)/(prod_lmt*beta_star))**0.5  &
     &          /Ref   )
		ARTk=min(ARTk,      &   !!!!vmulimit = 10000.   
     &         100000.0   * RmiuL* ARTO/pRo/cmua)

		V(6,I,J,K)=ARTk
		V(7,I,J,K)=ARTO

		Rmiu(I,J,K)=Cmua * ARTk/ARTO*pRo  

!	     else if(isst==2) then
!			  cmua=min(one,sqrt(cmu)*omgg/sqrt(0.5*(prod(i,j,k)+votg(i,j,k))))
		
       elseif (Kind_model == 2) then   !!!!!WD+2006
		    Wij=(DuDy-DvDx)/2.
	        Wik=(DuDz-DwDx)/2.
	        Wjk=(DvDz-DwDy)/2.

            Vort=SQRT(2.*(2.*Wij*Wij+2.*Wik*Wik+2.*Wjk*Wjk))
            
            Sii = DuDx-(DuDx+DvDy+DwDz)/3.
			Sjj = DvDy-(DuDx+DvDy+DwDz)/3.
			Skk = DwDz-(DuDx+DvDy+DwDz)/3.
			Sij = (DuDy+DvDx)/2.0
			Sik = (DuDz+DwDx)/2.0
			Sjk = (DvDz+DwDy)/2.0

	!           Vort=SQRT(2.*(2.*Wij*Wij+2.*Wik*Wik+2.*Wjk*Wjk))
		    Skl =SQRT(1.0*Sii*Sii+1.0*Sjj*Sjj+1.0*Skk*Skk      &
     &                  +2.0*Sij*Sij+2.0*Sik*Sik+2.0*Sjk*Sjk)
            
 !           Cmua = min(1.0,a1*ARTO/pRo/SQRT(Skl*Skl*0.5+Vort*Vort*0.5)*Ref)  
 

  !   limiting turbulence quantities
		ARTO=max(ARTO,      &
     &      (pRo*abs(prod)/(prod_lmt*beta_star))**0.5  &
     &          /Ref   )
		ARTk=min(ARTk,      &   !!!!vmulimit = 10000.   
     &         100000.0   * RmiuL* ARTO/pRo)

		V(6,I,J,K)=ARTk
		V(7,I,J,K)=ARTO
  !      Omg_=amax1(vOmega, 7./8. * SQRT(SF_kl/0.09) / Ref)
         
      
        ARTO =max(ARTO, 7./8. *pRO* SQRT( Skl* Skl*2./0.09) / Ref)

		Rmiu(I,J,K)=ARTk/ARTO*pRo      

	!     vmu(i,j,k)=min(vmu(i,j,k),(vmulimit+1)*vmul(i,j,k))
       end if 


!          if (IF_Transition ==1.and.Kind_model<4) then
!		else
            GOTO 220
!		endif



!		Den = V(1,I,J,K)    !by ydd
!		GMXX = V(8,I,J,K) + D(8,I,J,K)
!		GMXX = GMXX / Den
!
!	     IF(GMXX<=Rg_LMT_min) then
!		    GMXX=Rg_LMT_min
!		 endif
!	     IF(GMXX>=Rg_LMT_max) then
!		    GMXX=Rg_LMT_max
!	     endif
!	    V(8,I,J,K)=GMXX*Den
        
         



220	CONTINUE

	enddo
	enddo
	enddo

    
    ThisBlock%DRSa=DRSa
    ThisBlock%DUa=DUa
    ThisBlock%Ima=Ima
    ThisBlock%Jma=Jma
    ThisBlock%Kma=Kma
    
    ThisBlock%DRSk=DRSk
    ThisBlock%DUk=DUk
    ThisBlock%Imk=Imk
    ThisBlock%Jmk=Jmk
    ThisBlock%Kmk=Kmk

    ThisBlock%DRSo=DRSo
    ThisBlock%DUo=DUo
    ThisBlock%Imo=Imo
    ThisBlock%Jmo=Jmo
    ThisBlock%Kmo=Kmo

!    ThisBlock%DRSg=DRSg
!    ThisBlock%DUg=DUg
!    ThisBlock%Img=Img
!    ThisBlock%Jmg=Jmg
!    ThisBlock%Kmg=Kmg


END SUBROUTINE TIME_Update_V
