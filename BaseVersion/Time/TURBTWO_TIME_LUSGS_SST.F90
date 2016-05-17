!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE TURBTWO_TIME_LUSGS_SST
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

	!sources
	REAL::	source, sourceP, sourceN1, sourceNc, sourceCross  !sourceN=sourceN1*sourceNc
	REAL::	dudx, dudy, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz
	REAL::	dtkdx, dtkdy, dtkdz, dtodx, dtody, dtodz
	REAL::	Div, Votg
	REAL::	prod, tprod, tprodl1, CnegTk, CnegTO
	REAL::	Cdif, CDKO, arg1

	! DQ(m+1) = A*DeltaQ(m)   ,also used in U op
	REAL::    DRTkIm1, DRToIm1
	REAL::    DRTkJm1, DRToJm1
	REAL::    DRTkKm1, DRToKm1

	REAL:: ARo, ARTk, ARTO
	REAL:: AinRTk, AinRTO
	REAL:: DistN

	REAL:: roinv  !1./ro
	REAL:: pRo, pVx, pVy, pVz  !primitive vars
	REAL:: aSDx, aSDy, aSDz  ! SD in center point
	
	! Use F(3)-F(8) to save DVP&VP in I,J,K dir
	REAL:: DVP, VP  !V of DeltaQ & DeltaV
	REAL:: DifRTk, DifRTo !A(i-1)*DeltaQ(i-1)

	REAL:: Cmua, arg2
	REAL:: Sii,Sjj,Sij,Sik,Sjk, Skl
	REAL:: VV,PPP,TT  !,RXM2


!	REAL*8 prodMat(MI,MJ,MK)




    REAL:: AinRo,AinRx,AinRy,AinRz,AinRe
    REAL:: pE,pTT
    REAL:: AAc,UVW2,ALL,ubetac2

    REAL:: DRoIm1,DRvxIm1,DRvyIm1,DRvzIm1,DeReIm1
    REAL:: AinX,AinY,AinZ,AinP
    REAL:: DifRo,DifRvx,DifRvy,DifRvz,DifRe

    REAL:: RmiuL,GamaAB,BetaAB

    REAL:: AA,Mt

    
    REAL::W(3,3),S(3,3)
    REAL::WWSG,AK_OMG,FBETA,BETA_06,DLTCD
    INTEGER::II,JJ,KK

    REAL::Alf1,Alf2,Alf3
    parameter(Alf1=1.0,Alf2=0.4,Alf3=0.2)

    !!!!!!!!!!!DES

    REAL,EXTERNAL:: FSST,fd_function,fe_function
    REAL::fd,fe
    REAL::RmiuT,Fdes,Altddes,alt,Cdiff_DDES_S,alIDDE,Vomega,Den,Vk,aDst,Fsst_DES,Algrd
    REAL:: alIDDES,almax,almin,Cdes,dlt_IDDES
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!






    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1






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

		 do L=6,7
!		   Q(L,I,J,K) = -(Q(L,I,J,K) - D(L,I,J,K))   ! should be + or -?
	       D(L,I,J,K) = D(L,I,J,K) / Vol(I,J,K)     &
     &   -(3.0*V(L,I,J,K)-4.0*Wt0(L,I,J,K)+Wt1(L,I,J,K))*0.5/ &
     &	DT_LUSGS_Main   !stept
	!   should save for V and Wt1,  Vol, DTM right?  Dtm=step?
	    enddo
	enddo
	enddo
	enddo




if (kind_model ==1) then   !!!!SST

! source and implicit coefficent of negative source items
	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1
		Dudx = DqDxyz(1,I,J,K)
		Dudy = DqDxyz(2,I,J,K)
		Dudz = DqDxyz(3,I,J,K)
		Dvdx = DqDxyz(4,I,J,K)
		Dvdy = DqDxyz(5,I,J,K)
		Dvdz = DqDxyz(6,I,J,K)
		Dwdx = DqDxyz(7,I,J,K)
		Dwdy = DqDxyz(8,I,J,K)
		Dwdz = DqDxyz(9,I,J,K)

		DTkdx = DqDxyz(13,I,J,K)
		DTkdy = DqDxyz(14,I,J,K)
		DTkdz = DqDxyz(15,I,J,K)
		DTOdx = DqDxyz(16,I,J,K)
		DTOdy = DqDxyz(17,I,J,K)
		DTOdz = DqDxyz(18,I,J,K)

		!Prod
		Div  = DuDx + DvDy + DwDz
		Votg = (dudy-dvdx)**2+(dvdz-dwdy)**2+(dwdx-dudz)**2     !OmegaijOmegaij
	    prod = (2.0*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+      &
     &	 (dudy+dvdx)**2+(dvdz+dwdy)**2+(dwdx+dudz)**2-      &
     &	  2.0/3.0*Div*Div)      !prod=P/miuT=TaoijSij/miuT

!		prodmat(I,J,K) = prod

		cdif = dtkdx*dtodx+dtkdy*dtody+dtkdz*dtodz

		tprod = min(prod, votg) * Rmiu(I,J,K)/ Ref

	    !limiter only k-omega?
		ARo = V(1,I,J,K)
		ARTk = V(6,I,J,K)
		ARTO = V(7,I,J,K)
		tprodl1=Prod_lmt*Beta_star*  ARTk*ARTO/ARo *Ref
		tprod=min(tprod,tprodl1)


       TT = T(I,J,K)
       RmiuL=(1.0_8+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
	   RmiuT=Rmiu(i,j,k)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                                                    !
!  IDDES METHOD, 混合方法系数                                                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	Fdes=1._8
	alIDDES=1._8
        IF(Kind_Hybrid == 3 )Then
	          
			   aLt  =SQRT(ARTk/ARo)/(Beta_star*ARTO/ARo*Ref)
!!!!!!
!			   for Cross-Diffusion of DDES-Spalart
!
	           Cdiff_DDES_S=(2.*0.856/Vomega)/Ref*Cdif
	          
			   Den=ARo
			   vk=ARTk/ARo
	           vOmega=ARTO/ARo
           	   aDst  =  Dst(I,J,K)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!DES,DDES(F1,F2)!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
               Fsst_DES = FSST(Den,vK,vOmega,aDst,RmiuL,Cdiff_DDES_S)
	                         
	            aLgrd=((1.-arg1)*Cdes_epsi+arg1*Cdes_omeg)*aLeng(I,J,K)
	               
                Fdes =max(aLt*(1.-Fsst_DES)/aLgrd,1._8)
	               
			FunDES(I,J,K) = aLt*(1.-Fsst_DES)/aLgrd      !!如果大于1，则表示用DES方法
            
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            ! IDDES!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	        If (Kind_subHY == 5 )THEN
	              
				  aLmax =  aLeng(I,J,K)
				  aDst  =  Dst(I,J,K)
				  aLmin =  aSeng(I,J,K)

                  Cdes = (1.-arg1)*Ciddes_epsi +  arg1 * Ciddes_omeg
	              
				  dlt_IDDES=min(max(Cw*aLmax,Cw*aDst,aLmin),aLmax)

				  aLgrd=Cdes * dlt_IDDES
	              
	             
                  fd=0.
	              fe=0.
				  
				  fd =  fd_function (Den, RmiuT,  aDst,  aLmax,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz )

	              fe =  fe_function (Den,RmiuT,RMIUL,aDst, aLmax,dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz   ) 
	              
!	              fee(I,J,K) = fe_func
!	              fdd(I,J,K) = fd_func
!	              fbb(I,J,K) = fb
!                 fdr(I,J,K) = fd_real
!	              ardl(i,j,k) =rdl
!	              ardt(i,j,k)=rdt
!	               afl(i,j,k)=fl
!	               aft(i,j,k)=ft
                    
                    alIDDES= (1. + fe) * fd * aLt  + (1. - fd)* aLgrd
                    
	               ! alhyb(I,J,K)=alIDDES

	           FunDES(I,J,K) = aLt/alIDDES   !!如果大于1，则表示用DES方法
               
               end if


			

	        ENDIF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! compressilbe correction  add by dzw05 20130411  !!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        TT = T(I,J,K)
        RmiuL=(1.0+Csthlnd)/(TT*Csth1+Csthlnd)*(Csth1*TT)**1.5
		DistN = Dst(I,J,K)

	    CDKO=max(2.0*ARo*sigo2/(ARTo/ARo)*Cdif,1.0e-10)
		
          arg1=max(sqrt(ARTk/ARo)/Beta_Star/(ARTO/ARo)/DistN/Ref,500.0*RmiuL/DistN/DistN/ARTO/Ref/Ref)
          arg1=min(arg1, 4.0*sigo2*ARTk/CDKO/DistN/DistN)
          arg1=tanh(arg1**4) !fc1
         
         AA    =SQRT(1.4*PP(I,J,K)/ARo)     !by ydd
    !     AA    =SQRT(gam(I,J,K)*PP(I,J,K)/ARo)
         Mt=sqrt(ARTk/ARo)/AA
	    


SELECT CASE (IF_Comppress)
    CASE(0)

		!k eq source and coeff
		SourceP = tprod -2.0/3.0* min(Div, 0.0) * ARTk
        SourceP=SourceP*f_r1(i,j,k)     !added by ydd for rotation modification
    
	    SourceN1 = -2.0/3.0* max(Div, 0.0) - Beta_star*ARTO/ARo*Ref*Fdes
		SourceNc = ARTk
		Source = SourceP + SourceN1*SourceNc

        CnegTk = abs(SourceN1)

    if(If_fix_Transition)then
	  if(Xc(I,J,K)< Trans_location) Source=0.0_8
	endif

!!!!!!!!!!!!!!!!!!!!IDDES        
if(Kind_Hybrid == 3 .and. Kind_subHY == 5) then
	
	SourceN1 = -2.0_8/3.0_8* max(Div, 0.0_8) -sqrt(ARTK/ARo)/alIDDES

	SourceNc = ARTk
	
	Source = SourceP + SourceN1*SourceNc
	
	CnegTk = abs(-2.0_8/3.0_8* max(Div, 0.0_8) - 1.5*sqrt(ARTK/ARo)/alIDDES)

     if(If_fix_Transition)then
	  if(Xc(I,J,K)< Trans_location) Source=0.0_8
	  endif
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!k eq D
		D(6,I,J,K) = D(6,I,J,K) + Source


        !om eq source and coeff
		tprod = min(prod, votg) / Ref
	         
        GamaAB = Gama1*arg1 + Gama2*(1.0-arg1)
		BetaAB = Beta1*arg1 + Beta2*(1.0-arg1) 
        
        SourceP = GamaAB * ARo * tprod -2.0/3.0*GamaAB* min(Div, 0.0) * ARTO
        SourceP=SourceP*f_r1(i,j,k) !added by ydd for rotation modification
	    SourceN1 =-2.0/3.0*GamaAB* max(Div,0.0)-BetaAB*ARTO/ARo*Ref     !something maybe wrong, by ydd
		SourceNc = ARTo
		Source = SourceP + SourceN1*SourceNc

		SourceCross = 2.0*(1-arg1)*Aro*sigo2/(ARTo/ARo)*Cdif /Ref
		!limiter ????
		SourceCross = sign(1.0, SourceCross) *min(abs(SourceCross), abs(SourceP)*0.5)
		Source = Source + SourceCross

        !om eq D
		D(7,I,J,K) = D(7,I,J,K) + Source

		CnegTO = 2.0 * abs(SourceN1)  
    CASE(1)  !suzen copresss

        !k eq source and coeff
		SourceP = (1.-(1.-arg1)*Alf2*Mt*Mt)*tprod -2.0/3.0* min(Div, 0.0) * ARTk
        SourceP=SourceP*f_r1(i,j,k) !added by ydd for rotation modification
	    SourceN1 = -2.0/3.0* max(Div, 0.0) - (1.+(1.-arg1)*(ALF1-ALF3)*Mt*Mt)*Beta_star*ARTO/ARo*Ref*Fdes
		SourceNc = ARTk

        Source = SourceP + SourceN1*SourceNc
		CnegTk = abs(SourceN1)

     if(If_fix_Transition)then
	  if(Xc(I,J,K)< Trans_location) Source=0.0_8
	endif
  !!!!!!!!!!!!!!!!!!!!IDDES        
      if(Kind_Hybrid == 3 .and. Kind_subHY == 5) then
	
	SourceN1 = -2.0_8/3.0_8* max(Div, 0.0_8) -(1.+(1.-arg1)*(ALF1-ALF3)*Mt*Mt)*sqrt(ARTK/ARo)/alIDDES

	SourceNc = ARTk
	
	Source = SourceP + SourceN1*SourceNc
	
	CnegTk = abs(-2.0_8/3.0_8* max(Div, 0.0_8) - 1.5*sqrt(ARTK/ARo)/alIDDES)
	 if(If_fix_Transition)then
	  if(Xc(I,J,K)< Trans_location) Source=0.0_8
	endif
    
 endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		!k eq D
		D(6,I,J,K) = D(6,I,J,K) + Source
    
    
    !om eq source and coeff
		tprod = min(prod, votg) / Ref
	 

 		
        
        GamaAB = Gama1*arg1 + Gama2*(1.0-arg1)
		BetaAB = Beta1*arg1 + Beta2*(1.0-arg1) 
        
        SourceP = (GamaAB +(1.-arg1)*Alf2*Mt*Mt)* ARo* tprod -2.0/3.0*GamaAB* min(Div, 0.0) * ARTO
        SourceP=SourceP*f_r1(i,j,k) !added by ydd for rotation modification
	    SourceN1 =-2.0/3.0*GamaAB* max(Div,0.0)-(1.-(1.-arg1)*(ALF1-ALF3)*Mt*Mt)*BetaAB*ARTO/ARo*Ref
		SourceNc = ARTo
		Source = SourceP + SourceN1*SourceNc

		SourceCross = 2.0*(1-arg1)*Aro*sigo2/(ARTo/ARo)*Cdif /Ref
		!limiter ????
		SourceCross = sign(1.0, SourceCross) *min(abs(SourceCross), abs(SourceP)*0.5)
		Source = Source + SourceCross

        !om eq D
		D(7,I,J,K) = D(7,I,J,K) + Source
        
        CnegTO = 2.0 * abs(SourceN1)


END SELECT


!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		!diag Matrix
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		!diag for k eq
		F(1,I,J,K) = 1.+(Rds(1,I,J,K)+Rds(2,I,J,K)+Rds(3,I,J,K))    &
     &               *(stept/Vol(I,J,K))+   &
     &				(CnegTk+CnegTo)*stept
		!diag for om eq ??? why CnegTk+CnegTo???
		F(2,I,J,K) = F(1,I,J,K)

	enddo
	enddo
	enddo
    
    
    elseif (kind_model == 2 ) then    !!!!WD+ 2006
    ! source and implicit coefficent of negative source items
	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1
		Dudx = DqDxyz(1,I,J,K)
		Dudy = DqDxyz(2,I,J,K)
		Dudz = DqDxyz(3,I,J,K)
		Dvdx = DqDxyz(4,I,J,K)
		Dvdy = DqDxyz(5,I,J,K)
		Dvdz = DqDxyz(6,I,J,K)
		Dwdx = DqDxyz(7,I,J,K)
		Dwdy = DqDxyz(8,I,J,K)
		Dwdz = DqDxyz(9,I,J,K)

		DTkdx = DqDxyz(13,I,J,K)
		DTkdy = DqDxyz(14,I,J,K)
		DTkdz = DqDxyz(15,I,J,K)
		DTOdx = DqDxyz(16,I,J,K)
		DTOdy = DqDxyz(17,I,J,K)
		DTOdz = DqDxyz(18,I,J,K)

         W(1,1)=0.
		 W(2,2)=0.
		 W(3,3)=0.
	     W(1,2)=0.5*(dUdy - dVdx)
	     W(1,3)=0.5*(dUdz - dWdx)
	     W(2,3)=0.5*(dVdz - dWdy)
	     W(2,1)=-W(1,2)
	     W(3,1)=-W(1,3)
	     W(3,2)=-W(2,3)
!
!		 Sij = (dUi/dxj + dUj/dxi) /2.
!
	     S(1,1)=0.5*(dUdx + dUdx)
	     S(2,2)=0.5*(dVdy + dVdy)
	     S(3,3)=0.5*(dWdz + dWdz)
	     S(1,2)=0.5*(dUdy + dVdx)
	     S(1,3)=0.5*(dUdz + dWdx)
	     S(2,3)=0.5*(dVdz + dWdy)
	     S(2,1)=S(1,2)
	     S(3,1)=S(1,3)
	     S(3,2)=S(2,3)


		!Prod
		Div  = DuDx + DvDy + DwDz
		Votg = (dudy-dvdx)**2+(dvdz-dwdy)**2+(dwdx-dudz)**2
	    prod = (2.0*(dudx*dudx+dvdy*dvdy+dwdz*dwdz)+      &
     &	 (dudy+dvdx)**2+(dvdz+dwdy)**2+(dwdx+dudz)**2-      &
     &	  2.0/3.0*Div*Div)

!		prodmat(I,J,K) = prod


		tprod = min(prod, votg) * Rmiu(I,J,K)  / Ref

	    !limiter only k-omega?
		ARo = V(1,I,J,K)
		ARTk = V(6,I,J,K)
		ARTO = V(7,I,J,K)
		tprodl1=Prod_lmt*Beta_star*  ARTk*ARTO/ARo *Ref
		tprod=min(tprod,tprodl1)

	    
		!k eq source and coeff
		SourceP = tprod -2.0/3.0* min(Div, 0.0) * ARTk
	    SourceN1 = -2.0/3.0* max(Div, 0.0) - Beta_star*ARTO/ARo*Ref
		SourceNc = ARTk
		Source = SourceP + SourceN1*SourceNc
		CnegTk = abs(SourceN1)

		!k eq D
		D(6,I,J,K) = D(6,I,J,K) + Source
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

		!om eq source and coeff
		tprod = min(prod, votg) / Ref


        !!!!!!!!!!!!!!
        WWSg = 0.
	        do II=1,3
	        do JJ=1,3
	        do KK=1,3
		       WWSg = WWSg + W(II,JJ) * W(JJ,KK) * S(KK,II)
	        enddo
	        enddo
	        enddo
			aK_omg  = abs( WWSg / (0.09*ARTO/ARo*Ref)**3. )
			FBeta   = (1.+85.*aK_omg) / (1. + 100.*aK_omg) 
	
           Beta_06 = Beta1 * FBeta 

		cdif = dtkdx*dtodx+dtkdy*dtody+dtkdz*dtodz
			if (cdif <= 0.) dltCD =0.
			if (cdif >  0.) dltCD =0.125
		!!!!!!!!!!!!!!!!!!!
        
        SourceP = Gama1 * ARo * tprod -2.0/3.0*Gama1* min(Div, 0.0) * ARTO
	    SourceN1 =-2.0/3.0*Gama1* max(Div, 0.0) - Beta_06 *ARTO/ARo*Ref

		SourceNc = ARTo
		
        Source = SourceP + SourceN1*SourceNc

!!!!!!!!!!!!!!cross diffusion
        SourceCross = dltCD * ARo*ARo/ARTO /Ref * cdif
		Source = Source + SourceCross
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		CnegTO = 2.0 * abs(SourceN1)

		!om eq D
		D(7,I,J,K) = D(7,I,J,K) + Source
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		!diag Matrix
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		!diag for k eq
		F(1,I,J,K) = 1.+(Rds(1,I,J,K)+Rds(2,I,J,K)+Rds(3,I,J,K))    &
     &               *(stept/Vol(I,J,K))+   &
     &				(CnegTk+CnegTo)*stept
		!diag for om eq ??? why CnegTk+CnegTo???
		F(2,I,J,K) = F(1,I,J,K)

	enddo
	enddo
	enddo

    elseif(kind_model ==3) then !!!!!EARSM

    endif





!     implicit operators
!	temporarily use F as diag


!	L Op      
	do K=1,NK1
	do J=1,NJ1
	do I=1,NI1

		!stept = Dtm(I,J,K) 
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		cvs=stept/Vol(I,J,K) !beta

		!I dir
		IF (I==1) THEN !??? why
			DRTkIm1 =0.0
			DRToIm1 =0.0
		ELSE
			pRo = V(1,I-1,J,K)
			roinv = 1.0/pRo
			pVx = V(2,I-1,J,K)*roinv
			pVy = V(3,I-1,J,K)*roinv
			pVz = V(4,I-1,J,K)*roinv
			
		    aSDx = 0.5*(SD(1,1,I,J,K)+SD(1,1,I-1,J,K))
		    aSDy = 0.5*(SD(1,2,I,J,K)+SD(1,2,I-1,J,K))
		    aSDz = 0.5*(SD(1,3,I,J,K)+SD(1,3,I-1,J,K))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(3,I,J,K)

			DifRTk= D(6,I-1,J,K) *VP !+ V(6,I-1,J,K) *DVP
			DifRTO= D(7,I-1,J,K) *VP !+ V(7,I-1,J,K) *DVP

			RDSt=Rds(1,I-1,J,K)
			DRTkIm1 =0.5*(DifRTk + RDSt*D(6,I-1,J,K))
			DRTOIm1 =0.5*(DifRTo + RDSt*D(7,I-1,J,K))

		ENDIF



		!J dir
		IF (J==1) THEN !??? why
			DRTkJm1 =0.0
			DRToJm1 =0.0
		ELSE
			pRo = V(1,I,J-1,K)
			roinv = 1.0/pRo
			pVx = V(2,I,J-1,K)*roinv
			pVy = V(3,I,J-1,K)*roinv
			pVz = V(4,I,J-1,K)*roinv
			
		    aSDx = 0.5*(SD(2,1,I,J,K)+SD(2,1,I,J-1,K))
		    aSDy = 0.5*(SD(2,2,I,J,K)+SD(2,2,I,J-1,K))
		    aSDz = 0.5*(SD(2,3,I,J,K)+SD(2,3,I,J-1,K))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(4,I,J,K)

			DifRTk= D(6,I,J-1,K) *VP !+ V(6,I,J-1,K) *DVP
			DifRTO= D(7,I,J-1,K) *VP !+ V(7,I,J-1,K) *DVP

			RDSt=Rds(2,I,J-1,K)
			DRTkJm1 =0.5*(DifRTk + RDSt*D(6,I,J-1,K))
			DRTOJm1 =0.5*(DifRTo + RDSt*D(7,I,J-1,K))

		ENDIF



		!K dir
		IF (K==1) THEN !??? why
			DRTkKm1 =0.0
			DRToKm1 =0.0
		ELSE
			pRo = V(1,I,J,K-1)
			roinv = 1.0/pRo
			pVx = V(2,I,J,K-1)*roinv
			pVy = V(3,I,J,K-1)*roinv
			pVz = V(4,I,J,K-1)*roinv
			
		    aSDx = 0.5*(SD(3,1,I,J,K)+SD(3,1,I,J,K-1))
		    aSDy = 0.5*(SD(3,2,I,J,K)+SD(3,2,I,J,K-1))
		    aSDz = 0.5*(SD(3,3,I,J,K)+SD(3,3,I,J,K-1))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(5,I,J,K)

			DifRTk= D(6,I,J,K-1) *VP !+ V(6,I,J,K-1) *DVP
			DifRTO= D(7,I,J,K-1) *VP !+ V(7,I,J,K-1) *DVP

			RDSt=Rds(3,I,J,K-1)
			DRTkKm1 =0.5*(DifRTk + RDSt*D(6,I,J,K-1))
			DRTOKm1 =0.5*(DifRTo + RDSt*D(7,I,J,K-1))

		ENDIF

		! update -> for U Op
		!stept = Dtm(I,J,K) 
	  stept = 2.0*Dtm(I,J,K)/(3.0*Dtm(I,J,K)/DT_LUSGS_Main+2.0)
		D(6,I,J,K)=(D(6,I,J,K)*stept +      &
     &	  (DRTkIm1+DRTkJm1+DRTkKm1)*cvs) / F(1,I,J,K)
		D(7,I,J,K)=(D(7,I,J,K)*stept +      &
     &	  (DRTOIm1+DRTOJm1+DRTOKm1)*cvs) / F(2,I,J,K)

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
			DRTkIm1 =0.0
			DRToIm1 =0.0
		ELSE
			pRo = V(1,I+1,J,K)
			roinv = 1.0/pRo
			pVx = V(2,I+1,J,K)*roinv
			pVy = V(3,I+1,J,K)*roinv
			pVz = V(4,I+1,J,K)*roinv

		    aSDx = 0.5*(SD(1,1,I+2,J,K)+SD(1,1,I+1,J,K))
		    aSDy = 0.5*(SD(1,2,I+2,J,K)+SD(1,2,I+1,J,K))
		    aSDz = 0.5*(SD(1,3,I+2,J,K)+SD(1,3,I+1,J,K))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(6,I,J,K)

			DifRTk= D(6,I+1,J,K) *VP !+ V(6,I+1,J,K) *DVP
			DifRTO= D(7,I+1,J,K) *VP !+ V(7,I+1,J,K) *DVP

			RDSt=RDS(1,I+1,J,K)
			DRTkIm1 =0.5*(DifRTk - RDSt*D(6,I+1,J,K))
			DRTOIm1 =0.5*(DifRTo - RDSt*D(7,I+1,J,K))
		ENDIF


		!J dir
		IF (J==NJ1) THEN !??? why
			DRTkJm1 =0.0
			DRToJm1 =0.0
		ELSE
			pRo = V(1,I,J+1,K)
			roinv = 1.0/pRo
			pVx = V(2,I,J+1,K)*roinv
			pVy = V(3,I,J+1,K)*roinv
			pVz = V(4,I,J+1,K)*roinv

		    aSDx = 0.5*(SD(2,1,I,J+2,K)+SD(2,1,I,J+1,K))
		    aSDy = 0.5*(SD(2,2,I,J+2,K)+SD(2,2,I,J+1,K))
		    aSDz = 0.5*(SD(2,3,I,J+2,K)+SD(2,3,I,J+1,K))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(7,I,J,K)

			DifRTk= D(6,I,J+1,K) *VP !+ V(6,I,J+1,K) *DVP
			DifRTO= D(7,I,J+1,K) *VP !+ V(7,I,J+1,K) *DVP

			RDSt=RDS(2,I,J+1,K)
			DRTkJm1 =0.5*(DifRTk - RDSt*D(6,I,J+1,K))
			DRTOJm1 =0.5*(DifRTo - RDSt*D(7,I,J+1,K))
		ENDIF


		!K dir
		IF (K==NK1) THEN !??? why
			DRTkKm1 =0.0
			DRToKm1 =0.0
		ELSE
			pRo = V(1,I,J,K+1)
			roinv = 1.0/pRo
			pVx = V(2,I,J,K+1)*roinv
			pVy = V(3,I,J,K+1)*roinv
			pVz = V(4,I,J,K+1)*roinv

		    aSDx = 0.5*(SD(3,1,I,J,K+2)+SD(3,1,I,J,K+1))
		    aSDy = 0.5*(SD(3,2,I,J,K+2)+SD(3,2,I,J,K+1))
		    aSDz = 0.5*(SD(3,3,I,J,K+2)+SD(3,3,I,J,K+1))

			VP = pVx *aSDx +pVy *aSDy + pVz*aSDz
!			dVP= F(8,I,J,K)

			DifRTk= D(6,I,J,K+1) *VP !+ V(6,I,J,K+1) *DVP
			DifRTO= D(7,I,J,K+1) *VP !+ V(7,I,J,K+1) *DVP

			RDSt=RDS(3,I,J,K+1)
			DRTkKm1 =0.5*(DifRTk - RDSt*D(6,I,J,K+1))
			DRTOKm1 =0.5*(DifRTo - RDSt*D(7,I,J,K+1))
		ENDIF

		! update
		D(6,I,J,K)=D(6,I,J,K) -     &
     &	  (DRTkIm1+DRTkJm1+DRTkKm1)*cvs/F(1,I,J,K)
		D(7,I,J,K)=D(7,I,J,K) -     &
     &	  (DRTOIm1+DRTOJm1+DRTOKm1)*cvs/F(2,I,J,K)

	enddo
	enddo
	enddo








END SUBROUTINE TURBTWO_TIME_LUSGS_SST


