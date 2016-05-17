SUBROUTINE TIME_RK_V_Update_New
       USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    
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
    
    
    INTEGER:: I1,J1,K1
	INTEGER:: I2,J2,K2
    INTEGER:: NM
    REAL:: roinv  !1./ro
	REAL:: pRo, pVx, pVy, pVz, pE  !primitive vars
	REAL:: pRx, pRy, pRz, pRe, pPP,pTT,pTTv
    REAL:: pY1,pY2,pY3,pY4,pY5,pY6,pEv
    REAL:: cve1,cve2,cve3,cve4,cve5,cve6
    INTEGER:: II,JJ
    REAL:: RsT,Tvem1,Tvem,dEv,Eve_kn,Tm,Tm1,Tm2,Tcp,Ths1,Ths2,Ths,En_all
    REAL:: En(6)
    REAL:: pRo1,pRo2,pRo3,pRo4,pRo5,pRo6,pPn
    REAL:: aSDx, aSDy, aSDz  ! SD in center point
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

    
     Rg_LMT_min=0.0
	  Rg_LMT_max=0.999999999
     
    
     do K=1,NK1
	 do J=1,NJ1
	 do I=1,NI1
         
         pRo = V(1,I,J,K)
         roinv = 1.0/pRo
         pRx = V(2,I,J,K)
        pRy = V(3,I,J,K)
        pRz = V(4,I,J,K)
        pRe = V(5,I,J,K)
        
		pVx = V(2,I,J,K)*roinv
		pVy = V(3,I,J,K)*roinv
		pVz = V(4,I,J,K)*roinv
		pE  = V(5,I,J,K)*roinv
        
        DV=0.5*(pVx**2+pVy**2+pVz**2)
        !           gam(I,J,K)=1.4  by ydd
         
            VV=V(2,I,J,K)*V(2,I,J,K)+V(3,I,J,K)*V(3,I,J,K)+V(4,I,J,K)*V(4,I,J,K)
            !PP(I,J,K)=(gam(I,J,K)-1.0)*V(5,I,J,K)-0.5*(gam(I,J,K)-1.0)*VV/V(1,I,J,K)
            PP(I,J,K)=(1.4-1.0)*V(5,I,J,K)-0.5*(1.4-1.0)*VV/V(1,I,J,K)  ! by ydd
            T(I,J,K)=RXM2*PP(I,J,K)/V(1,I,J,K)
           
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
        
            VV=V(2,I,J,K)**2+V(3,I,J,K)**2+V(4,I,J,K)**2
            PPP=0.4*V(5,I,J,K)-0.2*VV/pRo       !by ydd
            !PPP=(gam(I,J,K)-1.0)*V(5,I,J,K)-0.5*(gam(I,J,K)-1.0)*VV/pRo
			!RXM2=Xm*Xm*1.4
			TT=RXM2*PPP/pRo
			RmiuL =(1.0+Csthlnd)/(TT+Csthlnd)*TT**1.5
        
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


			Cmua = min(1.0 , a1*ARTO/pRo/((Skl)*arg2+1.e-10)*Ref) !?
             
            
             
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



!		Den = V(1,I,J,K)        !by ydd
!		GMXX = V(8,I,J,K) + D(8,I,J,K)
!		GMXX = GMXX / Den

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
if(ISNAN(DRK)) then
        write(*,*)"k NaN,BLock=",ThisBlock%ID_Present_Blk,i,j,k
endif
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
        
!        DRSg=0.0      ! by ydd
!        DUg=0.0
!        DO K=1,NK1
!           DO J=1,NJ1
!              DO I=1,NI1
!                 DO L=1,ML-2
!	              Q(L,I,J,K)=DU_Total(L,I,J,K)
!               ENDDO
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
        
ENDSUBROUTINE



!SUBROUTINE SolveTem_New(pY1,pY2,pY3,pY4,pY5,pY6,Tm,pRo,pRe,DV,R_aveR,   Ths)
!    USE Global
!    IMPLICIT NONE
!    
!    REAL :: pY1,pY2,pY3,pY4,pY5,pY6,Tm,pRo,pRe,DV,R_aveR,Tm1
!    REAL :: Ths,En_all,Cp_all
!    INTEGER :: II
    
!    Tm1=Tm*Tinf
!    do II=1,6
!                Cps(II)=R/MW(II)*(CpA1(II)+CpA2(II)*Tm1+CpA3(II)*Tm1**2+CpA4(II)*Tm1**3+CpA5(II)*Tm1**4+CpA6(II)*Tm1**5+CpA7(II)*Tm1**6+CpA8(II)*Tm1**7)
!                Cps(II)=Cps(II)/(RXM2*R_avef)
!                Hs(II)=R/MW(II)*Tm1*(CpA1(II)+CpA2(II)*Tm1/2.0+CpA3(II)*Tm1**2/3.0+CpA4(II)*Tm1**3/4.0+CpA5(II)*Tm1**4/5.0+CpA6(II)*Tm1**5/6.0+CpA7(II)*Tm1**6/7.0+CpA8(II)*Tm1**7/8.0)+h0(II)
!                Hs(II)=Hs(II)/(RXM2*R_avef*Tinf)
!    enddo
!    
!    En_all=pY1*(Tm*Cps(1)-Hs(1))+pY2*(Tm*Cps(2)-Hs(2))+pY3*(Tm*Cps(3)-Hs(3))+pY4*(Tm*Cps(4)-Hs(4))+pY5*(Tm*Cps(5)-Hs(5))+pY6*(Tm*Cps(6)-Hs(6))
!    Cp_all=pY1*Cps(1)+pY2*Cps(2)+pY3*Cps(3)+pY4*Cps(4)+pY5*Cps(5)+pY6*Cps(6)
!     En_all=(Yf(1)*Hs(1)+Yf(2)*Hs(2)+Yf(3)*Hs(3)+Yf(4)*Hs(4)+Yf(5)*Hs(5)+Yf(6)*Hs(6))
!           
!        Ths=abs((En_all-DV+pRe/pRo)/(Cp_all-R_aveR/RXM2))
!        
!        continue
!    
!    
!END SUBROUTINE
