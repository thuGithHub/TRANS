SUBROUTINE CONVECT_AUSMPW_Plus
    USE Global
    IMPLICIT NONE
    
    INTEGER:: I,J,K,L,M
    INTEGER:: NI, NJ, NK
    INTEGER:: NI1, NJ1, NK1
    
    
    REAL:: ULo,URo,VLo,VRo,WLo,WRo,PLo,PRo,TLo,TRo
    REAL:: sav1,sav2,sav3,sav
    REAL:: RmiuT,AlagmR
    REAL:: ga,ga1
    REAL:: DroR,DrxR,DryR,DrzR,DreR
    REAL:: TTL
    REAL:: pGam
    
    NI=ThisBlock%NI
    NJ=ThisBlock%NJ
    NK=ThisBlock%NK
    NI1=ThisBlock%NI1
    NJ1=ThisBlock%NJ1
    NK1=ThisBlock%NK1
    
    
    
!   I direcition
    do k=1,NK1
    do j=1,NJ1
    do i=1,NI
!        VL(ML,-2:MI,-2:MJ,-2:MK)  ! U V W PP T k omega muT gm
       ULo=VL(1,1,i,j,k)
	   URo=VR(1,1,i,j,k)
            
         !!!!!!!!!!!!!!!!!
       VLo=VL(1,2,i,j,k)
	   VRo=VR(1,2,i,j,k)
         !!!!!!!!!!!!!!!!!
       WLo=VL(1,3,i,j,k)
	   WRo=VR(1,3,i,j,k)

       PLo=VL(1,4,i,j,k)
	   PRo=VR(1,4,i,j,k)
        

	   TLo=VL(1,5,i,j,k)
	   TRo=VR(1,5,i,j,k)
        
        
         
       
          
       
         !!!!!!!!!!!!!!!!!


          sav1=SD(1,1,I,J,K)  !aix(i,j,k)
          sav2=SD(1,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(1,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(1,I,J,K)   !aim(i,j,k)+small
		 RmiuT = VL(1,8,I,J,K)
        AlagmR = ALAGM(1,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I-1,J,K))
                
        CALL CONVECT_AUSMPW_Plus_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR)  !output
        
        
        D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I-1,J,K)=D(1,I-1,J,K)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I-1,J,K)=D(2,I-1,J,K)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I-1,J,K)=D(3,I-1,J,K)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I-1,J,K)=D(4,I-1,J,K)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I-1,J,K)=D(5,I-1,J,K)+dreR * sav
           IF(ISNAN(D(1,I,J,K)))then
                write(*,*) I,J,K
                write(*,*) KI, D(1,I,J,K)
                continue
        endif
        
          
    end do
    end do
    end do
    
!   J direction
    do k=1,NK1
    do j=1,NJ
    do i=1,NI1
       ULo=VL(2,1,i,j,k)
	   URo=VR(2,1,i,j,k)
         !!!!!!!!!!!!!!!!!
       VLo=VL(2,2,i,j,k)
	   VRo=VR(2,2,i,j,k)
         !!!!!!!!!!!!!!!!!
       WLo=VL(2,3,i,j,k)
	   WRo=VR(2,3,i,j,k)

       PLo=VL(2,4,i,j,k)
	   PRo=VR(2,4,i,j,k)

	   TLo=VL(2,5,i,j,k)
	   TRo=VR(2,5,i,j,k)
         !!!!!!!!!!!!!!!!!


          sav1=SD(2,1,I,J,K)  !aix(i,j,k)
          sav2=SD(2,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(2,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(2,I,J,K)   !aim(i,j,k)+small
		 RmiuT = VL(2,8,I,J,K)
        AlagmR = ALAGM(2,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I,J-1,K))  !by ydd
        
         CALL CONVECT_AUSMPW_Plus_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR)  !output
         D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I,J-1,K)=D(1,I,J-1,K)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I,J-1,K)=D(2,I,J-1,K)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I,J-1,K)=D(3,I,J-1,K)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I,J-1,K)=D(4,I,J-1,K)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I,J-1,K)=D(5,I,J-1,K)+dreR * sav
    end do
	end do
    end do
    
!   K direciton
    do k=1,NK
    do j=1,NJ1
    do i=1,NI1

         ULo=VL(3,1,i,j,k)
	   URo=VR(3,1,i,j,k)
         !!!!!!!!!!!!!!!!!
         VLo=VL(3,2,i,j,k)
	   VRo=VR(3,2,i,j,k)
         !!!!!!!!!!!!!!!!!
         WLo=VL(3,3,i,j,k)
	   WRo=VR(3,3,i,j,k)

         PLo=VL(3,4,i,j,k)
	   PRo=VR(3,4,i,j,k)

	   TLo=VL(3,5,i,j,k)
	   TRo=VR(3,5,i,j,k)
         !!!!!!!!!!!!!!!!!


          sav1=SD(3,1,I,J,K)  !aix(i,j,k)
          sav2=SD(3,2,I,J,K)  !aiy(i,j,k)
          sav3=SD(3,3,I,J,K)  !aiz(i,j,k)
          sav=Grad(3,I,J,K)   !aim(i,j,k)+small
		 RmiuT = VL(3,8,I,J,K)
        AlagmR = ALAGM(3,I,J,K)
        pGam=1.4    !0.5*(gam(I,J,K)+gam(I,J,K-1))
        CALL CONVECT_AUSMPW_Plus_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR)  !output
        
        D(1,I,J,K)=D(1,I,J,K)-droR * sav
          D(1,I,J,K-1)=D(1,I,J,K-1)+droR *sav

	    D(2,I,J,K)=D(2,I,J,K)-drxR * sav
          D(2,I,J,K-1)=D(2,I,J,K-1)+drxR * sav

		D(3,I,J,K)=D(3,I,J,K)-dryR * sav
          D(3,I,J,K-1)=D(3,I,J,K-1)+dryR * sav

		D(4,I,J,K)=D(4,I,J,K)-drzR * sav
          D(4,I,J,K-1)=D(4,I,J,K-1)+drzR * sav

		D(5,I,J,K)=D(5,I,J,K)-dreR * sav
          D(5,I,J,K-1)=D(5,I,J,K-1)+dreR * sav
          
          IF(ISNAN(D(1,I,J,K)))then
                     write(*,*) I,J,K
                     write(*,*) KI, D(1,I,J,K)
                     continue
          endif


	end do
	end do
	end do    
    
END SUBROUTINE CONVECT_AUSMPW_Plus


SUBROUTINE CONVECT_AUSMPW_Plus_Common(pGam,PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR     &   !input
            &                   ,droR,drxR,dryR,drzR,dreR)  !output

    Use Global
    IMPLICIT NONE

REAL:: PLo,PRo,TLo,TRo,ULo,URo,VLo,VRo,WLo,WRo,sav1,sav2,sav3,sav,RmiuT,AlagmR  !in
REAL:: pGam                     !in  gam
REAL:: droR,drxR,dryR,drzR,dreR !out
REAL:: rcpcv
REAL:: DENLo,DENRo
REAL:: ga,ga1
REAL:: VAMLo,VAMRo,HLo,HRo,acl,acr
REAL:: sav1n,sav2n,sav3n
REAL:: ULn,URn,VVL,VVR,Hnormal,astar_L,astar_R,Acc,MaL,MaR,aas_L,aas_R
REAL:: MaL_plus,MaR_plus,MaL_neg,MaR_neg,pL_plus,pL_neg,pR_plus,pR_neg
REAL:: pS,Mas,wP,fL,fR
REAL:: MMaL_plus,MMaR_neg

       rcpcv=PPF ! rcpcv=Cp-Cv=R
       DENLo=PLo/(rcpcv*TLo)
       DENRo=PRo/(rcpcv*TRo)
       
       sav1n=sav1/sav
       sav2n=sav2/sav
       sav3n=sav3/sav
          
	   ga=pGam
	   ga1=pGam-1.0
       VAMLo=ULo*ULo+VLo*VLo+WLo*WLo   !vaml=ul*ul+vl*vl+wl*wl
       VAMRo=URo*URo+VRo*VRo+WRo*WRo   !vamr=ur*ur+vr*vr+wr*wr
       HLo=PLo*ga/(DENLo*ga1)+0.5*VAMLo !hl=pl*ga/(rl*ga1)+0.5*vaml
       HRo=PRo*ga/(DENRo*ga1)+0.5*VAMRo !hr=pr*ga/(rr*ga1)+0.5*vamr
       acl=sqrt(ga*PLo/DENLo)
       acr=sqrt(ga*PRo/DENRo)
       
       
       !网格法向上投影速度
       ULn=(sav1n*ULo+sav2n*VLo+sav3n*WLo)
       URn=(sav1n*URo+sav2n*VRo+sav3n*WRo)
       !激波法线总焓
       VVL=(ULo*ULo+VLo*VLo+WLo*WLo)-ULn*ULn
       VVR=(URo*URo+VRo*VRo+WRo*WRo)-URn*URn
!       Hnormal=0.5*(HLo-0.5*VVL+HRo-0.5*VVR)
       
       astar_L=sqrt(2*ga1*HLo/(ga+1))
       astar_R=sqrt(2*ga1*HRo/(ga+1))
       aas_L=astar_L**2/max(sqrt(VAMLo),astar_L)
       aas_R=astar_R**2/max(sqrt(VAMRo),astar_R)
       Acc=min(aas_L,aas_R)
       
        !界面声速
!       IF((ULn+URn)>=0)then
!           Acc=astar_L*astar_L/max(abs(ULn),astar_L)
!       ELSE
!           Acc=astar_L*astar_L/max(abs(URn),astar_L)
!       ENDIF
       !特征马赫数
!       Acc=0.5*(acl+acr)
       MaL=ULn/Acc
       MaR=URn/Acc
       
       IF(abs(MaL)>1.0)then
           MaL_plus=0.5*(MaL+abs(MaL))
  !         MaL_neg=0.5*(MaL-abs(MaL))
           pL_plus=0.5*(MaL+abs(MaL))/(MaL+tiny)
  !         pL_neg=0.5*(MaL-abs(MaL))/MaL
       ELSE
           MaL_plus=0.25*(MaL+1)**2
   !        MaL_neg=(-1)*0.25*(MaL-1)**2
           pL_plus=0.25*(MaL+1)**2*(2-MaL)+3/16*MaL*(MaL*MaL-1)**2
    !       pL_neg=0.25*(MaL-1)**2*(2+MaL)-3/16*MaL*(MaL*MaL-1)**2
       ENDIF
       IF(abs(MaR)>1.0)then
  !         MaR_plus=0.5*(MaR+abs(MaR))
           MaR_neg=0.5*(MaR-abs(MaR))
   !        pR_plus=0.5*(MaR+abs(MaR))/MaR
           pR_neg=0.5*(MaR-abs(MaR))/(MaR+tiny)
       ELSE
 !          MaR_plus=0.25*(MaR+1)**2
           MaR_neg=(-1)*0.25*(MaR-1)**2
 !          pR_plus=0.25*(MaR+1)**2*(2-MaR)+3/16*MaR*(MaR*MaR-1)**2
           pR_neg=0.25*(MaR-1)**2*(2+MaR)-3/16*MaR*(MaR*MaR-1)**2
       ENDIF
       pS=pL_plus*PLo+pR_neg*PRo
       Mas=MaL_plus+MaR_neg
       
       wP=1.0-(min(PLo/PRo,PRo/PLo))**3
       IF(abs(MaL)<1)then
           fL=PLo/(pS+tiny)-1
       else
           fL=0
       ENDIF
       IF(abs(MaR)<1)then
           fR=PRo/(pS+tiny)-1
       else
           fR=0
       ENDIF
       IF(Mas>=0)then
           MMaL_plus=MaL_plus+MaR_neg*((1-wP)*(1+fR)-fL)
           MMaR_neg=MaR_neg*wP*(1+fR)
       else
           MMaL_plus=MaL_plus*wP*(1+fL)
           MMaR_neg=MaR_neg+MaL_plus*((1-wP)*(1+fL)-fR)
       endif
       
       droR=(-1)*MMaL_plus*Acc*DENLo-MMaR_neg*Acc*DENRo
       drxR=(-1)*MMaL_plus*Acc*DENLo*ULo-MMaR_neg*Acc*DENRo*URo-sav1n*(pL_plus*PLo+pR_neg*PRo)
       dryR=(-1)*MMaL_plus*Acc*DENLo*VLo-MMaR_neg*Acc*DENRo*VRo-sav2n*(pL_plus*PLo+pR_neg*PRo)
       drzR=(-1)*MMaL_plus*Acc*DENLo*WLo-MMaR_neg*Acc*DENRo*WRo-sav3n*(pL_plus*PLo+pR_neg*PRo)
       dreR=(-1)*MMaL_plus*Acc*DENLo*HLo-MMaR_neg*Acc*DENRo*HRo
    

END SUBROUTINE CONVECT_AUSMPW_Plus_Common


