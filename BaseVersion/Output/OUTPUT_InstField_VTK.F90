SUBROUTINE OUTPUT_InstField_VTK(filename)
    USE Global
    IMPLICIT NONE
    
    INTEGER*4:: I,J,K,iav,jav,kav
    real::count
    INTEGER*4:: NI, NJ, NK
    INTEGER*4:: NI1, NJ1, NK1

      CHARACTER(LEN=100):: OutputFile,filename
      CHARACTER(LEN=6):: cc

      REAL*4::Den,UU,VV,WW,Uabs,Vabs,Wabs,vK,Mach,Mabs,QQ
      REAL*4::Ts,Ps,P0rel,P0abs,T0rel,T0abs,asonic
      REAL*4::Ox,Oy,Oz
      REAL*4::DuDx,DuDy,DuDz
      REAL*4::DvDx,DvDy,DvDz
      REAL*4::DwDx,DwDy,DwDz
      REAL*4::W12,W13,W23
      REAL*4::S11,S22,S33,S12,S13,S23
      REAL*4::Vort,Skl,Qv
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      REAL*4::Den_av,UU_av,VV_av,WW_av,Uabs_av,Vabs_av,Wabs_av
      real*4::Ma_av,Mabs_av,Ps_av,Ts_av
      REAL*4::Rmiu_av,vK_av
      REAL*4::Ox_av,Oy_av,Oz_av,Vort_av,Qv_av
      REAL*4::angular
      real,allocatable::Var1(:,:,:),Var2(:,:,:),Var3(:,:,:)
      integer::FileID,iBlock,NPoints



!    call CornerPoints
!    do iBlock=1,Max_Block
!        write(*,*)"Processorid=",myid,"BlockID=",ThisBlock%ID_Present_Blk
!        CALL GLOBAL_SetPointers(iBlock)

        NI=ThisBlock%NI
        NJ=ThisBlock%NJ
        NK=ThisBlock%NK
        NI1=ThisBlock%NI1
        NJ1=ThisBlock%NJ1
        NK1=ThisBlock%NK1
        angular=omega(1)

!        OutputFile = trim(ThisBlock%FLPLT)
        Outputfile=trim(Filename)

        IF(Kind_hybrid==0 )THEN
                IF(KI >= KI_c)THEN
                !DO Is=0,999
                 !IF( KI == INT(KI_c+Is*Ks) )THEN
                        !WRITE(cc, '(I3.3)') Is
                               !WRITE(cc, '(I6.6)') KI
                        WRITE(cc, '(I6.6)') (KI-KI_c)/Ksss

                    !WRITE(ThisBlock%FLPLT, '(A,A,A,A)') trim(FilePathPrefix), 'resu/RASST_', cc,".PLT"
                        WRITE(OutputFile, '(A,A,A,A)'), trim(ThisBlock%FLPLT), '_C', cc, ".PLT"
               !ENDIF
              !ENDDO
                ENDIF
        ENDIF
    !WRITE(*,*) 'Output:', OutputFile

        allocate(Var1(NI,NJ,NK))
        allocate(Var2(NI,NJ,NK))
        allocate(Var3(NI,NJ,NK))

        FileID=100
        open(unit=FileID,file=trim(OutputFile))
        write(FileID,'(a)')"# vtk DataFile Version 3.0"
        write(FileID,'(a)')"VTK file for CFD!"
        write(FileID,'(a)')"ASCII"
        ! file header ok
        write(FileID,'(a)')"DATASET STRUCTURED_GRID"
        write(FileID,'(a,3i7)')"DIMENSIONS",NI,NJ,NK
        write(FileID,'(a,i7,a7)')"POINTS",NI*NJ*NK,"float"
        
        do k=1,NK
        do j=1,NJ
        do i=1,NI
                write(FileID,*)XX(i,j,k),YY(i,j,k),ZZ(i,j,k)       ! coordinate
        enddo
        enddo
        enddo
        
        NPoints=(NI-1)*(NJ-1)*(NK-1)
        write(FileID,'(a,i7)')"CELL_DATA",NPoints
        
        write(FileID,'(a)')"SCALARS density float"
        write(FileID,'(a)')"LOOKUP_TABLE default"

        do k=1,NK-1
        do j=1,NJ-1
        do i=1,NI-1                        
                write(FileID,*)V(1,i,j,k)*Roref
        enddo
        enddo
        enddo

        write(FileID,'(a,i7)')"FIELD variables",5

        write(FileID,'(a,2i7,a7)')"relative_velocity",3,NPoints,"float"  !relative and absolute velocity
        do k=1,NK-1
        do j=1,NJ-1
        do i=1,NI-1
                   Den=V(1,i,j,k)
                   UU=V(2,i,j,k)/Den
                   VV=V(3,i,j,k)/Den
                   WW=V(4,i,j,k)/Den
                Var1(i,j,k)=sqrt(UU**2.0+VV**2.0+WW**2.0)
                UU=UU*Vref
                VV=VV*Vref
                WW=WW*Vref
                write(FileID,*)UU,VV,WW
                
        enddo
        enddo
        enddo

        write(FileID,'(a,2i7,a7)')"absolute_velocity",3,NPoints,"float"  !relative and absolute velocity
        do k=1,NK-1
        do j=1,NJ-1
        do i=1,NI-1
                   Den=V(1,i,j,k)
                   Uabs=V(2,i,j,k)/Den
                   Vabs=V(3,i,j,k)/Den-angular*Zc(i,j,k)
                   Wabs=V(4,i,j,k)/Den+angular*Yc(i,j,k)
                Var2(i,j,k)=sqrt(Uabs**2.0+Vabs**2.0+Wabs**2.0)
                Uabs=Uabs*Vref
                Vabs=Vabs*Vref
                Wabs=Wabs*Vref
                write(FileID,*)Uabs,Vabs,Wabs
                
        enddo
        enddo
        enddo
        write(FileID,'(a,2i7,a7)')"Ma",2,NPoints,"float" !relative and absolute Mach number
        
        do k=1,NK-1       
        do j=1,NJ-1
        do i=1,NI-1
                        
                Ps=PP(i,j,k)
                Den=V(1,i,j,k)
                asonic=sqrt(1.4*Ps/Den)
                Mach=Var1(i,j,k)/asonic
                Mabs=Var2(i,j,k)/asonic
                
                Var1(i,j,k)=Mach
                Var2(i,j,k)=Mabs
                Var3(i,j,k)=Ps
                write(FileID,*)Mach,Mabs
        enddo
        enddo
        enddo

        write(FileID,'(a,2i7,a7)')"Pressure",3,NPoints,"float"  ! static, relative total and absolute total
        do k=1,NK-1
        do j=1,NJ-1
        do i=1,NI-1
                Ps=PP(i,j,k)*Pref
                Mach=Var1(i,j,k)
                Mabs=Var2(i,j,k)
                P0rel=Ps*(1.0+0.2*Mach*Mach)**3.5
                P0abs=Ps*(1.0+0.2*Mabs**2.0)**3.5
                write(FileID,*)Ps,P0rel,P0abs
        enddo
        enddo
        enddo

        write(FileID,'(a,2i7,a7)')"Temperature",3,NPoints,"float"! static, relative total and absolute total
        do k=1,NK-1
        do j=1,NJ-1
        do i=1,NI-1
                Ts=T(i,j,k)*Tinf
                Mach=Var1(i,j,k)
                Mabs=Var2(i,j,k)
                T0rel=Ts*(1.0+0.2*Mach*Mach)
                T0abs=Ts*(1.0+0.2*Mabs**2.0)
                write(FileID,*)Ts,T0rel,T0abs
        enddo
        enddo
        enddo

!        write(FileID,'(a,2i7,a7)')"Turbulence",2,NPoints,"float"! turbulence, miuT,k 
!        do k=1,NK-1
!        do j=1,NJ-1
!        do i=1,NI-1
!                Rmiu_av=Rmiu(i,j,k)
!                vk_av=V(6,i,j,k)/V(1,i,j,k)
!                write(FileID,*)Rmiu_av,vK_av
!        enddo
!        enddo
!        enddo
        
!        write(FileID,'(a,2i7,a7)')"Vorticity_and_Q",5,NPoints,"float" !components and modulus
!        do k=1,NK-1
!        do j=1,NJ-1
!        do i=1,NI-1
!                        Den=V(1,i,j,k)
!                        Dudx=DqDxyz(1,i,j,k)
!                        DuDy=DqDxyz(2,i,j,k)
!                        DuDz=DqDxyz(3,i,j,k)
!                        Dvdx=DqDxyz(4,i,j,k)
!                        DvDy=DqDxyz(5,i,j,k)
!                        DvDz=DqDxyz(6,i,j,k)
!                        Dwdx=DqDxyz(7,i,j,k)
!                        DwDy=DqDxyz(8,i,j,k)
!                        DwDz=DqDxyz(9,i,j,k)
        
!                        Ox=-DvDz+DwDy
!                        Oz=-DuDy+DvDx
!                        Oy= DuDz-DwDx
!                        Vort=0.5*sqrt(Ox*Ox+Oy*Oy+Oz*Oz)
!                        W12=Oz/2.0
!                        W13=Oy/2.0
!                        W23=Ox/2.0
!                        S11=DuDx
!                        S22=DvDy
!                        S33=DwDz
!                        S12=(DuDy+DvDx)/2.0
!                        S13=(DUDz+DwDx)/2.0
!                        S23=(DvDz+DwDy)/2.0
!                        Skl=sqrt(2.0*(S11*S11+S22*S22+S33*S33+2.0*S12*S12+2.0*S13*S13+2.0*S23*S23))
!                        Qv=Vort*Vort-Skl*Skl
!                write(FileID,*)Ox,Oy,Oz,Vort,Qv
!        enddo
!        enddo
!        enddo
        close(FileID)
       
        deallocate(Var1)
        deallocate(Var2)
        deallocate(Var3)
!    enddo
    return
end subroutine



