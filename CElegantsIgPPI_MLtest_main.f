      program HmIgMLTest
      implicit none

      integer P_num
      parameter (P_num=399)
      integer N_num
      parameter (N_num=1872)

      integer num_test
      parameter (num_test=10)
      integer num_bchmk
      parameter (num_bchmk=500)
      integer profile_length
      parameter (profile_length=100)
      integer num_CE
      parameter (num_CE=1081)

      integer mdim,ntrain,nclass,ntest

      parameter(mdim=profile_length,ntrain=num_bchmk,nclass=2,ntest=1)

      integer nsample1
      parameter (nsample1=num_bchmk)
      integer nsample2
      parameter (nsample2=1)

      integer L
      parameter (L=3)	  
	  
      real*8 x(mdim,ntrain),xts(mdim,ntest)
      integer cl(ntrain),clts(ntest)
      integer jests(ntest)
      real*8 qts(nclass,ntest)

      integer labelts,labeltr,nrnn
      integer nscale,nprot,imp,interact
      integer nsample,nrnodes,mimp,near,nprox
      integer ifprot,ifscale,iftest,mdim0,ntest0,nprot0,nscale0

      parameter(labelts=1,labeltr=1,imp=1,interact=0,
     &     nprox=1,nrnn=ntrain,nscale=0,nprot=0)
      parameter(
     &     nsample=(2-labeltr)*ntrain,
     &     nrnodes=2*nsample+1,
     &     mimp=imp*(mdim-1)+1,
     &     ifprot=nprot/(nprot-.1),
     &     ifscale=nscale/(nscale-.1),
     &     iftest=ntest/(ntest-.1),
     &     nprot0=(1-ifprot)+nprot,
     &     nscale0=(1-ifscale)+nscale,
     &     ntest0=(1-iftest)+ntest,
     &     mdim0=interact*(mdim-1)+1,
     &     near=nprox*(nsample-1)+1)

      real*8 input(nsample1,101),output(nsample1,101)
      real*8 input2(nsample2,101),output2(nsample2,101)

      integer N(L)
      real*8 testout(nsample2,101)
      character*1 rec1(100),rec2(100),rec3(100),rec4(100)
      real*8 w(L,101,101)

      integer n2, m, mval, y(2700),Num_ctilde,yts(2700)
      real*8 sigmasq, point(100,2700),ctil(100)
      integer miss(2700)

      real*8 Sensitivity,Specificity,Precision,Accuracy
      real*8 Sensitivity_SV,Specificity_SV,Precision_SV,Accuracy_SV
      real*8 Sensitivity_BP,Specificity_BP,Precision_BP,Accuracy_BP
      real*8 Sensitivity_RF,Specificity_RF,Precision_RF,Accuracy_RF
      integer TP,TN,FP,FN
      integer TP_SV,TN_SV,FP_SV,FN_SV
      integer TP_BP,TN_BP,FP_BP,FN_BP
      integer TP_RF,TN_RF,FP_RF,FN_RF
      integer Neg_Acc,Pos_Acc,Bind_Ind
	  
      integer i,j,k,i2,id,index

      integer bchmk_PN(num_bchmk)
      integer bchmk_num(num_bchmk)
      integer temp_num,temp_PN
      integer check_flag
      real*8 PN_flag
      integer sum_P,sum_N
      integer profile_P(P_num,profile_length)
      integer profile_N(N_num,profile_length)

      real*8 bchmk_input(profile_length,num_bchmk)
      integer bchmk_output(num_bchmk)

      integer CE_input(profile_length,num_CE)
      real*8 CEprofile_input(profile_length,num_CE)
      character*7 Gene(num_CE,2)	
      
      real*8 S(num_CE),R(num_CE),B(num_CE),C(num_CE)
      real*8 FinalScore(num_CE)

      real rand3
      double precision r3
      real rand4
      double precision r4
      real rand5
      double precision r5

      r3=5.0   
      r4=5.0
      r5=5.0         

cccccccccccccccccccccccccccccccccccccccccccccccccc

      open(unit=10,file='HSIPPI_Profile_positiveset.dat',status='old')
      do i=1,P_num
         read(10,*) 
         do k=1,profile_length
            read(10,1007) profile_P(i,k)
         enddo
         read(10,*) 
      enddo
      close(10)

      open(unit=10,file='HSIPPI_Profile_negativeset.dat',
     &     status='old')
      do i=1,N_num
         read(10,*) 
         do k=1,profile_length
            read(10,1007) profile_N(i,k)
         enddo
         read(10,*) 
      enddo
      close(10)

      open(unit=10,file='CElegantsIgPPI_Profile_Final.dat',status='old')
      do i=1,num_CE
         read(10,1006) Gene(i,1),Gene(i,2) 
         do k=1,profile_length
            read(10,1007) CE_input(k,i)
         enddo
         read(10,*) 
      enddo
      close(10)	  

      do i=1,num_CE
         do k=1,profile_length
            CEprofile_input(k,i)=real(CE_input(k,i))/100.0
         enddo
         S(i)=0
         R(i)=0
         B(i)=0
         C(i)=0
      enddo



      do i2=1,num_test

         do j=1,num_bchmk
            do k=1,profile_length
               bchmk_input(k,j)=0.0
            enddo
            bchmk_output(j)=9
         enddo

         do j=1,num_bchmk

            bchmk_PN(j)=9
            bchmk_num(j)=0



            PN_flag=rand3(r3)
            temp_num=0
            temp_PN=9
            if(PN_flag.gt.0.75)then

 1977          continue

               temp_PN=1
               temp_num=int(rand4(r4)*P_num)+1

               check_flag=0
               do k=1,j-1
                  if((temp_PN.eq.bchmk_PN(k)).and.
     &                 (temp_num.eq.bchmk_num(k)))then
                     check_flag=1
                  endif
               enddo
               
               if(check_flag.eq.1)then
                  goto 1977
               else
                  bchmk_PN(j)=temp_PN
                  bchmk_num(j)=temp_num
               endif

            elseif(PN_flag.le.0.75)then

 1978          continue

               temp_PN=0
               temp_num=int(rand5(r5)*N_num)+1

               check_flag=0
               do k=1,j-1
                  if((temp_PN.eq.bchmk_PN(k)).and.
     &                 (temp_num.eq.bchmk_num(k)))then
                     check_flag=1
                  endif
               enddo
               
               if(check_flag.eq.1)then
                  goto 1978
               else
                  bchmk_PN(j)=temp_PN
                  bchmk_num(j)=temp_num
               endif

            endif



c            print*,i,j,bchmk_PN(j),bchmk_num(j)

         enddo

         sum_P=0
         sum_N=0
         do j=1,num_bchmk
            if(bchmk_PN(j).eq.1)then
               sum_P=sum_P+1
            elseif(bchmk_PN(j).eq.0)then
               sum_N=sum_N+1
            endif
         enddo
c         print*,'index',i,sum_P,sum_N

         do j=1,num_bchmk
            if(bchmk_PN(j).eq.1)then
               do k=1,profile_length
                  bchmk_input(k,j)=real(profile_P(bchmk_num(j),k))/100.0
               enddo
               bchmk_output(j)=1
            elseif(bchmk_PN(j).eq.0)then
               do k=1,profile_length
                  bchmk_input(k,j)=real(profile_N(bchmk_num(j),k))/100.0
               enddo
               bchmk_output(j)=0
            endif

c            print*,i2,j,bchmk_input(1,j),bchmk_input(2,j),
c     &           bchmk_input(3,j),bchmk_output(j)
         enddo

		 
cccccccccccccccccccccccccccccccccccccccc
c  set the topology of the BP network
cccccccccccccccccccccccccccccccccccccccc

         N(1)=101
         N(2)=4
         N(3)=1

cccccccccccccccccccccccccccccccccccccccccccccc

         TP=0
         FP=0
         FN=0
         TN=0
         
         TP_SV=0
         FP_SV=0
         FN_SV=0
         TN_SV=0
         
         TP_BP=0
         FP_BP=0
         FN_BP=0
         TN_BP=0
         
         TP_RF=0
         FP_RF=0
         FN_RF=0
         TN_RF=0

         do id=1,num_CE

cccccccccc   SVM

            n2=profile_length
            m=num_bchmk
            mval=1
            
            sigmasq=1.0
            Num_ctilde=50
            
            do i=1,Num_ctilde
               ctil(i)=0.2*i
            enddo
            
            index=0
            do i=1,num_bchmk 
c              if(i.ne.id)then
                  index=index+1
                  do j=1,profile_length
                     point(j,index)=bchmk_input(j,i)
                  enddo
                  if(bchmk_output(i).eq.0)then
                     y(index)=-1
                  elseif(bchmk_output(i).eq.1)then
                     y(index)=1
                  endif
c              endif
            enddo
            
            do j=1,profile_length
               point(j,num_bchmk+1)=CEprofile_input(j,id)
            enddo
c            if(bchmk_output(id).eq.0)then
c               y(num_bchmk)=-1
c            elseif(bchmk_output(id).eq.1)then
               y(num_bchmk+1)=1
c            endif
            
            do i=1,2700
               miss(i)=0
            enddo
            
            do i=1,2700
               yts(i)=0
            enddo
            
            call svm(n2,m,mval,sigmasq,Num_ctilde,ctil,point,y,miss)
            
            if(miss(1).eq.0)then
               yts(1)=1
            elseif(miss(1).eq.1)then
               yts(1)=-1
c            elseif((y(num_bchmk).eq.1).AND.(miss(1).eq.0))then
c               yts(1)=1
c            elseif((y(num_bchmk).eq.1).AND.(miss(1).eq.1))then
c               yts(1)=-1
            endif
            
            do i=1,mval
c               print*,id,miss(i)
            enddo

cccccccccc   Random Forest

            index=0
            do i=1,num_bchmk 
c               if(i.ne.id)then
                  index=index+1
                  do j=1,profile_length
                     x(j,index)=bchmk_input(j,i)
                  enddo
                  if(bchmk_output(i).eq.0)then
                     cl(index)=1
                  elseif(bchmk_output(i).eq.1)then
                     cl(index)=2
                  endif
c               endif
            enddo
            
            do j=1,profile_length
               xts(j,1)=CEprofile_input(j,id)
            enddo
c            if(bchmk_output(id).eq.0)then
c               clts(1)=1
c            elseif(bchmk_output(id).eq.1)then
               clts(1)=1
c            endif
            
            do i=1,ntest
               jests(i)=0
               qts(1,i)=0
               qts(2,i)=0
            enddo
            
            call RandForest(mdim,ntrain,nclass,
     &           ntest,x,cl,xts,clts,
     &           jests,qts,
     &           labelts,labeltr,
     &           nrnn,interact,imp,nprot,nscale,nprox,
     &           nsample,nrnodes,mimp,near,
     &           ifprot,ifscale,iftest,mdim0,ntest0,nprot0,nscale0)
            
            do i=1,ntest
c              print*,id,jests(i),qts(1,i),qts(2,i),
c    &              bchmk_output(id)
            enddo		 
            
ccccccccccccccccccc   BPNN


            do i=1,nsample1
               do j=1,profile_length
                  input(i,j)=0
                  output(i,j)=0
               enddo
            enddo
            
            do i=1,nsample2
               do j=1,profile_length
                  input2(i,j)=0
                  output2(i,j)=0
               enddo
            enddo
            
            index=0
            do i=1,num_bchmk 
c               if(i.ne.id)then
                  index=index+1
                  do j=1,profile_length
                     input(index,j)=bchmk_input(j,i)
                  enddo
                  if(bchmk_output(i).eq.0)then
                     output(index,1)=-1
                  elseif(bchmk_output(i).eq.1)then
                     output(index,1)=1
                  endif
c               endif
            enddo
            
            do j=1,profile_length
               input2(1,j)=CEprofile_input(j,id)
            enddo
c            if(bchmk_output(id).eq.0)then
c               output2(1,1)=-1
c            elseif(bchmk_output(id).eq.1)then
               output2(1,1)=1
c            endif
            
cccccccccccccccccccccccccccccccccccccccccc
c  Initialize weight values
cccccccccccccccccccccccccccccccccccccccccc

            do i=2,L
               do j=1,N(i)
                  do k=1,N(i-1)
                     w(i,j,k)=(rand3(r3)-0.5)*0.2
                  enddo
               enddo
            enddo
            
ccccccccccccccccccccccccccccccccccccccc
c   Learning Process
ccccccccccccccccccccccccccccccccccccccc

            call BPlearn(input,output,N,w,L,nsample1)

ccccccccccccccccccccccccccccccccccccccc
c   Recall Process
ccccccccccccccccccccccccccccccccccccccc

            call BPrecall(input2,output2,testout,N,w,L,nsample2)


            do i=1,nsample2
c               print*,id,testout(i,1),bchmk_output(id)
            enddo

ccccccccccccccccccccccccccccccccccccc

            
ccccccccc   concensus result
            
            Neg_Acc=0
            Pos_Acc=0
            Bind_Ind=0
            
            if(yts(1).eq.-1)then
               S(id)=
     &              S(id)+0			   
               Neg_Acc=Neg_Acc+1
            else
               S(id)=
     &              S(id)+1			
               Pos_Acc=Pos_Acc+1
            endif
            
            if(jests(1).eq.1)then
               R(id)=
     &              R(id)+0			   			   
               Neg_Acc=Neg_Acc+1
            else
               R(id)=
     &              R(id)+1			   			   
               Pos_Acc=Pos_Acc+1
            endif
            
            if(testout(1,1).le.0)then
               B(id)=
     &              B(id)+0			   			   
               Neg_Acc=Neg_Acc+1
            else
               B(id)=
     &              B(id)+1			   
               Pos_Acc=Pos_Acc+1
            endif
            
            if(Neg_Acc.gt.Pos_Acc)then
               Bind_Ind=0
               C(id)=
     &              C(id)+0			   			   
            else
               Bind_Ind=1
               C(id)=
     &              C(id)+1		   			   
            endif

         enddo

ccccccccccccccccccccccccccccccc

ccccccccccccccccccccccccccccccccccccc

      enddo

      open(unit=10,file='CElegantIgPPI_MLtest_Result.dat',
     &     status='unknown')
      do id=1,num_CE
         S(id)=
     &        S(id)/real(num_test)
         R(id)=
     &        R(id)/real(num_test)
         B(id)=
     &        B(id)/real(num_test)
         C(id)=
     &        C(id)/real(num_test)

         FinalScore(id)=(S(id)+R(id)+B(id)+C(id))/4.0
         write(10,1008) "Predicted PPI for ",Gene(id,1)," and ",
     &        Gene(id,2)
     &        ,'SVM Score: ',S(id),
     &        'RF Socre: ',R(id),'BPNN Score: ',B(id),
     &        'Concensus Score: ',C(id),'Final Score: ',
     &        FinalScore(id)
      enddo
      close(10)
 1006 format(19x,A7,5x,A7)
 1007 format(6x,I5)
 1008 format(A18,A7,A5,A7,3x,A11,F8.3,3x,A10,F8.3,3x,A12,F8.3,3x,
     &     A15,F8.3,3x,A13,F8.3)

      stop
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand3(r3)
      double precision s,u,v,r3
      s=65536.0
      u=2053.0
      v=13849.0
      m=r3/s
      r3=r3-m*s
      r3=u*r3+v
      m=r3/s
      r3=r3-m*s
      rand3=r3/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc  
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand4(r4)
      double precision s,u,v,r4
      s=65536.0
      u=2053.0
      v=13849.0
      m=r4/s
      r4=r4-m*s
      r4=u*r4+v
      m=r4/s
      r4=r4-m*s
      rand4=r4/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                      
ccccccccccccccccccccccccccccccccccccccccccccccccc

      real  function rand5(r5)
      double precision s,u,v,r5
      s=65536.0
      u=2053.0
      v=13849.0
      m=r5/s
      r5=r5-m*s
      r5=u*r5+v
      m=r5/s
      r5=r5-m*s
      rand5=r5/s
      return
      end

ccccccccccccccccccccccccccccccccccccccccccccccccc                                  
