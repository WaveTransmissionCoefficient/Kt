close all
clc
clear

%-------------------------IMPORT DATA------------------------%

data=importdata('Database.csv');
database=data.data;
textdata=data.textdata;
textdata(1,:)=[];
TestID=textdata(:,2);
TestID=cell2table(TestID);


% load dataset
ds1=database(1:13,:); 

   
data = {ds1}; % cell array

%% General parameters 
for k = 1:numel(data) % indices  
    
   for n = 1:size(data{1,k})
     %  sop(k,n)=(2*pi*data{1,k}(n,11))./(9.81*(data{1,k}(n,13))^2); 
        L(k,n) = wavelength(data{1,k}(n,13),data{1,k}(n,24));      %wavelenght 
        sop(k,n)=data{1,k}(n,11)./L(k,n);                          %wave steepness 
        Ir(k,n) = (data{1,k}(n,48))^(-1)*(sop(k,n))^(-0.5);        %surf similarity parameter (Iribbaren)  

   end 
end 
   
%% Transmission Van der Meer (1991) (CIRIA/CUR)
% controllare il plot con Rc/Hmo, qualcosa non convince
for k = 1:numel(data) % indices  

for n = 1:size(data{1,k})
    RH(k,n)= data{1,k}(n,29)./data{1,k}(n,11); %Rc/Hmo 
 %calculation Kt= 0.46 - 0.3(Rc/Hm0)
    Kt0(k,n) = 0.46-0.3*(data{1,k}(n,29)./data{1,k}(n,11)); 
 %Range validity on Rc/Hmo  
    if data{1,k}(n,29)./data{1,k}(n,11)>-1.13 && data{1,k}(n,29)./data{1,k}(n,11)<1.2
          KtC(k,n) = Kt0(k,n);    
    elseif data{1,k}(n,29)./data{1,k}(n,11)>-2 && data{1,k}(n,29)./data{1,k}(n,11) <-1.13
          KtC(k,n)= 0.80;           
    elseif data{1,k}(n,29)./data{1,k}(n,11)>1.2 && RH(k,n)<2
          KtC(k,n)= 0.10;
     else 
          KtC(k,n)= nan;    
    end 
          
%     if KtC(k,n)==0.46
%           KtC(k,n) =nan;     
     end 
end 
    
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtC(k,n)==0 
        KtC(k,n)=nan;
end
end
end


%% Transmission Van der Meer and Daemen (1994) 
for k = 1:numel(data) % indices  
for n = 1:size(data{1,k})
%a parameter of Vdm-Daemen (a3)  
    if data{1,k}(n,11)./data{1,k}(n,35) > 1 && data{1,k}(n,11)./data{1,k}(n,35)< 6
        a3(k,n) = 0.031*(data{1,k}(n,11)./data{1,k}(n,35))-0.24;  
    else 
        a3(k,n)= nan;
    end
%b parameter of Vdm-Daemen (b3)
    if  data{1,k}(n,11)./data{1,k}(n,35) > 1 && data{1,k}(n,11)./data{1,k}(n,35)< 6 && sop(k,n) <0.05 && sop(k,n)> 0.01 
        b3(k,n) = -5.42*sop(k,n)+0.0323*data{1,k}(n,11)./data{1,k}(n,35)-0.0017*(data{1,k}(n,43)./data{1,k}(n,35))^(1.84)+0.51;   
    else 
        b3(k,n)= nan;
    end 
% calculation Kt (no range)   
    if data{1,k}(n,35)~=0
        KtVD(k,n)= b3(k,n)+a3(k,n)*data{1,k}(n,29)./data{1,k}(n,35);
    else
        KtVD(k,n)=nan;
    end
%range validity on Kt     
    if KtVD(k,n)>0.75
        KtVD(k,n)=0.75;
    elseif KtVD(k,n)< 0.075
        KtVD(k,n)= 0.075;
    end 
end
end 
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtVD(k,n)==0 
        KtVD(k,n)=nan;
end
end
end

%% Transmission Seabrook (2000)
%valid for submerged only 

for k=1:numel(data) %indices
    
for n = 1:size(data{1,k})      
        HB(k,n)= data{1,k}(n,11)./data{1,k}(n,43); %Hmo/B
        BRld(k,n)= (data{1,k}(n,44)*abs(data{1,k}(n,29)))./(L(k,n)*data{1,k}(n,35)); %BRc/LD50 (-Rc) 
        HRbd(k,n)= (data{1,k}(n,11)*abs(data{1,k}(n,29))./(data{1,k}(n,44)*data{1,k}(n,35))); %HRc/BD50 (-Rc)
        
   if BRld(k,n)<0 && BRld(k,n)>7.08
       BRld(k,n)=nan;
   end 
    
    if HRbd(k,n)<0 && HRbd(k,n)>2.14
       HRbd(k,n)= nan;
    end 
    
    if data{1,k}(n,29)<0 
        KtSea(k,n) = 1-(exp(-0.65*abs(data{1,k}(n,29))./data{1,k}(n,11)-1.09*HB(k,n))+0.047*BRld(k,n)-0.067*HRbd(k,n));
    else 
        KtSea(k,n) = nan;
    end
    if KtSea (k,n)>1
        KtSea (k,n)=nan;
    end 
end            
end

% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtSea(k,n)==0 
        KtSea(k,n)=nan;
end
end
end
%% Transmission Calabrese,Buccino, Vicinanza (2002) 

for k = 1:numel(data) % indices  
for n=1:size(data{1,k})  
%calculate b (b2) for Kt=a Rc/Bc + b
    if Ir(k,n)>=3 && Ir(k,n)<=5.2                                 %range Iribarren
        alfa(k,n)= 1-0.562*exp(-0.0507*Ir(k,n));
    else 
        alfa(k,n) =nan;
    end 
    if  data{1,k}(n,44)./data{1,k}(n,11)>=1.06 && data{1,k}(n,44)./data{1,k}(n,11)<=8.13 %range B/Hmo
        b2(k,n)= alfa(k,n)*exp(-0.0845*data{1,k}(n,43)./data{1,k}(n,11));
    else 
        b2(k,n) = nan;
    end 
%calculate a for Kt=a Rc/Bc + b
    if  data{1,k}(n,11)./data{1,k}(n,24)>=0.31 && data{1,k}(n,11)./data{1,k}(n,24)<=0.61 %range Hmo/h 
        beta(k,n)= 0.6957*data{1,k}(n,11)./data{1,k}(n,24)-0.7021;
    else 
        beta(k,n) = nan;
    end 

    if data{1,k}(n,44)./data{1,k}(n,11)>=1.06 && data{1,k}(n,44)./data{1,k}(n,11)<=8.13  %range B/Hmo
        a2(k,n) = beta(k,n)*exp(0.2568*data{1,k}(n,43)./data{1,k}(n,11));
    else
        a2(k,n)=nan;
    end 
%calcute Kt 
    if data{1,k}(n,29)./data{1,k}(n,44)>=-0.4 && data{1,k}(n,29)./data{1,k}(n,44)<=0.3 %range Rc/B
        KtCal(k,n)=a2(k,n)*data{1,k}(n,29)./data{1,k}(n,43)+b2(k,n);
    else 
        KtCal(k,n)=nan;
    end 
end
end 
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtCal(k,n)==0 
        KtCal(k,n)=nan;
end
end
end

%% Transmission Buccino (Conceptual approach 2007)
%for submerged only 
for k = 1:numel(data) % indices  
 for n=1:size(data{1,k})
         %parametro B/H*sop^0.5          
            C(k,n) = data{1,k}(n,44)./(data{1,k}(n,11)*L(k,n))^0.5; 
         %sol.1
             F1(k,n) = (1.18*(abs(data{1,k}(n,11)./data{1,k}(n,29)))^0.12+0.33*(abs(data{1,k}(n,11)./data{1,k}(n,29)))^1.5*C(k,n))^-1; 
         %sol.2
             F2(k,n) = (min(0.74,0.62*Ir(k,n)^0.17)-0.25*min(C(k,n),2.2))^2;
         %sol.3    
             F3(k,n) = (1.18*1.2^0.12+0.33*1.2^1.5*C(k,n))^-1;
        %interpolazione
            S1(k,n) = F2(k,n)+(F3(k,n)-F2(k,n))*(1/1.2)^(-1)*(abs(data{1,k}(n,11)./data{1,k}(n,29)))^-1;   
   %Kteffettivo Buccino
    if abs(data{1,k}(n,11)./data{1,k}(n,29))<1.2 
        KtB(k,n) = F1(k,n);
    else 
        KtB(k,n) = S1(k,n);
    end     
    
    if data{1,k}(n,29)>0
        KtB(k,n)=  nan;
    end
    
 end 
end
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtB(k,n)==0 
        KtB(k,n)=nan;
end
end
end


%% Transmission D'Agremond (1996)
 
for k=1:numel(data) %indices              
% calculate Kt 
 for n=1:size(data{1,k})
    if data{1,k}(n,35)==0
        KtDAG(k,n) = -0.4*data{1,k}(n,29)./data{1,k}(n,11)+0.8*(data{1,k}(n,44)./data{1,k}(n,11))^-0.31*(1-exp(-0.5*Ir(k,n))); 
    else
        KtDAG(k,n) = -0.4*data{1,k}(n,29)./data{1,k}(n,11)+0.64*(data{1,k}(n,44)./data{1,k}(n,11))^-0.31*(1-exp(-0.5*Ir(k,n))); 
    end
 end 
%range of Kt max/min
 for n=1:size(data{1,k})
    if KtDAG(k,n)>0.8
        KtDAG(k,n)=0.8;
    elseif KtDAG(k,n)<0.075
        KtDAG(k,n)=0.075;    
    end  
 end
 
end 
 % final check on zero-values
 for k = 1:numel(data) % indices     
 for n = 1:size(data{1,k})
     if KtDAG(k,n)==0 
         KtDAG(k,n)=nan;
 end
 end
 end

 
%% Briganti (DELOS 2003)
% general parameters 
for k=1:numel(data) %indices    

% calculate Kt 
 for n=1:size(data{1,k})
         if  (data{1,k}(n,43))./(data{1,k}(n,11))<=10  %B/Hmo<10 d'agremond
            KtDEL(k,n) =-0.4*data{1,k}(n,29)./data{1,k}(n,11)+0.8*(data{1,k}(n,44)./data{1,k}(n,11))^-0.31*(1-exp(-0.5*Ir(k,n)));
         elseif (data{1,k}(n,43))./(data{1,k}(n,11))<=10  && data{1,k}(n,35)==0 %SMOOTH 
             KtDEL(k,n) = -0.4*data{1,k}(n,29)./data{1,k}(n,11)+0.8*(data{1,k}(n,44)./data{1,k}(n,11))^-0.31*(1-exp(-0.5*Ir(k,n)));
        elseif data{1,k}(n,43)./data{1,k}(n,11)>10 %BRIGANTI B/Hmo >10 
          KtDEL(k,n) =  -0.35*data{1,k}(n,29)./data{1,k}(n,11)+0.51*(data{1,k}(n,44)./data{1,k}(n,11))^-0.65*(1-exp(-0.41*Ir(k,n)));
        else 
            KtDEL(k,n)= nan;
        end
%Kt max (delos)
        KtmaxBR(k,n) = -0.006*(data{1,k}(n,44)./data{1,k}(n,11))+0.93;
%kt effettivo
        if KtDEL(k,n)> KtmaxBR(k,n)
            KtDEL(k,n)= KtmaxBR(k,n);
        elseif KtDEL(k,n)< 0.05
            KtDEL(k,n)= 0.05;
        end 
 end
end 
 for k = 1:numel(data) % indices     
 for n = 1:size(data{1,k})
     if KtDEL(k,n)==0 
         KtDEL(k,n)=nan;
 end
 end
 end

    
%% Transmission Van der Meer - smooth and oblique wave attack (2005) 

for k = 1:numel(data) % indices        
 for n=1:size(data{1,k})
     %Calculate Kt without validity range 
%        if data{1,k}(n,35)== 0 
 KtV5(k,n) = -0.3*(data{1,k}(n,29)./data{1,k}(n,11))+0.75*(1-exp(-0.5*Ir(k,n)))*cos(data{1,k}(n,9))^2/3;
%          elseif 
%              KtV5(k,n) = -0.3*(data{1,k}(n,29)./data{1,k}(n,11))+0.75*(1-exp(-0.5*Ir(k,n)))*cos(data{1,k}(n,9))^2/3;
%         else 
%             KtV5(k,n)= nan;
 
        
      if  Ir(k,n)>3 && Ir(k,n)<1  && data{1,k}(n,44)./data{1,k}(n,11)>4 && data{1,k}(n,44)./data{1,k}(n,11)<1
           KtV5(k,n)=nan; 
      end 
   % Calculate final Kt  
        if KtV5(k,n) > 0.8
            KtV5(k,n)=0.8;
        elseif KtV5(k,n) < 0.075 
            KtV5(k,n)=0.075;
        end
 end 
end 
   
%  end
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtV5(k,n)==0 
        KtV5(k,n)=nan;
end
 end
 end
 
 
%% Transmission Goda & Ahrens (2008)
 for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
        if data{1,k}(n,29) > 0
            a(k,n) = 0.248*exp(-0.348*(log(data{1,k}(n,44)./data{1,k}(n,16))));
        elseif data{1,k}(n,29) == 0
            a(k,n) = 0.248*exp(-0.348*(log(((9*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./10)./data{1,k}(n,16))));
        else
            a(k,n) = 0.248*exp(-0.348*(log(((4*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./5)./data{1,k}(n,16))));
        end
 
        if data{1,k}(n,35) > 0
            a1(k,n) = max(0.2,a(k,n));
        else
            a1(k,n) = min(0.8,a(k,n));
        end
   
        if data{1,k}(n,35) == 0;
            F0(k,n) = 1;
        else
            F0(k,n) = max(0.5,min(1,(data{1,k}(n,11)./data{1,k}(n,35))));
        end

      KtoverGA(k,n) = max(0,1-exp(a1(k,n)*((data{1,k}(n,29)./data{1,k}(n,11))-F0(k,n))));
      KhGA(k,n) = min(1,(data{1,k}(n,24)+data{1,k}(n,29)/(data{1,k}(n,24)+data{1,k}(n,11))));
      L(k,n) = wavelength(data{1,k}(n,13),data{1,k}(n,24));
    
       if data{1,k}(n,29) > 0
           CGA(k,n) = 1.135*(data{1,k}(n,44)./data{1,k}(n,35))^0.65;
        elseif data{1,k}(n,29) == 0;
           CGA(k,n) = 1.135*(((9*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./10)./data{1,k}(n,35))^0.65;
        else
           CGA(k,n) = 1.135*(((4*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./5)./data{1,k}(n,35))^0.65;
       end  
       if data{1,k}(n,35) == 0;
            KthGA(k,n) = 0;
        else
            KthGA(k,n) = 1./(1+CGA(k,n)*(data{1,k}(n,11)/L(k,n))^0.5)^2;
       end  
       KtotGA(k,n) = min(1,((KtoverGA(k,n)^2 + KhGA(k,n)^2 * KthGA(k,n)^2)^0.5));
end
end
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtotGA(k,n)==0 
        KtotGA(k,n)=0.01;
end
end
end
%% Transmission Tomasicchio - D'Alessandro (2013) 
for k = 1:numel(data) % indices  
for n = 1:size(data{1,k})
        if data{1,k}(n,29) > 0
            a4(k,n) = 0.248*exp(-0.348*(log(data{1,k}(n,44)./data{1,k}(n,16))));
        elseif data{1,k}(n,29) == 0
            a4(k,n) = 0.248*exp(-0.348*(log(((9*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./10)./data{1,k}(n,16))));
        else
            a4(k,n) = 0.248*exp(-0.348*(log(((4*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./5)./data{1,k}(n,16))));
        end  
        if data{1,k}(n,35) > 0
            a5(k,n) = max(0.2,a4(k,n));
        else
            a5(k,n) = min(0.8,a4(k,n));
        end
        if data{1,k}(n,35) == 0;
            Rc0(k,n) = 1;
        else
            Rc0(k,n) = max(0.6,min(0.8,(data{1,k}(n,11)./data{1,k}(n,35))));
        end 
      KtoverTDA(k,n) = max(0,1-exp(a5(k,n)*((data{1,k}(n,29)./data{1,k}(n,11))-Rc0(k,n))));
      KhTDA(k,n) = min(0.8,(data{1,k}(n,24)+data{1,k}(n,29)/(data{1,k}(n,24)+data{1,k}(n,11))));   
      L(k,n) = wavelength(data{1,k}(n,13),data{1,k}(n,24)); 
        if data{1,k}(n,29) > 0
            CTDA(k,n) = 3.45*(data{1,k}(n,44)./data{1,k}(n,35))^0.65;
        elseif data{1,k}(n,29) == 0;
            CTDA(k,n) = 3.45*(((9*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./10)./data{1,k}(n,35))^0.65;
        else
            CTDA(k,n) = 3.45*(((4*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./5)./data{1,k}(n,35))^0.65;
        end
   
        if data{1,k}(n,35) == 0;
            KthTDA(k,n) = 0;
        else
            KthTDA(k,n) = 1./(1+CTDA(k,n)*(data{1,k}(n,11)/L(k,n))^0.5)^2;
        end  
       KtotTDA(k,n) = min(1,((KtoverTDA(k,n)^2 + KhTDA(k,n)^2 * KthTDA(k,n)^2)^0.5));   
end 
end
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtotTDA(k,n)==0 
        KtotTDA(k,n)=0.01;
end
end
end
%% Transmission Zhang 
%first for emerged, porous breakwaters, second is for submerged porous 
for k = 1:numel(data) % indices  
    for  n = 1:size(data{1,k})
         if data{1,k}(n,29)>= 0 & data{1,k}(n,35)~=0 & data{1,k}(n,29)./data{1,k}(n,11)<=1.0
             KtZ(k,n)=0.9*((1-1*data{1,k}(n,29)./data{1,k}(n,11))./(1+1*data{1,k}(n,29)./data{1,k}(n,11)))*exp(-0.18*(data{1,k}(n,44)./data{1,k}(n,11)))*(1-exp(-0.5*Ir(k,n)));   %emerged            
         elseif data{1,k}(n,29)<0 & data{1,k}(n,35)~=0
             KtZ(k,n)=0.50*((1-0.23*(data{1,k}(n,29)./data{1,k}(n,11)))./(1+0.23*(data{1,k}(n,29)./data{1,k}(n,11))))*exp(-0.18*(data{1,k}(n,44)./data{1,k}(n,11))); %submerged
         else 
             KtZ(k,n)=nan;
         end 
         if KtZ(k,n)>1
             KtZ(k,n)=nan;
         elseif KtZ(k,n)<0
             KtZ(k,n)=nan; 
         end 
    end 
end 
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtZ(k,n)==0 
        KtZ(k,n)=nan;
end
end
end
%% Transmission Sindhu 
%semi-empirical equation, Rc has to be negative 

for k = 1:numel(data) % indices       

for  n = 1:size(data{1,k})
    if data{1,k}(n,29)<0
        KtS(k,n)=(0.02*(-data{1,k}(n,29))./data{1,k}(n,44)+0.035*(data{1,k}(n,27)./data{1,k}(n,24)))*(data{1,k}(n,24)./data{1,k}(n,35)+0.25./sop(k,n)^0.5);
    else
        KtS(k,n)=nan;
    end           
    if KtS(k,n)<0
        KtS(k,n)= nan;
    end 
end
end 
% final check on zero-values
for k = 1:numel(data) % indices     
for n = 1:size(data{1,k})
    if KtS(k,n)==0 
        KtS(k,n)=nan;
end
end
end
%% Kt Kurdistani

for k = 1:numel(data) % indices  
  
   for n = 1:size(data{1,k})
       if data{1,k}(n,29) == 0
         Beff(k,n) = (9*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))))./10; 
       elseif data{1,k}(n,29)<0
          Beff(k,n) = (4*data{1,k}(n,44)+(data{1,k}(n,44)+((((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,48)))+(((data{1,k}(n,24))+(data{1,k}(n,29)))./(1./data{1,k}(n,49))))));
       else 
          Beff (k,n) = data{1,k}(n,44); 
       end 
      
   end 
   for n = 1:size(data{1,k})
       omega(k,n) = (1/2*pi)*tanh(2*pi*(data{1,k}(n,24))./L(k,n)); 
       fi(k,n) = ((data{1,k}(n,42))^0.5*data{1,k}(n,24)*Beff(k,n))./(data{1,k}(n,44)*data{1,k}(n,11));
       Leff(k,n)= L(k,n)./Beff(k,n);        
   end 

   for n = 1:size(data{1,k})
       if data{1,k}(n,35)~=0  &&  data{1,k}(n,29)<0       
          KtKU(k,n)= 0.923+(0.576.*log(0.428*(1+data{1,k}(n,48))^(0.042).*(1+(abs(data{1,k}(n,29))./data{1,k}(n,11)))^0.75.*(Beff(k,n)./data{1,k}(n,35))^(0.125).*omega(k,n)^(0.413).*Leff(k,n)^(0.39).*fi(k,n)^(-0.18)));
       else 
           KtKU(k,n)=nan;
       end 
           
   end 
end 

