close all
clc
clear

%-------------------------IMPORT DATA------------------------%

data=importdata('Updated_Database.csv');
database=data.data;
textdata=data.textdata;
textdata(1,:)=[];
TestID=textdata(:,2);
TestID=cell2table(TestID);

% | TestID                                  | TestNumber range 
%-------------------------------------------------------------
% | ds1 Seelig (1980)- smooth               | 1 - 13
% | ds2 Seelig (1980)- rubble mound         | 14 - 82
% | ds3 Allsop (1983                        | 83 - 103
% | ds4 Daemrich&Kahle (1985)- smooth       | 104 - 250
% | ds5 Daemrich&Kahle (1985)- smooth       | 250 - 446
% | ds6 Powell&Allsop (1985)                | 447 - 488
% | ds7 Delft M2090 (1985) - smooth         | 489 - 526
% | ds8 Delft M2090 (1985) -rubble mound    | 489 - 526
% | ds9 Ahrens (1987)                       | 527 - 727
% | ds10 Van der Meer (1988)                | 728 - 758
% | ds11 Delft H524 (1990)                  | 759 - 772
% | ds12 Daemen (1991)                      | 773 - 825
% | ds13 Delft H1872 (1994)                 | 826 - 864
% | ds14 Delft H2061 (1994)                 | 865 - 896
% | ds15 Delft H2014 (1994)                 | 897 - 907
% | ds16 Delft H1974 (1994)                 | 908 - 917
% | ds17 TU Delft (1997)                    | 918 - 1054
% | ds18 Taveira Pinto (1987)               | 1055 - 1607
% | ds19 Seebrook&Hall (1998)               | 1608 - 2239
% | ds20 Zannutigh (2000)                   | 2240 - 2295
% | ds21 Van der Meer (2000)                | 2296 - 2323
% | ds22 UCA (2001)                         | 2324 - 2376
% | ds23 Daemrich,Mai,Ohle (2001)           | 2377 - 2476
% | ds24 Kimura (2002)                      | 2477 - 2566
% | ds25 Aquareef (2002)                    | 2567 - 3629
% | ds26 UPC (2002)                         | 3630 - 3649
% | ds27 Wang (2002) - rubble mound         | 3650 - 3733
% | ds28 Wang (2002) - smooth               | 3734 - 3817
% | ds29 Melito&Melby (2002)                | 3818 - 3939
% | ds30 GWK (2002)                         | 3940 - 3984
% | ds31 Delft H4087  (2002)                | 3985 - 4004
% | ds32 Delft H4171 (2003)                 | 4005 - 4013
% | ds33 Ruol and Faedo (2004)              | 4014 - 4024 
% | ds34 Mori and Cappietti  (2005)         | 4025 - 4081
% | ds35 Kubowicz-Grajewska (2017)          | 4082 - 4129
% | ds36 Mahmoudof (2021)                   | 4130 - 4144

% load dataset
ds1=database(1:13,:); 
ds2=database(14:82,:);
ds3=database(83:103,:);
ds4=database(104:250,:);  
ds5=database(251:446,:); 
ds6=database(447:488,:); 
ds7=database(489:495,:);
ds8=database(496:526,:);
ds9=database(527:727,:); 
ds10=database(728:758,:); 
ds11=database(759:772,:); 
ds12=database(773:825,:); 
ds13=database(826:864,:); 
ds14=database(865:896,:); 
ds15=database(897:907,:); 
ds16=database(908:917,:);
ds17=database(918:1054,:); 
ds18=database(1055:1606,:); 
ds19=database(1607:2239,:);
ds20=database(2240:2295,:); 
ds21=database(2296:2323,:); 
ds22=database(2324:2376,:); 
ds23=database(2377:2476,:); 
ds24=database(2477:2566,:);
ds25=database(2567:3629,:);
ds26=database(3630:3649,:);
ds27=database(3650:3733,:); 
ds28=database(3734:3817,:);  
ds29=database(3818:3939,:);
ds30=database(3940:3984,:);
ds31=database(3985:4004,:);
ds32=database(4005:4013,:);
ds33=database(4014:4024,:); 
ds34=database(4025:4081,:); 
ds35=database(4082:4129,:); 
ds36=database(4130:4144,:); 

   
data = {ds1,ds2,ds3,ds4,ds5,ds6,ds7,ds8,ds9,ds10,ds11,ds12,ds13,ds14,ds15,ds16,ds17,ds18,ds19,ds20,ds21,ds22,ds23,ds24,ds25,ds26,ds27,ds28,ds29,ds30,ds31,ds32,ds33,ds34,ds35,ds36}; % cell array

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
 
%% Variables for plot 
for k = 1:numel(data) % indices  
 K1= KtC;      %ciria/vandermeer
 K2= KtVD;     % vdm-daemen
 K3= KtSea;    %Seabrook and Hall
 K4= KtCal;    %calabrese
 K5= KtB;      %buccino conceptual
 K6= KtDAG;    %d'agremond 
 K7= KtDEL;    %Briganti (DELOS)
 K8= KtV5;     %vdm smooth
 K9= KtotGA;   %goda ahrens
 K10= KtotTDA; %tom dales
 K11= KtZ;     %zhang
 K12=KtS;      %sindhu
 K13=KtKU;     %kurdistani
 
end 

%% Plot Ktobs-Ktcalc 

fig1=figure
x0=[0 1];
y0=[0 1];
t=tiledlayout('flow', 'TileSpacing','compact');
Marker= ['o','+','*','.','x','s','d','X','^','.','>','<','p','h','H','D','o','.','*','.','x','s','d','X','.','v','>','<','p','h','H','D','o','>','<','p'];
color= ["#FF0000";"#00FF00";"#034030";"#00FFFF";"#FF00FF";"#AF88FF";"#0072BD";"#D95319";"#EDB120";"#7E2F8E";"#77AC30";"#4DBEEE";"#A2142F";"#FF0000";"#00FF00";"#E2702A";"#00FFFF";"#2D8525";"#FF00FF";"#0072BD";"#D95319";"#EDB120";"#7E2F8E";"#77AC30";"#041291";"#A2142F";"#FF0000";"#00FF00";"#E2702A";"#40033F";"#FF00FF";"#FFC09D";"#0072BD";"#D95319";"#FF3D2A";"#7E2F8E"];
  %  "#808000";"#FF3D2A";"#800080"; "#1AAC8B";"#F19EEF";"#AF88FF";"#28685B";"#CDC108";"#CD7A08";"#CDCF54";"#EB7D2D";"#7E8288";"#FF9E12";"#96BDC7";"#CFA800";"#00532D";"#0FD1D5";"#638000";"#9B009D";"#FFC09D";"#CDC108";"#CD7A08";"#CDCF54";"#EB7D2D"]
    

for z=1:13
   nexttile  
   v=num2str(z);
   Kappa=eval("K"+z);

for dsi=1:k
      plot(data{1,dsi}(:,21),Kappa(dsi,1:size(data{1,dsi})), Marker(dsi),'MarkerEdgeColor',color(dsi));
      hold on 
      grid on 
      grid minor
end
plot(x0,y0,'k--')

    
if z==1
    title('Van der Meer (1990)', 'Fontsize',11);
elseif z==2
    title ('Van der Meer - Daemen (1994)', 'Fontsize',11);
elseif z==3
    title ('Seabrook and Hall (2000)', 'Fontsize',11);           
elseif z==4
    title ('Calabrese (2002)', 'Fontsize',11);
elseif z==5
    title ('Conceptual Buccino (2007)', 'Fontsize',11);
elseif z==6
    title ('D''Agremond (1996)', 'Fontsize',11);
elseif z==7
    title ('Briganti (DELOS 2003)', 'Fontsize', 11);
elseif z==8
    title ('Van der Meer - smooth/oblique (2005)', 'Fontsize', 11);
elseif z==9
    title ('Goda & Ahrens (2008)', 'Fontsize',11);
elseif z==10
    title ('Tomasicchio & D''Alessandro (2013)', 'Fontsize',11);
elseif z==11
    title ('Zhang (2014)', 'Fontsize',11);
elseif z==12
    title ('Sindhu (2015)', 'Fontsize',11);
elseif z==13
    title ('Kurdistani (2021)', 'Fontsize',11);
end
           
axis equal
xlim([0 1]);
ylim([0 1]);
end 
  
%Common characteristic tiledlayout
title(t,'Kt')
xlabel(t,{'kt observed'},'FontSize',14);
ylabel(t,{'kt calculated'},'FontSize',14);
axis equal
legend('Seelig (1980)- smooth','Seelig (1980)- rubble mound','Allsop (1983)',...
    'Daemrich and Kahle (1985)- smooth','Daemrich and Kahle (1985)- rubble mound',...
    'Powell and Allsop (1985)', 'Delft M2090 (1985) - smooth','Delft M2090 (1985) - rubble mound',...
    'Ahrens (1987)','Van der Meer (1988)','Delft H524 (1990)','Daemen (1991)','Delft H1872 (1994)',...
    'Delft H2061 (1994)','Delft H2014 (1994)','Delft H1974 (1994)','TU Delft (1997)','Taveira Pinto (1987)',...
    'Seebrook and Hall (1998)','Zannutigh (2000)','Van der Meer (2000)', 'UCA (2001)','Daemrich,Mai,Ohle (2001)',...
    'Kimura (2002)','Aquareef (2002)','UPC (2002)','Wang (2002)- rubble mound','Wang (2002)- smooth',...
    'Melito and Melby (2002)','GWK (2002)','Delft H4087  (2002)','Delft H4171 (2003)','Ruol and Faedo (2004)',...
    'Mori and Cappietti  (2005)','Kubowicz and Grajewska (2017)','Mahmoudof (2021)','FontSize',13,'Interpreter','latex'); 
lgd=legend 
lgd.Layout.Tile = 'east';
hold on

set(gcf,'position',[400,400,2000,2000])






%% write results

for j=1:size(K1,1)  
Kt1_C(:,j)= transpose(K1(j,:));
Kt2_VD(:,j)= transpose(K2(j,:));
Kt3_Sea(:,j)= transpose(K3(j,:));
Kt4_Cal(:,j)= transpose(K4(j,:));
Kt5_B(:,j)= transpose(K5(j,:));
Kt6_DAG(:,j)= transpose(K6(j,:));
Kt7_DEL(:,j)= transpose(K7(j,:));
Kt8_V5(:,j)= transpose(K8(j,:));
Kt9_GA(:,j)= transpose(K9(j,:));
Kt10_TDA(:,j)= transpose(K10(j,:));
Kt11_Z(:,j)= transpose(K11(j,:));
Kt12_S(:,j)= transpose(K12(j,:));
Kt13_KU(:,j)= transpose(K13(j,:));
end

%%

Kt1_C = Kt1_C(:); 
Kt1_C= nonzeros(Kt1_C);
Kt2_VD = Kt2_VD(:); 
Kt2_VD= nonzeros(Kt2_VD);
Kt3_Sea=Kt3_Sea(:);
Kt3_Sea=nonzeros(Kt3_Sea);
Kt4_Cal = Kt4_Cal(:); 
Kt4_Cal= nonzeros(Kt4_Cal);
Kt5_B = Kt5_B(:);
Kt5_B= nonzeros(Kt5_B);
Kt6_DAG = Kt6_DAG(:); 
Kt6_DAG= nonzeros(Kt6_DAG);
Kt7_DEL = Kt7_DEL(:); 
Kt7_DEL= nonzeros(Kt7_DEL);
Kt8_V5 = Kt8_V5(:); 
Kt8_V5= nonzeros(Kt8_V5);
Kt9_GA = Kt9_GA(:);
Kt9_GA= nonzeros(Kt9_GA);
Kt10_TDA = Kt10_TDA(:);
Kt10_TDA= nonzeros(Kt10_TDA);
Kt11_Z = Kt11_Z(:);
Kt11_Z= nonzeros(Kt11_Z);
Kt12_S = Kt12_S(:);
Kt12_S= nonzeros(Kt12_S);
Kt13_KU = Kt13_KU(:);
Kt13_KU= nonzeros(Kt13_KU);


csv=cat(2,TestID,array2table(database(:,21)),array2table(Kt1_C),array2table(Kt2_VD),array2table(Kt3_Sea), array2table(Kt4_Cal),array2table(Kt5_B),array2table(Kt6_DAG),array2table(Kt7_DEL),array2table(Kt8_V5),array2table(Kt9_GA),array2table(Kt10_TDA),array2table(Kt11_Z),array2table(Kt12_S),array2table(Kt13_KU));

writetable(csv,'myDataFile.csv');

