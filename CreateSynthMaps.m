function [sigArray , t_Synth,X,Y,T] =CreateSynthMaps()

elecArray = reshape(1:100,10,10);
elecSpacing=5;
fs=30;
%% 1 . Create synthetic signal
elec = creaetSynthSig(fs);

%% 2. Assemble electrode array
x = 1:elecSpacing:elecSpacing*size(elecArray,1);
y = 1:elecSpacing:elecSpacing*size(elecArray,2);
[X,Y]=meshgrid(x,y);


%Create synthetic function --- T1, T2, T3 

% Linear activation 
a =(sqrt(2*60^2)/5)/(2*60);
T1 = ((-a.*X.^1) + (-a.*Y)) + 5;  %linear wave

% Source activations - Elliptical wavefronts 
a=0.02;   
b=0.1;   
a=0.02;b=0.04; 
kr=0.025;  
ftr=-30;  
X_off=25;
Y_off=25;

T2 = ftr*sin(-kr.*sqrt(a.*(X-X_off).^2+b.*(Y-Y_off).^2)) + 5; %contourf(T);colorbar

% Clahsing wavefronts 
a=0.05;   
a1=0.05;a1=a;
b=0.1;   
kr=0.025;  
ftr=-25;  
X_off1=0;
X_off2=60;
Y_off1=0;
Y_off2=60;
% T = ftr*sin(-kr.*sqrt(a.*X.^2+b.*Y.^2)); 
 
TLr = ftr*sin(-kr.*sqrt(a1.*(X+X_off1).^2+b.*(Y+Y_off1).^2));   TL=TLr; 
TUp = ftr*sin(-kr.*sqrt(a.*(X-X_off2).^2+b.*(Y-Y_off2).^2))-3;    TU=TUp;

TUp(TU>TL)=TLr(TU>TL);

T3 = TUp+5;

% Choice of synthetic signal 
T = T3;

t_Synth = (0:(max(max(T))*fs)+(600*fs)-1)*(1/fs);

%Place signals in apporiate leads
sigArray = zeros(length(T(:)),length(t_Synth));                 %zero sig
sin(2*pi*t_Synth*(12/60));

[~,pos] = min(gradient(elec));

%% 3. Defining patterns


refract = [0:40:580];
TT_pattrn = {T2,T2,T2};
PatternSwtch = [5 5 5];
n = 1;
for i = 1:length(PatternSwtch)
    for j = 1:PatternSwtch(i)
        TT{n} = TT_pattrn{i}+refract(n);
        n=n+1;
    end
    
end



for k = 1:length(refract)
    T=TT{k};
    
    %Placing signals in the the electrode position in time
    for i = 1:size(elecArray,1)
        for j = 1:size(elecArray,2)
            elecNo = elecArray(i,j);
            timeVal = T(i,j);
            fprintf('k %d, i %d, j %d \n',k,i,j)
            sigArray(elecNo,round(timeVal*fs)-pos:round(timeVal*fs)-pos+length(elec)-1)  ...
                = sigArray(elecNo,round(timeVal*fs)-pos:round(timeVal*fs)-pos+length(elec)-1)+ elec;
        end
    end

end


end



%% Function to create a synthetic signal 

% input -fs - sampling freq
% output - elec - electrode according to samp freq

function elec = creaetSynthSig(fs)

%load data
H = load ('sampleData.mat');
Sig = H.P_est2;

Sig = Sig - Sig(1);
fsSig = 30;
t = (0:length(Sig)-1) * (1/fsSig);

%Interpolate to 0 at edge
interpL = 0.5*fsSig - 5;
t_vec = (0:length(Sig)-1+interpL) * (1/fsSig);

tInt = [t t_vec(end)];
SigInt = [Sig' 0];
Sig = spline(tInt,SigInt,t_vec);

%Resample date accordingly
[p,q] = rat(fs/fsSig,0.0001);
elec = resample(Sig,p,q);

t25=tukeywin(length(elec),0.25);
elec =elec .*t25';

end





