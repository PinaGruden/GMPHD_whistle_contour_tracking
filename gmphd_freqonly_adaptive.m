function [Xk_m,XTag]=gmphd_freqonly_adaptive(Zset,parameters,models)
%gmphd_freqonly_adaptive is a function that performs Gaussian Mixture 
%Probability Hypothesis Density filtering on a set of noisy measurements  
%for the purpose of frequency contour extraction. Measurements consist of
%frequency information only. States consist of frequency and rate of cahnge
%of frequency information.

%Inputs:
% 1) Zset = cell array of measurements (frequency spectral peaks), each 
%cell contains collection of measurements at a particular time step
% 2) parameters = stucture containing GMPHD parameters:
% parameters.pdet = probability of detection (used in Update calculation)
% parameters.psurv = probability of survival
% parameters.clutter = Clutter intensity
% parameters.wth = weight threshold used in State estimation
% parameters.Jmax = max number of Gaussian componenets used in pruning
% parameters.U = pruning threshold used in Pruning & merging
% parameters.Tr = truncation threshold used in Pruning & merging
% 3) models = stucture containing system and measurement models parameters:
% models.F = system (dynamic) matrix (state transition matrix)
% models.Q = covariance matrix for the system noise
% models.H = measurement matrix
% models.R = covariance matrix for the measurement noise
% models.birthpdf = a matrix containing a pdf of birth objects

%Outputs:
% Xk_m = State estimates (means of the Gaussian components) 
% XTag = tags of the components

% Pina Gruden, Institute of Sound and Vibration Research, April 2015

%-------------------------------------------------------------------------
% For detailed explanation see: Gruden, P. and White, P. (2016). 
%-------------------------------------------------------------------------

T= size(Zset,2); %number of time steps 

%get the filter parameters
wth=parameters.wth; 
pdet = parameters.pdet; 
psurv = parameters.psurv;  
clutter=parameters.clutter;
Jmax=parameters.Jmax; 
U=parameters.U; 
Tr=parameters.Tr;

%get the system, observation and birth model parameters
F = models.F; % state transition matrix
Q = models.Q; % noise covariance matrix for the system noise 
H = models.H; % measurement matrix
R = models.R; % noise covariance matrix for the measurement noise
I = [1. 0.; 0. 1.]; % identity matrix
birthpdf=models.birthpdf;

%pre-allocate
Jk=cell(1,T);weights=Jk;means=Jk;cov=Jk;Tags=Jk;Xk_w=Jk;Xk_m=Jk;XTag=Jk;

for k=1:T

%% Initialization
Zk=Zset{k}; %measurement set at time k- needs to be xdim x Nmeasurements (each measurement separate column)
%Initialize existing targets from previous step parameters
if k==1 %for the first time step
    Jk_1=randi(10); %initial number of Gaussian components- draws a random value between 1 and 10
    wk_1=zeros(1,Jk_1);
    mk_1=zeros(2,Jk_1);
    Pk_1=cell(1,Jk_1);
    Tag=cell(1,Jk_1);
    for n=1:Jk_1
        wk_1(n)=1/Jk_1; % initial weights
        m=(30000-2000).*rand + 2000; %draw a frequency from uniform distribution between 2-30kHz
        mk_1(:,n)= [m;0]; % initial means
        Pk_1{n}=Q;%or do Pk_1=Q; % initial covariances
        Tag{n}=[sprintf('%d', k),'.',sprintf('%d', n)];%Assign unique tag to each Gaussian component.
    end
else
    Jk_1 = Jk{k-1}; %number of Gaussian components at previous step
    wk_1= weights{k-1}; %weights from previous step
    mk_1= means{k-1}; %means from previous step
    Pk_1=cov{1,k-1}; %covariances from previous step- Pk_1 is a cell array of Jk x 1 where
    %every cell corresponds to covariance for a particular Gaussian
    Tag=Tags{k-1};
end

% Initialize birth parameters

if ~isempty(Zk) 
    J_birth = numel(Zk);
    P_birth=cell(1,J_birth);
    w_birth=zeros(1,J_birth);
    for n=1:J_birth
        P_birth{n} = Q;
        v=birthpdf(1,:)-Zk(1,n);
        [~,l]=min(abs(v));%index of the frequency that is closest to the measurement
        val=birthpdf(2,l);
        w_birth(n)=val/J_birth;
    end
    
    Qb = zeros(2,2,J_birth);
    Mub = zeros(J_birth,2);
    mixb = ones(1,J_birth)./J_birth;
    for i = 1:J_birth
        if Zk(1,i)<5000
            mb = [6000 0];
        else
            mb = [Zk(1,i) 0];
        end
        Mub(i,:) =mb;
        Qb(:,:,i) = abs(.01*diag([Zk(1,i) Zk(1,i)]));
    end
    objb = gmdistribution(Mub,Qb,mixb);
    m_birth = random(objb,J_birth)' ;
    
else
    m_birth=[];w_birth=[];P_birth={};
    J_birth = [];
end

%% Prediction
% Using Kalman prediction equations

% Prediction for newborn targets (newborns do not change in this step)
T_birth=cell(1,J_birth);
for j=1:J_birth %number of Gaussian components in birth model
    T_birth{j}=[sprintf('%d', k),'.',sprintf('%d', j)];
end
% Prediction for existing targets (predict each Gaussian component) 
for j=1:Jk_1 %number of Gaussian components of existing targets
    wk_1(j)= psurv*wk_1(j);
    mk_1(:,j)= F*mk_1(:,j);
    Pk_1{j}= F*Pk_1{j}*F'+ Q; %a priori estimate of P
    %tags stay the same as previous time step
end 


%% Update
%Using Kalman update equations
%Update exsisitng and new born separately!(as per Ristic et al. (2012)).
%if there are no measurements there will be no newborn targets either

if ~isempty(Zk) %if there are measurements to update with
    %---------------Construct Update Variables-----------------------
    Hmk_1=H*mk_1;
    S=zeros(1,Jk_1);
    for j=1:Jk_1,S(j)=H*Pk_1{j}*H'+R;end
    Hm_birth=H*m_birth;
    Sb=zeros(1,J_birth);
    for j=1:J_birth,Sb(j)=H*P_birth{j}*H'+R;end
    qk=zeros(length(Zk),Jk_1);
    X=Zk';
    for j=1:Jk_1,qk(:,j)= mynormpdf(X,repmat(Hmk_1(:,j), size(X)),repmat(sqrt(S(j)),size(X)));end
    qkn=zeros(length(Zk),J_birth);
    for j=1:J_birth,qkn(:,j)= mynormpdf(X,repmat(Hm_birth(:,j),size(X)),repmat(sqrt(Sb(j)),size(X)));end
    
    %pre-allocate
    wk=zeros(1,Jk_1+Jk_1*size(Zk,2)+J_birth*size(Zk,2));
    mk=zeros(2,Jk_1+Jk_1*size(Zk,2)+J_birth*size(Zk,2));
    Pk=cell(1,Jk_1+Jk_1*size(Zk,2)+J_birth*size(Zk,2));
    Tag=[Tag,cell(1,Jk_1*size(Zk,2)+J_birth*size(Zk,2))];
    
    %---------------------EXSISTING targets--------------------------
    %first we assume we did not detect all targets- update weights,means, cov:
    for j=1:Jk_1
        wk(j)=(1-pdet)*wk_1(j); %1-pdet= probability of missed detection
        mk(:,j)=mk_1(:,j);
        Pk{j}=Pk_1{j};
    end
    
    indx=0;
    for n=1:size(Zk,2) %iterate through measurements
        for j=1:Jk_1
            indx= n*Jk_1+j;
            wk(indx)= pdet*wk_1(j)*qk(n,j)/(clutter+(w_birth*qkn(n,:)')+(pdet*(wk_1*qk(n,:)')));
            y = Zk(:,n) - (Hmk_1(:,j)); %error in prediction
            K = Pk_1{j}*H'/(S(j)); %Kalman gain (works the same as P*H'*inv(H*P*H'+R) but faster)
            mk(:,indx) = mk_1(:,j) + (K*y); %parameter estimates
            Pk{indx}=(I-(K*H))*Pk_1{j}; %Error covariance matrix
            Tag{indx}= Tag{j}; %assign same tag to each updated component as its associated predicted term
        end
    end
    
    %---------------------NEWBORN targets--------------------------
    for n=1:size(Zk,2) %iterate through measurements
        for j=1:J_birth
            indxn= indx+j;
            wk(indxn)= w_birth(j)*qkn(n,j)/(clutter+(w_birth*qkn(n,:)')+(pdet*(wk_1*qk(n,:)')));
            y = Zk(:,n) - (Hm_birth(:,j)); %error in prediction
            K = P_birth{j}*H'/(Sb(j)); %Kalman gain (works the same as P*H'*inv(H*P*H'+R) but faster)
            mk(:,indxn) = m_birth(:,j) + (K*y); %parameter estimates
            Pk{indxn}=(I-(K*H))*P_birth{j}; %Error covariance matrix
            Tag{indxn}= T_birth{j}; %assign same tag to each updated component as its associated predicted term
        end
        indx=indxn;
    end
       
else %There are no measurements to update with
   %-----------------------EXSISTING targets-----------------------
    %we assume we did not detect all targets
     wk=(1-pdet).*wk_1; %1-pdet= probability of missed detection
     mk=mk_1;
     Pk=Pk_1;
    %Tags stay the same
    % no newborn
end

%% Pruning and Merging
Indx = find(wk >= Tr); % Pruning- find targets with high enough weights

%define new weights:
wk_new = wk(Indx).*(sum(wk)/sum(wk(Indx)));
mk_new = mk(:,Indx);
Pk_new = Pk(Indx);
Tag_new = Tag(Indx);
In=1:length(wk_new);

l=0; 
%preallocate (to the max possible size, delete empty spaces later)
wk_l=zeros(1,length(In)); mk_l=zeros(2,length(In));Pk_l=cell(1,length(In));
Tag_l=Pk_l;
while isempty(In) == 0
    l=l+1;
    [~,j] = max(wk_new(In)); %find biggest of the pruned weights
    %now find indices of the means of Gaussian components that are within
    %distance U from the Gaussian component with highest weight
    %Use Mahalanobis distance
    difs=mk_new(:,In)-repmat(mk_new(:,In(j)),size(In));
    mahal=diag(difs'*(Pk_new{In(j)}\difs))';
    L=  In(mahal<=U);%get indices
    
    %calculate new weights,means and covariances for merged components
    wk_l(l)=sum(wk_new(L));
    mk_l(:,l)=(wk_new(L)*mk_new(:,L)')/wk_l(l); %the wk*mk' does the summation
    P_l=zeros(2,2);
    for i=1:length(L)
        P_l = P_l + wk_new(L(i))*(Pk_new{L(i)}+(mk_l(:,l) - mk_new(:,L(i)))*(mk_l(:,l) - mk_new(:,L(i)))');
    end
    Pk_l{l}= P_l/wk_l(l);
    
    Tag_l{l}=Tag_new{In(j)};
    
    % remove elements of L from I (I:= I\L)
    matches=zeros(size(L));
    for ix=1:length(L),matches(ix)=find(L(ix)==In);end
    In(matches)=[];
end
%now delete empty spaces
wk_l=wk_l(1:l);mk_l=mk_l(:,1:l);Pk_l=Pk_l(1:l);Tag_l=Tag_l(1:l);

%if number of Gaussians exceeds maximum allowed number of Gaussians- Jmax:
if numel(wk_l) > Jmax
    [~,ranking] = sort(wk_l,'descend'); %gives indices of ranking in descending order
    indx=ranking(1:Jmax); %indices that we want to keep 
    wk_l = wk_l(indx); 
    mk_l=mk_l(:,indx);
    Pk_l= Pk_l(indx);
    Tag_l=Tag_l(indx);
end

%save final Jk, wk, mk and Pk for this time step k
Jk{k} =  numel(wk_l); %number of Gaussians at time k (after pruning)
weights{k} = wk_l;
means{k} = mk_l;
cov{1,k} = Pk_l; 
Tags{k}=Tag_l;

%% State estimation 
%take Gaussians with weights above wth
if k>1
    i = find(wk_l > wth);
    w = wk_l(i);
    m = mk_l(:,i);
    Tagest = Tag_l(i);
   
    %select the branch with highest weight
    if numel(unique(Tagest)) ~= numel(Tagest)
        [idx,gridx]=group2cell(1:numel(Tagest),Tagest);
        Xk_wn=zeros(1,size(idx,1));Xk_mn=zeros(2,size(idx,1));XTagn=cell(1,size(idx,1));
        for n=1:size(idx,1)
            if numel(idx{n})>1
                [val,pos]=max(w(idx{n}));
                Xk_wn(n)=val;
                Xk_mn(:,n)=m(:,idx{n}(pos));
                XTagn{n}=Tagest{gridx(n)};
            else
                Xk_wn(n)=w(idx{n});
                Xk_mn(:,n)=m(:,idx{n});
                XTagn{n}=Tagest{gridx(n)};
            end
        end
        Xk_w{k} = Xk_wn;
        Xk_m{k} = Xk_mn;
        XTag{k} = XTagn;
    else
        Xk_w{k} = w;
        Xk_m{k} = m;
        XTag{k} = Tagest;
    end
          
else
    continue
end

end


end

function y = mynormpdf(x,mu,sigma)

    y = exp(-0.5 * ((x - mu)./sigma).^2) ./ (sqrt(2*pi) .* sigma);

end