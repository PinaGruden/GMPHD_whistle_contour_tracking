function [DT,Target,ind] = tracktarget(XTag,Xk_m,dt,tl)

%Inputs:
%XTag - all estimated tags from GMPHD filter,
%Xk_m - all estimated target means from GMPHD filter,
%tl - minimum length of target track (minumum number of time steps)

%Output
%Target - list of all targets identified by PHD filter. It is a struct of 
%means,timesteps where individual target occured and targets old label.
%ind- indices of all targets within Target list that are longer than
%specified minimum length tl

tagfull=XTag;
targmeanfull=Xk_m;

targetlist=[];

for k=1:size(tagfull,2)
    tagk=tagfull{k}; %tags at time step k
    targmeank=targmeanfull{k};
    if ~isempty(tagk)
        for n=1:numel(tagk)
            tagcmpr= strcmp(tagk(n),targetlist);
            if sum(tagcmpr)>0 %new target matches one of the old
                indx=find(tagcmpr==1); %index of the target that the tagk belongs to
                Target(indx).freq=[Target(indx).freq;targmeank(1,n)];
                Target(indx).time=[Target(indx).time; k*dt];
                Target(indx).label=[Target(indx).label;tagk(n)];
            else %new target does not match any of the old, it is new
                targetlist=[targetlist,tagk(n)];
                indx=length(targetlist);
                Target(indx).freq=targmeank(1,n);
                Target(indx).time= k*dt;
                Target(indx).label=tagk(n);
            end
        end
    else
        continue
    end
end

%find targets of certain length:
count=1;
for l=1:size(Target,2)
    if numel(Target(l).time)>tl
ind(count)=l;
count=count+1;
    else
        continue
    end
end

DT=Target(ind);
end

