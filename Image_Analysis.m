clc
clear
close all


% A
seq=dir('../Mass9');
seq=seq(5:end);
for i=1:numel(seq)
    clear x z
    dirname = ['../Mass9','/',seq(i).name];
    dirlist = dir([dirname '/*.jpg']);
    for c = 1:numel(dirlist)
        im = imread([dirname '/' dirlist(c).name]);
        imagesc(im)
        colormap ('gray')
        axis equal
        tmp = ginput(1);
        if numel(tmp) == 2
            x(c) = tmp(1);
            z(c) = tmp(2);
        else
            x(c) = NaN;
            z(c) = NaN;
        end
    end
    save(seq(i).name,'x','z')
end

%% B
% Save particle positions as Before Bouncing
namesBF={'BFflat1','BFflat1soft','BFflat2','BFflat2soft','BF' ...
    'slope21','BFslope21soft'...
    ,'BFslope212','BFslope212soft','BFslope33','BFslope33soft','B' ...
    'Fslope332','BFslope332soft'}';
% Save particle positions as After Bouncing
namesAF={'AFflat1','AFflat1soft','AFflat2','AFflat2soft','AF' ...
    'slope21','AFslope21soft'...
    'AFslope212','AFslope212soft','AFslope33','AFslope33soft','A' ...
    'Fslope332','AFslope332soft'}';

matstr=dir('*mat');
for q=1:numel(matstr)
    t=load(matstr(q).name);
    for i=5:numel(t.z) % i from 5 because I wanted to get rif of first
        % particles which are in same position in some falls
        if t.z(i)>t.z(i+1)
            sBF.(namesBF{q})=[t.x(1:i);t.z(1:i)];
            sAF.(namesAF{q})=[t.x(i+1:end);t.z(i+1:end)];
            break
        end
    end
end



% % TO TEST SEPARATED VARIABLES
plot( sBF.(namesBF{1})(1,:),sBF.(namesBF{1})(2,:),'r+')
hold on
plot( sAF.(namesAF{1})(1,:),sAF.(namesAF{1})(2,:),'b+')
axis equal
set(gca,'ydir','reverse')




%% C,D,E

d=zeros(2,1);

for i=1:12
    
    firstPart=numel(sBF.(namesBF{i})(1,:));
    secondPart=numel(sAF.(namesAF{i})(1,:));
    A=zeros(firstPart+secondPart,5);
    
    % First Half of Matrix
    A(1:firstPart,1)=1;
    A(1:firstPart,2)=0;
    A(1:firstPart,3)=1:firstPart;
    A(1:firstPart,4)=0;
    A(1:firstPart,5)=0.5*(1:firstPart).^2;
    
    % Second Half of Matrix
    A(firstPart+1:firstPart+secondPart,1)=0;
    A(firstPart+1:firstPart+secondPart,2)=1;
    
    A(firstPart+1:firstPart+secondPart,3)=0;
    A(firstPart+1:firstPart+secondPart,4)=firstPart+1:firstPart+ ...
        secondPart;
    A(firstPart+1:firstPart+secondPart,5)=0.5*[firstPart+1:first...
        Part+secondPart].^2;
    
    % Locations
    xinits= [sBF.(namesBF{i})(1,:)';sAF.(namesAF{i})(1,:)'];
    zinits= [sBF.(namesBF{i})(2,:)';sAF.(namesAF{i})(2,:)'];
    
    % New Matrix by combined vectors
    px=A\xinits;
    pz=A\zinits;
    
    
    seq=dir;
    seq=seq(5:2:end);
    
    clear im
    dirname = seq(i).name;
    dirlist = dir(['./' dirname '/*.jpg']);
    for c=1:numel(dirlist)
        im(c,:,:)=imread([dirname '/' dirlist(c).name]);
    end
    ims(i,:,:)=squeeze(min(im));
    figure(i)
    imagesc(squeeze(ims(i,:,:)))
    colormap ('gray')
    hold on
    xfall=   @(t) px(1)+px(3).* t + 0.5 * px(5)*t.^2;
    zfall=   @(t) pz(1)+pz(3).* t + 0.5 * pz(5)*t.^2;
    
    xbounce= @(t) px(2)+px(4).* t + 0.5 * px(5)*t.^2;
    zbounce= @(t) pz(2)+pz(4).* t + 0.5 * pz(5)*t.^2;
    
    % Distance
    
    xdist  =   @(t) xfall(t(1))-xbounce(t(2));
    zdist  =   @(t) zfall(t(1))-zbounce(t(2));
    findist=   @(t) sqrt(xdist(t)^2+zdist(t)^2);
    
    if     i==any([1 2])
        
        f=[10 10];
        d      =  fminsearch(findist,f);
    elseif i==4
        f=[10 10];
        d      =  fminsearch(findist,f);
    elseif any(i== [5  6  8  10 11 12])
        f=[15 15];
        d      =  fminsearch(findist,f);
    elseif any(i== [ 3 7 9])
        f=[20 20];
        d      =  fminsearch(findist,f);
    end
    
    % Velocities
    
    % Fall
    
    VxBF = px(3)+px(5)*d(1);
    VzBF = pz(3)+pz(5)*d(1);
    % After Rebound
    VxAF = px(4)+px(5)*d(2);
    VzAF = pz(4)+pz(5)*d(2);
    
    % Rt and Rn calculations
    if any(i== [9 10 11 12])
        diff1=(530-257);
        diff2=(360-534);
    elseif any(i==[1 2 3])
        diff1=(465-230);
        diff2=0;
    elseif (i== 4)
        diff1=(460-230);
        diff2=0;
        
    elseif any(i==[5 6 7 8])
        diff1=(530-260);
        diff2=(430-527);
        
    end
    tvec=[diff1;diff2];
    %abs(diff(numbersX));abs(diff(numbersZ))];%(467-228);(532-53
    
    n=[tvec(2);-tvec(1)];
    Rt=dot([VxAF;VzAF],tvec)/dot([VxBF,VzBF],tvec);
    Rn=-dot([VxAF;VzAF],n)/dot([VxBF,VzBF],n);
    
    
    Table(i,:)= [ {seq(i).name} ,Rt ,Rn];
    
    plot(sBF.(namesBF{i})(1,:),sBF.(namesBF{i})(2,:),'b+','LineWidth',1)
    plot(sAF.(namesAF{i})(1,:),sAF.(namesAF{i})(2,:),'b+','LineWidth',1)
    fplot(xfall,zfall,[2,25],'r','LineWidth',2)
    fplot(xbounce,zbounce,[0 ,45],'g','LineWidth',3)
    plot(xfall(d(1)),zfall(d(1)),'y*','LineWidth',3)
    
    FallName = {'flat_sphere1'...
        ;'flat_sphere1_soft';'flat_sphere2'...
        ;'flat_sphere2_soft';...
        'slope21_sphere1';'slope21_sphere1_soft';'slope21_sphere2';  'slope21_sphere2_soft';'slope33_sphere1';'slope33_sphere1_so',...
        'ft';...
        'slope33_sphere2';'slope33_sphere2_soft'};
    title(FallName{i});
    xlabel('x[px]')
    ylabel('z[px]')
    legend('particle center','particle center','FallTrace','Rebound Trace','intersection')
    set(gca,'ydir','reverse')
    %axis equal
end


%
Results=cell2mat(Table(:,2:3));
Rt=Results(:,1);
Rn=Results(:,2);
table(FallName,Rt,Rn)

%% INTERPERATIONS


% C

% xo= initial position on x direction
% zo= initial position on z direction
% Vx= initial velocity on x direction
% Vz= initial velocity on z direction
% gx= gravational constant on x direction
% gz= gravational constant on z direction
% Vx Vz xo xz will differ between fall and after bouncing
% movement.

% gx and gz will not change.

%% D

% Most remarkable thing in table is that there are 2
% negative results one for Rt and one for Rn which is not reasonable.
% The reason could be that the ball drop for fall part was not done 
% properly which caused negative values as well as such a high values 
%(over 1) in Rt coloumn. The rest of values look okay. And the difference
% between the fall in which we used a scarf and used no scarf is that Rn 
% values are kind of lower for the ones with scarf since it reduces 
% the effect of impact.
